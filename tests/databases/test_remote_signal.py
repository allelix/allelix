# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for fetch_remote_text + head_request_headers freshness primitives."""

from __future__ import annotations

import http.server
import socket
import threading
from typing import ClassVar

import pytest

from allelix.databases.manager import (
    USER_AGENT,
    fetch_remote_text,
    head_request_headers,
)


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class _MD5Handler(http.server.BaseHTTPRequestHandler):
    """Serves an NCBI-style .md5 body to GET; records UA."""

    captured_ua: ClassVar[str] = ""
    body: ClassVar[bytes] = b"abcdef0123456789  clinvar.vcf.gz\n"

    def do_GET(self) -> None:
        type(self).captured_ua = self.headers.get("User-Agent", "")
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.body)))
        self.end_headers()
        self.wfile.write(self.body)

    def log_message(self, *args, **kwargs) -> None:
        pass


class _HeadHandler(http.server.BaseHTTPRequestHandler):
    """Responds to HEAD with configurable Last-Modified / ETag."""

    last_modified: ClassVar[str] = "Wed, 21 Oct 2025 07:28:00 GMT"
    etag: ClassVar[str | None] = None

    def do_HEAD(self) -> None:
        self.send_response(200)
        if self.last_modified:
            self.send_header("Last-Modified", self.last_modified)
        if self.etag:
            self.send_header("ETag", self.etag)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_GET(self) -> None:
        self.do_HEAD()

    def log_message(self, *args, **kwargs) -> None:
        pass


def _serve(handler):
    port = _free_port()
    server = http.server.HTTPServer(("127.0.0.1", port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, thread, f"http://127.0.0.1:{port}/x"


@pytest.fixture
def md5_server():
    server, thread, url = _serve(_MD5Handler)
    try:
        yield url
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


@pytest.fixture
def head_server():
    server, thread, url = _serve(_HeadHandler)
    try:
        yield url
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


class TestFetchRemoteText:
    def test_returns_body(self, md5_server: str):
        body = fetch_remote_text(md5_server)
        assert body is not None
        assert "clinvar.vcf.gz" in body

    def test_uses_allelix_user_agent(self, md5_server: str):
        fetch_remote_text(md5_server)
        assert _MD5Handler.captured_ua == USER_AGENT

    def test_returns_none_on_connection_refused(self):
        assert fetch_remote_text(f"http://127.0.0.1:{_free_port()}/missing") is None

    def test_returns_none_on_invalid_url(self):
        assert fetch_remote_text("not-a-valid-url") is None


class TestHeadRequestHeaders:
    def test_returns_last_modified(self, head_server: str):
        _HeadHandler.last_modified = "Mon, 01 Jan 2026 00:00:00 GMT"
        _HeadHandler.etag = None
        headers = head_request_headers(head_server)
        assert headers is not None
        assert headers.get("Last-Modified") == "Mon, 01 Jan 2026 00:00:00 GMT"

    def test_returns_etag_when_present(self, head_server: str):
        _HeadHandler.last_modified = ""
        _HeadHandler.etag = '"abc123"'
        headers = head_request_headers(head_server)
        assert headers is not None
        assert headers.get("ETag") == '"abc123"'

    def test_returns_none_on_connection_refused(self):
        assert head_request_headers(f"http://127.0.0.1:{_free_port()}/x") is None

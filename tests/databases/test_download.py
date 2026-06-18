# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for atomic, timeout-bounded downloads (M-2)."""

from __future__ import annotations

import http.server
import socket
import threading
import time
from typing import TYPE_CHECKING, ClassVar

import pytest

from allelix.databases import manager
from allelix.databases.manager import USER_AGENT, download, verify_file_hash

if TYPE_CHECKING:
    from pathlib import Path


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class _CapturingHandler(http.server.BaseHTTPRequestHandler):
    """Records request headers; serves a small fixed payload."""

    captured_headers: ClassVar[dict[str, str]] = {}
    payload: ClassVar[bytes] = b"hello-allelix-download-test"

    def do_GET(self) -> None:
        type(self).captured_headers = dict(self.headers.items())
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, *args, **kwargs) -> None:  # silence test output
        pass


@pytest.fixture
def http_server():
    port = _free_port()
    server = http.server.HTTPServer(("127.0.0.1", port), _CapturingHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{port}/data"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


class TestDownload:
    def test_writes_payload_to_dest(self, http_server: str, tmp_path: Path):
        dest = tmp_path / "out.bin"
        download(http_server, dest)
        assert dest.read_bytes() == _CapturingHandler.payload

    def test_sends_allelix_user_agent(self, http_server: str, tmp_path: Path):
        dest = tmp_path / "out.bin"
        download(http_server, dest)
        ua = _CapturingHandler.captured_headers.get("User-Agent", "")
        assert ua == USER_AGENT
        assert "allelix" in ua

    def test_no_part_file_left_after_success(self, http_server: str, tmp_path: Path):
        dest = tmp_path / "out.bin"
        download(http_server, dest)
        assert not (tmp_path / "out.bin.part").exists()

    def test_part_file_cleaned_up_on_failure(self, tmp_path: Path):
        dest = tmp_path / "out.bin"
        # Bad URL — connection refused. Should clean up the .part file.
        with pytest.raises(OSError):
            download(f"http://127.0.0.1:{_free_port()}/", dest)
        assert not (tmp_path / "out.bin.part").exists()
        assert not dest.exists()


class _TruncatingHandler(http.server.BaseHTTPRequestHandler):
    """Lies about Content-Length to simulate a truncated upstream.

    #79 missing-branch coverage: the canonical _CapturingHandler always
    matches header bytes to delivered bytes, so the size-mismatch check
    at databases/manager.py:117-124 never fires under tests. This handler
    advertises a longer payload than it actually delivers — exercising
    the OSError("Download truncated") branch.
    """

    advertised: ClassVar[int] = 1000  # claim 1000 bytes
    actual_payload: ClassVar[bytes] = b"partial-allelix-payload"  # 23 bytes

    def do_GET(self) -> None:
        self.send_response(200)
        self.send_header("Content-Length", str(self.advertised))
        self.end_headers()
        # Deliver fewer bytes than advertised, then close — simulates a
        # mid-transfer connection drop after a CDN sent the header.
        self.wfile.write(self.actual_payload)

    def log_message(self, *args, **kwargs) -> None:
        pass


@pytest.fixture
def truncating_server():
    port = _free_port()
    server = http.server.HTTPServer(("127.0.0.1", port), _TruncatingHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{port}/data"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


class TestDownloadTruncation:
    """#79 missing-branch coverage: server lies about Content-Length, the
    download function must detect the size mismatch and refuse the cache."""

    def test_truncated_response_raises_download_truncated(
        self, truncating_server: str, tmp_path: Path
    ) -> None:
        dest = tmp_path / "truncated.bin"
        with pytest.raises(OSError, match="Download truncated"):
            download(truncating_server, dest)

    def test_truncated_response_reports_both_sizes(
        self, truncating_server: str, tmp_path: Path
    ) -> None:
        """The error message should carry the expected and actual byte counts
        so the operator can tell what was advertised vs delivered."""
        dest = tmp_path / "truncated.bin"
        with pytest.raises(OSError) as exc:
            download(truncating_server, dest)
        msg = str(exc.value)
        assert "1,000" in msg or "1000" in msg  # advertised
        assert "23" in msg  # actual delivered

    def test_truncated_response_leaves_no_part_file(
        self, truncating_server: str, tmp_path: Path
    ) -> None:
        """Truncation must NOT leave a partial .part file or a half-written dest."""
        dest = tmp_path / "truncated.bin"
        with pytest.raises(OSError):
            download(truncating_server, dest)
        assert not (tmp_path / "truncated.bin.part").exists()
        assert not dest.exists()


class TestVerifyFileHash:
    """Integrity verification for downloaded database files."""

    def test_correct_hash_passes(self, tmp_path: Path):
        payload = b"allelix-integrity-test-payload"
        f = tmp_path / "good.bin"
        f.write_bytes(payload)
        import hashlib

        expected = hashlib.sha256(payload).hexdigest()
        verify_file_hash(f, "sha256", expected)
        assert f.exists()

    def test_flipped_byte_fails(self, tmp_path: Path):
        payload = bytearray(b"allelix-integrity-test-payload")
        f = tmp_path / "corrupt.bin"
        f.write_bytes(bytes(payload))
        import hashlib

        expected = hashlib.sha256(bytes(payload)).hexdigest()
        payload[0] ^= 0xFF
        f.write_bytes(bytes(payload))
        with pytest.raises(OSError, match="Integrity check failed"):
            verify_file_hash(f, "sha256", expected)
        assert not f.exists()

    def test_md5_verification(self, tmp_path: Path):
        payload = b"clinvar-vcf-test-data"
        f = tmp_path / "clinvar.vcf.gz"
        f.write_bytes(payload)
        import hashlib

        expected = hashlib.md5(payload).hexdigest()
        verify_file_hash(f, "md5", expected)
        assert f.exists()

    def test_md5_mismatch_deletes_file(self, tmp_path: Path):
        f = tmp_path / "clinvar.vcf.gz"
        f.write_bytes(b"real data")
        with pytest.raises(OSError, match="Integrity check failed"):
            verify_file_hash(f, "md5", "0" * 32)
        assert not f.exists()


class _SlowHandler(http.server.BaseHTTPRequestHandler):
    """Server that sleeps before responding, to trigger a download timeout."""

    sleep_seconds: ClassVar[float] = 3.0

    def do_GET(self) -> None:
        time.sleep(self.sleep_seconds)
        self.send_response(200)
        self.end_headers()

    def log_message(self, *args, **kwargs) -> None:
        pass


class TestDownloadTimeout:
    """m-timeout: a hung server must trip the configured timeout, not block forever."""

    def test_slow_response_triggers_timeout(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr(manager, "DOWNLOAD_TIMEOUT_SECONDS", 0.5)
        port = _free_port()
        server = http.server.HTTPServer(("127.0.0.1", port), _SlowHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with pytest.raises((TimeoutError, OSError)):
                download(f"http://127.0.0.1:{port}/", tmp_path / "out.bin")
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=2)
        # Cleanup invariant from M-2 still holds under timeout.
        assert not (tmp_path / "out.bin.part").exists()

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Automatic pre-analysis refresh follows the same HTTP evidence as db update."""

from __future__ import annotations

import io
import os
import time

import pytest
from rich.console import Console

from allelix.cli import _helpers
from tests.test_http_freshness_fields import (
    ETAG,
    LM,
    NEW_ETAG,
    NEW_LM,
    _install,
    _legacy,
    _payload_state,
)
from tests.test_http_freshness_fields import remote as remote


def _age(data, source):
    stale = time.time() - _helpers._STALENESS_SECONDS - 3600
    path = data / (source + ".sqlite")
    os.utime(path, (stale, stale))


def _refresh(data, remote, monkeypatch):
    stream = io.StringIO()
    monkeypatch.setattr(_helpers, "get_annotators", lambda folder, **kwargs: [remote.cls(folder)])
    monkeypatch.setattr(_helpers, "console", Console(file=stream, width=160, color_system=None))
    _helpers._maybe_refresh_databases(data)
    return " ".join(stream.getvalue().lower().split())


@pytest.mark.parametrize(
    "case",
    [
        "lm-to-both",
        "etag-to-lm-disjoint",
        "both-to-lm",
        "etag-precedes-changed-lm",
        "ETag",
        "Etag",
        "etag",
    ],
)
def test_stale_cache_does_not_download_for_header_presentation_changes(
    tmp_path, remote, monkeypatch, case
):
    if case in ("both-to-lm", "etag-precedes-changed-lm"):
        data = _install(tmp_path, remote)
    else:
        data = _legacy(tmp_path, remote, "lm" if case == "lm-to-both" else "etag")
    if case in ("etag-to-lm-disjoint", "both-to-lm"):
        remote.headers = {"Last-Modified": LM}
    elif case == "etag-precedes-changed-lm":
        remote.headers = {"ETag": ETAG, "Last-Modified": NEW_LM}
    elif case in ("ETag", "Etag", "etag"):
        remote.headers = {case: ETAG, "last-modified": LM}
    _age(data, remote.source)
    before = _payload_state(data, remote.source)
    output = _refresh(data, remote, monkeypatch)
    assert remote.heads == 1, "A stale cache must actually probe the remote headers"
    assert remote.downloads == remote.cpic_fetches == 0, output
    assert _payload_state(data, remote.source) == before
    assert "updating " not in output
    if case != "etag-to-lm-disjoint":
        assert "when online" not in output and "can't verify" not in output, output


@pytest.mark.parametrize("kind", ["etag", "lm"])
def test_stale_cache_with_changed_shared_field_really_refreshes(
    tmp_path, remote, monkeypatch, kind
):
    data = _legacy(tmp_path, remote, kind)
    remote.headers = {"ETag": NEW_ETAG} if kind == "etag" else {"Last-Modified": NEW_LM}
    remote.revision = 1
    _age(data, remote.source)
    before = _payload_state(data, remote.source)
    output = _refresh(data, remote, monkeypatch)
    assert remote.downloads == 1, output
    assert _payload_state(data, remote.source) != before
    assert "updated" in output, output
    # A subsequent stale check reads the replacement's stored evidence.
    remote.reset()
    _age(data, remote.source)
    output = _refresh(data, remote, monkeypatch)
    assert remote.heads == 1 and remote.downloads == 0, output


def test_fresh_cache_does_not_probe_or_download(tmp_path, remote, monkeypatch):
    data = _legacy(tmp_path, remote, "etag")
    before = _payload_state(data, remote.source)
    _refresh(data, remote, monkeypatch)
    assert remote.heads == remote.downloads == remote.cpic_fetches == 0
    assert _payload_state(data, remote.source) == before

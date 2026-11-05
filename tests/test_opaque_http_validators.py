# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A quoted HTTP entity tag is opaque, even when it resembles stamp fields."""

from __future__ import annotations

import pytest

from tests.test_automatic_http_refresh import _age, _refresh
from tests.test_http_freshness_fields import LM, NEW_LM, _install, _payload_state, _update
from tests.test_http_freshness_fields import remote as remote


def _run(data, remote, monkeypatch, path):
    if path == "automatic":
        _age(data, remote.source)
        return _refresh(data, remote, monkeypatch)
    return _update(data, remote)


@pytest.mark.parametrize("path", ["db-update", "automatic"])
@pytest.mark.parametrize("weak", ["", "W/"])
@pytest.mark.parametrize(
    ("old", "new"),
    [
        ('"release|lm:first"', '"release|lm:second"'),
        ('"first|etag:tail"', '"second|etag:tail"'),
        ('"release|cpic:first"', '"release|cpic:second"'),
    ],
)
def test_changed_opaque_etag_refreshes_real_content(
    tmp_path, remote, monkeypatch, path, weak, old, new
):
    # RFC 9110 section 8.8.3 permits pipe and colon inside the opaque tag.
    remote.headers = {"ETag": weak + old, "Last-Modified": LM}
    data = _install(tmp_path, remote)
    before = _payload_state(data, remote.source)
    remote.headers = {"etag": weak + new, "last-modified": LM}
    remote.revision = 1
    output = _run(data, remote, monkeypatch, path)
    assert remote.downloads == 1, output
    assert _payload_state(data, remote.source) != before
    remote.reset()
    _run(data, remote, monkeypatch, path)
    assert remote.downloads == remote.cpic_fetches == 0


@pytest.mark.parametrize("path", ["db-update", "automatic"])
@pytest.mark.parametrize("weak", ["", "W/"])
def test_opaque_etag_does_not_invent_a_last_modified_field(
    tmp_path, remote, monkeypatch, path, weak
):
    remote.headers = {"ETag": weak + '"release|lm:opaque-data"'}
    data = _install(tmp_path, remote)
    before = _payload_state(data, remote.source)
    # The real responses share no HTTP field. Text inside ETag is not one.
    remote.headers = {"Last-Modified": LM}
    output = _run(data, remote, monkeypatch, path)
    assert remote.downloads == remote.cpic_fetches == 0, output
    assert _payload_state(data, remote.source) == before


@pytest.mark.parametrize("path", ["db-update", "automatic"])
@pytest.mark.parametrize("weak", ["", "W/"])
def test_unchanged_opaque_etag_still_takes_precedence_over_modified_date(
    tmp_path, remote, monkeypatch, path, weak
):
    tag = weak + '"release|lm:data|etag:data|cpic:data"'
    remote.headers = {"ETag": tag, "Last-Modified": LM}
    data = _install(tmp_path, remote)
    before = _payload_state(data, remote.source)
    remote.headers = {"etag": tag, "last-modified": NEW_LM}
    output = _run(data, remote, monkeypatch, path)
    assert remote.downloads == remote.cpic_fetches == 0, output
    assert _payload_state(data, remote.source) == before


@pytest.mark.parametrize("remote", ["pharmgkb"], indirect=True)
@pytest.mark.parametrize("path", ["db-update", "automatic"])
@pytest.mark.parametrize("weak", ["", "W/"])
def test_missing_cpic_does_not_truncate_the_opaque_etag(tmp_path, remote, monkeypatch, path, weak):
    remote.headers = {"ETag": weak + '"release|cpic:opaque-data"', "Last-Modified": LM}
    remote.cpic = None
    monkeypatch.setattr(remote.module, "fetch_cpic_allele_functions", lambda: {})
    data = _install(tmp_path, remote)
    before = _payload_state(data, remote.source)
    output = _run(data, remote, monkeypatch, path)
    assert remote.downloads == 0, output
    assert _payload_state(data, remote.source) == before

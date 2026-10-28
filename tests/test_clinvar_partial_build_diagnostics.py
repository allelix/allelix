# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Diagnostics retain insertion progress when a source stream fails mid-build."""

import pytest
from test_clinvar_build_diagnostics import _log_build, _metric, _old_cache, _sources

from allelix.databases import manager


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("compressed", [False, True])
def test_partial_insertion_is_logged_before_rollback(tmp_path, monkeypatch, build, compressed):
    vs, ss = _sources(tmp_path, build, compressed)
    cache = tmp_path / "clinvar.sqlite"
    before = _old_cache(cache)
    log = tmp_path / "build.log"
    original = manager.iter_clinvar_tsv_records

    def interrupted(*args, **kwargs):
        records = original(*args, **kwargs)
        try:
            yield next(records)
            raise RuntimeError("synthetic source interruption after first insert")
        finally:
            records.close()

    monkeypatch.setattr(manager, "iter_clinvar_tsv_records", interrupted)
    monkeypatch.setattr(manager, "INSERT_BATCH_SIZE", 1)

    def call():
        with pytest.raises(RuntimeError, match="synthetic source interruption"):
            manager.load_clinvar_tsv(vs, ss, cache, build)

    _log_build(log, call)
    assert cache.read_bytes() == before
    assert not cache.with_name(cache.name + ".partial").exists()
    _metric(log.read_text(), ("records yielded",), 1)
    _metric(log.read_text(), ("rows inserted",), 1)

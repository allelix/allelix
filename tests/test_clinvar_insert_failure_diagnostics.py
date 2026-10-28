# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An interrupted SQL batch retains the count of writes already performed."""

import sqlite3

import pytest
from test_clinvar_build_diagnostics import _log_build, _metric, _old_cache, _sources

from allelix.databases import manager


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("compressed", [False, True])
def test_sql_batch_failure_keeps_partial_insert_count(tmp_path, monkeypatch, build, compressed):
    vs, ss = _sources(tmp_path, build, compressed)
    cache = tmp_path / "clinvar.sqlite"
    before = _old_cache(cache)
    log = tmp_path / "build.log"
    connect = sqlite3.connect

    class InterruptedConnection(sqlite3.Connection):
        def executemany(self, sql, parameters):
            if sql == manager._CLINVAR_INSERT:
                rows = list(parameters)
                super().executemany(sql, rows[:1])
                raise sqlite3.OperationalError("synthetic interrupted insert batch")
            return super().executemany(sql, parameters)

    def interrupted_connect(*args, **kwargs):
        return connect(*args, **kwargs, factory=InterruptedConnection)

    monkeypatch.setattr(manager.sqlite3, "connect", interrupted_connect)

    def call():
        with pytest.raises(sqlite3.OperationalError, match="synthetic interrupted insert"):
            manager.load_clinvar_tsv(vs, ss, cache, build)

    _log_build(log, call)
    assert cache.read_bytes() == before
    assert not cache.with_name(cache.name + ".partial").exists()
    _metric(log.read_text(), ("records yielded",), 2)
    _metric(log.read_text(), ("rows inserted",), 1)

"""A late recovery refusal must disclose the ClinVar build already published."""

import os
import sqlite3
from pathlib import Path

import pytest
from test_cache_publication_journal_class import _old_state, _read
from test_clinvar_integrity import BUILDS, _network, _old_caches

from allelix.annotators import clinvar
from allelix.databases import manager


def test_second_build_recovery_refusal_reports_the_actual_partial_update(tmp_path, monkeypatch):
    data = tmp_path / "data"
    _old_caches(data, dict.fromkeys(BUILDS, 2))
    _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    first = data / clinvar.clinvar_db_filename("GRCh37")
    second = data / clinvar.clinvar_db_filename("GRCh38")
    old_second_inode = second.stat().st_ino
    original_replace, original_connect = os.replace, sqlite3.connect
    injected = []

    def replace(src, dst, *args, **kwargs):
        original_replace(src, dst, *args, **kwargs)
        if Path(dst) == first:
            # Both pre-checks already passed. Another user can leave recovery
            # bytes before the later build's individual publication begins.
            _old_state(second, "hot")
            injected.append(True)

    def connect(path, *args, **kwargs):
        if injected and str(path) == str(second):
            raise sqlite3.OperationalError("synthetic late recovery refusal")
        return original_connect(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(os, "replace", replace)
        patch.setattr(sqlite3, "connect", connect)
        with clinvar.ClinVarAnnotator(data) as annotator, pytest.raises(RuntimeError) as error:
            annotator.setup()

    assert injected, "The first new build actually reached publication"
    first_info = manager.get_database_info(first, "clinvar.GRCh37")
    assert first_info is not None and first_info["record_count"] == 3
    assert second.stat().st_ino == old_second_inode
    assert "publication_old_only" in _read(second)
    second_info = manager.get_database_info(second, "clinvar.GRCh38")
    assert second_info is not None and second_info["record_count"] == 2
    assert "already updated: GRCh37" in str(error.value), (
        "A recovery refusal after one successful publication must not hide that update"
    )
    assert "GRCh38" in str(error.value)
    assert "old cache kept" in str(error.value)
    assert "synthetic late recovery refusal" in str(error.value)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A one-build ClinVar install is freshness-checked and shown like a two-build one.

`db update --build grch37` (or grch38) is a supported install: ClinVar is
then ready for input on that build. Its cache must be probed and refreshed
when stale before an analysis, and `db status` must show it ready.

The first test is reviewer C's (outside review of 2.6.0 RC4, finding R4-03),
unchanged but for line wrapping. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import io
import os
import sqlite3
import time
import urllib.request
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.annotators import clinvar as clinvar_mod
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.cli import _helpers, main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    """No test may reach the network."""

    def refuse(*args, **kwargs):
        raise AssertionError("network use in an offline regression test")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _clinvar_cache(path: Path, build: str, signal: str) -> None:
    """A ready, current ClinVar cache for ``build`` with one row."""
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.execute(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES ('rs1801133', '1', 11856378, 'G', 'A', 'Pathogenic', 'x', 'MTHFR', "
            "'criteria provided, single submitter', 1)"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://', '20260101', '2026-01-01', 1, ?, ?)",
            (f"clinvar.{build}", signal, f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def test_r03_single_build_clinvar_is_freshness_checked_before_analysis(tmp_path, monkeypatch):
    """README Privacy / ADR-0012: the switched-on sources get a freshness check before
    analysis. `db update --build grch37` is a supported install (ClinVar is then
    ready for GRCh37 input), so its stale cache must be probed and refreshed too."""
    db = tmp_path / "clinvar.GRCh37.sqlite"
    _clinvar_cache(db, "GRCh37", "md5:" + "0" * 32)
    stale = time.time() - _helpers._STALENESS_SECONDS - 3600
    os.utime(db, (stale, stale))
    assert ClinVarAnnotator(tmp_path, builds=("GRCh37",)).is_ready()

    probes: list[str] = []
    setups: list[tuple[str, ...]] = []
    monkeypatch.setattr(clinvar_mod, "_fetch_checksum", lambda url: probes.append(url) or "f" * 32)
    monkeypatch.setattr(ClinVarAnnotator, "setup", lambda self: setups.append(self._builds))
    monkeypatch.setattr(_helpers, "console", Console(file=io.StringIO()))

    _helpers._maybe_refresh_databases(tmp_path)

    assert probes, "the stale GRCh37-only ClinVar cache was never freshness-checked"
    assert setups and "GRCh37" in setups[0], "a changed ClinVar release was not refreshed"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_stale_single_build_is_refreshed_for_that_build_only(tmp_path, monkeypatch, build):
    db = tmp_path / f"clinvar.{build}.sqlite"
    _clinvar_cache(db, build, "md5:" + "0" * 32)
    stale = time.time() - _helpers._STALENESS_SECONDS - 3600
    os.utime(db, (stale, stale))
    setups: list[tuple[str, ...]] = []
    monkeypatch.setattr(clinvar_mod, "_fetch_checksum", lambda url: "f" * 32)
    monkeypatch.setattr(ClinVarAnnotator, "setup", lambda self: setups.append(self._builds))
    monkeypatch.setattr(_helpers, "console", Console(file=io.StringIO()))

    _helpers._maybe_refresh_databases(tmp_path)

    assert setups == [(build,)]


def test_offline_stale_single_build_gets_the_age_notice(tmp_path, monkeypatch):
    db = tmp_path / "clinvar.GRCh37.sqlite"
    _clinvar_cache(db, "GRCh37", "md5:" + "0" * 32)
    stale = time.time() - _helpers._STALENESS_SECONDS - 3 * 86400
    os.utime(db, (stale, stale))
    monkeypatch.setattr(clinvar_mod, "_fetch_checksum", lambda url: None)
    out = io.StringIO()
    monkeypatch.setattr(_helpers, "console", Console(file=out, width=200))

    _helpers._maybe_refresh_databases(tmp_path)

    assert "ClinVar database is 10 days old" in out.getvalue()


def test_both_builds_are_still_checked_together(tmp_path, monkeypatch):
    for build in ("GRCh37", "GRCh38"):
        db = tmp_path / f"clinvar.{build}.sqlite"
        _clinvar_cache(db, build, "md5:" + "0" * 32)
        stale = time.time() - _helpers._STALENESS_SECONDS - 3600
        os.utime(db, (stale, stale))
    setups: list[tuple[str, ...]] = []
    monkeypatch.setattr(clinvar_mod, "_fetch_checksum", lambda url: "f" * 32)
    monkeypatch.setattr(ClinVarAnnotator, "setup", lambda self: setups.append(self._builds))
    monkeypatch.setattr(_helpers, "console", Console(file=io.StringIO()))

    _helpers._maybe_refresh_databases(tmp_path)

    assert setups == [("GRCh37", "GRCh38")]


def _clinvar_row(output: str) -> str:
    return next(line for line in output.splitlines() if "ClinVar" in line)


def test_db_status_shows_a_single_build_install_ready(tmp_path):
    _clinvar_cache(tmp_path / "clinvar.GRCh37.sqlite", "GRCh37", "md5:" + "0" * 32)

    result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(tmp_path)])

    assert result.exit_code == 0, result.output
    row = _clinvar_row(result.output)
    assert "yes" in row
    assert "no" not in row.replace("ClinVar", "")


def test_db_status_shows_no_clinvar_not_ready(tmp_path):
    result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(tmp_path)])

    assert result.exit_code == 0, result.output
    assert "no" in _clinvar_row(result.output)

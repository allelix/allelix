# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Outside review of 2.6.1 RC1 (reviewer A): the findings fixed for RC2.

The tests, their helpers and their fixtures are reviewer A's, unchanged (from
~/reports/261-outside-review-rc1/A/out/test_review_rc1.py); only the imports are
narrowed to what these tests use. Finding 1's three tests (ClinVar, GWAS
Catalog and ClinPGx publishing beside an old journal) are not here: that code
is for 2.7. Offline and synthetic only.
"""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import shutil
import sqlite3
import urllib.request
from email.message import Message
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.databases import manager
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

_REAL_URLOPEN = urllib.request.urlopen


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    """No test may reach the network."""

    def refuse(*args, **kwargs):
        raise AssertionError("network use in an offline review test")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)


def _flat(text: str) -> str:
    return " ".join(text.split())


def _gnomad_like(path: Path, value: str) -> None:
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.execute(
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT, "
            "version TEXT, downloaded_at TEXT, record_count INTEGER, remote_signal TEXT, "
            "local_version_tag TEXT)"
        )
        conn.execute("CREATE TABLE gnomad_af (v TEXT)")
        conn.execute("INSERT INTO gnomad_af VALUES (?)", (value,))
        conn.execute("INSERT INTO database_versions VALUES ('gnomad','x','1','2026',1,NULL,NULL)")
        conn.commit()


def test_prebuilt_install_is_not_refused_over_a_journal_that_holds_nothing(tmp_path):
    """New in 2.6.1 (loader_utils._settle). A zero-length -journal (a write cut off
    by power loss before SQLite wrote the journal header) is not hot: SQLite
    ignores it and the old cache reads fine. _settle opens the old file, SQLite
    leaves the journal in place, and every install of that cache is refused with
    UnrecoveredJournalError, whose advice ("close any program using it and run db
    update again") can never succeed. 2.6.0 installed the new cache.
    Applies to every prebuilt cache (gnomAD, AlphaMissense, CADD, SNPedia, rsID
    identity): a schema bump then leaves the source not ready for good.
    """
    from allelix.databases import loader_utils
    from allelix.databases.cache_state import connect_readonly

    db = tmp_path / "gnomad.sqlite"
    _gnomad_like(db, "old")
    (tmp_path / "gnomad.sqlite-journal").write_bytes(b"")
    with contextlib.closing(connect_readonly(db)) as conn:  # the old cache reads fine
        assert conn.execute("SELECT v FROM gnomad_af").fetchall() == [("old",)]
    source = tmp_path / "source.sqlite"
    _gnomad_like(source, "new")
    gz = tmp_path / "gnomad.sqlite.gz"
    with open(source, "rb") as raw, gzip.open(gz, "wb") as packed:
        shutil.copyfileobj(raw, packed)

    loader_utils.install_prebuilt_gz_cache(gz, db, "gnomad", source_url="test://")

    with contextlib.closing(connect_readonly(db)) as conn:
        assert conn.execute("SELECT v FROM gnomad_af").fetchall() == [("new",)]


# --------------------------------------------------------------------------
# Finding 3: db status hides a damaged second ClinVar build (new in 2.6.1).
# --------------------------------------------------------------------------


def _ready_clinvar_cache(path: Path, build: str) -> None:
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
            "VALUES (?, 'test://', '20260101', '2026-01-01', 1, 'md5:x', ?)",
            (f"clinvar.{build}", f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def test_db_status_does_not_call_clinvar_ready_with_a_damaged_build(tmp_path, monkeypatch):
    """2.6.1 builds db status's ClinVar from `installed_clinvar_builds`, which
    returns the builds whose caches are READY. A two-build install whose GRCh38
    cache is damaged (or left outdated) therefore shows ClinVar "yes", and the
    damaged file is not mentioned anywhere in `db status`. 2.6.0 showed "no".
    GRCh38 input then runs with ClinVar skipped.
    """
    monkeypatch.setenv("COLUMNS", "200")
    _ready_clinvar_cache(tmp_path / "clinvar.GRCh37.sqlite", "GRCh37")
    (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"\0" * 4096)

    result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(tmp_path)])

    assert result.exit_code == 0, result.output
    row = next(line for line in result.output.splitlines() if "ClinVar" in line)
    assert " no " in row, f"db status reports ClinVar ready with a damaged GRCh38 cache: {row}"


# --------------------------------------------------------------------------
# Finding 4: db clean deletes a user directory shaped like a staging directory.
# --------------------------------------------------------------------------


def test_db_clean_keeps_a_users_hidden_backup_directory(tmp_path):
    """CHANGELOG: "A directory is removed only when it is one of allelix's staging
    directories"; ADR-0056: "Every other entry is kept". The test is only the
    shape `.<source>-<8 of [a-z0-9_]>`, so a user's dated backup
    `.clinvar-20261101/` (holding a copy of the caches) is removed with
    shutil.rmtree, without --force. RC1's own test keeps `.gnomad-mine`.
    """
    (tmp_path / "gnomad.sqlite").write_bytes(b"x")
    backup = tmp_path / ".clinvar-20261101"
    backup.mkdir()
    (backup / "clinvar.GRCh38.sqlite").write_bytes(b"my backup")

    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])

    assert result.exit_code == 0, result.output
    assert (backup / "clinvar.GRCh38.sqlite").exists(), "db clean deleted a user's directory"


# --------------------------------------------------------------------------
# Finding 5: db clean keeps the legacy clinvar.sqlite and says allelix did not
# create it (regression against 2.6.0).
# --------------------------------------------------------------------------


def test_db_clean_removes_the_legacy_single_clinvar_cache(tmp_path):
    """CHANGELOG v0.8 migration: "v0.7.x users have a single clinvar.sqlite cache
    ... The legacy clinvar.sqlite can be deleted". It is allelix's file. RC1's
    allow-list has only clinvar.GRCh37/GRCh38.sqlite, so db clean keeps it and
    prints "Kept 1 item(s) allelix did not create: clinvar.sqlite". 2.6.0
    removed it.
    """
    (tmp_path / "clinvar.sqlite").write_bytes(b"x" * 1024)
    (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"x")

    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])

    assert result.exit_code == 0, result.output
    assert "did not create: clinvar.sqlite" not in _flat(result.output)
    assert not (tmp_path / "clinvar.sqlite").exists()


# --------------------------------------------------------------------------
# Finding 6 (low): a CPIC table reply that passes the completeness check but is
# not a list of row objects fails the whole ClinPGx update.
# --------------------------------------------------------------------------


class _Reply(io.BytesIO):
    def __init__(self, body: object, content_range: str) -> None:
        super().__init__(json.dumps(body).encode())
        self.headers = Message()
        self.headers["Content-Range"] = content_range

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def test_cpic_table_of_non_objects_takes_the_cpic_unavailable_path(monkeypatch):
    """ADR-0020 amendment: anything but a whole table raises CpicIncompleteError and
    "ClinPGx setup now catches every probe error there (PROBE_ERRORS)". A reply
    `["x"]` with `Content-Range: 0-0/1` passes `_check_whole_table` (a list, range
    and total agree) and then raises AttributeError in fetch_cpic_allele_functions,
    which PROBE_ERRORS does not catch: the ClinPGx update fails instead of being
    built without CPIC rows.
    """
    from allelix.databases import cpic_loader

    monkeypatch.setattr(cpic_loader, "CPIC_RETRY_BACKOFF_SECONDS", (0.0, 0.0, 0.0))
    monkeypatch.setattr(
        urllib.request, "urlopen", lambda request, timeout=None: _Reply(["x"], "0-0/1")
    )
    probe_errors = getattr(manager, "PROBE_ERRORS", (OSError, ValueError))

    try:
        cpic_loader.fetch_cpic_allele_functions()
    except probe_errors:
        return  # the ClinPGx setup's fallback path catches it
    except Exception as exc:  # the point of the test
        pytest.fail(f"{type(exc).__name__} escapes the CPIC-unavailable path: {exc}")

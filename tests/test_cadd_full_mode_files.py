# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""CADD full mode without its tabix files, on every surface.

Full mode reads only CADD's own files (`whole_genome_SNVs.tsv.gz` and its
`.tbi`), which the user downloads from CADD and places in the data
directory; `db update` cannot fetch them and full mode never reads the
SQLite cache. So `db update --cadd` downloads nothing for full mode and
reports the missing files as a failure; `analyze` names them and does not
send the user to `db update`; `db status` shows CADD not ready.

The first test is reviewer C's (outside review of 2.6.0 RC4, finding
R4-06), unchanged but for line wrapping. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import gzip
import sqlite3
import urllib.request
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators import cadd as cadd_mod
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli import _helpers, main
from allelix.cli import db as db_cli
from allelix.databases import manager

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    """No test may reach the network."""

    def refuse(*args, **kwargs):
        raise AssertionError("network use in an offline regression test")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    monkeypatch.setattr(manager, "open_url", refuse)


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_r06_cadd_full_mode_does_not_redownload_the_sqlite_cache(tmp_path, monkeypatch):
    """ADR-0030: a code-driven cache is downloaded once. With options.cadd_full
    true and the tabix files not yet in place, each `db update --cadd` downloads
    and installs the ~5.8 GB SQLite cache again (full mode never reads it) and
    prints "✓ cadd ready"."""
    raw = tmp_path / "cadd-raw.sqlite"
    with contextlib.closing(sqlite3.connect(raw)) as conn:
        conn.execute(
            "CREATE TABLE cadd_scores (chrom TEXT, pos INT, ref TEXT, alt TEXT, phred REAL)"
        )
        conn.execute("INSERT INTO cadd_scores VALUES ('1', 1, 'A', 'G', 20.0)")
        conn.execute(
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
            "version TEXT, downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
            "remote_signal TEXT, local_version_tag TEXT)"
        )
        conn.execute("INSERT INTO database_versions VALUES ('cadd','x','v1.7','now',1,NULL,NULL)")
        conn.commit()
    payload = gzip.compress(raw.read_bytes())
    downloads: list[str] = []

    def fake_download(url, dest):
        downloads.append(url)
        dest.write_bytes(payload)

    monkeypatch.setattr(cadd_mod, "download", fake_download)
    monkeypatch.setattr(cadd_mod, "verify_file_hash", lambda *a, **k: None)
    monkeypatch.setattr(
        db_cli,
        "get_annotators",
        lambda folder, **kw: [CaddAnnotator(folder, full_mode=kw.get("cadd_full", False))],
    )
    data = tmp_path / "data"
    data.mkdir()
    runner = CliRunner()
    assert (
        runner.invoke(
            main, ["config", "set", "options.cadd_full", "true", "--data-dir", str(data)]
        ).exit_code
        == 0
    )

    first = runner.invoke(main, ["db", "update", "--cadd", "--data-dir", str(data)], input="y\n")
    after_first = len(downloads)
    second = runner.invoke(main, ["db", "update", "--cadd", "--data-dir", str(data)], input="y\n")

    # The first update may install the SQLite cache, or say the full mode files
    # are missing; either way the second has nothing new to download.
    assert after_first <= 1, first.output
    assert len(downloads) == after_first, (
        f"the second db update downloaded the CADD cache again:\n{second.output}"
    )


@pytest.fixture
def full_mode(tmp_path, monkeypatch):
    """A data directory with CADD full mode on and only CADD managed by `db update`."""
    downloads: list[str] = []

    def fake_download(url, dest):
        downloads.append(url)
        raise AssertionError("full mode must not download the SQLite cache")

    monkeypatch.setattr(cadd_mod, "download", fake_download)
    monkeypatch.setattr(
        db_cli,
        "get_annotators",
        lambda folder, **kw: [CaddAnnotator(folder, full_mode=kw.get("cadd_full", False))],
    )
    data = tmp_path / "data"
    data.mkdir()
    runner = CliRunner()
    for key in ("options.cadd_full", "sources.cadd"):
        assert (
            runner.invoke(main, ["config", "set", key, "true", "--data-dir", str(data)]).exit_code
            == 0
        )
    return data, downloads


def _place_full_mode_files(data: Path) -> None:
    for name in ("whole_genome_SNVs.tsv.gz", "whole_genome_SNVs.tsv.gz.tbi"):
        (data / name).write_bytes(b"x")


def test_db_update_reports_missing_full_mode_files_as_a_failure(full_mode):
    data, downloads = full_mode

    result = CliRunner().invoke(
        main, ["db", "update", "--cadd", "--data-dir", str(data)], input="y\n"
    )

    text = _flat(result.output)
    assert downloads == []
    assert result.exit_code == 1, result.output
    assert "whole_genome_SNVs.tsv.gz" in text
    assert "cannot fetch" in text
    assert "cadd ready" not in text


def test_db_update_with_full_mode_files_in_place_downloads_nothing(full_mode):
    data, downloads = full_mode
    _place_full_mode_files(data)

    result = CliRunner().invoke(
        main, ["db", "update", "--cadd", "--data-dir", str(data)], input="y\n"
    )

    assert downloads == []
    assert result.exit_code == 0, result.output


def test_force_does_not_download_the_cache_full_mode_never_reads(full_mode):
    data, downloads = full_mode
    _place_full_mode_files(data)

    result = CliRunner().invoke(
        main, ["db", "update", "--cadd", "--force", "--data-dir", str(data)], input="y\n"
    )

    assert downloads == []
    assert result.exit_code == 0, result.output


def test_analyze_line_names_the_files_and_not_db_update(tmp_path):
    cadd = CaddAnnotator(tmp_path, full_mode=True)

    line = _helpers.unready_line(cadd, cadd.cache_status())

    assert "whole_genome_SNVs.tsv.gz" in line
    assert "cadd.gs.washington.edu" in line
    assert "Run `allelix db update`" not in line


def test_analyze_line_names_a_missing_index(tmp_path):
    (tmp_path / "whole_genome_SNVs.tsv.gz").write_bytes(b"x")
    cadd = CaddAnnotator(tmp_path, full_mode=True)

    line = _helpers.unready_line(cadd, cadd.cache_status())

    assert "whole_genome_SNVs.tsv.gz.tbi" in line
    assert "Run `allelix db update`" not in line


@pytest.mark.parametrize(
    "annotator",
    [
        lambda d: CaddAnnotator(d),
        lambda d: GnomadAnnotator(d),
    ],
    ids=["cadd-cache-mode", "gnomad"],
)
def test_other_missing_caches_still_point_to_db_update(tmp_path, annotator):
    source = annotator(tmp_path)

    line = _helpers.unready_line(source, source.cache_status())

    assert line.endswith("not installed. Run `allelix db update`.")


def test_db_status_shows_full_mode_without_files_not_ready(full_mode):
    data, _ = full_mode

    result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])

    assert result.exit_code == 0, result.output
    row = next(line for line in result.output.splitlines() if "CADD" in line)
    assert "no" in row
    assert "yes" not in row

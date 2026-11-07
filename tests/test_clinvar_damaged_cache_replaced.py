# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A ClinVar cache file that is not a database is replaced by `db update`.

`analyze` calls such a file damaged and says to run `db update`; `db update`
refused it, because ClinVar's relative floor (ADR-0053) needs the existing
cache's record count. A file that is not an SQLite database has no baseline
to protect, so it is replaced like a missing cache (the absolute floor still
applies). A readable database with a missing or invalid count is still
refused, as before.

The first test and its fixtures are reviewer C's (outside review of 2.6.0
RC4, finding R4-05), unchanged but for line wrapping. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import os
import sqlite3
import urllib.request
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators import clinvar as clinvar_mod
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.cli import _helpers, main
from allelix.databases import manager
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.cache_state import CacheState
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

_REAL_URLOPEN = urllib.request.urlopen


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


def _clinvar_sources(folder: Path, build: str) -> dict[str, bytes]:
    """Synthetic variant_summary / submission_summary pair giving two records."""
    vs_header = [manager._VS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(43)]
    vs_header[0] = "#AlleleID"
    vs_rows = []
    for identifier, pos in ((101, "1000"), (102, "2000")):
        row = [""] * 43
        for index, value in {
            0: str(identifier),
            1: "single nucleotide variant",
            2: "NM_SYNTHETIC:c.1A>G",
            4: "GENE",
            6: "Pathogenic",
            9: str(900252900 + identifier),
            16: build,
            18: "1",
            24: "criteria provided, single submitter",
            30: str(identifier),
            31: pos,
            32: "A",
            33: "G",
        }.items():
            row[index] = value
        vs_rows.append("\t".join(row))
    ss_header = [manager._SS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(16)]
    ss_header[0] = "#VariationID"
    ss_rows = []
    for number, identifier in enumerate((101, 102), 1):
        row = [""] * 16
        for index, value in {
            0: str(identifier),
            1: "Pathogenic",
            5: "na:Synthetic condition",
            6: "criteria provided, single submitter",
            10: f"SCV9002529{number:02}.1",
            15: "yes",
        }.items():
            row[index] = value
        ss_rows.append("\t".join(row))
    out = {}
    for name, header, rows in (
        ("variant_summary.txt.gz", vs_header, vs_rows),
        ("submission_summary.txt.gz", ss_header, ss_rows),
    ):
        text = "\t".join(header) + "\n" + "\n".join(rows) + "\n"
        out[name] = gzip.compress(text.encode(), mtime=0)
    return out


def test_r05_damaged_clinvar_cache_is_replaced_by_db_update(tmp_path, monkeypatch):
    """ADR-0048: a corrupt cache is reported as "damaged ... Run `allelix db update`",
    which replaces only the damaged source. For ClinVar, db update refuses forever
    because the damaged file has no record count (with or without --force)."""
    data = tmp_path / "data"
    data.mkdir()
    damaged = data / "clinvar.GRCh37.sqlite"
    damaged.write_bytes(b"\x00" * 4096)
    status = ClinVarAnnotator(data, builds=("GRCh37",)).cache_status()
    assert status.state is CacheState.CORRUPT
    assert "Run `allelix db update`" in _helpers.unready_line(ClinVarAnnotator(data), status)

    payloads = _clinvar_sources(tmp_path, "GRCh37")
    digests = {name: hashlib.md5(body).hexdigest() for name, body in payloads.items()}
    monkeypatch.setattr(
        clinvar_mod, "_fetch_checksum", lambda url: digests[url.rsplit("/", 1)[-1]]
    )
    monkeypatch.setattr(
        clinvar_mod,
        "download",
        lambda url, dest: dest.write_bytes(payloads[url.rsplit("/", 1)[-1]]),
    )
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", {"GRCh37": 1, "GRCh38": 1})

    result = CliRunner().invoke(
        main, ["db", "update", "--data-dir", str(data), "--only", "clinvar", "--build", "grch37"]
    )

    assert result.exit_code == 0, result.output
    assert ClinVarAnnotator(data, builds=("GRCh37",)).is_ready()


@pytest.fixture
def offline_sources(tmp_path, monkeypatch):
    payloads = _clinvar_sources(tmp_path, "GRCh37")
    digests = {name: hashlib.md5(body).hexdigest() for name, body in payloads.items()}
    monkeypatch.setattr(
        clinvar_mod, "_fetch_checksum", lambda url: digests[url.rsplit("/", 1)[-1]]
    )
    monkeypatch.setattr(
        clinvar_mod,
        "download",
        lambda url, dest: dest.write_bytes(payloads[url.rsplit("/", 1)[-1]]),
    )
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", {"GRCh37": 1, "GRCh38": 1})
    data = tmp_path / "data"
    data.mkdir()
    return data


def _update(data, *args):
    return CliRunner().invoke(
        main, ["db", "update", "--data-dir", str(data), "--only", "clinvar", *args]
    )


@pytest.mark.parametrize(
    "content",
    [b"", b"not a database at all", b"SQLite format 3\x00" + b"\xff" * 200],
    ids=["empty", "text", "bad-header"],
)
def test_a_file_that_is_not_a_database_is_replaced(offline_sources, content):
    data = offline_sources
    (data / "clinvar.GRCh37.sqlite").write_bytes(content)

    result = _update(data, "--build", "grch37")

    assert result.exit_code == 0, result.output
    assert ClinVarAnnotator(data, builds=("GRCh37",)).is_ready()
    assert "not an SQLite database" in " ".join(result.output.split())


def test_a_readable_cache_with_no_count_is_still_refused(offline_sources):
    data = offline_sources
    path = data / "clinvar.GRCh37.sqlite"
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.commit()
    before = path.read_bytes()

    result = _update(data, "--build", "grch37")

    assert result.exit_code == 1
    assert "record count is missing or invalid" in " ".join(result.output.split())
    assert path.read_bytes() == before


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads any file")
def test_an_unreadable_cache_is_still_refused(offline_sources):
    data = offline_sources
    path = data / "clinvar.GRCh37.sqlite"
    path.write_bytes(b"not a database at all")
    path.chmod(0)
    try:
        result = _update(data, "--build", "grch37")
    finally:
        path.chmod(0o644)

    assert result.exit_code == 1
    assert path.read_bytes() == b"not a database at all"


def test_with_both_builds_a_damaged_file_does_not_block_the_other(
    offline_sources, monkeypatch, tmp_path
):
    data = offline_sources
    grch38 = _clinvar_sources(tmp_path, "GRCh38")
    grch37 = _clinvar_sources(tmp_path, "GRCh37")

    # One variant_summary carries both builds' rows, as NCBI's does.
    def both(name):
        a = gzip.decompress(grch37[name]).decode().splitlines()
        b = gzip.decompress(grch38[name]).decode().splitlines()
        rows = a[1:] + (b[1:] if name.startswith("variant") else [])
        return gzip.compress(("\n".join([a[0], *rows]) + "\n").encode(), mtime=0)

    payloads = {name: both(name) for name in grch37}
    digests = {name: hashlib.md5(body).hexdigest() for name, body in payloads.items()}
    monkeypatch.setattr(
        clinvar_mod, "_fetch_checksum", lambda url: digests[url.rsplit("/", 1)[-1]]
    )
    monkeypatch.setattr(
        clinvar_mod,
        "download",
        lambda url, dest: dest.write_bytes(payloads[url.rsplit("/", 1)[-1]]),
    )
    (data / "clinvar.GRCh37.sqlite").write_bytes(b"\x00" * 4096)
    _clinvar_cache(data / "clinvar.GRCh38.sqlite", "GRCh38", "md5:" + "0" * 32)

    result = _update(data, "--build", "both", "--force")

    assert result.exit_code == 0, result.output
    for build in ("GRCh37", "GRCh38"):
        assert ClinVarAnnotator(data, builds=(build,)).is_ready(), build

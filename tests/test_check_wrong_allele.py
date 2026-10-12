# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Smoke tests for scripts/check_wrong_allele.py on synthetic caches."""

from __future__ import annotations

import contextlib
import json
import sqlite3
import sys
from pathlib import Path

import pytest

from allelix.databases.schema import (
    ALPHAMISSENSE_SCHEMA,
    CADD_SCHEMA,
    CLINVAR_SCHEMA,
    GNOMAD_SCHEMA,
    RSID_IDENTITY_SCHEMA,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from check_wrong_allele import main

ALL_SOURCES = ["clinvar", "rsid_identity", "gnomad", "alphamissense", "cadd"]


def _db(path: Path, schema: str, sql: str, rows: list[tuple]) -> None:
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(schema)
        if rows:
            conn.executemany(sql, rows)
        conn.commit()


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    """rs1 at GRCh37 1:90 / GRCh38 1:100 G>A (another allele, G>T, has other
    values); rs9 only inside a ';'-joined gnomAD rsID field at 1:200 C>T."""
    d = tmp_path / "data"
    d.mkdir()
    _db(
        d / "rsid_identity.sqlite",
        RSID_IDENTITY_SCHEMA,
        "INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)",
        [("rs1", "1", 90, "G", "A", 100, "G", "A", 0)],
    )
    _db(
        d / "gnomad.sqlite",
        GNOMAD_SCHEMA,
        "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) VALUES (?,?,?,?,?,?)",
        [
            ("1", 100, "G", "A", "rs1", 0.25),
            ("1", 100, "G", "T", "rs1", 0.01),
            ("1", 200, "C", "T", "rs8;rs9", 0.5),
        ],
    )
    _db(
        d / "alphamissense.sqlite",
        ALPHAMISSENSE_SCHEMA,
        "INSERT INTO alphamissense_scores (chrom, pos, ref, alt, am_pathogenicity, am_class) "
        "VALUES (?,?,?,?,?,?)",
        [("1", 100, "G", "A", 0.9, "likely_pathogenic")],
    )
    _db(
        d / "cadd.sqlite",
        CADD_SCHEMA,
        "INSERT INTO cadd_scores VALUES (?,?,?,?,?)",
        [("1", 100, "G", "A", 25.0), ("1", 100, "G", "T", 3.0)],
    )
    for build in ("GRCh37", "GRCh38"):
        _db(d / f"clinvar.{build}.sqlite", CLINVAR_SCHEMA, "", [])
    return d


def _report(tmp_path: Path, build: str, annotations: list[dict], sources=None) -> Path:
    path = tmp_path / "report.json"
    path.write_text(
        json.dumps(
            {
                "input": {"build": build},
                "annotators": [{"name": s} for s in (sources or ALL_SOURCES)],
                "annotations": annotations,
            }
        )
    )
    return path


def _ann(rsid: str = "rs1", alt: str | None = "A", **values) -> dict:
    return {"rsid": rsid, "source": "clinvar", "alt": alt, "genotype_match": "GA", **values}


def _input(tmp_path: Path, rows: str = "") -> Path:
    path = tmp_path / "input.vcf"
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + rows)
    return path


def _run(report: Path, data: Path, inp: Path, *extra: str) -> int:
    return main([str(report), "--input", str(inp), "--data-dir", str(data), *extra])


DIRECT = {"allele_frequency": 0.25, "am_pathogenicity": 0.9, "cadd_phred": 25.0}


@pytest.mark.parametrize("build", ["GRCh38", "GRCh37"])
def test_allele_direct_values_pass(tmp_path: Path, data_dir: Path, build: str) -> None:
    """GRCh37 keys are translated through the identity table by locus and alleles."""
    report = _report(tmp_path, build, [_ann(**DIRECT)])
    assert _run(report, data_dir, _input(tmp_path)) == 0


def test_alt_less_row_uses_the_carried_allele(tmp_path: Path, data_dir: Path) -> None:
    report = _report(tmp_path, "GRCh38", [_ann(alt=None, **DIRECT)])
    assert _run(report, data_dir, _input(tmp_path)) == 0


def test_another_alleles_value_fails(tmp_path: Path, data_dir: Path, capsys) -> None:
    report = _report(tmp_path, "GRCh38", [_ann(allele_frequency=0.01, cadd_phred=3.0)])
    assert _run(report, data_dir, _input(tmp_path)) == 1
    assert "0/1 allele-direct, 1 non-direct" in capsys.readouterr().out


def test_value_matching_nothing_fails(tmp_path: Path, data_dir: Path, capsys) -> None:
    report = _report(tmp_path, "GRCh38", [_ann(**{**DIRECT, "allele_frequency": 0.7})])
    assert _run(report, data_dir, _input(tmp_path)) == 1
    assert "1 mismatch" in capsys.readouterr().out


def test_rsid_no_cache_places_needs_the_input(tmp_path: Path, data_dir: Path, capsys) -> None:
    """rs9 sits only inside a ';'-joined gnomAD field: the input's ID column places it."""
    report = _report(tmp_path, "GRCh38", [_ann("rs9", "T", allele_frequency=0.5)])
    assert _run(report, data_dir, _input(tmp_path)) == 1
    assert "1 unresolved" in capsys.readouterr().out
    placed = _input(tmp_path, "1\t200\trs9\tC\tT\t50\tPASS\t.\n")
    assert _run(report, data_dir, placed) == 0


def test_input_is_required(tmp_path: Path, data_dir: Path) -> None:
    report = _report(tmp_path, "GRCh38", [_ann(**DIRECT)])
    with pytest.raises(SystemExit) as exc:
        main([str(report), "--data-dir", str(data_dir)])
    assert exc.value.code == 2


def test_only_provenance_caches_are_required(tmp_path: Path, data_dir: Path, capsys) -> None:
    (data_dir / "cadd.sqlite").unlink()
    without_cadd = [s for s in ALL_SOURCES if s != "cadd"]
    values = {k: v for k, v in DIRECT.items() if k != "cadd_phred"}
    report = _report(tmp_path, "GRCh38", [_ann(**values)], without_cadd)
    assert _run(report, data_dir, _input(tmp_path)) == 0
    report = _report(tmp_path, "GRCh38", [_ann(**values)])
    assert _run(report, data_dir, _input(tmp_path)) == 2
    assert "cadd is in the provenance" in capsys.readouterr().err


def test_value_from_a_source_not_in_provenance_fails(tmp_path: Path, data_dir: Path) -> None:
    without_cadd = [s for s in ALL_SOURCES if s != "cadd"]
    report = _report(tmp_path, "GRCh38", [_ann(**DIRECT)], without_cadd)
    assert _run(report, data_dir, _input(tmp_path)) == 1


def test_no_enriched_row(tmp_path: Path, data_dir: Path) -> None:
    report = _report(tmp_path, "GRCh38", [_ann()])
    assert _run(report, data_dir, _input(tmp_path)) == 1
    assert _run(report, data_dir, _input(tmp_path), "--allow-empty") == 0

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Weak outranked pathogenic assertions cannot outrank the aggregate in a report."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import load_clinvar_tsv
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.models import Variant
from tests.databases.test_clinvar_tsv_loader import _SS_HEADER, _VS_HEADER, _ss_row, _vs_row

RSID = "rs900260346"
WEAK = "no assertion criteria provided"
STRONG = "criteria provided, single submitter"


def _cache(tmp_path, build, aggregate, sig="Pathogenic", *, outside=1, review=WEAK, fresh=False):
    data = tmp_path / "data"
    data.mkdir()
    db = data / clinvar_db_filename(build)
    if fresh:
        vs = tmp_path / "variants.tsv"
        ss = tmp_path / "submissions.tsv"
        vs.write_text(
            _VS_HEADER
            + "\n"
            + _vs_row(
                gene="MOCK",
                rs=RSID[2:],
                assembly=build,
                chrom="1",
                position_vcf="12345",
                ref="A",
                alt="G",
                aggregate=aggregate,
                name="NM_MOCK:c.1A>G",
            )
            + "\n"
        )
        ss.write_text(
            _SS_HEADER
            + "\n"
            + _ss_row(
                significance=aggregate.split(";")[0],
                review_status=STRONG,
                reported_phenotype="na:Mock condition",
                scv="SCV900260346.1",
            )
            + "\n"
            + _ss_row(
                significance=sig,
                review_status=review,
                contributes="no",
                reported_phenotype="na:Mock condition",
                scv="SCV900260347.1",
            )
            + "\n"
        )
        load_clinvar_tsv(vs, ss, db, build, record_name=clinvar_record_name(build))
    else:
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            conn.execute(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
                "submitters_same, submitters_total, outside_aggregate) "
                "VALUES (?, '1', 12345, 'A', 'G', ?, 'Mock condition', "
                "'MOCK', ?, 100, ?, 1, 1, ?)",
                (RSID, sig, review, aggregate, outside),
            )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
                "record_count, local_version_tag) VALUES (?, 'test://mock', 'synthetic', "
                "'2001-01-01', 1, ?)",
                (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    return data


def _rows(data, build, sig):
    with ClinVarAnnotator(data, builds=(build,), include_benign=True) as ann:
        rows = ann.annotate(Variant(RSID, "1", 12345, "A", "G", ref="A", build=build))
    label = "clinvar_" + sig.lower().replace(" ", "_")
    return [row for row in rows if row.significance == label]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "aggregate", ["Likely risk allele", "other", "Likely risk allele; drug response"]
)
@pytest.mark.parametrize(
    "sig", ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
)
@pytest.mark.parametrize("fresh", [False, True], ids=["cache-row", "real-loader"])
def test_outranked_weak_pathogenic_row_is_kept_but_capped(tmp_path, build, aggregate, sig, fresh):
    data = _cache(tmp_path, build, aggregate, sig, fresh=fresh)
    rows = _rows(data, build, sig)
    assert len(rows) == 1, "The ruled fix caps this finding; it must not silently disappear"
    row = rows[0]
    assert row.magnitude == 4.0
    assert row.aggregate_capped is True
    assert row.outside_aggregate is True
    assert row.clinvar_aggregate == aggregate
    assert (row.submitters_same, row.submitters_total) == (1, 1)


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "aggregate,expected",
    [
        ("Conflicting classifications of pathogenicity", 4.0),
        ("Benign", 1.0),
        ("Likely benign", 2.0),
        ("Benign/Likely benign", 1.5),
        ("Pathogenic", 9.0),
        ("Likely pathogenic", 9.0),
        ("drug response", 9.0),
        ("risk factor", 9.0),
    ],
)
def test_other_aggregate_policies_stay_as_they_are(tmp_path, build, aggregate, expected):
    data = _cache(tmp_path, build, aggregate)
    (row,) = _rows(data, build, "Pathogenic")
    assert row.magnitude == expected


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("aggregate", ["Likely risk allele", "other"])
@pytest.mark.parametrize(
    "sig,expected",
    [
        ("drug response", 6.5),
        ("risk factor", 6.0),
        ("protective", 5.0),
    ],
)
def test_non_disease_rows_do_not_acquire_the_pathogenic_cap(
    tmp_path, build, aggregate, sig, expected
):
    data = _cache(tmp_path, build, aggregate, sig)
    (row,) = _rows(data, build, sig)
    assert row.magnitude == expected
    assert row.aggregate_capped is False


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("aggregate", ["Likely risk allele", "other"])
def test_contributing_pathogenic_row_is_outside_this_change(tmp_path, build, aggregate):
    data = _cache(tmp_path, build, aggregate, outside=0)
    (row,) = _rows(data, build, "Pathogenic")
    assert row.magnitude == 9.0
    assert row.outside_aggregate is False


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("aggregate", ["Likely risk allele", "other"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_the_actual_report_shows_the_scored_down_finding(tmp_path, build, aggregate, surface):
    data = _cache(tmp_path, build, aggregate)
    source = tmp_path / "mock.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMOCK\n"
        f"1\t12345\t{RSID}\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    )
    out = tmp_path / f"report.{surface}"
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build.lower(),
        "--min-magnitude",
        "0",
    ]
    if surface != "terminal":
        args += ["--output", str(out)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    if surface == "json":
        payload = json.loads(out.read_text())
        (row,) = [a for a in payload["annotations"] if a["rsid"] == RSID]
        assert row["magnitude"] == 4.0
        assert row["aggregate_capped"] is True
        assert payload["clinvar"]["pathogenic_rows_scored_down"] == 1
        assert any(a["rsid"] == RSID for a in payload["clinvar"]["scored_down_variants"])
    else:
        text = result.output if surface == "terminal" else out.read_text()
        assert RSID in text
        assert "scored down" in text.lower()
        assert aggregate.lower() in text.lower()

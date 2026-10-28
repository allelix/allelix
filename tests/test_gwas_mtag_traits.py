# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Roll up MTAG suffixes from actual GWAS annotations without losing distinct findings."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GWAS_SCHEMA
from allelix.models import Variant
from allelix.reports._pipeline import rollup_gwas_duplicates

RSID = "rs900252701"
TRAIT = "Synthetic vascular trait"
SUFFIXES = (" (MTAG)", " (mtag)", " (MtAg)", "  ( MTAG )  ", "(MTAG)", " \t( mTaG\t) ")


def _annotations(tmp_path, traits, build, *, batch=True, rsid=RSID):
    """traits is (trait, risk allele, study); a G/C carrier supports both risk alleles."""
    data = tmp_path / "data"
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.executescript(GWAS_SCHEMA)
        for trait, allele, study in traits:
            conn.execute(
                "INSERT INTO gwas_associations (rsid, risk_allele, trait, p_value, or_beta, "
                "ci_text, gene, study_accession, pubmed_id, trait_category) VALUES "
                "(?, ?, ?, 1e-30, 1.5, '[1.2-1.8]', 'MTHFR', ?, '900252701', 'disease')",
                (rsid, allele, trait, study),
            )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('gwas', 'test://synthetic', 'synthetic', "
            "'2001-01-01', ?, ?)",
            (len(traits), f"cv:{_CATEGORIZER_VERSION}"),
        )
        conn.commit()
    write_identity_db(
        [(rsid, "1", 1000, "A", alt, 2000, "A", alt, 0) for alt in ("G", "C")],
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://identity",
        downloaded_at="2001-01-01",
    )
    variant = Variant(
        rsid,
        "1",
        1000 if build == "GRCh37" else 2000,
        "G",
        "C",
        build=build,
        ref="A",
        ref_from_input=True,
    )
    with (
        RsidIdentityAnnotator(data) as identity,
        GWASCatalogAnnotator(data, forward_site_provider=identity.forward_sites) as annotator,
    ):
        assert annotator.is_ready()
        rows = list(annotator.batch_annotate([variant])) if batch else annotator.annotate(variant)
    assert len(rows) == len(traits), rows
    assert {row.trait for row in rows} == {trait.strip() for trait, _, _ in traits}
    return rows, data, variant


@pytest.mark.parametrize("suffix", SUFFIXES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("same_study", [False, True])
def test_plain_twin_wins_with_real_trait_suffix(tmp_path, suffix, build, batch, same_study):
    rows, _, _ = _annotations(
        tmp_path,
        [
            (TRAIT, "G", "GCST900252701"),
            (TRAIT + suffix, "G", "GCST900252701" if same_study else "GCST900252702"),
        ],
        build,
        batch=batch,
    )
    plain = next(row for row in rows if row.trait == TRAIT)
    for ordered in (rows, rows[::-1]):
        assert rollup_gwas_duplicates(ordered) == [plain]


@pytest.mark.parametrize("suffix", SUFFIXES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_mtag_only_trait_is_retained(tmp_path, suffix, build):
    rows, _, _ = _annotations(tmp_path, [(TRAIT + suffix, "G", "GCST900252701")], build)
    assert rollup_gwas_duplicates(rows) == rows


@pytest.mark.parametrize("difference", ["trait", "allele", "internal-mtag"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_other_differences_are_not_mtag_twins(tmp_path, difference, build):
    other = {
        "trait": ("Another synthetic trait (MTAG)", "G", "GCST900252701"),
        "allele": (TRAIT + " (MTAG)", "C", "GCST900252701"),
        "internal-mtag": (TRAIT + " (MTAG) subtype", "G", "GCST900252701"),
    }[difference]
    rows, _, _ = _annotations(tmp_path, [(TRAIT, "G", "GCST900252701"), other], build)
    for ordered in (rows, rows[::-1]):
        assert len(rollup_gwas_duplicates(ordered)) == 2


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_different_rsids_are_not_mtag_twins(tmp_path, build):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    plain, _, _ = _annotations(first, [(TRAIT, "G", "GCST900252701")], build)
    mtag, _, _ = _annotations(
        second, [(TRAIT + " (MTAG)", "G", "GCST900252701")], build, rsid="rs900252702"
    )
    for ordered in (plain + mtag, mtag + plain):
        result = rollup_gwas_duplicates(ordered)
        assert len(result) == 2
        assert all(row in result for row in ordered)


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("report_format", ["json", "html", "terminal", "diff"])
@pytest.mark.parametrize("same_study", [False, True])
def test_cli_reports_roll_up_actual_mtag_twins(tmp_path, build, report_format, same_study):
    annotations, data, variant = _annotations(
        tmp_path,
        [
            (TRAIT, "G", "GCST900252701"),
            (TRAIT + " (MTAG)", "G", "GCST900252701" if same_study else "GCST900252702"),
        ],
        build,
    )
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{variant.position}\t{RSID}\tA\tG,C\t.\tPASS\t.\tGT\t1/2\n"
    )
    output = tmp_path / f"report.{report_format}"
    args = [
        "analyze",
        str(source),
        "--build",
        build.lower(),
        "--data-dir",
        str(data),
        "--no-update",
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
    ]
    if report_format != "terminal":
        args += [
            "--output",
            str(output),
            "--report-format",
            "json" if report_format == "diff" else report_format,
        ]
    if report_format == "diff":
        from allelix.models import annotation_to_public_dict
        from allelix.reports.json_report import SCHEMA_VERSION

        baseline = tmp_path / "previous.json"
        plain = next(row for row in annotations if row.trait == TRAIT)
        baseline.write_text(
            json.dumps(
                {
                    "schema_version": SCHEMA_VERSION,
                    "generated_at": "2001-01-01",
                    "input": {"file": str(source)},
                    "annotations": [annotation_to_public_dict(plain)],
                }
            )
        )
        args += ["--diff", str(baseline)]
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
    assert result.exit_code == 0, result.output
    if report_format == "terminal":
        assert TRAIT in result.output
        assert "(MTAG)" not in result.output
    elif report_format == "diff":
        diff = json.loads(output.read_text())["diff"]
        assert diff["new"] == [] and diff["removed"] == [] and diff["changed"] == [], diff
    elif report_format == "json":
        rows = [
            row for row in json.loads(output.read_text())["annotations"] if row["source"] == "gwas"
        ]
        assert len(rows) == 1, rows
        assert rows[0]["trait"] == TRAIT
    else:
        # The suffix belongs to neither an HTML finding nor its sidebar JSON payload.
        assert TRAIT in output.read_text()
        assert TRAIT + " (MTAG)" not in output.read_text()

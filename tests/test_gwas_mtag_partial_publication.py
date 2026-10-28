# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Partial publication metadata never splits an MTAG twin from its plain row."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.reports._pipeline import rollup_gwas_duplicates
from tests.test_gwas_mtag_traits import RSID, TRAIT, _annotations


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["rollup", "json", "html", "terminal"])
@pytest.mark.parametrize("missing", ["plain", "mtag", "both"])
@pytest.mark.parametrize("study_kind", ["same", "different", "none"])
def test_mtag_twins_roll_up_with_partial_publication_metadata(
    tmp_path, build, surface, missing, study_kind
):
    first, second = {
        "same": ("GCST900252701", "GCST900252701"),
        "different": ("GCST900252701", "GCST900252702"),
        "none": (None, None),
    }[study_kind]
    _, data, variant = _annotations(
        tmp_path,
        [(TRAIT, "G", first), (TRAIT + " (MTAG)", "G", second)],
        build,
    )
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as connection:
        for trait in (TRAIT, TRAIT + " (MTAG)"):
            if missing == "both" or (trait == TRAIT) == (missing == "plain"):
                connection.execute(
                    "UPDATE gwas_associations SET pubmed_id = NULL WHERE trait = ?", (trait,)
                )
        connection.commit()
    if surface == "rollup":
        with (
            RsidIdentityAnnotator(data) as identity,
            GWASCatalogAnnotator(data, forward_site_provider=identity.forward_sites) as annotator,
        ):
            rows = annotator.annotate(variant)
        assert len(rows) == 2
        for row in rows:
            accession = first if row.trait == TRAIT else second
            if accession is not None:
                assert f"gwas:{accession}" in row.references
            else:
                assert not any(ref.startswith("gwas:") for ref in row.references)
        plain = next(row for row in rows if row.trait == TRAIT)
        for ordered in (rows, rows[::-1]):
            assert rollup_gwas_duplicates(ordered) == [plain]
        return
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{variant.position}\t{RSID}\tA\tG,C\t.\tPASS\t.\tGT\t1/2\n"
    )
    output = tmp_path / f"report.{surface}"
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
    if surface != "terminal":
        args += ["--output", str(output), "--report-format", surface]
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    if surface == "json":
        rows = [
            row for row in json.loads(output.read_text())["annotations"] if row["source"] == "gwas"
        ]
        assert len(rows) == 1, rows
        assert rows[0]["trait"] == TRAIT
    else:
        text = result.output if surface == "terminal" else output.read_text()
        assert TRAIT in text
        mtag_visible = "(MTAG)" in text
        assert not mtag_visible, (surface, missing)

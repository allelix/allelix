# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A hidden PheCode sibling cannot erase a visible or scored-down diff identity."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_diff_annotator_identity import RSID, _args, _cache, _input, _visible_text

from allelix.cli import main


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("output_kind", ["json", "html", "terminal"])
@pytest.mark.parametrize("change", ["unchanged", "newly-visible", "scored-down"])
def test_unfiltered_rollup_does_not_erase_a_diff_identity(tmp_path, build, output_kind, change):
    """Use real GWAS annotation, scoring, filtering, rollup and CLI rendering."""
    data = tmp_path / "data"
    source = tmp_path / "sample.vcf"
    _cache(data, "gwas", build)
    _input(source, build)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute(
            "UPDATE gwas_associations SET trait='Synthetic parent (PheCode 250)', "
            "p_value=1e-30, or_beta=3.5"
        )
        conn.execute(
            "INSERT INTO gwas_associations "
            "(rsid,risk_allele,trait,p_value,or_beta,ci_text,gene,study_accession,"
            "pubmed_id,trait_category) VALUES (?, 'G','Synthetic child (PheCode 250.2)',"
            "1e-40,1.5,'[1.2-1.8]','GENE_ORIGINAL','GCST900252302','900252302','disease')",
            (RSID,),
        )
        conn.execute("UPDATE database_versions SET record_count=2")
        conn.commit()

    args = _args("analyze", source, data, build)
    baseline = tmp_path / "baseline.json"
    runner = CliRunner(env={"COLUMNS": "240", "NO_COLOR": "1"})
    result = runner.invoke(
        main,
        [
            *args,
            "--min-magnitude",
            "9" if change == "newly-visible" else "8",
            "--output",
            str(baseline),
        ],
    )
    assert result.exit_code == 0, result.output
    old = json.loads(baseline.read_text())["annotations"]
    assert [(a["condition"], a["magnitude"]) for a in old] == (
        [] if change == "newly-visible" else [("Synthetic parent", 8.0)]
    )
    if change == "scored-down":
        with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
            conn.execute("UPDATE gwas_associations SET or_beta=1.5 WHERE p_value=1e-30")
            conn.commit()

    output = tmp_path / f"current.{output_kind}"
    current_args = [*args, "--min-magnitude", "8", "--diff", str(baseline)]
    if output_kind != "terminal":
        current_args += ["--output", str(output)]
    result = runner.invoke(main, current_args)
    assert result.exit_code == 0, result.output
    if output_kind == "json":
        payload = json.loads(output.read_text())
        assert [(a["condition"], a["magnitude"]) for a in payload["annotations"]] == (
            [] if change == "scored-down" else [("Synthetic parent", 8.0)]
        )
        diff = payload["diff"]
        expected = {
            "unchanged": (0, 0, 0),
            "newly-visible": (0, 1, 0),
            "scored-down": (1, 0, 0),
        }[change]
        assert tuple(len(diff[key]) for key in ("changed", "new", "removed")) == expected
        if change == "scored-down":
            row = diff["changed"][0]
            assert row["condition"] == "Synthetic parent"
            assert (row["previous_magnitude"], row["magnitude"]) == (8.0, 7.0)
            assert row["below_display_floor"] is True
    else:
        text = (
            _visible_text(output.read_text())
            if output_kind == "html"
            else " ".join(result.output.split())
        )
        expected = {
            "unchanged": "No changes since previous report.",
            "newly-visible": "1 new",
            "scored-down": "1 changed",
        }[change]
        assert expected in text and "1 removed" not in text, text
        if change == "scored-down":
            assert "8.0" in text and "7.0" in text
            assert "below display floor" in text or "below the display floor" in text

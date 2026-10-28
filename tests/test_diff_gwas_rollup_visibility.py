# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Rollup suppression must not be described as failure of a display floor."""

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
def test_rollup_hidden_finding_is_not_below_a_zero_display_floor(tmp_path, build, output_kind):
    """Both GWAS rows clear every floor; only the PheCode rollup hides the parent."""
    data, source = tmp_path / "data", tmp_path / "sample.vcf"
    _cache(data, "gwas", build)
    _input(source, build)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute(
            "UPDATE gwas_associations SET trait='Synthetic parent (PheCode 250)', "
            "p_value=1e-30, or_beta=3.5"
        )
        conn.commit()
    args = _args("analyze", source, data, build)  # Global and source floors are all zero.
    runner = CliRunner(env={"COLUMNS": "240", "NO_COLOR": "1"})
    baseline = tmp_path / "baseline.json"
    first = runner.invoke(main, [*args, "--output", str(baseline)])
    assert first.exit_code == 0, first.output
    old = json.loads(baseline.read_text())["annotations"]
    assert [(a["condition"], a["magnitude"]) for a in old] == [("Synthetic parent", 8.0)]
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        # Changing p-value changes the real generated prose, but not magnitude.
        conn.execute("UPDATE gwas_associations SET p_value=1e-29")
        conn.execute(
            "INSERT INTO gwas_associations "
            "(rsid,risk_allele,trait,p_value,or_beta,ci_text,gene,study_accession,"
            "pubmed_id,trait_category) VALUES (?, 'G','Synthetic child (PheCode 250.2)',"
            "1e-40,1.5,'[1.2-1.8]','GENE_ORIGINAL','GCST900252302','900252302','disease')",
            (RSID,),
        )
        conn.execute("UPDATE database_versions SET record_count=2")
        conn.commit()
    output = tmp_path / f"current.{output_kind}"
    current_args = [*args, "--diff", str(baseline)]
    if output_kind != "terminal":
        current_args += ["--output", str(output)]
    current = runner.invoke(main, current_args)
    assert current.exit_code == 0, current.output
    if output_kind == "json":
        payload = json.loads(output.read_text())
        assert [(a["condition"], a["magnitude"]) for a in payload["annotations"]] == [
            ("Synthetic child", 7.0)
        ]
        diff = payload["diff"]
        assert (len(diff["changed"]), len(diff["new"]), len(diff["removed"])) == (1, 1, 0)
        row = diff["changed"][0]
        assert row["condition"] == "Synthetic parent"
        assert (row["previous_magnitude"], row["magnitude"]) == (8.0, 8.0)
        assert row.get("below_display_floor", False) is False
    else:
        text = (
            _visible_text(output.read_text())
            if output_kind == "html"
            else " ".join(current.output.split())
        )
        assert "1 changed" in text and "1 new" in text and "1 removed" not in text
        assert "below display floor" not in text and "below the display floor" not in text

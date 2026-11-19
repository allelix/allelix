# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The HTML report's ClinVar review stars follow ClinVar's own rubric."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import _review_stars, render_html

# Every aggregate and per-submission review status ClinVar defines, and its stars.
RUBRIC = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the single variant": 0,
    "no classification for the individual variant": 0,
    "no interpretation for the single variant": 0,
    "flagged submission": 0,
}


@pytest.mark.parametrize("status,stars", RUBRIC.items())
@pytest.mark.parametrize(
    "spell",
    [str, lambda s: s.replace(" ", "_"), str.capitalize],
    ids=["clinvar", "underscored", "capitalized"],
)
def test_every_status_gets_clinvar_stars(status, stars, spell):
    assert _review_stars(spell(status)) == "★" * stars + "☆" * (4 - stars)


@pytest.mark.parametrize(
    "status", ["", None, "reviewed by an expert", "criteria provided", "multiple submitters"]
)
def test_unknown_status_gets_no_stars(status):
    assert _review_stars(status) == ""


def _report_rows(tmp_path: Path, review_status: str) -> list[dict]:
    ann = Annotation(
        source="clinvar",
        rsid="rs4646437",
        significance="clinvar_drug_response",
        category="clinical",
        magnitude=6.0,
        description="ClinVar classifies this allele as drug response",
        attribution="ClinVar",
        genotype_match="CT",
        gene="CYP3A4",
        condition="tacrolimus response",
        review_status=review_status,
    )
    result = AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="STARS",
        build="GRCh37",
        total_variants=1,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101")],
        annotations=[ann],
    )
    out = tmp_path / "report.html"
    render_html(result, output_path=out)
    body = out.read_text()
    start = body.index(">", body.index('id="variant-data"')) + 1
    data = json.loads(body[start : body.index("</script>", start)].replace("<\\/", "</"))
    return [a for entry in data for a in entry["annotations"] if a["source"] == "clinvar"]


@pytest.mark.parametrize("status", ["reviewed by an expert", "criteria provided"])
def test_unknown_status_row_shows_its_text_without_stars(tmp_path: Path, status):
    (row,) = _report_rows(tmp_path, status)
    assert row["reviewStatus"] == status
    assert row["reviewStars"] == status


def test_expert_panel_row_shows_three_stars(tmp_path: Path):
    (row,) = _report_rows(tmp_path, "reviewed by expert panel")
    assert row["reviewStatus"] == "reviewed by expert panel"
    assert row["reviewStars"] == "★★★☆"

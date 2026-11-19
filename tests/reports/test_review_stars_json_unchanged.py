# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The review-star fix changes the HTML report only: JSON stays byte-identical."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest

import allelix.reports.html as html_report
import allelix.reports.json_report as json_report
from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from tests.reports.test_review_stars import RUBRIC


def _review_stars_270(status: int | str | None) -> str:
    """The 2.7.0 renderer, which guessed stars from the status text."""
    if status is None or status == "" or status == "—":
        return ""
    if isinstance(status, str):
        star_count = status.count("_") if "_" in status else 0
        if "expert" in status.lower() or "practice" in status.lower():
            star_count = 4
        elif "multiple" in status.lower():
            star_count = 3
        elif "single" in status.lower():
            star_count = 1
        elif "no_assertion" in status.lower() or "no assertion" in status.lower():
            star_count = 0
    else:
        star_count = int(status)
    star_count = max(0, min(star_count, 4))
    return "★" * star_count + "☆" * (4 - star_count)


STATUSES = [
    *RUBRIC,
    *(s.replace(" ", "_") for s in RUBRIC),
    "reviewed by an expert",
    "criteria provided",
    "",
]


def _result() -> AnalysisResult:
    annotations = [
        Annotation(
            source="clinvar",
            rsid=f"rs90060700{i:02d}",
            significance="clinvar_pathogenic",
            category="clinical",
            magnitude=9.0,
            description="ClinVar classifies this allele as Pathogenic",
            attribution="ClinVar",
            genotype_match="AG",
            gene="SYNGENE",
            condition="synthetic condition",
            review_status=status,
        )
        for i, status in enumerate(STATUSES)
    ]
    return AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="STARS_JSON",
        build="GRCh37",
        total_variants=len(annotations),
        skipped_count=0,
        annotators_used=[("clinvar", "20260101")],
        annotations=annotations,
    )


class _FixedClock(datetime):
    @classmethod
    def now(cls, tz=None):  # type: ignore[override]
        return datetime(2026, 11, 13, 12, 0, 0, tzinfo=UTC)


@pytest.mark.parametrize("min_magnitude", [0.0, 5.0])
def test_json_report_is_byte_identical_with_and_without_the_star_fix(
    tmp_path, monkeypatch, min_magnitude
):
    monkeypatch.setattr(json_report, "datetime", _FixedClock)
    outputs = {}
    pages = {}
    for name, stars in (("2.7.0", _review_stars_270), ("fixed", html_report._review_stars)):
        monkeypatch.setattr(html_report, "_review_stars", stars)
        out = tmp_path / f"{name}.json"
        render_json(_result(), output_path=out, min_magnitude=min_magnitude)
        outputs[name] = out.read_bytes()
        page = tmp_path / f"{name}.html"
        render_html(_result(), output_path=page)
        pages[name] = page.read_text()
    assert outputs["2.7.0"] == outputs["fixed"]
    # The switch is real: the HTML pages differ where the stars do.
    assert pages["2.7.0"] != pages["fixed"]

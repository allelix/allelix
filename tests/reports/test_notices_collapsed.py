# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The top of the HTML report: action warnings visible, other notices collapsed.

The short disclaimer and the warnings that change what a reader should do
stay outside ``<details>``. Every other notice is inside "Technical
details", word for word as before, so nothing goes silent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from allelix.parsers.vcf import VcfQualityStats
from allelix.reports import REGULATORY_NOTICE, vcf_quality_lines
from allelix.reports._pipeline import AnalysisResult, BuildDiagnostics
from allelix.reports._sources import SKIPPED, ReportWarning, SourceRecord
from allelix.reports.association_text import strand_abstention_lines
from allelix.reports.clinvar_text import common_allele_line, scored_down_line
from allelix.reports.html import (
    _DETAILS_SUMMARY_PLAIN,
    _DETAILS_UNCHECKED,
    _SHORT_DISCLAIMER,
    _escape,
    render_html,
)
from allelix.utils.allele import STRAND_ABSTAIN_REASONS

_WARNINGS = [
    ReportWarning(
        "clinvar_skipped", "ClinVar was not used, so this report has no ClinVar findings."
    ),
    ReportWarning("truncated_input", "The input file may be cut short."),
    ReportWarning("forced_format", "The file was read with a forced --format."),
    ReportWarning("multiple_homozygous_plp", "3 variants are called homozygous for an allele."),
]


def _result(**fields) -> AnalysisResult:
    result = AnalysisResult(
        file_path=Path("genotype.vcf"),
        parser_name="vcf",
        parser_display_name="VCF",
        sample_id="S1",
        build="GRCh38",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101")],
    )
    for name, value in fields.items():
        setattr(result, name, value)
    return result


def _full() -> AnalysisResult:
    """Every notice the top of the report can show, at once."""
    return _result(
        report_warnings=list(_WARNINGS),
        build_diagnostics=BuildDiagnostics(
            effective_build="GRCh38",
            source="seeds",
            file_label="GRCh37",
            detected_build="GRCh38",
            matched_count=5,
            inspected_count=5,
            position_confident=True,
            label_overridden=True,
        ),
        sources=[
            SourceRecord("gwas", "GWAS Catalog", SKIPPED, "missing", message="GWAS: missing."),
            SourceRecord("cadd", "CADD", SKIPPED, "disabled"),
        ],
        input_method="low-pass-wgs-imputed",
        vcf_quality=VcfQualityStats(filtered_rows=7, filtered_by_label={"LowQual": 7}),
        clinvar_scored_down=12,
        clinvar_unrecognised=2,
        clinvar_common_allele=1,
        strand_abstentions={"gwas": {STRAND_ABSTAIN_REASONS[0]: 4}},
    )


def _split(page: str) -> tuple[str, str]:
    """(the page without the details section, the details section)."""
    start = page.index('<details class="tech-details">')
    end = page.index("</details>", start) + len("</details>")
    return page[:start] + page[end:], page[start:end]


def _render(result: AnalysisResult, tmp_path: Path, **kwargs) -> str:
    out = tmp_path / "r.html"
    render_html(result, output_path=out, **kwargs)
    return out.read_text(encoding="utf-8")


@pytest.fixture
def page(tmp_path: Path) -> str:
    return _render(_full(), tmp_path, high_value_no_calls=["rs6025 (F5)"])


def test_one_details_section(page: str) -> None:
    assert page.count('<details class="tech-details">') == 1
    assert "<summary>Technical details</summary>" in page


def test_action_warnings_are_visible(page: str) -> None:
    outside, inside = _split(page)
    for warning in _WARNINGS:
        assert f"<strong>Warning.</strong> {warning.message}" in outside
        assert warning.message not in inside
    assert "Build mismatch" in outside
    assert "<strong>Sources not used.</strong><ul><li>GWAS: missing.</li>" in outside
    assert "Sources not used" not in inside
    assert "GWAS: missing." not in inside
    assert "<strong>High-value no-calls.</strong>" in outside
    assert "rs6025 (F5)" in outside


def test_short_disclaimer_visible_full_notice_inside(page: str) -> None:
    outside, inside = _split(page)
    assert len(_SHORT_DISCLAIMER.split(". ")) == 2
    assert _SHORT_DISCLAIMER in outside
    assert REGULATORY_NOTICE not in outside
    assert f"<strong>Informational only.</strong> {REGULATORY_NOTICE}" in inside


def test_counts_and_notes_are_inside(page: str) -> None:
    outside, inside = _split(page)
    texts = [
        *vcf_quality_lines(_full().vcf_quality),
        scored_down_line(12),
        common_allele_line(1),
        "2 ClinVar row(s) not shown: their classification is not one Allelix scores yet.",
        *strand_abstention_lines(_full().strand_abstentions),
        "Not run by choice: CADD (turned off in config).",
    ]
    texts = [_escape(text) for text in texts]
    for text in texts:
        assert text and text in inside, text
        assert text not in outside, text


def test_imputed_input_caution_is_visible(page: str) -> None:
    outside, inside = _split(page)
    assert "<strong>Input method.</strong>" in outside
    assert "confirm any finding by direct genotyping" in outside
    assert "<strong>Input method.</strong>" not in inside


def test_measured_input_method_is_inside(tmp_path: Path) -> None:
    page = _render(_result(input_method="wgs-30x"), tmp_path)
    outside, inside = _split(page)
    assert "<strong>Input method.</strong> Wgs-30x.</div>" in inside
    assert "<strong>Input method.</strong>" not in outside


def test_summary_line_sits_above_the_details(page: str) -> None:
    line = (
        '<p class="details-summary">Some calls could not be checked. '
        "13 ClinVar rows were scored down. 2 ClinVar rows were not shown. "
        "Details below.</p>"
    )
    assert line in page
    assert page.index(line) < page.index('<details class="tech-details">')


def test_plain_summary_without_counts(tmp_path: Path) -> None:
    page = _render(_result(), tmp_path)
    outside, inside = _split(page)
    assert _DETAILS_SUMMARY_PLAIN in outside
    assert _DETAILS_UNCHECKED not in page
    assert "scored down" not in page
    assert REGULATORY_NOTICE in inside
    assert 'class="notice-warn"' not in outside


def test_details_closed_by_default(page: str) -> None:
    assert '<details class="tech-details" open' not in page


def test_notice_reader_uses_utf8_with_a_cp1252_platform_default(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    read_text = Path.read_text

    def platform_reader(path, encoding=None, errors=None):
        return read_text(path, encoding=encoding or "cp1252", errors=errors)

    monkeypatch.setattr(Path, "read_text", platform_reader)
    page = _render(_full(), tmp_path)
    test_short_disclaimer_visible_full_notice_inside(page)
    plain = _render(_result(), tmp_path)
    assert REGULATORY_NOTICE in _split(plain)[1]

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The top of the HTML report says a source failure once, and the line above
"Technical details" names only what the section holds."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from allelix.parsers.base import ArrayInputStats
from allelix.parsers.vcf import VcfQualityStats
from allelix.reports import input_left_unchecked, vcf_quality_lines
from allelix.reports._pipeline import AnalysisResult
from allelix.reports._sources import (
    READY_NOT_USED,
    SKIPPED,
    USED,
    SourceRecord,
    clinvar_skipped_warning,
)
from allelix.reports.html import (
    _DETAILS_SUMMARY_PLAIN,
    _DETAILS_UNCHECKED,
    _escape,
    render_html,
)
from allelix.utils.allele import STRAND_ABSTAIN_REASONS

_CLINVAR_OUTDATED = SourceRecord(
    "clinvar",
    "ClinVar",
    SKIPPED,
    "outdated",
    message="ClinVar: out of date (built for iv:5, needs iv:9). Run `allelix db update`.",
)
_GWAS_MISSING = SourceRecord("gwas", "GWAS Catalog", SKIPPED, "missing", message="GWAS: missing.")


def _result(**fields) -> AnalysisResult:
    result = AnalysisResult(
        file_path=Path("genotype.vcf"),
        parser_name="vcf",
        parser_display_name="VCF",
        sample_id="S1",
        build="GRCh38",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("gwas", "20260101")],
    )
    for name, value in fields.items():
        setattr(result, name, value)
    return result


def _render(result: AnalysisResult, tmp_path: Path) -> str:
    out = tmp_path / "r.html"
    render_html(result, output_path=out)
    return out.read_text()


def _split(page: str) -> tuple[str, str]:
    """(the page without the details section, the details section)."""
    start = page.index('<details class="tech-details">')
    end = page.index("</details>", start) + len("</details>")
    return page[:start] + page[end:], page[start:end]


def _summary(page: str) -> str:
    start = page.index('<p class="details-summary">') + len('<p class="details-summary">')
    return page[start : page.index("</p>", start)]


def _clinvar_not_used(sources: list[SourceRecord]) -> AnalysisResult:
    warning = clinvar_skipped_warning(sources)
    assert warning is not None
    return _result(sources=sources, report_warnings=[warning])


# --- The ClinVar-not-used message, once ----------------------------------------


def test_clinvar_failure_is_said_once(tmp_path: Path) -> None:
    page = _render(_clinvar_not_used([_CLINVAR_OUTDATED]), tmp_path)
    outside, _inside = _split(page)
    message = _escape(_CLINVAR_OUTDATED.message)
    assert page.count(message) == 1
    assert "<strong>Warning.</strong> ClinVar was not used" in outside
    assert message in outside
    # Its only entry is in the warning, so there is no "Sources not used" box.
    assert "Sources not used" not in page


def test_other_failed_sources_keep_their_box(tmp_path: Path) -> None:
    page = _render(_clinvar_not_used([_CLINVAR_OUTDATED, _GWAS_MISSING]), tmp_path)
    outside, _inside = _split(page)
    assert page.count(_escape(_CLINVAR_OUTDATED.message)) == 1
    assert "<strong>Sources not used.</strong><ul><li>GWAS: missing.</li></ul>" in outside


def test_failed_clinvar_without_the_warning_stays_in_the_box(tmp_path: Path) -> None:
    # The box leaves a source out only when a warning on the page names it.
    page = _render(_result(sources=[_CLINVAR_OUTDATED]), tmp_path)
    outside, _inside = _split(page)
    assert f"<li>{_escape(_CLINVAR_OUTDATED.message)}</li>" in outside
    assert page.count(_escape(_CLINVAR_OUTDATED.message)) == 1


def test_only_failed_sources_are_listed(tmp_path: Path) -> None:
    sources = [
        SourceRecord("clinvar", "ClinVar", USED, version="2026-10-10"),
        SourceRecord("pharmgkb", "ClinPGx", READY_NOT_USED),
        SourceRecord("cadd", "CADD", SKIPPED, "disabled"),
        _GWAS_MISSING,
    ]
    page = _render(_result(sources=sources), tmp_path)
    outside, _inside = _split(page)
    start = outside.index("<strong>Sources not used.</strong>")
    box = outside[start : outside.index("</div>", start)]
    assert box == "<strong>Sources not used.</strong><ul><li>GWAS: missing.</li></ul>"


# --- The summary line ---------------------------------------------------------


def test_only_general_notes_give_the_neutral_line(tmp_path: Path) -> None:
    # As on the outdated-ClinVar mock page: a plain VCF, nothing scored down.
    plain = _clinvar_not_used([_CLINVAR_OUTDATED])
    plain.vcf_quality = VcfQualityStats()
    page = _render(plain, tmp_path)
    _outside, inside = _split(page)
    assert _escape(vcf_quality_lines(plain.vcf_quality)[0]) in inside
    assert _summary(page) == _DETAILS_SUMMARY_PLAIN
    assert "scored down" not in _summary(page)
    assert _DETAILS_UNCHECKED not in page


@pytest.mark.parametrize(
    ("fields", "expected"),
    [
        (
            {"vcf_quality": VcfQualityStats(is_gvcf=True, low_quality_calls=3)},
            "Some calls could not be checked. Details below.",
        ),
        ({"clinvar_scored_down": 1}, "1 ClinVar row was scored down. Details below."),
        ({"clinvar_common_allele": 2}, "2 ClinVar rows were scored down. Details below."),
        (
            {"clinvar_scored_down": 28, "clinvar_common_allele": 1},
            "29 ClinVar rows were scored down. Details below.",
        ),
        ({"clinvar_unrecognised": 1}, "1 ClinVar row was not shown. Details below."),
        (
            {"strand_abstentions": {"snpedia": {STRAND_ABSTAIN_REASONS[0]: 5}}},
            "Some calls could not be checked. Details below.",
        ),
        (
            {"probe_classes": {"i_probes_at_plp_loci": 4}},
            "Some calls could not be checked. Details below.",
        ),
        (
            {"array_quality": ArrayInputStats(di_calls=2)},
            "Some calls could not be checked. Details below.",
        ),
        (
            {
                "vcf_quality": VcfQualityStats(low_quality_calls=1),
                "clinvar_scored_down": 28,
                "clinvar_common_allele": 1,
                "strand_abstentions": {"gwas": {STRAND_ABSTAIN_REASONS[0]: 4}},
            },
            "Some calls could not be checked. 29 ClinVar rows were scored down. Details below.",
        ),
    ],
)
def test_summary_names_what_the_section_holds(tmp_path: Path, fields, expected) -> None:
    page = _render(_result(**fields), tmp_path)
    assert _summary(page) == expected


@pytest.mark.parametrize(
    "stats",
    [
        VcfQualityStats(),  # no reference blocks: a note on the file
        VcfQualityStats(is_gvcf=True, filtered_rows=3, kept_filtered=True),
        VcfQualityStats(is_gvcf=True, trimmed_matches=2),
        VcfQualityStats(is_gvcf=True, reference_carried_calls=2),
    ],
)
def test_descriptive_notes_alone_give_the_neutral_line(tmp_path: Path, stats) -> None:
    page = _render(_result(vcf_quality=stats), tmp_path)
    _outside, inside = _split(page)
    assert "Input quality" in inside
    assert _summary(page) == _DETAILS_SUMMARY_PLAIN


# Every count field of VcfQualityStats is either a call not checked or a
# note on the file. A new field fails here until it is put in one group.
_UNCHECKED_FIELDS = {
    "filtered_rows",
    "low_quality_calls",
    "spanning_deletion_calls",
    "half_calls",
    "clinvar_ref_mismatches",
    "reference_source_conflicts",
    "reference_ploidy_conflicts",
    "reference_carried_not_in_file",
    "malformed_rows",
    "no_genotype_rows",
    "unusable_genotype_calls",
    "symbolic_allele_calls",
    "symbolic_rows",
}
_DESCRIPTIVE_FIELDS = {
    "trimmed_matches",
    "reference_carried_calls",
    "genotype_records",
    "hom_ref_calls",
    "non_reference_calls",
}


def _int_fields() -> set[str]:
    return {
        f.name
        for f in dataclasses.fields(VcfQualityStats)
        if f.type in ("int", int) and not f.name.startswith("min_")
    }


def test_every_count_field_is_classified() -> None:
    assert _int_fields() == _UNCHECKED_FIELDS | _DESCRIPTIVE_FIELDS


@pytest.mark.parametrize("name", sorted(_UNCHECKED_FIELDS))
def test_unchecked_field(name: str) -> None:
    stats = dataclasses.replace(VcfQualityStats(is_gvcf=True), **{name: 1})
    assert input_left_unchecked(stats, None, None)


@pytest.mark.parametrize("name", sorted(_DESCRIPTIVE_FIELDS))
def test_descriptive_field(name: str) -> None:
    stats = dataclasses.replace(VcfQualityStats(is_gvcf=True), **{name: 1})
    assert not input_left_unchecked(stats, None, None)


def test_kept_filtered_rows_are_checked() -> None:
    stats = VcfQualityStats(is_gvcf=True, filtered_rows=3, kept_filtered=True)
    assert not input_left_unchecked(stats, None, None)

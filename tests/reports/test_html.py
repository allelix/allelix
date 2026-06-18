# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the HTML report renderer."""

from __future__ import annotations

import json
import re
from pathlib import Path

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult, BuildDiagnostics
from allelix.reports.html import render_html


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="MHG_HTML",
        build="GRCh37",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101")],
        annotations=annotations,
    )


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "ClinVar classifies this allele as Pathogenic",
        "attribution": "ClinVar",
        "genotype_match": "AG",
        "gene": "MTHFR",
        "condition": "MTHFR deficiency",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


# ---- Existing tests (updated for new structure) ----


class TestPanelCoverageHtml:
    """GH #75: HTML report grows a Panel Coverage section near the top
    when --filter-file supplied a rsid panel."""

    def _result_with_panel(
        self,
        annotations: list[Annotation],
        *,
        panel_rsids: frozenset[str],
        genotyped_panel_rsids: frozenset[str],
        panel_annotated_rsids: frozenset[str],
    ) -> AnalysisResult:
        base = _result(annotations)
        base.panel_rsids = panel_rsids
        base.genotyped_panel_rsids = genotyped_panel_rsids
        base.panel_annotated_rsids = panel_annotated_rsids
        return base

    def test_no_panel_omits_section(self, tmp_path: Path):
        out = tmp_path / "r.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Panel coverage" not in body

    def test_section_includes_counts_and_missing(self, tmp_path: Path):
        out = tmp_path / "r.html"
        result = self._result_with_panel(
            [_ann(rsid="rs1801133")],
            panel_rsids=frozenset({"rs1801133", "rs4680", "rs999"}),
            genotyped_panel_rsids=frozenset({"rs1801133", "rs4680"}),
            panel_annotated_rsids=frozenset({"rs1801133"}),
        )
        render_html(result, output_path=out)
        body = out.read_text()
        assert "Panel coverage" in body
        assert "2/3 requested rsIDs genotyped" in body
        # State 3: missing rsID surfaced with chip note.
        assert "rs999" in body
        assert "not genotyped on this chip" in body
        # State 2: genotyped-but-no-findings rsID surfaced.
        assert "rs4680" in body
        assert "no annotations matched" in body

    def test_panel_gt_column_renders_when_genotypes_supplied(self, tmp_path: Path):
        """GH #120: when `result.panel_genotypes` is supplied (the
        rsid → "A/G" map run_analysis collected during the batch
        flush), the HTML Panel Coverage table grows a GT column.
        Genotyped rsIDs show their literal call; state-3 rows show
        em dash; no-call rows surface their dash genotype (-/-, A/-)
        so the no-call is visually obvious even though the v6
        accounting still folds it into no_findings.
        """
        result = self._result_with_panel(
            [_ann(rsid="rs1801133")],
            panel_rsids=frozenset({"rs1801133", "rs4680", "rs_nocall", "rs999"}),
            genotyped_panel_rsids=frozenset({"rs1801133", "rs4680", "rs_nocall"}),
            panel_annotated_rsids=frozenset({"rs1801133"}),
        )
        result.panel_genotypes = {
            "rs1801133": "A/G",
            "rs4680": "C/T",
            "rs_nocall": "-/-",
        }
        out = tmp_path / "r.html"
        render_html(result, output_path=out)
        body = out.read_text()
        # GT column header present.
        assert "<th>GT</th>" in body
        # State 2 row: rs4680 with its genotype rendered.
        assert "C/T" in body
        # No-call row: literal -/- visible to the user.
        assert "rs_nocall" in body
        assert "-/-" in body
        # State 3 row: em dash in the GT cell.
        assert "rs999" in body
        # Em dash specifically — easier to grep without false-positives.
        assert "<td>—</td>" in body

    def test_panel_gt_column_omitted_when_no_genotype_map(self, tmp_path: Path):
        """Back-compat: a result without panel_genotypes (e.g. older
        callers, manual AnalysisResult construction) renders the
        original two-column shape with no GT column."""
        result = self._result_with_panel(
            [_ann(rsid="rs1801133")],
            panel_rsids=frozenset({"rs1801133", "rs4680", "rs999"}),
            genotyped_panel_rsids=frozenset({"rs1801133", "rs4680"}),
            panel_annotated_rsids=frozenset({"rs1801133"}),
        )
        # panel_genotypes deliberately left unset (None).
        out = tmp_path / "r.html"
        render_html(result, output_path=out)
        body = out.read_text()
        assert "Panel coverage" in body
        # GT column header NOT present.
        assert "<th>GT</th>" not in body


class TestRenderHtml:
    def test_writes_self_contained_file(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert body.startswith("<!DOCTYPE html>")
        assert "<style>" in body
        assert '<link rel="stylesheet"' not in body
        assert "<script src=" not in body

    def test_disclaimer_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Informational only" in body
        assert "not medical advice" in body

    def test_attribution_visible(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "ClinVar" in body
        assert "rs1801133" in body
        assert "MTHFR" in body

    def test_html_escapes_user_supplied_strings(self, tmp_path: Path):
        evil = _ann(condition='<script>alert("xss")</script>')
        out = tmp_path / "report.html"
        render_html(_result([evil]), output_path=out)
        body = out.read_text()
        html_part = body.split('<script id="variant-data"')[0]
        assert "<script>alert" not in html_part
        assert "&lt;script&gt;" in body
        assert "</script>" not in body.split('<script id="variant-data"')[1].split("</script>")[0]

    def test_empty_annotations(self, tmp_path: Path):
        out = tmp_path / "report.html"
        count = render_html(_result([]), output_path=out)
        assert count == 0
        body = out.read_text()
        assert "No annotations" in body

    def test_min_magnitude_filter(self, tmp_path: Path):
        anns = [_ann(rsid="lo", magnitude=2.0), _ann(rsid="hi", magnitude=8.0)]
        out = tmp_path / "report.html"
        count = render_html(_result(anns), output_path=out, min_magnitude=5.0)
        assert count == 1
        body = out.read_text()
        assert "hi" in body
        assert ">lo<" not in body

    def test_genes_filter(self, tmp_path: Path):
        anns = [_ann(rsid="m", gene="MTHFR"), _ann(rsid="b", gene="BRCA1")]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out, genes={"MTHFR"})
        body = out.read_text()
        assert "MTHFR" in body
        assert ">BRCA1<" not in body


class TestEducationSection:
    def test_education_section_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "Reading This Report" in body

    def test_education_pseudogene_warning(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert "cross-hybridize" in out.read_text()

    def test_education_carrier_vs_affected(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert "Carrier vs. affected" in out.read_text()

    def test_education_confirmatory_testing(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert "Confirmatory testing" in out.read_text()


class TestBuildMismatchBanner:
    def test_no_banner_without_mismatch(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh37",
            effective_build="GRCh37",
            override=False,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        assert "Build mismatch" not in out.read_text()

    def test_banner_on_mismatch(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh38",
            effective_build="GRCh38",
            override=False,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "Build mismatch" in body
        assert "GRCh37" in body
        assert "GRCh38" in body
        assert "notice-warn" in body

    def test_no_banner_with_override(self, tmp_path: Path):
        r = _result([_ann()])
        r.build_diagnostics = BuildDiagnostics(
            header_build="GRCh37",
            detected_build="GRCh38",
            effective_build="GRCh38",
            override=True,
            matched_count=5,
            inspected_count=5,
        )
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        assert "Build mismatch" not in out.read_text()

    def test_no_banner_without_diagnostics(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert "Build mismatch" not in out.read_text()


class TestLicenseAttributions:
    def _result_with_annotators(self, annotators: list[tuple[str, str | None]]) -> AnalysisResult:
        return AnalysisResult(
            file_path=Path("genotype.txt"),
            parser_name="myhappygenes",
            parser_display_name="MyHappyGenes (Tempus)",
            sample_id="MHG_LICENSE",
            build="GRCh37",
            total_variants=10,
            skipped_count=0,
            annotators_used=annotators,
            annotations=[_ann()],
        )

    def test_pharmgkb_attribution_present(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("pharmgkb", "2026-01")])
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "ClinPGx" in body
        assert "CC BY-SA 4.0" in body

    def test_snpedia_attribution_present(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("snpedia", None)])
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "SNPedia" in body
        assert "CC BY-NC-SA 3.0 US" in body

    def test_gnomad_attribution_present(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("gnomad", "4.1")])
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "ODbL" in body
        assert "gnomAD" in body

    def test_gwas_attribution_present(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("gwas", "2026-01")])
        out = tmp_path / "report.html"
        render_html(r, output_path=out)
        body = out.read_text()
        assert "NHGRI-EBI GWAS Catalog" in body


class TestReputeClassification:
    def test_pathogenic_is_bad(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann(significance="clinvar_pathogenic")]), output_path=out)
        body = out.read_text()
        assert "badge-bad" in body
        assert "pill-bad" in body

    def test_benign_is_good(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann(significance="clinvar_benign")]), output_path=out)
        body = out.read_text()
        assert "badge-good" in body
        assert "pill-good" in body

    def test_vus_is_neutral(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(
            _result([_ann(significance="clinvar_uncertain_significance")]),
            output_path=out,
        )
        body = out.read_text()
        assert "badge-neutral" in body
        assert "pill-neutral" in body

    def test_snpedia_bad(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann(significance="snpedia_bad")]), output_path=out)
        assert "badge-bad" in out.read_text()

    def test_snpedia_good(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann(significance="snpedia_good")]), output_path=out)
        assert "badge-good" in out.read_text()


# ---- New tests for v1.8.0 redesign ----


class TestHtmlFiveColumns:
    def test_html_five_columns(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        thead = body.split("<thead>")[1].split("</thead>")[0]
        th_count = thead.count("<th")
        assert th_count == 5


class TestHtmlNoHorizontalScroll:
    def test_no_overflow_x(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert "overflow-x" not in out.read_text()


class TestHtmlVariantGrouping:
    def test_two_annotations_same_variant_one_row(self, tmp_path: Path):
        anns = [
            _ann(source="clinvar", attribution="ClinVar", magnitude=9.0),
            _ann(
                source="snpedia", attribution="SNPedia", magnitude=4.0, significance="snpedia_bad"
            ),
        ]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out)
        body = out.read_text()
        tbody = body.split("<tbody>")[1].split("</tbody>")[0]
        tr_count = tbody.count("<tr ")
        assert tr_count == 1

    def test_different_variants_different_rows(self, tmp_path: Path):
        anns = [
            _ann(rsid="rs1", genotype_match="AG"),
            _ann(rsid="rs2", genotype_match="CC"),
        ]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out)
        body = out.read_text()
        tbody = body.split("<tbody>")[1].split("</tbody>")[0]
        tr_count = tbody.count("<tr ")
        assert tr_count == 2


class TestHtmlDataAttributes:
    def test_data_attributes_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "data-row-id=" in body
        assert "data-magnitude=" in body
        assert "data-gene=" in body
        assert "data-genotype=" in body
        assert "data-repute=" in body
        assert "data-search-text=" in body


class TestHtmlVariantDataJson:
    def _extract_json(self, html_body: str) -> list:
        start = html_body.index('id="variant-data"')
        json_start = html_body.index(">", start) + 1
        json_end = html_body.index("</script>", json_start)
        raw = html_body[json_start:json_end].replace("<\\/", "</")
        return json.loads(raw)

    def test_variant_data_json(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        data = self._extract_json(out.read_text())
        assert isinstance(data, list)
        assert len(data) == 1
        entry = data[0]
        assert "rsid" in entry
        assert "gene" in entry
        assert "genotype" in entry
        assert "zygosity" in entry
        assert "annotations" in entry
        assert len(entry["annotations"]) >= 1

    def test_enrichment_fields_numeric(self, tmp_path: Path):
        ann = _ann(
            allele_frequency=0.23,
            am_pathogenicity=0.057,
            am_class="likely_benign",
            cadd_phred=38.0,
        )
        out = tmp_path / "report.html"
        render_html(_result([ann]), output_path=out)
        data = self._extract_json(out.read_text())
        v = data[0]
        assert v["allele_frequency"] == 0.23
        assert isinstance(v["allele_frequency"], float)
        assert v["am_pathogenicity"] == 0.057
        assert isinstance(v["am_pathogenicity"], float)
        assert v["am_class"] == "likely_benign"
        assert v["cadd_phred"] == 38.0
        assert isinstance(v["cadd_phred"], float)
        a = v["annotations"][0]
        assert "allele_frequency" not in a
        assert "am_pathogenicity" not in a
        assert "cadd_phred" not in a


class TestHtmlSearchTextAllSources:
    def test_search_text_contains_all_sources(self, tmp_path: Path):
        anns = [
            _ann(source="clinvar", attribution="ClinVar"),
            _ann(
                source="snpedia", attribution="SNPedia", magnitude=4.0, significance="snpedia_bad"
            ),
        ]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out)
        body = out.read_text()
        # Find data-search-text attribute value
        idx = body.index("data-search-text=")
        search_text = body[idx : idx + 500]
        assert "clinvar" in search_text
        assert "snpedia" in search_text


class TestHtmlSidebarExists:
    def test_sidebar_elements_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert 'id="detail-panel"' in body
        assert 'id="backdrop"' in body


class TestHtmlViewportMeta:
    def test_viewport_meta_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert '<meta name="viewport"' in out.read_text()


class TestHtmlSelfContained:
    def test_no_external_resources(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert '<link rel="stylesheet"' not in body
        assert "<script src=" not in body


class TestHtmlMagnitudeBadge:
    def test_badge_class_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        assert 'class="badge' in out.read_text()


class TestHtmlDefaultSortDesc:
    def test_first_row_has_highest_magnitude(self, tmp_path: Path):
        anns = [
            _ann(rsid="rs_low", magnitude=2.0, genotype_match="AA"),
            _ann(rsid="rs_high", magnitude=8.0, genotype_match="CC"),
        ]
        out = tmp_path / "report.html"
        render_html(_result(anns), output_path=out)
        body = out.read_text()
        tbody = body.split("<tbody>")[1].split("</tbody>")[0]
        first_tr_idx = tbody.index("<tr ")
        second_tr_idx = tbody.index("<tr ", first_tr_idx + 1)
        first_row = tbody[first_tr_idx:second_tr_idx]
        assert 'data-magnitude="8.0"' in first_row


class TestHtmlDarkMode:
    def test_css_custom_properties_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert "prefers-color-scheme: dark" in body
        assert "[data-theme=" in body

    def test_toggle_button_present(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        assert 'id="theme-toggle"' in body
        assert "theme-toggle" in body

    def test_no_hardcoded_colors_in_component_css(self, tmp_path: Path):
        out = tmp_path / "report.html"
        render_html(_result([_ann()]), output_path=out)
        body = out.read_text()
        style_start = body.index("<style>") + len("<style>")
        style_end = body.index("</style>")
        css = body[style_start:style_end]
        for line in css.split("\n"):
            stripped = line.strip()
            if not stripped or stripped.startswith(":root") or stripped.startswith("--"):
                continue
            if "@media (prefers-color-scheme" in stripped:
                continue
            if "[data-theme" in stripped:
                continue
            if "badge-bad" in stripped or "badge-good" in stripped:
                continue
            if "badge-neutral" in stripped:
                continue
            if "pill-bad" in stripped or "pill-good" in stripped:
                continue
            if "card-bad" in stripped or "card-good" in stripped:
                continue
            if ".filter-btn.active" in stripped:
                continue
            if "#search:focus" in stripped:
                continue
            if "box-shadow: inset" in stripped:
                continue
            if "border-color: #1976d2" in stripped:
                continue
            assert "color: #" not in stripped or "var(--" in stripped, (
                f"Hardcoded color found: {stripped}"
            )


class TestHtmlEscaping:
    def test_escaping_in_table_and_data_attrs(self, tmp_path: Path):
        ann = _ann(condition='x < y & "z"')
        out = tmp_path / "report.html"
        render_html(_result([ann]), output_path=out)
        body = out.read_text()
        assert "x &lt; y &amp; &quot;z&quot;" in body
        assert 'x < y & "z"' not in body.split("<script")[0]


def _extract_json(html_body: str) -> list:
    """Extract the variant-data JSON blob from rendered HTML."""
    start = html_body.index('id="variant-data"')
    json_start = html_body.index(">", start) + 1
    json_end = html_body.index("</script>", json_start)
    raw = html_body[json_start:json_end].replace("<\\/", "</")
    return json.loads(raw)


def _extract_js(html_body: str) -> str:
    """Extract the main DOMContentLoaded script block."""
    marker = 'document.addEventListener("DOMContentLoaded"'
    start = html_body.index(marker)
    block_start = html_body.rindex("<script>", 0, start)
    block_end = html_body.index("</script>", start)
    return html_body[block_start:block_end]


def _rich_report(tmp_path: Path) -> str:
    """Render a report with all enrichment fields populated."""
    ann = _ann(
        allele_frequency=0.254,
        am_pathogenicity=0.320,
        am_class="ambiguous",
        cadd_phred=0.1,
        review_status="criteria_provided,_single_submitter",
        references=["clinvar:allele/100001"],
    )
    out = tmp_path / "report.html"
    render_html(_result([ann]), output_path=out)
    return out.read_text()


class TestJsHtmlContract:
    """Verify the JS selectors and field references match what the HTML emits.

    These tests guard against the class of bugs where Python changes a DOM id,
    data attribute, or JSON field name but the JS still references the old one.
    """

    def test_variant_table_id_unique(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        assert body.count('id="variant-table"') == 1

    def test_js_element_ids_exist_in_html(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        js = _extract_js(body)
        ids = re.findall(r'getElementById\("([^"]+)"\)', js)
        assert len(ids) >= 5
        for eid in ids:
            assert f'id="{eid}"' in body, f"JS references #{eid} but HTML has no such id"

    def test_variant_table_tbody_selector(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        assert 'querySelector("#variant-table tbody")' in body
        assert 'id="variant-table"' in body
        start = body.index('id="variant-table"')
        table_region = body[start : start + 2000]
        assert "<tbody>" in table_region

    def test_data_sort_attrs_match_dataset_keys(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        sort_keys = set(re.findall(r'data-sort="([^"]+)"', body))
        assert sort_keys
        js = _extract_js(body)
        for key in sort_keys:
            assert f"dataset.{key}" in js or "dataset[key]" in js, (
                f'data-sort="{key}" in HTML but JS never reads dataset.{key}'
            )

    def test_row_data_attrs_match_js_references(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        js = _extract_js(body)
        js_dataset_refs = set(re.findall(r"\.dataset\.([a-zA-Z]+)", js))
        html_before_script = body.split("<script>")[0]
        for attr in js_dataset_refs:
            kebab = re.sub(r"([A-Z])", r"-\1", attr).lower()
            assert f"data-{kebab}=" in html_before_script, (
                f"JS reads dataset.{attr} but no data-{kebab} attribute in HTML"
            )

    def test_variant_json_keys_match_js_reads(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        data = _extract_json(body)
        js = _extract_js(body)
        v_fields = set(re.findall(r"v\.([a-z_]+)", js)) - {"locale"}
        v = data[0]
        for field in v_fields:
            assert field in v, f"JS reads v.{field} but _build_variant_data does not emit it"

    def test_annotation_json_keys_match_js_reads(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        data = _extract_json(body)
        js = _extract_js(body)
        a_fields = set(re.findall(r"\ba\.([a-z_A-Z]+)\b", js))
        a_fields -= {"com", "length", "source", "dataset", "append"}
        a_fields -= {"classList", "scrollIntoView", "querySelector"}
        a_fields -= {"preventDefault", "stopPropagation", "textContent"}
        ann = data[0]["annotations"][0]
        for field in a_fields:
            assert field in ann, f"JS reads a.{field} but annotation dict does not contain it"

    def test_enrichment_on_variant_not_annotation(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        data = _extract_json(body)
        v = data[0]
        assert "allele_frequency" in v
        assert "am_pathogenicity" in v
        assert "cadd_phred" in v
        ann = v["annotations"][0]
        assert "allele_frequency" not in ann
        assert "am_pathogenicity" not in ann
        assert "cadd_phred" not in ann

    def test_filter_btn_data_filter_attrs(self, tmp_path: Path):
        body = _rich_report(tmp_path)
        js = _extract_js(body)
        html_filters = set(re.findall(r'data-filter="([^"]+)"', body))
        assert {"all", "bad", "good", "neutral"} <= html_filters
        for f in html_filters:
            assert f'"{f}"' in js

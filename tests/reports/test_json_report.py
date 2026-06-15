# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the JSON report renderer."""

from __future__ import annotations

import json
from pathlib import Path

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.json_report import REGULATORY_NOTICE, SCHEMA_VERSION, render_json


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("genotype.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="MHG_X",
        build="GRCh37",
        total_variants=10,
        skipped_count=0,
        annotators_used=[("clinvar", "20260101"), ("pharmgkb", "2026-05-11")],
        annotations=annotations,
    )


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "x",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "MTHFR",
        "condition": "MTHFR deficiency",
        "references": ["clinvar:allele/100001"],
    }
    defaults.update(overrides)
    return Annotation(**defaults)


class TestRenderJson:
    def test_writes_valid_json(self, tmp_path: Path):
        out = tmp_path / "report.json"
        count = render_json(_result([_ann()]), output_path=out)
        assert count == 1
        payload = json.loads(out.read_text())
        assert payload["schema_version"] == SCHEMA_VERSION
        assert payload["input"]["sample_id"] == "MHG_X"
        assert payload["input"]["total_variants"] == 10
        assert payload["regulatory_notice"] == REGULATORY_NOTICE

    def test_annotators_listed_with_versions(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann()]), output_path=out)
        payload = json.loads(out.read_text())
        assert {"name": "clinvar", "version": "20260101"} in payload["annotators"]
        assert {"name": "pharmgkb", "version": "2026-05-11"} in payload["annotators"]

    def test_annotation_round_trip(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann(rsid="rs99")]), output_path=out)
        payload = json.loads(out.read_text())
        ann = payload["annotations"][0]
        assert ann["rsid"] == "rs99"
        assert ann["attribution"] == "ClinVar"  # ADR-0003: every row attributed
        assert ann["significance"].startswith("clinvar_")

    def test_filters_recorded_in_payload(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(
            _result([_ann()]),
            output_path=out,
            min_magnitude=5.0,
            category="clinical",
            genes={"MTHFR"},
        )
        payload = json.loads(out.read_text())
        assert payload["filters"]["min_magnitude"] == 5.0
        assert payload["filters"]["category"] == "clinical"
        assert payload["filters"]["genes"] == ["MTHFR"]

    def test_min_magnitude_filters_annotations(self, tmp_path: Path):
        anns = [_ann(rsid="lo", magnitude=2.0), _ann(rsid="hi", magnitude=8.0)]
        out = tmp_path / "r.json"
        count = render_json(_result(anns), output_path=out, min_magnitude=5.0)
        assert count == 1
        payload = json.loads(out.read_text())
        assert [a["rsid"] for a in payload["annotations"]] == ["hi"]

    def test_empty_annotations(self, tmp_path: Path):
        out = tmp_path / "r.json"
        count = render_json(_result([]), output_path=out)
        assert count == 0
        payload = json.loads(out.read_text())
        assert payload["annotations"] == []


class TestZygosityField:
    def test_zygosity_in_json_heterozygous(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann(genotype_match="A/G")]), output_path=out)
        payload = json.loads(out.read_text())
        assert payload["annotations"][0]["zygosity"] == "Heterozygous"

    def test_zygosity_in_json_homozygous(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann(genotype_match="A/A")]), output_path=out)
        payload = json.loads(out.read_text())
        assert payload["annotations"][0]["zygosity"] == "Homozygous"

    def test_zygosity_in_json_no_call(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann(genotype_match="A/-")]), output_path=out)
        payload = json.loads(out.read_text())
        assert payload["annotations"][0]["zygosity"] == "No Call"


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

    def test_pharmgkb_attribution_in_json(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("pharmgkb", "2026-01")])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        pharmgkb = [a for a in attrs if a["source"] == "ClinPGx"]
        assert len(pharmgkb) == 1
        assert pharmgkb[0]["license"] == "CC-BY-SA-4.0"

    def test_v3_attribution_shape_pinned(self, tmp_path: Path):
        """Pin the v3 license_attributions shape: source_url + license_url as separate keys."""
        r = self._result_with_annotators([("clinvar", "20260101"), ("pharmgkb", "2026-01")])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        pharmgkb = next(a for a in attrs if a["source"] == "ClinPGx")
        assert pharmgkb["source_url"] == "https://www.clinpgx.org"
        assert pharmgkb["license_url"] == "https://creativecommons.org/licenses/by-sa/4.0/"
        assert pharmgkb["license"] == "CC-BY-SA-4.0"

    def test_clinvar_only_still_has_attribution(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101")])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        assert len(attrs) == 1
        assert attrs[0]["source"] == "ClinVar"

    def test_gwas_attribution_in_json(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("gwas", "2026-01")])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        gwas_attr = [a for a in attrs if a["source"] == "GWAS Catalog"]
        assert len(gwas_attr) == 1
        assert gwas_attr[0]["license"] == "custom-embl-ebi"

    def test_snpedia_attribution_in_json(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("snpedia", None)])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        snpedia = [a for a in attrs if a["source"] == "SNPedia"]
        assert len(snpedia) == 1
        assert snpedia[0]["license"] == "CC-BY-NC-SA-3.0-US"

    def test_both_attributions_in_json(self, tmp_path: Path):
        r = self._result_with_annotators(
            [("clinvar", "20260101"), ("pharmgkb", "2026-01"), ("snpedia", None)]
        )
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        sources = {a["source"] for a in attrs}
        assert "ClinPGx" in sources
        assert "SNPedia" in sources

    def test_alphamissense_attribution_in_json(self, tmp_path: Path):
        r = self._result_with_annotators([("clinvar", "20260101"), ("alphamissense", "2023.2")])
        out = tmp_path / "r.json"
        render_json(r, output_path=out)
        payload = json.loads(out.read_text())
        attrs = payload["license_attributions"]
        am_attr = [a for a in attrs if a["source"] == "AlphaMissense"]
        assert len(am_attr) == 1
        assert am_attr[0]["license"] == "CC-BY-4.0"
        assert "citation" in am_attr[0]

    def test_am_fields_in_annotation_output(self, tmp_path: Path):
        out = tmp_path / "r.json"
        a = _ann(am_pathogenicity=0.95, am_class="likely_pathogenic")
        render_json(_result([a]), output_path=out)
        payload = json.loads(out.read_text())
        ann = payload["annotations"][0]
        assert ann["am_pathogenicity"] == 0.95
        assert ann["am_class"] == "likely_pathogenic"

    def test_am_fields_default_when_absent(self, tmp_path: Path):
        out = tmp_path / "r.json"
        render_json(_result([_ann()]), output_path=out)
        payload = json.loads(out.read_text())
        ann = payload["annotations"][0]
        assert ann["am_pathogenicity"] is None
        assert ann["am_class"] == ""

    def test_pharmgkb_am_caveat(self, tmp_path: Path):
        out = tmp_path / "r.json"
        a = _ann(source="pharmgkb", am_pathogenicity=0.85, am_class="likely_pathogenic")
        render_json(_result([a]), output_path=out)
        payload = json.loads(out.read_text())
        ann = payload["annotations"][0]
        assert ann["am_caveat"] == "protein structure impact only"

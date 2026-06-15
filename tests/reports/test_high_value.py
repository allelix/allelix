# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for high-value SNP no-call detection."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from allelix.models import Variant
from allelix.reports.high_value import (
    HighValueSNP,
    format_warnings,
    load_high_value_snps,
    scan_no_calls,
)


def _v(rsid: str, *, no_call: bool = False) -> Variant:
    a1 = "-" if no_call else "A"
    a2 = "-" if no_call else "G"
    return Variant(rsid=rsid, chromosome="1", position=1, allele1=a1, allele2=a2)


class TestLoadHighValueSnps:
    def test_loads_builtin_file(self) -> None:
        snps = load_high_value_snps()
        assert len(snps) >= 12
        assert "rs429358" in snps
        assert snps["rs429358"].gene == "APOE"
        assert snps["rs429358"].cluster == "APOE"

    def test_user_file_overrides_builtin(self, tmp_path: Path) -> None:
        custom = tmp_path / "custom.yaml"
        custom.write_text(
            "- rsid: rs429358\n  gene: CUSTOM_GENE\n  cluster: CUSTOM\n  note: overridden\n",
            encoding="utf-8",
        )
        snps = load_high_value_snps(extra_paths=[custom])
        assert snps["rs429358"].gene == "CUSTOM_GENE"
        assert len(snps) >= 12

    def test_malformed_yaml_raises_clear_error(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.yaml"
        bad.write_text("not: valid: yaml: {{{{", encoding="utf-8")
        with pytest.raises(ValueError, match="Failed to read"):
            load_high_value_snps(extra_paths=[bad])

    def test_missing_rsid_raises_clear_error(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.yaml"
        bad.write_text("- gene: BRCA1\n  note: missing rsid\n", encoding="utf-8")
        with pytest.raises(ValueError, match="missing required 'rsid' field"):
            load_high_value_snps(extra_paths=[bad])

    def test_non_list_yaml_raises_clear_error(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.yaml"
        bad.write_text("rsid: rs1\ngene: BRCA1\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Expected a YAML list"):
            load_high_value_snps(extra_paths=[bad])

    def test_user_file_adds_new_entry(self, tmp_path: Path) -> None:
        custom = tmp_path / "custom.yaml"
        custom.write_text(
            "- rsid: rs9999999\n  gene: NEWGENE\n  note: custom snp\n",
            encoding="utf-8",
        )
        snps = load_high_value_snps(extra_paths=[custom])
        assert "rs9999999" in snps
        assert snps["rs9999999"].gene == "NEWGENE"


class TestScanNoCallsDefaultLoad:
    def test_default_loads_builtin_high_value(self) -> None:
        variants = [_v("rs429358", no_call=True)]
        warnings = scan_no_calls(variants)
        assert len(warnings) == 1
        assert warnings[0].snp.rsid == "rs429358"


class TestScanNoCalls:
    @pytest.fixture()
    def high_value(self) -> dict[str, HighValueSNP]:
        return {
            "rs429358": HighValueSNP("rs429358", "APOE", "APOE", "APOE SNP 1"),
            "rs7412": HighValueSNP("rs7412", "APOE", "APOE", "APOE SNP 2"),
            "rs4680": HighValueSNP("rs4680", "COMT", "", "COMT activity"),
        }

    def test_no_call_on_high_value_snp_flagged(self, high_value: dict[str, HighValueSNP]) -> None:
        variants = [_v("rs429358", no_call=True), _v("rs7412")]
        warnings = scan_no_calls(variants, high_value)
        assert len(warnings) == 1
        assert warnings[0].snp.rsid == "rs429358"

    def test_called_snps_not_flagged(self, high_value: dict[str, HighValueSNP]) -> None:
        variants = [_v("rs429358"), _v("rs7412"), _v("rs4680")]
        warnings = scan_no_calls(variants, high_value)
        assert warnings == []

    def test_apoe_cluster_incomplete(self, high_value: dict[str, HighValueSNP]) -> None:
        """rs429358 no-call + rs7412 called → cluster_incomplete=True."""
        variants = [_v("rs429358", no_call=True), _v("rs7412")]
        warnings = scan_no_calls(variants, high_value)
        assert len(warnings) == 1
        assert warnings[0].cluster_incomplete is True

    def test_both_apoe_no_call_not_incomplete(self, high_value: dict[str, HighValueSNP]) -> None:
        """Both APOE SNPs are no-calls — cluster_incomplete is False
        because no member was called.
        """
        variants = [_v("rs429358", no_call=True), _v("rs7412", no_call=True)]
        warnings = scan_no_calls(variants, high_value)
        assert len(warnings) == 2
        for w in warnings:
            assert w.cluster_incomplete is False

    def test_non_cluster_no_call(self, high_value: dict[str, HighValueSNP]) -> None:
        variants = [_v("rs4680", no_call=True)]
        warnings = scan_no_calls(variants, high_value)
        assert len(warnings) == 1
        assert warnings[0].snp.gene == "COMT"
        assert warnings[0].cluster_incomplete is False

    def test_empty_variants(self, high_value: dict[str, HighValueSNP]) -> None:
        warnings = scan_no_calls([], high_value)
        assert warnings == []

    def test_non_high_value_no_call_ignored(self, high_value: dict[str, HighValueSNP]) -> None:
        variants = [_v("rs9999999", no_call=True)]
        warnings = scan_no_calls(variants, high_value)
        assert warnings == []


class TestFormatWarnings:
    def test_simple_warning(self) -> None:
        snp = HighValueSNP("rs4680", "COMT", "", "COMT activity")
        lines = format_warnings(
            [
                __import__("allelix.reports.high_value", fromlist=["NoCallWarning"]).NoCallWarning(
                    snp=snp, cluster_incomplete=False
                )
            ]
        )
        assert len(lines) == 1
        assert "rs4680" in lines[0]
        assert "COMT" in lines[0]

    def test_cluster_incomplete_warning(self) -> None:
        from allelix.reports.high_value import NoCallWarning

        snp = HighValueSNP("rs429358", "APOE", "APOE", "APOE SNP 1")
        lines = format_warnings([NoCallWarning(snp=snp, cluster_incomplete=True)])
        assert "APOE genotype cannot be determined" in lines[0]

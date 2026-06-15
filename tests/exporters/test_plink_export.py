# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for PLINK1 binary format exporter."""

from __future__ import annotations

import shutil

import pytest

from allelix.exporters.plink import _orient_genotype, export_plink
from allelix.models import Variant


def _v(rsid="rs1", chrom="1", pos=100, a1="A", a2="A", build="GRCh37"):
    return Variant(rsid, chrom, pos, a1, a2, build=build)


class TestExportProducesThreeFiles:
    def test_export_produces_three_files(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v()]), prefix, "GRCh37")
        assert (tmp_path / "out.bed").exists()
        assert (tmp_path / "out.bim").exists()
        assert (tmp_path / "out.fam").exists()


class TestBedMagicHeader:
    def test_bed_magic_header(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v()]), prefix, "GRCh37")
        data = (tmp_path / "out.bed").read_bytes()
        assert data[:3] == bytes([0x6C, 0x1B, 0x01])


class TestBimFormat:
    def test_bim_format(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v("rs99", "5", 12345, "G", "G")]), prefix, "GRCh37")
        lines = (tmp_path / "out.bim").read_text().strip().split("\n")
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert len(parts) == 6
        assert parts[0] == "5"
        assert parts[1] == "rs99"
        assert parts[2] == "0"
        assert parts[3] == "12345"
        assert parts[4] == "G"
        assert parts[5] == "0"


class TestBimChromosomeCodes:
    @pytest.mark.parametrize(
        ("chrom", "expected"),
        [("X", "23"), ("Y", "24"), ("MT", "26"), ("1", "1"), ("22", "22")],
    )
    def test_bim_chromosome_codes(self, tmp_path, chrom, expected):
        prefix = tmp_path / "out"
        export_plink(iter([_v(chrom=chrom)]), prefix, "GRCh37")
        line = (tmp_path / "out.bim").read_text().strip()
        assert line.startswith(f"{expected}\t")


class TestFamSingleSample:
    def test_fam_single_sample(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v()]), prefix, "GRCh37")
        lines = (tmp_path / "out.fam").read_text().strip().split("\n")
        assert len(lines) == 1
        assert lines[0] == "0\tSAMPLE\t0\t0\t0\t-9"


class TestHomRefWithRefAlt:
    def test_hom_ref_with_ref_alt(self, tmp_path):
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        export_plink(iter([_v("rs1", a1="A", a2="A")]), prefix, "GRCh37", ref_alt_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x00


class TestHomAltWithRefAlt:
    def test_hom_alt_with_ref_alt(self, tmp_path):
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        export_plink(iter([_v("rs1", a1="G", a2="G")]), prefix, "GRCh37", ref_alt_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x03


class TestHomFallbackNoRefAlt:
    def test_hom_fallback_no_ref_alt(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v("rs1", a1="A", a2="A")]), prefix, "GRCh37")
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "0"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x00


class TestHetWithRefAlt:
    def test_het_with_ref_alt(self, tmp_path):
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        export_plink(iter([_v("rs1", a1="A", a2="G")]), prefix, "GRCh37", ref_alt_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x02


class TestHetFallbackNoRefAlt:
    def test_het_fallback_no_ref_alt(self, tmp_path):
        prefix = tmp_path / "out"
        export_plink(iter([_v("rs1", a1="A", a2="G")]), prefix, "GRCh37")
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x02


class TestHetMinusStrand:
    def test_het_minus_strand(self, tmp_path):
        """User alleles C/T on minus strand; gnomAD says (A, G) forward."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        export_plink(iter([_v("rs1", a1="C", a2="T")]), prefix, "GRCh37", ref_alt_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x02


class TestCrossSubstitutionCollisionFallsBack:
    def test_cross_substitution_collision_falls_back(self, tmp_path):
        """User het A/G, gnomAD (A, C): alleles don't fit, must fall back."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "C")}
        written, _, _, _ = export_plink(
            iter([_v("rs1", a1="A", a2="G")]), prefix, "GRCh37", ref_alt_map
        )
        assert written == 1
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "G"
        assert parts[5] != "C"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x02


class TestNoCallsSkipped:
    def test_no_calls_skipped(self, tmp_path):
        prefix = tmp_path / "out"
        variants = [_v("rs1", a1="A", a2="A"), _v("rs2", a1="-", a2="A")]
        written, skipped, _, _ = export_plink(iter(variants), prefix, "GRCh37")
        assert written == 1
        assert skipped == 1
        bim_lines = (tmp_path / "out.bim").read_text().strip().split("\n")
        assert len(bim_lines) == 1
        assert "rs2" not in bim_lines[0]


class TestBimBedAlignmentWithNocalls:
    def test_bim_bed_alignment_with_nocalls(self, tmp_path):
        prefix = tmp_path / "out"
        variants = [
            _v("rs1", a1="A", a2="G"),
            _v("rs2", a1="-", a2="-"),
            _v("rs3", a1="C", a2="C"),
            _v("rs4", a1="A", a2="-"),
            _v("rs5", a1="T", a2="T"),
        ]
        export_plink(iter(variants), prefix, "GRCh37")
        bim_lines = (tmp_path / "out.bim").read_text().strip().split("\n")
        bed_data = (tmp_path / "out.bed").read_bytes()
        assert len(bim_lines) == len(bed_data) - 3


class TestMonomorphicCount:
    def test_monomorphic_count(self, tmp_path):
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        variants = [
            _v("rs1", a1="A", a2="A"),
            _v("rs2", a1="C", a2="C"),
            _v("rs3", a1="T", a2="T"),
        ]
        written, _, _, mono = export_plink(iter(variants), prefix, "GRCh37", ref_alt_map)
        assert written == 3
        assert mono == 2


class TestRoundtripWithPlink:
    @pytest.mark.integration
    def test_roundtrip_with_plink(self, tmp_path):
        if shutil.which("plink2") is None:
            pytest.skip("plink2 not installed")
        import subprocess

        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G"), "rs2": ("C", "T")}
        variants = [
            _v("rs1", a1="A", a2="G"),
            _v("rs2", a1="C", a2="C"),
        ]
        export_plink(iter(variants), prefix, "GRCh37", ref_alt_map)

        result = subprocess.run(
            ["plink2", "--bfile", str(prefix), "--freq", "--out", str(tmp_path / "check")],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0


class TestOrientGenotypeIndelReturnsNone:
    def test_orient_genotype_indel_returns_none(self):
        """Multi-base alleles (indels) must fall back, not crash or miscode."""
        assert _orient_genotype("A", "AT", "A", "AT") is None
        assert _orient_genotype("AT", "A", "A", "AT") is None
        assert _orient_genotype("AT", "AT", "A", "AT") is None


class TestOrientGenotypePalindromicReturnsNone:
    def test_orient_genotype_palindromic_returns_none(self):
        """Palindromic ref/alt (A/T, C/G) can't determine strand — must return None."""
        assert _orient_genotype("A", "A", "A", "T") is None
        assert _orient_genotype("A", "T", "A", "T") is None
        assert _orient_genotype("C", "G", "C", "G") is None
        assert _orient_genotype("G", "G", "C", "G") is None


class TestPalindromicRefAltFallsBack:
    def test_palindromic_ref_alt_falls_back(self, tmp_path):
        """Hom variant with palindromic ref/alt in map should use fallback coding."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "T")}
        written, _, _, mono = export_plink(
            iter([_v("rs1", a1="A", a2="A")]), prefix, "GRCh37", ref_alt_map
        )
        assert written == 1
        assert mono == 1
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "A"
        assert parts[5] == "0"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x00


class TestEmptyVariantList:
    def test_empty_variant_list(self, tmp_path):
        """Empty input produces valid PLINK files with zero variants."""
        prefix = tmp_path / "out"
        written, skipped, indel_skip, mono = export_plink(iter([]), prefix, "GRCh37")
        assert written == 0
        assert skipped == 0
        assert indel_skip == 0
        assert mono == 0
        assert (tmp_path / "out.fam").exists()
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed == bytes([0x6C, 0x1B, 0x01])
        assert (tmp_path / "out.bim").read_text() == ""


class TestMultiAllelicStrandCollision:
    def test_forward_match_preferred_over_complement(self, tmp_path):
        """At ref=G alts=[A,T], user T/T must pick (G,T) not (G,A).

        (G,A) would match via complement (T->A), producing hom-A — wrong.
        The CLI two-pass loop should prefer the forward-matching coord.
        This test verifies the exporter produces the correct output when
        given the right ref_alt_map entry.
        """
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("G", "T")}
        export_plink(iter([_v("rs1", a1="T", a2="T")]), prefix, "GRCh37", ref_alt_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[4] == "G"
        assert parts[5] == "T"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x03

    def test_complement_collision_would_miscode(self, tmp_path):
        """Prove that picking the wrong coord (G,A) miscodes the genotype.

        User carries T/T forward. If ref_alt_map wrongly says (G,A),
        _orient_genotype maps via complement: T->A, yielding hom-alt for
        the A allele — the user gets exported as A/A instead of T/T.
        """
        prefix = tmp_path / "out"
        wrong_map = {"rs1": ("G", "A")}
        export_plink(iter([_v("rs1", a1="T", a2="T")]), prefix, "GRCh37", wrong_map)
        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[5] == "A"
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x03


class TestIndelsSkipped:
    def test_indels_skipped(self, tmp_path):
        """Multi-character alleles (indels) are skipped, not written to BIM."""
        prefix = tmp_path / "out"
        variants = [
            _v("rs1", a1="A", a2="A"),
            _v("rs2", a1="A", a2="AT"),
            _v("rs3", a1="GAC", a2="G"),
        ]
        written, skipped, indel_skip, _ = export_plink(iter(variants), prefix, "GRCh37")
        assert written == 1
        assert skipped == 0
        assert indel_skip == 2
        bim_lines = (tmp_path / "out.bim").read_text().strip().split("\n")
        assert len(bim_lines) == 1
        assert "rs1" in bim_lines[0]


class TestExporterPreservesInputOrder:
    def test_exporter_preserves_input_order(self, tmp_path):
        """Exporter writes variants in iterator order, does not sort."""
        prefix = tmp_path / "out"
        variants = [
            _v("rs1", "2", 200, "G", "G"),
            _v("rs2", "1", 100, "A", "A"),
        ]
        export_plink(iter(variants), prefix, "GRCh37")
        lines = (tmp_path / "out.bim").read_text().strip().split("\n")
        assert lines[0].split("\t")[1] == "rs1"
        assert lines[1].split("\t")[1] == "rs2"


class TestAllNoCalls:
    def test_all_no_calls(self, tmp_path):
        """All no-call input produces magic-only .bed and empty .bim."""
        prefix = tmp_path / "out"
        variants = [_v("rs1", a1="-", a2="-"), _v("rs2", a1="A", a2="-")]
        written, skipped, _, _ = export_plink(iter(variants), prefix, "GRCh37")
        assert written == 0
        assert skipped == 2
        bed = (tmp_path / "out.bed").read_bytes()
        assert len(bed) == 3

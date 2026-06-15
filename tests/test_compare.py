# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for genotype file comparison with strand awareness."""

from __future__ import annotations

from allelix.compare import compare_variants
from allelix.models import Variant


def _v(rsid: str, a1: str, a2: str, *, chrom: str = "1", build: str = "GRCh38") -> Variant:
    return Variant(rsid=rsid, chromosome=chrom, position=1, allele1=a1, allele2=a2, build=build)


class TestConcordance:
    def test_identical_genotypes_concordant(self) -> None:
        v1 = [_v("rs1", "A", "G")]
        v2 = [_v("rs1", "A", "G")]
        r = compare_variants(v1, v2)
        assert r.concordant == 1
        assert r.discordant == 0

    def test_allele_order_independent(self) -> None:
        """G/A and A/G are the same genotype."""
        v1 = [_v("rs1", "G", "A")]
        v2 = [_v("rs1", "A", "G")]
        r = compare_variants(v1, v2)
        assert r.concordant == 1

    def test_discordant_genotypes(self) -> None:
        v1 = [_v("rs1", "A", "A")]
        v2 = [_v("rs1", "G", "G")]
        r = compare_variants(v1, v2)
        assert r.discordant == 1
        assert r.concordant == 0


class TestStrandFlip:
    def test_het_complement_match(self) -> None:
        """A/G matches T/C (complement) — 4 distinct alleles, not palindromic."""
        v1 = [_v("rs1", "A", "G")]
        v2 = [_v("rs1", "T", "C")]
        r = compare_variants(v1, v2)
        assert r.strand_flip_match == 1
        assert r.discordant == 0

    def test_het_complement_reversed(self) -> None:
        """C/A matches G/T (complement, allele order independent)."""
        v1 = [_v("rs1", "C", "A")]
        v2 = [_v("rs1", "G", "T")]
        r = compare_variants(v1, v2)
        assert r.strand_flip_match == 1


class TestStrandAmbiguous:
    def test_at_het_reordered_is_concordant(self) -> None:
        """A/T vs T/A — same genotype (order independent), concordant."""
        v1 = [_v("rs1", "A", "T")]
        v2 = [_v("rs1", "T", "A")]
        r = compare_variants(v1, v2)
        assert r.concordant == 1
        assert r.strand_ambiguous == 0

    def test_at_homozygous_disagreement_flagged(self) -> None:
        """A/A vs T/T on a palindromic SNP — ambiguous, not discordant."""
        v1 = [_v("rs1", "A", "A")]
        v2 = [_v("rs1", "T", "T")]
        r = compare_variants(v1, v2)
        assert r.strand_ambiguous == 1
        assert r.discordant == 0

    def test_cg_het_concordant(self) -> None:
        """C/G vs G/C — same genotype (order independent), concordant."""
        v1 = [_v("rs1", "C", "G")]
        v2 = [_v("rs1", "G", "C")]
        r = compare_variants(v1, v2)
        assert r.concordant == 1

    def test_cg_homozygous_ambiguous(self) -> None:
        """C/C vs G/G on a CG SNP — strand ambiguous."""
        v1 = [_v("rs1", "C", "C")]
        v2 = [_v("rs1", "G", "G")]
        r = compare_variants(v1, v2)
        assert r.strand_ambiguous == 1
        assert r.discordant == 0

    def test_hom_vs_het_palindromic_ambiguous(self) -> None:
        """A/A vs T/A — palindromic pair across files, strand ambiguous."""
        v1 = [_v("rs1", "A", "A")]
        v2 = [_v("rs1", "T", "A")]
        r = compare_variants(v1, v2)
        assert r.strand_ambiguous == 1
        assert r.discordant == 0


class TestNoCalls:
    def test_no_call_in_file1(self) -> None:
        v1 = [_v("rs1", "-", "A")]
        v2 = [_v("rs1", "A", "G")]
        r = compare_variants(v1, v2)
        assert r.no_call == 1
        assert r.concordant == 0

    def test_no_call_in_both(self) -> None:
        v1 = [_v("rs1", "-", "-")]
        v2 = [_v("rs1", "-", "A")]
        r = compare_variants(v1, v2)
        assert r.no_call == 1


class TestCoverage:
    def test_file1_only_and_file2_only(self) -> None:
        v1 = [_v("rs1", "A", "A"), _v("rs2", "G", "G")]
        v2 = [_v("rs1", "A", "A"), _v("rs3", "C", "C")]
        r = compare_variants(v1, v2)
        assert r.shared == 1
        assert r.file1_only == 1
        assert r.file2_only == 1

    def test_totals(self) -> None:
        v1 = [_v("rs1", "A", "G"), _v("rs2", "C", "T"), _v("rs3", "A", "A")]
        v2 = [_v("rs1", "A", "G"), _v("rs4", "T", "T")]
        r = compare_variants(v1, v2)
        assert r.file1_total == 3
        assert r.file2_total == 2


class TestEmptyInputs:
    def test_both_empty(self) -> None:
        r = compare_variants([], [])
        assert r.file1_total == 0
        assert r.file2_total == 0
        assert r.shared == 0
        assert r.build1 == "unknown"

    def test_one_empty(self) -> None:
        v1 = [_v("rs1", "A", "G")]
        r = compare_variants(v1, [])
        assert r.file1_total == 1
        assert r.file2_total == 0
        assert r.shared == 0
        assert r.file1_only == 1


class TestBuildWarning:
    def test_different_builds_reported(self) -> None:
        v1 = [_v("rs1", "A", "G", build="GRCh37")]
        v2 = [_v("rs1", "A", "G", build="GRCh38")]
        r = compare_variants(v1, v2)
        assert r.build1 == "GRCh37"
        assert r.build2 == "GRCh38"


class TestChromosomeCounts:
    def test_per_chromosome_breakdown(self) -> None:
        v1 = [_v("rs1", "A", "G", chrom="1"), _v("rs2", "A", "C", chrom="2")]
        v2 = [_v("rs1", "A", "G", chrom="1"), _v("rs2", "T", "G", chrom="2")]
        r = compare_variants(v1, v2)
        assert r.chromosome_counts["1"]["concordant"] == 1
        assert r.chromosome_counts["2"]["strand_flip_match"] == 1

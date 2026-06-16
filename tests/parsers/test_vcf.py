# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the VCF / gVCF parser."""

from __future__ import annotations

import gzip
from typing import TYPE_CHECKING

import pytest

from allelix.models import NO_CALL_MARKER
from allelix.parsers.vcf import (
    MultiSampleError,
    SampleNotFoundError,
    VcfParser,
    format_sample_list,
)

if TYPE_CHECKING:
    from pathlib import Path

FIXTURES = "tests/fixtures"
MOCK_VCF = f"{FIXTURES}/mock_vcf.vcf"
MOCK_GVCF = f"{FIXTURES}/mock_gvcf.g.vcf"
MOCK_MULTISAMPLE = f"{FIXTURES}/mock_multisample.vcf"


def _path(rel: str) -> Path:
    from pathlib import Path

    return Path(rel)


# ── Auto-detection ──────────────────────────────────────────────


class TestCanParse:
    def test_recognizes_plain_vcf(self) -> None:
        assert VcfParser().can_parse(_path(MOCK_VCF))

    def test_recognizes_gvcf(self) -> None:
        assert VcfParser().can_parse(_path(MOCK_GVCF))

    def test_recognizes_multisample_vcf(self) -> None:
        assert VcfParser().can_parse(_path(MOCK_MULTISAMPLE))

    def test_rejects_non_vcf(self, tmp_path: Path) -> None:
        f = tmp_path / "not_a_vcf.txt"
        f.write_text("# Just some text\nrs123\t1\t100\tA\tG\n")
        assert not VcfParser().can_parse(f)

    def test_rejects_empty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.vcf"
        f.write_text("")
        assert not VcfParser().can_parse(f)

    def test_recognizes_gzipped_vcf(self, tmp_path: Path) -> None:
        gz = tmp_path / "small.vcf.gz"
        with gzip.open(gz, "wt") as h:
            h.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        assert VcfParser().can_parse(gz)


# ── vcf_type detection ──────────────────────────────────────────


class TestVcfType:
    def test_plain_vcf_detected(self) -> None:
        assert VcfParser().vcf_type(_path(MOCK_VCF)) == "plain"

    def test_gvcf_detected_via_non_ref_alt_header(self) -> None:
        """##ALT=<ID=NON_REF,...> is the strong gVCF marker."""
        assert VcfParser().vcf_type(_path(MOCK_GVCF)) == "gvcf"

    def test_gvcf_detected_via_end_tag_in_data(self, tmp_path: Path) -> None:
        """END= in INFO is the fallback gVCF signal when header lacks NON_REF."""
        f = tmp_path / "endtag.g.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "1\t10\t.\tA\tG\t.\t.\tEND=20\tGT\t0/0\n"
        )
        assert VcfParser().vcf_type(f) == "gvcf"


# ── list_samples ────────────────────────────────────────────────


class TestListSamples:
    def test_single_sample(self) -> None:
        assert VcfParser().list_samples(_path(MOCK_VCF)) == ["NA12878"]

    def test_multi_sample(self) -> None:
        assert VcfParser().list_samples(_path(MOCK_MULTISAMPLE)) == [
            "SAMPLE_A",
            "SAMPLE_B",
            "SAMPLE_C",
        ]


# ── get_metadata ────────────────────────────────────────────────


class TestGetMetadata:
    def test_extracts_build_from_contig_assembly(self) -> None:
        meta = VcfParser().get_metadata(_path(MOCK_VCF))
        assert meta["build"] == "GRCh38"
        assert meta["format"] == "vcf"
        assert meta["sample_id"] == "NA12878"

    def test_multisample_unbound_returns_empty_sample_id(self) -> None:
        """Multi-sample without a --sample binding leaves sample_id empty.

        Doesn't raise — get_metadata is cheap and side-effect-free.
        The CLI raises on parse() instead, when actually consuming data.
        """
        meta = VcfParser().get_metadata(_path(MOCK_MULTISAMPLE))
        assert meta["sample_id"] == ""

    def test_multisample_bound_returns_chosen_sample(self) -> None:
        meta = VcfParser(sample="SAMPLE_B").get_metadata(_path(MOCK_MULTISAMPLE))
        assert meta["sample_id"] == "SAMPLE_B"

    def test_build_assembly_strips_hg19_alias(self, tmp_path: Path) -> None:
        f = tmp_path / "hg19.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=hg19>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        )
        assert VcfParser().get_metadata(f)["build"] == "GRCh37"

    def test_chr_prefix_observed_chr1(self, tmp_path: Path) -> None:
        """GH #38: chr1 contig signals GRCh38 convention."""
        f = tmp_path / "chr1.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=248956422>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        )
        assert VcfParser().get_metadata(f)["chr_prefix_observed"] is True

    def test_chr_prefix_observed_chr_only_after_chr1(self, tmp_path: Path) -> None:
        """GH #38 widening: per-chromosome VCFs that omit chr1 still get
        detected (the previous regex only matched chr1 / chrX)."""
        for cid in ("chr2", "chr22", "chr10", "chrY", "chrM", "chrMT"):
            f = tmp_path / f"{cid}.vcf"
            f.write_text(
                f"##fileformat=VCFv4.2\n"
                f"##contig=<ID={cid},length=200000000>\n"
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            )
            assert VcfParser().get_metadata(f)["chr_prefix_observed"] is True, (
                f"chr-prefix not detected for {cid}"
            )

    def test_chr_prefix_not_observed_for_alt_contigs(self, tmp_path: Path) -> None:
        """GH #38: alt contigs and decoys must NOT flip the signal — they
        don't disambiguate the build the same way as standard
        chromosomes do."""
        f = tmp_path / "decoys.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=hs37d5,length=35477943>\n"
            "##contig=<ID=GL000207.1,length=4262>\n"
            "##contig=<ID=NC_007605,length=171823>\n"
            "##contig=<ID=chr1_KI270706v1_random,length=175055>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        )
        # chr1_KI270706v1_random has an `ID=chr1_...` prefix but is NOT
        # a standard chromosome — the terminator `[,>]` in the regex
        # prevents the match.
        assert VcfParser().get_metadata(f)["chr_prefix_observed"] is False

    def test_chr_prefix_not_observed_bare_chromosomes(self, tmp_path: Path) -> None:
        """GH #38: bare numeric / X / MT chromosomes do NOT signal chr-prefix."""
        f = tmp_path / "bare.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621>\n"
            "##contig=<ID=22,length=51304566>\n"
            "##contig=<ID=X,length=155270560>\n"
            "##contig=<ID=MT,length=16569>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        )
        assert VcfParser().get_metadata(f)["chr_prefix_observed"] is False


# ── parse: single-sample plain VCF ─────────────────────────────


class TestParsePlainVcf:
    def test_yields_carrier_variants(self) -> None:
        variants = list(VcfParser().parse(_path(MOCK_VCF)))
        rsids = [v.rsid for v in variants]
        # Carriers: rs1801133 (0/1), rs1801131 (1/1), rs4680 (1/1),
        # rs7412 (0/1), rs900000001 (1/2), rs900000002 (0/1),
        # rs900000003 (1/1), rs900000004 (1 haploid), rs900000006 (0/1 indel),
        # rs900000007 (0/1)
        assert "rs1801133" in rsids
        assert "rs4680" in rsids
        assert "rs1801131" in rsids
        assert "rs7412" in rsids

    def test_skips_homozygous_reference(self) -> None:
        """User carries no ALT → nothing for annotators to bind to."""
        variants = list(VcfParser().parse(_path(MOCK_VCF)))
        # rs429358 (T) is 0/0 in the fixture
        assert "rs429358" not in [v.rsid for v in variants]
        # rs900000005 (A) is 0/0 in the fixture
        assert "rs900000005" not in [v.rsid for v in variants]

    def test_skips_alt_dot(self) -> None:
        """ALT='.' is a reference-only record. Always skipped."""
        variants = list(VcfParser().parse(_path(MOCK_VCF)))
        # The mock has rs900000005 with ALT=. (and GT=0/0 — already
        # hom-ref so doubly skipped). Test passes via absence.
        assert all(v.rsid != "rs900000005" for v in variants)

    def test_normalizes_chromosome(self) -> None:
        """chr-prefixed and bare chromosomes both reach canonical form."""
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        # Bare '1' stays '1'
        assert by_rsid["rs1801133"].chromosome == "1"
        # 'chr1' becomes '1'
        assert by_rsid["rs900000002"].chromosome == "1"
        # 'chrX' becomes 'X'
        assert by_rsid["rs900000003"].chromosome == "X"
        # 'chrM' becomes 'MT'
        assert by_rsid["rs900000004"].chromosome == "MT"

    def test_haploid_call_duplicated(self) -> None:
        """GT='1' (haploid, MT/Y) → diploid (C, C) for downstream uniformity."""
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        # rs900000004 on chrM: REF=T, ALT=C, GT=1 → C/C
        v = by_rsid["rs900000004"]
        assert v.allele1 == "C"
        assert v.allele2 == "C"

    def test_multi_allelic_gt_1_2(self) -> None:
        """GT='1/2' resolves to (first ALT, second ALT)."""
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        # rs900000001: REF=A, ALT=G,C, GT=1/2 → G/C
        v = by_rsid["rs900000001"]
        assert {v.allele1, v.allele2} == {"G", "C"}

    def test_no_call_genotype_yielded_as_no_call(self) -> None:
        """./. GT yields a Variant with NO_CALL_MARKER alleles.

        Matches the array-parser convention: parsers yield no-call
        variants and the annotator pipeline checks ``is_no_call``.
        Keeping the variant in the stream means high-value-SNP no-call
        detection can pick it up downstream (R-2).
        """
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        v = by_rsid["rs121918506"]
        assert v.allele1 == NO_CALL_MARKER
        assert v.allele2 == NO_CALL_MARKER
        assert v.is_no_call

    def test_indel_yields_correctly(self) -> None:
        """REF=ATTC, ALT=A, GT=0/1 → (ATTC, A)."""
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        v = by_rsid["rs900000006"]
        assert {v.allele1, v.allele2} == {"ATTC", "A"}

    def test_empty_rsid_from_dot(self) -> None:
        """ID='.' yields a Variant with empty rsid (uncatalogued site)."""
        variants = list(VcfParser().parse(_path(MOCK_VCF)))
        empty = [v for v in variants if v.rsid == ""]
        # One variant in mock_vcf has ID='.' (the 1\t400000 line)
        assert len(empty) == 1
        assert empty[0].chromosome == "1"
        assert empty[0].position == 400000

    def test_semicolon_separated_ids_picks_rs_prefix(self) -> None:
        """ID='rs900000007;COSV12345' uses the rs-prefixed one."""
        by_rsid = {v.rsid: v for v in VcfParser().parse(_path(MOCK_VCF))}
        assert "rs900000007" in by_rsid
        assert "COSV12345" not in by_rsid

    def test_uppercase_rs_normalized_to_lowercase(self, tmp_path: Path) -> None:
        """ID='RS1801133' normalizes to 'rs1801133'.

        Annotation SQL lookups all assume the lowercase 'rs' convention.
        Returning the original-case ID would silently miss every
        annotation hit on these variants.
        """
        f = tmp_path / "uppercase_rs.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t100\tRS1801133\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
            "1\t200\tRs4680\tA\tG\t100\tPASS\t.\tGT\t0/1\n"
            "1\t300\trs429358;RS7412\tT\tC\t100\tPASS\t.\tGT\t0/1\n"
        )
        variants = list(VcfParser().parse(f))
        rsids = [v.rsid for v in variants]
        assert "rs1801133" in rsids
        assert "RS1801133" not in rsids
        assert "rs4680" in rsids
        # When multiple semicolon-separated IDs exist, the first rs-prefixed
        # one wins — and is lowercased.
        assert "rs429358" in rsids

    def test_non_rs_id_preserves_case(self, tmp_path: Path) -> None:
        """COSMIC IDs (COSV*) preserve case — they have their own conventions."""
        f = tmp_path / "cosmic.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t100\tCOSV12345\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
        )
        variants = list(VcfParser().parse(f))
        assert variants[0].rsid == "COSV12345"


class TestEdgeCases:
    """Tests for the awkward corners: half-no-calls, missing GT, sites-only."""

    def test_half_no_call_yields_asymmetric_variant(self, tmp_path: Path) -> None:
        """GT='0/.' → variant with one called allele and one no-call.

        Matches the array-parser convention: an asymmetric no-call is
        not the same as a hom-ref, not the same as a full no-call. The
        carrier-rule check (ADR-0007) treats it as a no-call because
        is_no_call returns True if either allele is the marker.
        """
        f = tmp_path / "half_no_call.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t100\trs1\tG\tA\t100\tPASS\t.\tGT\t0/.\n"
            "1\t200\trs2\tG\tA\t100\tPASS\t.\tGT\t./0\n"
        )
        variants = list(VcfParser().parse(f))
        assert len(variants) == 2
        for v in variants:
            assert v.is_no_call
            # One allele called as REF, the other no-call
            alleles = {v.allele1, v.allele2}
            assert "G" in alleles
            assert NO_CALL_MARKER in alleles

    def test_format_without_gt_subfield_skipped(self, tmp_path: Path) -> None:
        """FORMAT field without GT (e.g., just 'DP') → record skipped.

        The parser can't yield a variant without genotype information,
        but it shouldn't crash either. The line is silently dropped.
        """
        f = tmp_path / "no_gt.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t100\trs1\tG\tA\t100\tPASS\t.\tDP\t30\n"
            "1\t200\trs2\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
        )
        variants = list(VcfParser().parse(f))
        # rs1 (no GT) is dropped; rs2 still yields normally
        rsids = [v.rsid for v in variants]
        assert "rs1" not in rsids
        assert "rs2" in rsids

    def test_sites_only_vcf_yields_nothing(self, tmp_path: Path) -> None:
        """Sites-only VCF (no FORMAT/sample columns) → no variants yielded.

        Allelix annotates on the carrier rule (ADR-0007); without
        genotype information there's nothing to bind. The parser
        doesn't crash on the missing columns.
        """
        f = tmp_path / "sites_only.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\trs1\tG\tA\t100\tPASS\t.\n"
            "1\t200\trs2\tG\tA\t100\tPASS\t.\n"
        )
        variants = list(VcfParser().parse(f))
        assert variants == []

    def test_malformed_line_short_columns_skipped(self, tmp_path: Path) -> None:
        """Lines with fewer than 8 columns are logged and skipped."""
        f = tmp_path / "short.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t100\trs_bad\tG\tA\n"  # 5 columns, malformed
            "1\t200\trs_good\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
        )
        variants = list(VcfParser().parse(f))
        rsids = [v.rsid for v in variants]
        assert "rs_bad" not in rsids
        assert "rs_good" in rsids

    def test_malformed_line_non_numeric_position_skipped(self, tmp_path: Path) -> None:
        """Non-numeric POS column → line skipped."""
        f = tmp_path / "bad_pos.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\tNOT_A_NUMBER\trs_bad\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
            "1\t200\trs_good\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
        )
        variants = list(VcfParser().parse(f))
        rsids = [v.rsid for v in variants]
        assert "rs_bad" not in rsids
        assert "rs_good" in rsids


class TestPhasedGenotypes:
    """GT delimiter '|' (phased) vs '/' (unphased) — semantically equivalent
    for Allelix's carrier-rule annotation. Phasing tells you which
    chromosome each allele lives on, which the annotation pipeline doesn't
    use. The parser must accept both delimiters.
    """

    def test_phased_heterozygote(self, tmp_path: Path) -> None:
        """GT='0|1' resolves the same as '0/1'."""
        f = tmp_path / "phased.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t1000\trs1801133\tG\tA\t100\tPASS\t.\tGT\t0|1\n"
        )
        variants = list(VcfParser().parse(f))
        assert len(variants) == 1
        v = variants[0]
        assert v.rsid == "rs1801133"
        assert {v.allele1, v.allele2} == {"G", "A"}

    def test_phased_homozygous_alt(self, tmp_path: Path) -> None:
        """GT='1|1' resolves the same as '1/1'."""
        f = tmp_path / "phased_hom.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t1000\trs1801133\tG\tA\t100\tPASS\t.\tGT\t1|1\n"
        )
        variants = list(VcfParser().parse(f))
        assert len(variants) == 1
        v = variants[0]
        assert v.allele1 == "A"
        assert v.allele2 == "A"

    def test_phased_multi_allelic(self, tmp_path: Path) -> None:
        """GT='1|2' on a multi-allelic site resolves to (first ALT, second ALT)."""
        f = tmp_path / "phased_multi.vcf"
        f.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
            "1\t1000\trs900000001\tA\tG,C\t100\tPASS\t.\tGT\t1|2\n"
        )
        variants = list(VcfParser().parse(f))
        assert len(variants) == 1
        v = variants[0]
        assert {v.allele1, v.allele2} == {"G", "C"}

    def test_phased_and_unphased_yield_identical_variants(self, tmp_path: Path) -> None:
        """Phasing is metadata; the carrier-rule-relevant alleles are identical."""
        f_phased = tmp_path / "phased.vcf"
        f_unphased = tmp_path / "unphased.vcf"
        header = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100,assembly=GRCh38>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n"
        )
        f_phased.write_text(header + "1\t1000\trs1801133\tG\tA\t100\tPASS\t.\tGT\t0|1\n")
        f_unphased.write_text(header + "1\t1000\trs1801133\tG\tA\t100\tPASS\t.\tGT\t0/1\n")
        v_phased = list(VcfParser().parse(f_phased))
        v_unphased = list(VcfParser().parse(f_unphased))
        assert v_phased == v_unphased


# ── parse: gVCF reference-block skipping ───────────────────────


class TestParseGvcf:
    def test_reference_blocks_are_skipped(self) -> None:
        """gVCF reference blocks (<NON_REF> + END) are NOT yielded as variants.

        The brief's original test phrasing said "yielded as homozygous
        reference" but the corrected v2.0 behavior is SKIP entirely —
        reference blocks match nothing in any annotation database and
        would burn pipeline cycles for zero hits. The
        tested-vs-untested distinction is reserved for v2.1+.
        """
        variants = list(VcfParser().parse(_path(MOCK_GVCF)))
        # No variant should have ALT containing < or be the placeholder
        # reference block (REF-only).
        for v in variants:
            assert "<" not in v.allele1
            assert "<" not in v.allele2

    def test_gvcf_reference_block_skipped_even_when_spanning_known_rsid(self) -> None:
        """Reference block spanning a position with a known rsid is still skipped.

        The mock has block 1:11796322-11854475 immediately after the
        rs1801133 variant call (which IS yielded). The block itself,
        even though it spans positions that have rsIDs in dbSNP,
        produces no Variant.
        """
        variants = list(VcfParser().parse(_path(MOCK_GVCF)))
        rsids = [v.rsid for v in variants]
        # The variant line for rs1801133 IS yielded:
        assert "rs1801133" in rsids
        # But the reference block lines (which have ID='.') produce no
        # variants.
        empty_rsid_variants = [v for v in variants if not v.rsid]
        assert empty_rsid_variants == []

    def test_gvcf_variant_lines_with_non_ref_alt_yield_correctly(self) -> None:
        """ALT='A,<NON_REF>' on a variant line resolves GT=0/1 to (REF, A).

        The <NON_REF> is the second ALT (index 2 in 1-based GT), so
        GT=0/1 picks REF and the first ALT. <NON_REF> itself should
        never appear in the resulting Variant.
        """
        variants = list(VcfParser().parse(_path(MOCK_GVCF)))
        by_rsid = {v.rsid: v for v in variants}
        # rs1801133: REF=G, ALT=A,<NON_REF>, GT=0/1 → G/A
        v = by_rsid["rs1801133"]
        assert {v.allele1, v.allele2} == {"G", "A"}


# ── parse: multi-sample handling ───────────────────────────────


class TestMultiSample:
    def test_raises_without_sample_binding(self) -> None:
        with pytest.raises(MultiSampleError) as exc_info:
            list(VcfParser().parse(_path(MOCK_MULTISAMPLE)))
        msg = str(exc_info.value)
        assert "SAMPLE_A" in msg
        assert "SAMPLE_B" in msg
        assert "SAMPLE_C" in msg

    def test_raises_on_unknown_sample(self) -> None:
        with pytest.raises(SampleNotFoundError) as exc_info:
            list(VcfParser(sample="NOT_A_SAMPLE").parse(_path(MOCK_MULTISAMPLE)))
        assert "NOT_A_SAMPLE" in str(exc_info.value)

    def test_picks_sample_a_correctly(self) -> None:
        """SAMPLE_A: rs1801133 0/1 (carrier), rs4680 0/0 (skip), rs1801131 0/1."""
        variants = list(VcfParser(sample="SAMPLE_A").parse(_path(MOCK_MULTISAMPLE)))
        rsids = {v.rsid for v in variants}
        assert "rs1801133" in rsids
        assert "rs4680" not in rsids
        assert "rs1801131" in rsids

    def test_picks_sample_b_correctly(self) -> None:
        """SAMPLE_B: rs1801133 1/1, rs4680 0/1, rs1801131 0/0 (skip)."""
        variants = list(VcfParser(sample="SAMPLE_B").parse(_path(MOCK_MULTISAMPLE)))
        rsids = {v.rsid for v in variants}
        assert "rs1801133" in rsids
        assert "rs4680" in rsids
        assert "rs1801131" not in rsids

    def test_picks_sample_c_correctly(self) -> None:
        """SAMPLE_C: rs1801133 0/0 (skip), rs4680 1/1, rs1801131 1/1."""
        variants = list(VcfParser(sample="SAMPLE_C").parse(_path(MOCK_MULTISAMPLE)))
        rsids = {v.rsid for v in variants}
        assert "rs1801133" not in rsids
        assert "rs4680" in rsids
        assert "rs1801131" in rsids


# ── parse: gzipped VCF ─────────────────────────────────────────


class TestGzippedVcf:
    def test_parse_round_trip_through_gz(self, tmp_path: Path) -> None:
        """gzip-compressed VCF yields the same variants as uncompressed."""
        uncompressed = list(VcfParser().parse(_path(MOCK_VCF)))
        gz = tmp_path / "mock_vcf.vcf.gz"
        with open(MOCK_VCF, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        compressed = list(VcfParser().parse(gz))
        assert uncompressed == compressed

    def test_detect_gzip_by_magic_not_extension(self, tmp_path: Path) -> None:
        """File without .gz extension but with gzip magic bytes still parses.

        Sniffs ``1f 8b`` magic, not the extension — robust to misnaming.
        """
        # Note: the parser detects gzip by magic bytes, regardless of
        # the file extension.
        misnamed = tmp_path / "misnamed.vcf"
        with open(MOCK_VCF, "rb") as src:
            content = src.read()
        with gzip.open(misnamed, "wb") as dst:
            dst.write(content)
        variants = list(VcfParser().parse(misnamed))
        assert len(variants) > 0


# ── Registry integration ───────────────────────────────────────


class TestRegistryRegistration:
    def test_vcf_parser_in_registry(self) -> None:
        """VcfParser is registered in parsers/__init__.py:PARSERS."""
        from allelix.parsers import PARSERS, get_parser_by_name

        names = [p.name for p in PARSERS]
        assert "vcf" in names
        # Lookup by name works
        parser = get_parser_by_name("vcf")
        assert isinstance(parser, VcfParser)

    def test_auto_detect_picks_vcf_parser(self) -> None:
        """detect_parser() resolves mock_vcf.vcf to VcfParser."""
        from allelix.parsers import detect_parser

        parser = detect_parser(_path(MOCK_VCF))
        assert isinstance(parser, VcfParser)


# ── NO_CALL_MARKER export check ────────────────────────────────


class TestNoCallEncoding:
    def test_no_call_marker_constant(self) -> None:
        """NO_CALL_MARKER stays at '-' across the codebase."""
        assert NO_CALL_MARKER == "-"


class TestFormatSampleList:
    """1000 Genomes-scale sample lists are truncated for human-readable errors."""

    def test_under_threshold_shows_all(self) -> None:
        assert format_sample_list(["A", "B", "C"]) == "A, B, C"

    def test_at_threshold_shows_all(self) -> None:
        samples = [f"S{i}" for i in range(10)]
        assert format_sample_list(samples, max_shown=10) == ", ".join(samples)

    def test_over_threshold_truncates(self) -> None:
        samples = [f"S{i}" for i in range(15)]
        result = format_sample_list(samples, max_shown=10)
        assert "S0" in result and "S9" in result
        assert "S10" not in result and "S14" not in result
        assert "... and 5 more" in result

    def test_thousand_genomes_scale(self) -> None:
        """3,202 samples render to under 200 chars instead of 30KB."""
        samples = [f"HG{i:05d}" for i in range(3202)]
        result = format_sample_list(samples)
        assert result.startswith("HG00000, HG00001")
        assert "... and 3192 more" in result
        assert len(result) < 200

    def test_empty_list(self) -> None:
        assert format_sample_list([]) == ""

    def test_message_in_sample_error_uses_helper(self) -> None:
        """The error path actually routes through the helper end-to-end."""
        # The committed mock_multisample.vcf has 3 samples — under threshold.
        # Construct a parser with a sample binding that won't match, exercising
        # the SampleNotFoundError path with a small list.
        with pytest.raises(SampleNotFoundError) as exc_info:
            list(VcfParser(sample="NOPE").parse(_path(MOCK_MULTISAMPLE)))
        # All three sample names appear (3 ≤ default max_shown=10)
        for s in ("SAMPLE_A", "SAMPLE_B", "SAMPLE_C"):
            assert s in str(exc_info.value)
        assert "... and" not in str(exc_info.value)


class TestCanonicalRsid:
    """Parser-layer normalization of the VCF ID column."""

    def test_dot_returns_empty(self) -> None:
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid(".") == ""

    def test_empty_returns_empty(self) -> None:
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid("") == ""

    def test_rs_lowercased(self) -> None:
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid("RS1801133") == "rs1801133"
        assert _canonical_rsid("rs1801133") == "rs1801133"

    def test_cosmic_preserved(self) -> None:
        """Non-rs external-DB IDs pass through unchanged."""
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid("COSV12345") == "COSV12345"

    def test_positional_synthetic_returns_empty(self) -> None:
        """1000 Genomes positional IDs (`22:10519265:CA:C`) → empty.

        The variant is then routed through ClinVar position-keyed
        resolution rather than carrying a meaningless string in its
        rsid field. GH #8.
        """
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid("22:10519265:CA:C") == ""
        assert _canonical_rsid("1:11856378:G:A") == ""

    def test_rs_in_semicolon_list_wins(self) -> None:
        """If the ID column has multiple ;-separated values, first rs wins."""
        from allelix.parsers.vcf import _canonical_rsid

        assert _canonical_rsid("COSV1;rs1801133") == "rs1801133"
        assert _canonical_rsid("rs1801133;RCV000123") == "rs1801133"


class TestValidateSample:
    """validate_sample() pre-flight without consuming variant data."""

    def test_single_sample_does_not_raise(self) -> None:
        VcfParser().validate_sample(_path(MOCK_VCF))

    def test_multi_sample_unbound_raises(self) -> None:
        with pytest.raises(MultiSampleError):
            VcfParser().validate_sample(_path(MOCK_MULTISAMPLE))

    def test_multi_sample_with_valid_binding_does_not_raise(self) -> None:
        VcfParser(sample="SAMPLE_A").validate_sample(_path(MOCK_MULTISAMPLE))

    def test_unknown_sample_raises(self) -> None:
        with pytest.raises(SampleNotFoundError):
            VcfParser(sample="NOPE").validate_sample(_path(MOCK_MULTISAMPLE))

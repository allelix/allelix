# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for PLINK1 binary format exporter."""

from __future__ import annotations

import shutil

import pytest

from allelix.exporters.plink import EmptyPlinkSetError, _fits_site, export_plink
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
        """C/T against gnomAD's forward (A, G) is not read as its complement: missing."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "G")}
        counts = export_plink(iter([_v("rs1", a1="C", a2="T")]), prefix, "GRCh37", ref_alt_map)
        assert counts.allele_mismatches == 1
        parts = (tmp_path / "out.bim").read_text().strip().split("\t")
        assert parts[4:6] == ["A", "G"]
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x01


class TestCrossSubstitutionCollision:
    def test_cross_substitution_collision_is_missing(self, tmp_path):
        """User het A/G, gnomAD (A, C): G is neither allele, so the call is missing."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "C")}
        counts = export_plink(iter([_v("rs1", a1="A", a2="G")]), prefix, "GRCh37", ref_alt_map)
        assert counts.written == 1
        assert counts.allele_mismatches == 1
        parts = (tmp_path / "out.bim").read_text().strip().split("\t")
        assert parts[4:6] == ["A", "C"]
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x01


class TestNoCallsSkipped:
    def test_no_calls_skipped(self, tmp_path):
        prefix = tmp_path / "out"
        variants = [_v("rs1", a1="A", a2="A"), _v("rs2", a1="-", a2="A")]
        written, skipped, _, _, _ = export_plink(iter(variants), prefix, "GRCh37")
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
        written, _, _, mono, _ = export_plink(iter(variants), prefix, "GRCh37", ref_alt_map)
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


class TestFitsSite:
    def test_indel_alleles_do_not_fit(self):
        """Multi-base alleles (indels) never fit a site."""
        assert not _fits_site("A", "AT", "A", "AT")
        assert not _fits_site("AT", "A", "A", "AT")
        assert not _fits_site("AT", "AT", "A", "AT")

    def test_alleles_fit_only_as_written(self):
        """No complement reading: T/T does not fit a C>A site."""
        assert _fits_site("C", "A", "C", "A")
        assert _fits_site("A", "A", "C", "A")
        assert not _fits_site("G", "T", "C", "A")
        assert not _fits_site("T", "T", "C", "A")

    def test_palindromic_alleles_fit_as_written(self):
        assert _fits_site("A", "A", "A", "T")
        assert _fits_site("C", "G", "C", "G")
        assert not _fits_site("C", "C", "A", "T")


class TestPalindromicRefAltFallsBack:
    def test_palindromic_ref_alt_falls_back(self, tmp_path):
        """Hom variant with palindromic ref/alt in map should use fallback coding."""
        prefix = tmp_path / "out"
        ref_alt_map = {"rs1": ("A", "T")}
        written, _, _, mono, _ = export_plink(
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
        """Empty input is refused: no PLINK set with zero variants is written."""
        prefix = tmp_path / "out"
        with pytest.raises(EmptyPlinkSetError):
            export_plink(iter([]), prefix, "GRCh37")
        assert list(tmp_path.iterdir()) == []


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

    def test_a_call_that_fits_neither_allele_is_missing(self, tmp_path):
        """T/T at a G>A site is not read as its complement (A/A): it is missing, counted."""
        prefix = tmp_path / "out"
        wrong_map = {"rs1": ("G", "A")}
        counts = export_plink(iter([_v("rs1", a1="T", a2="T")]), prefix, "GRCh37", wrong_map)
        assert counts.written == 1
        assert counts.allele_mismatches == 1
        parts = (tmp_path / "out.bim").read_text().strip().split("\t")
        assert parts[4:6] == ["G", "A"]
        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x01


class TestNoComplementOnExport:
    def test_hom_call_at_a_site_it_does_not_fit_is_missing(self, tmp_path):
        """The issue's example: A/A at a C>T site was written as T/T."""
        prefix = tmp_path / "out"
        counts = export_plink(
            iter([_v("rs1", a1="A", a2="A")]), prefix, "GRCh37", {"rs1": ("C", "T")}
        )
        assert counts.allele_mismatches == 1
        assert counts.monomorphic == 0
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x01

    def test_het_call_at_a_site_it_does_not_fit_is_missing(self, tmp_path):
        prefix = tmp_path / "out"
        counts = export_plink(
            iter([_v("rs1", a1="G", a2="A")]), prefix, "GRCh37", {"rs1": ("C", "T")}
        )
        assert counts.allele_mismatches == 1
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x01

    def test_palindromic_site_with_a_call_that_does_not_fit_is_missing(self, tmp_path):
        prefix = tmp_path / "out"
        counts = export_plink(
            iter([_v("rs1", a1="C", a2="C")]), prefix, "GRCh37", {"rs1": ("A", "T")}
        )
        assert counts.allele_mismatches == 1
        assert counts.monomorphic == 0
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x01

    def test_calls_that_fit_are_not_counted(self, tmp_path):
        prefix = tmp_path / "out"
        counts = export_plink(
            iter([_v("rs1", a1="C", a2="T"), _v("rs2", pos=200, a1="T", a2="T")]),
            prefix,
            "GRCh37",
            {"rs1": ("C", "T"), "rs2": ("C", "T")},
        )
        assert counts.allele_mismatches == 0
        assert (tmp_path / "out.bed").read_bytes()[3:] == bytes([0x02, 0x03])


class TestIndelsSkipped:
    def test_indels_skipped(self, tmp_path):
        """Multi-character alleles (indels) are skipped, not written to BIM."""
        prefix = tmp_path / "out"
        variants = [
            _v("rs1", a1="A", a2="A"),
            _v("rs2", a1="A", a2="AT"),
            _v("rs3", a1="GAC", a2="G"),
        ]
        written, skipped, indel_skip, _, _ = export_plink(iter(variants), prefix, "GRCh37")
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
        """All no-call input is refused, with the skipped count; no files remain."""
        prefix = tmp_path / "out"
        variants = [_v("rs1", a1="-", a2="-"), _v("rs2", a1="A", a2="-")]
        with pytest.raises(EmptyPlinkSetError) as exc:
            export_plink(iter(variants), prefix, "GRCh37")
        assert (exc.value.no_calls, exc.value.indels) == (2, 0)
        assert list(tmp_path.iterdir()) == []


class TestMultiAllelicPrefersAPlainRecord:
    """A call that fits a palindromic and a plain record takes the plain one, in either order."""

    @pytest.mark.parametrize(
        "coords",
        [
            [("1", 100, "A", "T"), ("1", 100, "A", "C")],
            [("1", 100, "A", "C"), ("1", 100, "A", "T")],
        ],
    )
    def test_hom_ref_call_is_coded_against_the_plain_record(self, tmp_path, monkeypatch, coords):
        from allelix.exporters.plink import resolve_ref_alt_via_gnomad

        class _FakeGnomad:
            def __init__(self, *a, **kw):
                pass

            def is_ready(self):
                return True

            def bulk_resolve_coordinates(self, rsids):
                return {"rs1": coords}

            def close(self):
                pass

        monkeypatch.setattr("allelix.annotators.gnomad.GnomadAnnotator", _FakeGnomad)
        call = _v("rs1", a1="A", a2="A")
        ref_alt_map = resolve_ref_alt_via_gnomad({"rs1": call}, tmp_path)
        assert ref_alt_map == {"rs1": ("A", "C")}
        counts = export_plink(iter([call]), tmp_path / "out", "GRCh37", ref_alt_map)
        assert counts.monomorphic == 0
        assert (tmp_path / "out.bim").read_text().split("\t")[4:6] == ["A", "C\n"]
        assert (tmp_path / "out.bed").read_bytes()[3] == 0x00

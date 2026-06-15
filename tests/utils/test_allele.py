# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for strand-flip / complement / ambiguity helpers."""

from __future__ import annotations

from allelix.utils.allele import complement, flip_genotype, is_strand_ambiguous, resolve_strand


class TestComplement:
    def test_single_bases(self):
        assert complement("A") == "T"
        assert complement("T") == "A"
        assert complement("C") == "G"
        assert complement("G") == "C"

    def test_no_call_unchanged(self):
        assert complement("-") == "-"
        assert complement("") == ""

    def test_unknown_letter_unchanged(self):
        assert complement("N") == "N"

    def test_multibase_indel_reverses_and_complements(self):
        # CTT (forward) → AAG (reverse complement)
        assert complement("CTT") == "AAG"
        assert complement("AAG") == "CTT"


class TestFlipGenotype:
    def test_diploid_flip(self):
        assert flip_genotype("C", "T") == ("G", "A")

    def test_no_call_preserved(self):
        assert flip_genotype("-", "A") == ("-", "T")


class TestIsStrandAmbiguous:
    def test_at_pair_is_ambiguous(self):
        assert is_strand_ambiguous("A", "T")
        assert is_strand_ambiguous("T", "A")

    def test_cg_pair_is_ambiguous(self):
        assert is_strand_ambiguous("C", "G")
        assert is_strand_ambiguous("G", "C")

    def test_normal_pair_not_ambiguous(self):
        assert not is_strand_ambiguous("A", "G")
        assert not is_strand_ambiguous("C", "T")

    def test_indel_not_ambiguous(self):
        assert not is_strand_ambiguous("CTT", "C")

    def test_unknown_letter_not_ambiguous(self):
        assert not is_strand_ambiguous("A", "N")


class TestResolveStrand:
    def test_forward_match_ref(self):
        assert resolve_strand("A", "A", "G") == "A"

    def test_forward_match_alt(self):
        assert resolve_strand("G", "A", "G") == "G"

    def test_minus_strand_complement(self):
        assert resolve_strand("T", "A", "G") == "A"

    def test_minus_strand_complement_alt(self):
        assert resolve_strand("C", "A", "G") == "G"

    def test_palindromic_site_direct_match_returns_allele(self):
        assert resolve_strand("T", "A", "T") == "T"
        assert resolve_strand("A", "A", "T") == "A"
        assert resolve_strand("G", "C", "G") == "G"
        assert resolve_strand("C", "C", "G") == "C"

    def test_palindromic_site_no_external_allele_match(self):
        assert resolve_strand("G", "A", "T") is None

    def test_indel_passes_through(self):
        assert resolve_strand("AC", "A", "AC") == "AC"

    def test_no_match_returns_none(self):
        assert resolve_strand("A", "C", "G") is None

    def test_non_acgt_returns_none(self):
        assert resolve_strand("N", "A", "G") is None

    def test_palindromic_complement_match_returns_none(self, monkeypatch):
        """Palindromic guard rejects complement-resolved allele at A/T and C/G sites.

        In production, palindromic sites are always caught by the direct-match
        check (complement(ref)=alt for these pairs). This test exercises the
        defense-in-depth guard on line 72 by injecting a synthetic complement
        mapping that bypasses the direct-match path.
        """
        from allelix.utils import allele

        patched = dict(allele._COMPLEMENT)
        patched["X"] = "A"
        monkeypatch.setattr(allele, "_COMPLEMENT", patched)
        assert resolve_strand("X", "A", "T") is None

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
    """GH #18: complement-resolution is intentionally NOT performed.

    At multi-allelic sites, the complement of the user's true forward
    allele can coincidentally equal a different alt at the same
    position, so a complement fallback can stamp a wrong-allele CADD
    score. The function returns None for any allele that isn't directly
    in ``{ref, alt}``; minus-strand handling is deferred (ADR-0010).
    """

    def test_forward_match_ref(self):
        assert resolve_strand("A", "A", "G") == "A"

    def test_forward_match_alt(self):
        assert resolve_strand("G", "A", "G") == "G"

    def test_complement_no_longer_resolves(self):
        # Was: resolve_strand("T", "A", "G") == "A" — minus-strand fallback.
        # Now returns None; the caller skips enrichment.
        assert resolve_strand("T", "A", "G") is None
        assert resolve_strand("C", "A", "G") is None

    def test_audit_reproduction_no_false_complement(self):
        # GH #18 reproduction: user "A" at a C→T site.
        # Old code returned "T" (complement) and the consumer stamped a
        # CADD score for the C→T transition onto an annotation describing
        # the user's "A" carrier. Now returns None and enrichment is
        # skipped.
        assert resolve_strand("A", "C", "T") is None

    def test_palindromic_direct_match_returns_allele(self):
        assert resolve_strand("T", "A", "T") == "T"
        assert resolve_strand("A", "A", "T") == "A"
        assert resolve_strand("G", "C", "G") == "G"
        assert resolve_strand("C", "C", "G") == "C"

    def test_indel_passes_through(self):
        assert resolve_strand("AC", "A", "AC") == "AC"

    def test_no_match_returns_none(self):
        assert resolve_strand("A", "C", "G") is None

    def test_non_acgt_returns_none(self):
        assert resolve_strand("N", "A", "G") is None

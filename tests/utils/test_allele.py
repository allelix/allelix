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


class TestStrandAwareCarrierMatch:
    """ADR-0035 PR 4: per-row carrier match for ClinVar-style {ref, alt} rows."""

    def test_direct_alt_match(self):
        from allelix.utils.allele import strand_aware_carrier_match

        # User A/T, ClinVar A>T → carrier of T.
        assert strand_aware_carrier_match("A", "A", "T", "A", "T") is True

    def test_direct_no_match_without_variant_ref(self):
        """Conservative: no variant_ref context → direct match only, no flip."""
        from allelix.utils.allele import strand_aware_carrier_match

        # User T/T, ClinVar A>C → no carrier; without ref context, no flip.
        assert strand_aware_carrier_match(None, "T", "T", "A", "C") is False

    def test_strand_flip_match_when_variant_ref_is_complement(self):
        """User reverse-strand T/T, forward ref is A (complement of T); ClinVar A>G.
        complement(G)=C; C not in {T,T} so no flip-match for this row.
        Now ClinVar A>T: complement(T)=A; A not in {T,T}. Hmm rs5742904 case is C/T.
        Use the canonical rs5742904: forward C/T (het, carrier of T) vs coding G/A.
        Coding-side variant.ref = G (complement of forward C). User reads A,G.
        """
        from allelix.utils.allele import strand_aware_carrier_match

        # rs5742904 forward C/T: ref=C, alt=T. variant_ref=C; user {C,T}; direct.
        assert strand_aware_carrier_match("C", "C", "T", "C", "T") is True
        # rs5742904 coding G/A: variant_ref=G (=complement of C); user {G,A}.
        # complement(T)=A; A in {G,A} → strand-flip carrier match. CORRECT.
        assert strand_aware_carrier_match("G", "G", "A", "C", "T") is True

    def test_palindromic_site_skips_strand_flip(self):
        """A/T site: complement(T)=A is the other ref/alt — strand-flip ambiguous."""
        from allelix.utils.allele import strand_aware_carrier_match

        # User A/A (forward), ClinVar A>T at a palindromic site: direct ref only,
        # so user isn't a carrier of T. Direct check fails; flip skipped.
        assert strand_aware_carrier_match("A", "A", "A", "A", "T") is False
        # User T/T reverse (variant.ref=T, forward ref=A). Direct: T in {T,T}; True.
        # But the test point is palindrome: we should NOT use flip to pull from
        # T/T → complement(T)=A → match A. Direct already fires here so the
        # function returns True. Use a case where direct fails to test the flip
        # guard: user "A/A" hom-ref forward but variant_ref=T (claimed reverse).
        # complement(T)=A; A != T (source_ref) → no flip path applies cleanly.
        assert strand_aware_carrier_match("T", "A", "A", "A", "T") is False

    def test_multi_allelic_safety_variant_ref_disagrees_no_flip(self):
        """variant_ref doesn't match source_ref or its complement → no flip fires."""
        from allelix.utils.allele import strand_aware_carrier_match

        # User T/T, variant_ref=C (disagrees with A and complement(A)=T).
        # The audit-reproduction shape from #18: complement(user)=A coincides
        # with source_ref but the user isn't reverse-stranded — abstain.
        assert strand_aware_carrier_match("C", "T", "T", "A", "G") is False

    def test_forward_orientation_blocks_complement_path(self):
        """variant_ref == source_ref → user is forward; direct must be the only try."""
        from allelix.utils.allele import strand_aware_carrier_match

        # User T/T forward at C>T site (correct carrier path).
        assert strand_aware_carrier_match("C", "T", "T", "C", "T") is True
        # User G/G forward at C>T (not a carrier, no flip allowed since ref matches).
        assert strand_aware_carrier_match("C", "G", "G", "C", "T") is False


class TestStrandAwareGenotypeMatch:
    """ADR-0035 PR 4: per-diploid carrier match for ClinPGx / SNPedia."""

    def test_direct_match(self):
        from allelix.utils.allele import strand_aware_genotype_match

        # User A/G, source "AG" → match.
        assert strand_aware_genotype_match("A", "G", "AG", "A") is True
        # Order independence (helper sorts).
        assert strand_aware_genotype_match("G", "A", "AG", "A") is True

    def test_strand_flip_match(self):
        from allelix.utils.allele import strand_aware_genotype_match

        # User C/T reading on reverse strand; forward source ref is A, so
        # variant.ref (the user's claimed REF) is complement(A) = T.
        # User's complement diploid = {complement(C), complement(T)} = {G, A}
        # sorted = "AG" — matches source. variant.ref ('T') is NOT in source
        # 'AG' AND complement('T') = 'A' IS in source → orientation confirmed
        # reverse, palindrome guard clears → strand-flip carrier.
        assert strand_aware_genotype_match("C", "T", "AG", "T") is True

    def test_palindromic_genotype_skipped(self):
        from allelix.utils.allele import strand_aware_genotype_match

        # Source "AT" is palindromic; strand-flip never fires.
        # User T/A on reverse would map to {A,T} forward — also "AT".
        # The function must NOT claim a match via flip path (ambiguous).
        # Direct match: A/T → "AT" matches source; direct path fires (True).
        assert strand_aware_genotype_match("A", "T", "AT", "A") is True
        # Reverse path: user T/A with variant_ref=T (complement of forward A).
        # complement({T,A}) = {A,T} = "AT"; would otherwise match, but
        # palindrome guard blocks the flip path. Direct gives "AT" too, so True.
        # Use a case where direct fails: user G/C, source "AT", variant_ref=T.
        # Direct: GC != AT. Complement: {C,G}→{C,G}="CG"!="AT". No match.
        assert strand_aware_genotype_match("G", "C", "AT", "T") is False

    def test_no_variant_ref_skips_flip(self):
        from allelix.utils.allele import strand_aware_genotype_match

        # No reference context → direct match only.
        assert strand_aware_genotype_match("C", "T", "AG", None) is False

    def test_variant_ref_in_source_blocks_flip(self):
        """variant_ref ∈ source_geno → forward orientation; flip is wrong."""
        from allelix.utils.allele import strand_aware_genotype_match

        # User C/T forward at site source "AG". variant_ref=A is in "AG".
        # Direct: CT != AG. Complement: {G,A} = "AG" matches, but variant_ref=A
        # is in source — user is forward, complement match is coincidence.
        assert strand_aware_genotype_match("C", "T", "AG", "A") is False


class TestDeriveAltFromDiploid:
    """ADR-0035 PR 2: alt derivation for SNPedia / ClinPGx matched genotypes."""

    def test_none_ref_returns_empty(self):
        """Array data prior to PR 4's ref population has no reference context."""
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid(None, "A", "G") == ""

    def test_heterozygous_picks_non_ref(self):
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid("A", "A", "G") == "G"
        assert derive_alt_from_diploid("A", "G", "A") == "G"

    def test_homozygous_alt_returns_alt(self):
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid("A", "G", "G") == "G"

    def test_homozygous_ref_returns_empty(self):
        """Defensive: hom-ref should be filtered earlier (ADR-0023) but if it
        reaches the helper, return "" rather than guess.
        """
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid("A", "A", "A") == ""

    def test_multi_allelic_het_neither_equals_ref_returns_empty(self):
        """Rare: matched row's alleles disagree with REF on both sides.
        Conservative: return "" rather than guess the user's actual carried alt.
        """
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid("A", "C", "T") == ""

    def test_indel_ref_heterozygous(self):
        from allelix.utils.allele import derive_alt_from_diploid

        assert derive_alt_from_diploid("ATG", "ATG", "A") == "A"
        assert derive_alt_from_diploid("ATG", "A", "ATG") == "A"

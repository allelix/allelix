# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Strand flipping, complement logic, and ambiguous-SNP detection.

A SNP read on the reverse strand has its alleles complemented (A↔T, C↔G).
Two databases reporting the "same" variant on opposite strands will list
opposite allele letters. For most SNPs this is unambiguous and reversible.
For A/T and C/G SNPs (palindromic), the complement equals the alternative —
so a strand-flip is undetectable from sequence alone and is best handled by
extra information (allele frequency, surrounding context).

ADR-0010 documents the design.
"""

from __future__ import annotations

from allelix.models import NO_CALL_MARKER

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# A/T and C/G SNPs are palindromic; their complement equals the alternative,
# so strand orientation cannot be inferred from the alleles alone.
_AMBIGUOUS_PAIRS: frozenset[frozenset[str]] = frozenset(
    {frozenset({"A", "T"}), frozenset({"C", "G"})}
)


def complement(allele: str) -> str:
    """Return the reverse-complement of a single allele string.

    A → T, T → A, C → G, G → C. The no-call marker `-` and any unrecognized
    character are returned unchanged. Handles indels (multi-base alleles) by
    complementing each base in reverse order.
    """
    if allele == NO_CALL_MARKER or not allele:
        return allele
    if len(allele) == 1:
        return _COMPLEMENT.get(allele, allele)
    return "".join(_COMPLEMENT.get(b, b) for b in reversed(allele))


def flip_genotype(allele1: str, allele2: str) -> tuple[str, str]:
    """Return both alleles complemented (the reverse-strand reading)."""
    return complement(allele1), complement(allele2)


def resolve_strand(user_allele: str, gnomad_ref: str, gnomad_alt: str) -> str | None:
    """Return reference-forward allele, or None if not directly present.

    Maps an array-reported allele to its reference-forward equivalent
    using gnomAD's ref/alt as the ground truth. If the user allele
    matches ref or alt directly, it's already forward and is returned.

    GH #18: previously, when the user allele did NOT directly match
    ref or alt, this function fell back to a complement check
    (``complement(user_allele) ∈ {ref, alt}`` → return the complement)
    on the assumption that the array was read on the minus strand. At
    multi-allelic sites that fallback is unsound — the complement of
    the user's true forward allele can coincidentally equal a different
    alt at the same position, returning an allele the user does not
    carry. The downstream CADD enrichment path then stamps a wrong-
    allele score onto the annotation. The complement fallback is
    removed: if the user allele is not in ``{ref, alt}``, return None
    and let the caller skip enrichment rather than risk a coincidental
    match. (ADR-0010 already defers proper strand handling; this makes
    the deferral explicit at the resolution layer.)

    Only operates on single-base alleles. Multi-base alleles (indels)
    pass through as-is — array indels are rare and not minus-strand
    reported.
    """
    if len(user_allele) != 1:
        return user_allele
    if user_allele in (gnomad_ref, gnomad_alt):
        return user_allele
    return None


def is_strand_ambiguous(ref: str, alt: str) -> bool:
    """True if (ref, alt) is an A/T or C/G pair — strand cannot be inferred.

    Multi-base indels and any allele containing a no-call or unknown letter
    are reported as not ambiguous (they have other ways to disambiguate).
    """
    if len(ref) != 1 or len(alt) != 1:
        return False
    if ref not in _COMPLEMENT or alt not in _COMPLEMENT:
        return False
    return frozenset({ref, alt}) in _AMBIGUOUS_PAIRS

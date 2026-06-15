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


_PALINDROMIC = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


def resolve_strand(user_allele: str, gnomad_ref: str, gnomad_alt: str) -> str | None:
    """Return reference-forward allele, or None if ambiguous.

    Maps an array-reported allele to its reference-forward equivalent
    using gnomAD's ref/alt as the ground truth. If the user allele
    matches ref or alt directly, it's already forward. If the
    complement matches, the array was minus-strand. Palindromic SNPs
    (A/T, C/G ref/alt pairs) cannot be resolved and return None.

    Only operates on single-base alleles. Multi-base alleles (indels)
    pass through as-is — array indels are rare and not minus-strand
    reported.
    """
    if len(user_allele) != 1:
        return user_allele
    if user_allele in (gnomad_ref, gnomad_alt):
        return user_allele
    comp = _COMPLEMENT.get(user_allele)
    if comp is None:
        return None
    if comp in (gnomad_ref, gnomad_alt):
        if (gnomad_ref, gnomad_alt) in _PALINDROMIC:
            return None
        return comp
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

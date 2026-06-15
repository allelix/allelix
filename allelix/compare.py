# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Compare two genotype files for coverage and concordance.

Classifies each shared rsID into one of five buckets: concordant,
strand-flip match, discordant, strand-ambiguous, or no-call. Uses
:func:`allelix.utils.allele.complement` and
:func:`allelix.utils.allele.is_strand_ambiguous` for strand awareness
so that complementary genotypes are not reported as disagreements.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from allelix.utils.allele import complement, is_strand_ambiguous

if TYPE_CHECKING:
    from allelix.models import Variant


@dataclass(frozen=True)
class CompareResult:
    """Result of comparing two genotype files."""

    file1_total: int
    file2_total: int
    shared: int
    file1_only: int
    file2_only: int
    concordant: int
    strand_flip_match: int
    discordant: int
    strand_ambiguous: int
    no_call: int
    build1: str
    build2: str
    chromosome_counts: dict[str, Counter[str]] = field(default_factory=dict)


def _is_ambiguous_across(v1: Variant, v2: Variant) -> bool:
    """True if the disagreement could be a strand flip on a palindromic SNP.

    Checks whether all alleles observed across both files form an A/T or
    C/G pair — in which case a strand flip is indistinguishable from a
    genuine genotype difference.
    """
    all_alleles = {v1.allele1, v1.allele2, v2.allele1, v2.allele2}
    if len(all_alleles) != 2:
        return False
    a, b = all_alleles
    return is_strand_ambiguous(a, b)


def _genotype_set(v: Variant) -> frozenset[str]:
    """Unordered allele pair for comparison."""
    return frozenset((v.allele1, v.allele2))


def _flipped_set(v: Variant) -> frozenset[str]:
    """Complement of the allele pair."""
    return frozenset((complement(v.allele1), complement(v.allele2)))


def compare_variants(
    variants1: list[Variant],
    variants2: list[Variant],
    *,
    build1: str = "",
    build2: str = "",
) -> CompareResult:
    """Compare two lists of variants by rsID.

    For each shared rsID, classifies the pair as concordant,
    strand-flip match, discordant, strand-ambiguous, or no-call.

    Args:
        variants1: Parsed variants from the first file.
        variants2: Parsed variants from the second file.
        build1: Detected build for file 1. Falls back to first variant's
            build field if empty, or ``"unknown"`` if the list is empty.
        build2: Detected build for file 2. Same fallback logic.
    """
    map1: dict[str, Variant] = {v.rsid: v for v in variants1}
    map2: dict[str, Variant] = {v.rsid: v for v in variants2}

    shared_rsids = set(map1) & set(map2)
    file1_only = len(map1) - len(shared_rsids)
    file2_only = len(map2) - len(shared_rsids)

    concordant = 0
    strand_flip_match = 0
    discordant = 0
    strand_ambiguous = 0
    no_call = 0

    chrom_counts: dict[str, Counter[str]] = {}

    if not build1:
        build1 = variants1[0].build if variants1 else "unknown"
    if not build2:
        build2 = variants2[0].build if variants2 else "unknown"

    for rsid in sorted(shared_rsids):
        v1 = map1[rsid]
        v2 = map2[rsid]

        chrom = v1.chromosome
        if chrom not in chrom_counts:
            chrom_counts[chrom] = Counter()

        if v1.is_no_call or v2.is_no_call:
            no_call += 1
            chrom_counts[chrom]["no_call"] += 1
            continue

        g1 = _genotype_set(v1)
        g2 = _genotype_set(v2)

        if g1 == g2:
            concordant += 1
            chrom_counts[chrom]["concordant"] += 1
        elif _is_ambiguous_across(v1, v2):
            strand_ambiguous += 1
            chrom_counts[chrom]["strand_ambiguous"] += 1
        elif _flipped_set(v1) == g2:
            strand_flip_match += 1
            chrom_counts[chrom]["strand_flip_match"] += 1
        else:
            discordant += 1
            chrom_counts[chrom]["discordant"] += 1

    return CompareResult(
        file1_total=len(map1),
        file2_total=len(map2),
        shared=len(shared_rsids),
        file1_only=file1_only,
        file2_only=file2_only,
        concordant=concordant,
        strand_flip_match=strand_flip_match,
        discordant=discordant,
        strand_ambiguous=strand_ambiguous,
        no_call=no_call,
        build1=build1,
        build2=build2,
        chromosome_counts=chrom_counts,
    )

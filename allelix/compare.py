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

from allelix.parsers._helpers import normalize_chromosome, summary_chromosome
from allelix.parsers.vcf import _REFERENCE_BLOCK_ALTS
from allelix.utils.allele import complement, is_strand_ambiguous
from allelix.utils.build_detect import SUPPORTED_OVERRIDES

if TYPE_CHECKING:
    from allelix.models import Variant


@dataclass(frozen=True)
class CompareResult:
    """Comparison coverage and genotype categories.

    Input named/unnamed counts include repeated parsed calls. Totals and
    eligible named/unnamed counts describe unique comparison keys; omitted
    unnamed counts describe parsed observations lacking eligible identity.
    """

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
    file1_input_named: int = 0
    file2_input_named: int = 0
    file1_input_unnamed: int = 0
    file2_input_unnamed: int = 0
    file1_named: int = 0
    file2_named: int = 0
    file1_unnamed: int = 0
    file2_unnamed: int = 0
    file1_unnamed_omitted: int = 0
    file2_unnamed_omitted: int = 0
    shared_named: int = 0
    shared_unnamed: int = 0
    # Chip calls not compared, by reason (NO_EXPLAINING_RECORD,
    # SEVERAL_EXPLAINING_RECORDS); not counted as file-only either.
    chip_not_compared: Counter[str] = field(default_factory=Counter)


ComparisonKey = (
    tuple[str, str]
    | tuple[str, str, str, tuple[str, ...]]
    | tuple[str, str, int, str, tuple[str, ...]]
)
# Why a chip call was not compared with a VCF file's records under its rsID
# (ADR-0054): none of them, or more than one, explains its genotype.
NO_EXPLAINING_RECORD = "no_explaining_record"
SEVERAL_EXPLAINING_RECORDS = "several_explaining_records"


@dataclass(frozen=True)
class _Observation:
    """A call in comparison form: its alleles trimmed like its record's REF/ALT."""

    variant: Variant
    allele1: str
    allele2: str
    # The record's trimmed REF and ALTs; empty for a call without them (a chip).
    declared: tuple[str, ...] = ()


def _trimmed(variant: Variant) -> tuple[int, str, tuple[str, ...], dict[str, str]] | None:
    """The record's position, REF and ALTs without shared bases, and its allele map.

    REF and every ALT lose the bases they all share at the end, then at the
    start, while each keeps one base (ADR-0047, ADR-0061; no left-alignment).
    A dropped leading base moves the position on by one. A gVCF's
    reference-block ALT (``<NON_REF>``, ``<*>``) names no allele and is left
    out, so ``A G,<NON_REF>`` is the record ``A G``. None when the record has
    no literal input REF and ALT.
    """
    if not (
        variant.ref_from_input
        and variant.ref not in (None, "", ".", "-")
        and variant.input_alt
        and all(alt not in ("", ".") for alt in variant.input_alt.split(","))
    ):
        return None
    assert variant.ref is not None and variant.input_alt is not None
    alts = [a for a in variant.input_alt.upper().split(",") if a not in _REFERENCE_BLOCK_ALTS]
    if not alts:
        return None
    alleles = [variant.ref.upper(), *alts]
    start, end, position = 0, 0, variant.position
    while (
        all(len(a) - start - end > 1 for a in alleles)
        and len({a[len(a) - end - 1] for a in alleles}) == 1
    ):
        end += 1
    while all(len(a) - start - end > 1 for a in alleles) and len({a[start] for a in alleles}) == 1:
        start += 1
        position += 1
    trimmed = [a[start : len(a) - end] for a in alleles]
    return position, trimmed[0], tuple(trimmed[1:]), dict(zip(alleles, trimmed, strict=True))


def _comparison_map(
    variants: list[Variant], *, positional: bool
) -> tuple[dict[ComparisonKey, _Observation], int]:
    """Keep the last observation per key, counting ineligible unnamed rows.

    A named record with REF and ALT is keyed by its rsID and trimmed alleles,
    so the records of a split multi-allelic site keep their own identities; a
    named call without them (a chip) by its rsID. An unnamed record is keyed
    by its trimmed locus and alleles (ADR-0054).
    """
    found: dict[ComparisonKey, _Observation] = {}
    omitted = 0
    for variant in variants:
        trimmed = _trimmed(variant)
        if trimmed is None:
            if variant.rsid:
                found[("id", variant.rsid)] = _Observation(
                    variant, variant.allele1, variant.allele2
                )
            else:
                omitted += 1
            continue
        position, ref, alts, mapping = trimmed
        observation = _Observation(
            variant,
            mapping.get(variant.allele1.upper(), variant.allele1),
            mapping.get(variant.allele2.upper(), variant.allele2),
            (ref, *alts),
        )
        if variant.rsid:
            found[("id", variant.rsid, ref, alts)] = observation
        elif positional:
            found[("locus", normalize_chromosome(variant.chromosome), position, ref, alts)] = (
                observation
            )
        else:
            omitted += 1
    return found, omitted


def _explains(chip: _Observation, record: _Observation) -> bool:
    """True when the record's declared alleles hold the chip's call, on either strand."""
    if chip.variant.is_no_call:
        return True
    declared = set(record.declared)
    call = {chip.allele1, chip.allele2}
    return call <= declared or {complement(a) for a in call} <= declared


def _chip_partners(
    mine: dict[ComparisonKey, _Observation],
    theirs: dict[ComparisonKey, _Observation],
    taken: set[ComparisonKey],
) -> tuple[list[tuple[ComparisonKey, ComparisonKey]], Counter[str]]:
    """Pair each chip call with the one record under its rsID that explains it.

    The partner is chosen by the record's declared alleles, never by whether
    the genotypes agree. None or several explaining records leave the call
    uncompared, counted by reason.
    """
    pairs: list[tuple[ComparisonKey, ComparisonKey]] = []
    reasons: Counter[str] = Counter()
    # Records with alleles, by rsID: built once, so each chip call is one lookup.
    records: dict[str, list[ComparisonKey]] = {}
    for other in theirs:
        if other[0] == "id" and len(other) == 4:
            records.setdefault(str(other[1]), []).append(other)
    for key, chip in mine.items():
        if len(key) != 2 or key in theirs:
            continue
        candidates = [other for other in records.get(str(key[1]), ()) if other not in taken]
        if not candidates:
            continue
        explaining = [other for other in candidates if _explains(chip, theirs[other])]
        if len(explaining) == 1:
            pairs.append((key, explaining[0]))
            taken.add(explaining[0])
        else:
            reasons[SEVERAL_EXPLAINING_RECORDS if explaining else NO_EXPLAINING_RECORD] += 1
    return pairs, reasons


def _is_ambiguous_across(o1: _Observation, o2: _Observation) -> bool:
    """True if the disagreement could be a strand flip on a palindromic SNP.

    Checks whether all alleles observed across both files form an A/T or
    C/G pair — in which case a strand flip is indistinguishable from a
    genuine genotype difference.
    """
    all_alleles = {o1.allele1, o1.allele2, o2.allele1, o2.allele2}
    if len(all_alleles) != 2:
        return False
    a, b = all_alleles
    return is_strand_ambiguous(a, b)


def _genotype_set(o: _Observation) -> frozenset[str]:
    """Unordered allele pair for comparison."""
    return frozenset((o.allele1, o.allele2))


def _flipped_set(o: _Observation) -> frozenset[str]:
    """Complement of the allele pair."""
    return frozenset((complement(o.allele1), complement(o.allele2)))


def compare_variants(
    variants1: list[Variant],
    variants2: list[Variant],
    *,
    build1: str = "",
    build2: str = "",
) -> CompareResult:
    """Compare named records and eligible literal unnamed records (ADR-0054).

    A named record with REF and ALT is matched by its rsID and trimmed
    alleles, whatever the order of the records. An unnamed record also needs
    the same known effective build. A chip call (no REF/ALT) is matched with
    the other file's call under its rsID, or with the one record there whose
    alleles explain it. For each pair, classifies the call as concordant,
    strand-flip match, discordant, strand-ambiguous, or no-call.

    Args:
        variants1: Parsed variants from the first file.
        variants2: Parsed variants from the second file.
        build1: Detected build for file 1. Falls back to first variant's
            build field if empty, or ``"unknown"`` if the list is empty.
        build2: Detected build for file 2. Same fallback logic.
    """
    if not build1:
        build1 = variants1[0].build if variants1 else "unknown"
    if not build2:
        build2 = variants2[0].build if variants2 else "unknown"
    positional = build1 == build2 and build1 in SUPPORTED_OVERRIDES
    map1, omitted1 = _comparison_map(variants1, positional=positional)
    map2, omitted2 = _comparison_map(variants2, positional=positional)

    shared_keys = set(map1) & set(map2)
    pairs: list[tuple[ComparisonKey, ComparisonKey]] = [
        (k, k) for k in sorted(shared_keys, key=repr)
    ]
    taken1, taken2 = set(shared_keys), set(shared_keys)
    chip_pairs1, reasons1 = _chip_partners(map1, map2, taken2)
    chip_pairs2, reasons2 = _chip_partners(map2, map1, taken1)
    pairs += chip_pairs1 + [(mine, chip) for chip, mine in chip_pairs2]
    not_compared1, not_compared2 = sum(reasons1.values()), sum(reasons2.values())

    concordant = 0
    strand_flip_match = 0
    discordant = 0
    strand_ambiguous = 0
    no_call = 0

    chrom_counts: dict[str, Counter[str]] = {}

    for key1, key2 in pairs:
        o1 = map1[key1]
        o2 = map2[key2]

        # A standard chromosome under its lookup name, so 23 in one file and X
        # in the other count together; any other contig as the reader gave it.
        chrom = summary_chromosome(o1.variant.chromosome)
        if chrom not in chrom_counts:
            chrom_counts[chrom] = Counter()

        if o1.variant.is_no_call or o2.variant.is_no_call:
            no_call += 1
            chrom_counts[chrom]["no_call"] += 1
            continue

        g1 = _genotype_set(o1)
        g2 = _genotype_set(o2)

        if g1 == g2:
            concordant += 1
            chrom_counts[chrom]["concordant"] += 1
        elif _is_ambiguous_across(o1, o2):
            strand_ambiguous += 1
            chrom_counts[chrom]["strand_ambiguous"] += 1
        elif _flipped_set(o1) == g2:
            strand_flip_match += 1
            chrom_counts[chrom]["strand_flip_match"] += 1
        else:
            discordant += 1
            chrom_counts[chrom]["discordant"] += 1

    return CompareResult(
        file1_total=len(map1),
        file2_total=len(map2),
        shared=len(pairs),
        file1_only=len(map1) - len(pairs) - not_compared1,
        file2_only=len(map2) - len(pairs) - not_compared2,
        concordant=concordant,
        strand_flip_match=strand_flip_match,
        discordant=discordant,
        strand_ambiguous=strand_ambiguous,
        no_call=no_call,
        build1=build1,
        build2=build2,
        chromosome_counts=chrom_counts,
        file1_input_named=sum(bool(v.rsid) for v in variants1),
        file2_input_named=sum(bool(v.rsid) for v in variants2),
        file1_input_unnamed=sum(not v.rsid for v in variants1),
        file2_input_unnamed=sum(not v.rsid for v in variants2),
        file1_named=sum(key[0] == "id" for key in map1),
        file2_named=sum(key[0] == "id" for key in map2),
        file1_unnamed=sum(key[0] == "locus" for key in map1),
        file2_unnamed=sum(key[0] == "locus" for key in map2),
        file1_unnamed_omitted=omitted1,
        file2_unnamed_omitted=omitted2,
        shared_named=sum(key[0] == "id" for key, _ in pairs),
        shared_unnamed=sum(key[0] == "locus" for key, _ in pairs),
        chip_not_compared=reasons1 + reasons2,
    )

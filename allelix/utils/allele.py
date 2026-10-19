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

from typing import TYPE_CHECKING

from allelix.models import NO_CALL_MARKER

if TYPE_CHECKING:
    from collections.abc import Iterable

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

    previously, when the user allele did NOT directly match
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


def strand_aware_carrier_match(
    variant_ref: str | None,
    allele1: str,
    allele2: str,
    source_ref: str,
    source_alt: str,
) -> bool:
    """Return True if the user carries ``source_alt`` directly or via strand-flip.

    Multi-allelic safety (since v2.0.1): only fires strand-flip when the
    ``variant_ref`` context confirms the user is reverse-stranded relative to the
    source's ref. Without that context — ``variant_ref`` is None, or it agrees
    with ``source_ref`` (forward) — strand-flip never fires and the carrier
    match falls back to direct membership only. Skips palindromic ``(ref, alt)``
    pairs (A/T, C/G) where strand cannot be inferred from sequence alone
    (ADR-0010 documents the limitation).

    Args:
        variant_ref: User's claimed REF at the variant locus. ``None`` means
            the parser couldn't supply one (array data prior to gnomAD-backed
            population, or the source didn't carry a REF column).
        allele1: First user allele (upper-cased at ``Variant`` construction).
        allele2: Second user allele.
        source_ref: REF from the annotator's source row (forward strand).
        source_alt: ALT from the annotator's source row (forward strand).

    Returns:
        True when the user carries ``source_alt`` on either strand. False
        otherwise — including the abstain cases (palindromic site, no
        ``variant_ref`` context with no direct match, ref disagreement that
        isn't a clean complement).
    """
    # Defensive case normalization: the Variant model upper-cases its own
    # fields on construction, but `source_ref` / `source_alt` come from
    # external annotator queries (ClinVar VCF, ClinPGx TSV) and a future
    # source variant or a soft-masked-reference VCF could leak lowercase.
    # A lowercase letter would silently false-negative every comparison
    # below — normalize at the entry point so the function is robust.
    source_ref = source_ref.upper()
    source_alt = source_alt.upper()
    if variant_ref is not None:
        variant_ref = variant_ref.upper()

    user = {allele1, allele2}
    if source_alt in user:
        return True
    if variant_ref is None:
        return False
    if is_strand_ambiguous(source_ref, source_alt):
        return False
    if variant_ref == source_ref:
        return False
    if variant_ref == complement(source_ref):
        return complement(source_alt) in user
    return False


def carried_in_source_orientation(allele1: str, allele2: str, source_alt: str) -> tuple[str, str]:
    """Return the user's alleles in the orientation of ``source_alt``.

    Called after ``strand_aware_carrier_match`` accepted a row: either the
    user carries ``source_alt`` directly (forward reading, alleles returned
    unchanged) or the match fired through the strand-flip path, in which
    case the complemented pair is the user's call on the source's strand.
    Feeds ``Annotation.carried_alleles`` for the zero-carried-copies
    invariant.
    """
    if source_alt in (allele1, allele2):
        return allele1, allele2
    return flip_genotype(allele1, allele2)


def trim_allele_pair(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    """The VCF-style minimal form of one REF/ALT pair (ADR-0047).

    Drops the bases REF and ALT share at the end, then at the start, while
    both keep at least one base, so an indel keeps its anchor base. A
    leading base dropped moves the position on by one. ``GTT>GT`` becomes
    ``GT>G`` at the same position; ``CA>TA`` becomes ``C>T``; ``ACG>ATG``
    becomes ``C>T`` one base on. No left-alignment: that needs the
    reference sequence.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, position = ref[1:], alt[1:], position + 1
    return position, ref, alt


def resolve_rsid_at_locus(
    variant_ref: str | None,
    allele1: str,
    allele2: str,
    rows: Iterable[tuple[str, str, str]],
) -> tuple[str, tuple[str, str] | None] | None:
    """Pick the rsID for an rsID-less call from the ``(ref, alt, rsid)`` rows at its locus.

    Shared by every position resolver (ClinVar, gnomAD) so they apply one
    rule. A row is consistent with the call when its REF differs from its
    ALT (a REF == ALT row names the reference base itself, which no
    non-reference call explains), its REF equals ``variant_ref``
    when that is known, and the user carries its ALT. The call resolves
    when those rows:

    - explain every allele the user carries under exactly one REF (the
      REF itself or one of that REF's ALTs). This picks the right anchor
      between co-located indels (``CA>C`` vs ``CAG>C``) even when
      ``variant_ref`` is unknown, and lets a ``1/2`` call carrying two
      alts resolve when both alts belong to the same rsID, and
    - name exactly one rsID. Duplicate rows (per-SCV submissions repeat
      ``(ref, alt, rsid)``) are one rsID, not an ambiguity.

    Returns ``(rsid, (ref, alt))`` when the user carries exactly one alt,
    ``(rsid, None)`` when they carry two, or ``None`` to abstain.
    """
    user = {allele1, allele2}
    by_ref: dict[str, set[tuple[str, str]]] = {}
    for ref, alt, rsid in rows:
        if rsid and ref != alt and alt in user and (variant_ref is None or ref == variant_ref):
            by_ref.setdefault(ref, set()).add((alt, rsid))
    explaining = [
        (ref, entries)
        for ref, entries in by_ref.items()
        if user <= {alt for alt, _rsid in entries} | {ref}
    ]
    if len(explaining) != 1:
        return None
    ref, entries = explaining[0]
    rsids = {rsid for _alt, rsid in entries}
    if len(rsids) != 1:
        return None
    alts = {alt for alt, _rsid in entries}
    pair = (ref, next(iter(alts))) if len(alts) == 1 else None
    return rsids.pop(), pair


def strand_aware_genotype_match(
    allele1: str,
    allele2: str,
    source_geno: str,
    variant_ref: str | None,
    *,
    site_palindromic: bool = False,
) -> bool:
    """Per-genotype carrier match with strand-flip support for ClinPGx / SNPedia.

    PharmGKB / ClinPGx and SNPedia store the user's full diploid (sorted SNV
    pair like ``"AG"``) rather than per-allele rows. Direct equality remains
    the primary match. Strand-flip fires only when ``variant_ref`` confirms
    reverse orientation against the source diploid AND the source diploid is
    not palindromic (A/T or C/G).

    Args:
        allele1: First user allele (upper-cased single base for SNVs).
        allele2: Second user allele.
        source_geno: Source row's sorted diploid string (e.g. ``"AG"``).
        variant_ref: User's claimed REF at this locus. ``None`` means no
            reference context; strand-flip is skipped (direct match only).
        site_palindromic: The identity table flags the site palindromic
            (its full allele set equals its own complement, e.g. a
            multi-allelic site carrying all four bases). The strand cannot
            be told from sequence there, so the complement path is skipped;
            a direct match still counts.

    Returns:
        True when the user matches ``source_geno`` directly or via a
        complement reading. False otherwise.
    """
    # Defensive case normalization (see note on strand_aware_carrier_match):
    # source_geno comes from a third-party SQLite cache; lowercase would
    # silently false-negative the sorted-pair equality below.
    source_geno = source_geno.upper()
    if variant_ref is not None:
        variant_ref = variant_ref.upper()

    if len(allele1) != 1 or len(allele2) != 1 or len(source_geno) != 2:
        return False
    user_normalized = "".join(sorted((allele1, allele2)))
    if user_normalized == source_geno:
        return True
    if variant_ref is None:
        return False
    if site_palindromic:
        return False
    user_complement_pair = sorted((complement(allele1), complement(allele2)))
    if "".join(user_complement_pair) != source_geno:
        return False
    if variant_ref in source_geno:
        return False
    if complement(variant_ref) not in source_geno:
        return False
    # Site-level palindrome guard: include variant_ref so homozygous
    # genotypes at palindromic sites are correctly caught.  set("GG") =
    # {"G"} missed the C/G palindrome; {"G"} | {"C"} = {"C","G"} catches
    # it.
    site_alleles = set(source_geno) | {variant_ref}
    return site_alleles != {"A", "T"} and site_alleles != {"C", "G"}


def derive_alt_from_diploid(ref: str | None, allele1: str, allele2: str) -> str:
    """Derive the user's carried alt from a matched diploid given a known REF.

    For per-genotype annotators (SNPedia, ClinPGx) the matched row describes
    the user's full diploid pair; the alt the user carries is whichever side
    of the pair does not equal REF. Returns the empty string when alt cannot
    be cleanly identified:

    - ``ref`` is None — no reference context yet (e.g., array data prior to
      strand-aware ref population, ADR-0035).
    - Both alleles equal ref — hom-ref; should be filtered earlier
      (ADR-0023 for ClinPGx; SNPedia rarely matches hom-ref) but defensive.
    - Neither allele equals ref — multi-allelic het where the matched row's
      alleles disagree with the reference (rare for SNPs); conservatively
      skip rather than guess.

    Conservative-on-uncertainty mirrors the v2.0.1 suppress-half: when
    the alt can't be cleanly identified, return "" rather than a guess that
    could stamp the wrong allele's gnomAD / AlphaMissense / CADD value.
    """
    if ref is None:
        return ""
    if allele1 == ref and allele2 == ref:
        return ""
    if allele1 == ref:
        return allele2
    if allele2 == ref:
        return allele1
    if allele1 == allele2:
        return allele1
    return ""


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


# Strand readings of a source's alleles against a site's forward alleles
# (ADR-0024): the two readings, then the three reasons to abstain.
STRAND_DIRECT = "direct"
STRAND_COMPLEMENT = "complement"
STRAND_NOT_IN_IDENTITY = "not_in_identity_table"
STRAND_PALINDROMIC = "palindromic"
STRAND_NO_FIT = "no_fit"
STRAND_ABSTAIN_REASONS: tuple[str, ...] = (
    STRAND_NOT_IN_IDENTITY,
    STRAND_PALINDROMIC,
    STRAND_NO_FIT,
)
_BASES = frozenset("ACGT")


def is_single_base(allele: str) -> bool:
    """True for one of A, C, G, T."""
    return allele in _BASES


def forward_strand_reading(alleles: Iterable[str], forward: frozenset[str]) -> str:
    """How a source's single-base ``alleles`` read on the forward strand.

    ``forward`` is the site's REF plus ALTs on the input's build (the rsID
    identity table). The abstain-first rule (ADR-0024):

    1. Not an A/T or C/G site, and every allele is forward: ``direct``,
       even where the complements are forward too.
    2. Otherwise, every complement is forward: ``complement``.
    3. An A/T or C/G site is ``palindromic``; anything else is ``no_fit``.
    """
    if forward in _AMBIGUOUS_PAIRS:
        return STRAND_PALINDROMIC
    wanted = set(alleles)
    if not wanted or not wanted <= _BASES:
        return STRAND_NO_FIT
    if wanted <= forward:
        return STRAND_DIRECT
    if {complement(a) for a in wanted} <= forward:
        return STRAND_COMPLEMENT
    return STRAND_NO_FIT

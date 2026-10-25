# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""rsID-less resolution: one step for every command.

Variant callers emit ``ID=.`` (GATK, DeepVariant, DRAGEN), positional
pseudo-IDs (``22:10519265:CA:C``) or vendor IDs. Everything keyed by rsID
(the annotators, the high-value list, PLINK allele coding, extract) needs
the rsID first. This module is the only place that recovers it, so
``analyze``, ``stats``, ``extract`` and ``export plink`` all see the same
result.

``resolve_batch`` runs, in order, on variants whose ``build`` is already
the effective build:

1. ``Variant.ref`` from the identity table, for rsID-carrying calls only
   (ADR-0040). This runs first so ClinVar's REF check can confirm array
   calls. A vendor probe without an rsID never gets a REF from position
   alone, so it stays under ClinVar's no-REF, no-rsID guard.
2. ClinVar position resolution (selected build).
3. Identity-table position resolution (GRCh37 and GRCh38).
4. gnomAD position resolution (GRCh38): the full dbSNP rsID universe for
   pharmacogenomic / GWAS / SNPedia rsIDs ClinVar does not hold.
5. ``Variant.ref`` from gnomAD's rsID coordinates (GRCh38), accepted only
   when the coordinate is the variant's own locus.

Every resolver applies ``resolve_rsid_at_locus``: exactly one rsID whose
rows explain the user's alleles, abstain otherwise. A resolved variant
keeps its pre-resolution ID in ``Variant.original_rsid`` (debug only).
Variants no resolver can place pass through unchanged.

Build detection cannot use this step (resolution needs the build); it
votes on rsID-less records by position instead (``utils.build_detect``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.annotators.base import variant_locus
from allelix.utils.allele import resolve_rsid_at_locus
from allelix.utils.build_detect import BUILD_GRCH37, BUILD_GRCH38

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.models import Variant

# Resolved rsID -> the (chrom, pos, ref, alt) it was resolved at, for the
# annotated-VCF writer's ALLELIX_ORIGINAL_ID and position-keyed lookups.
ResolvedCoords = dict[str, tuple[str, int, str, str]]

_IDENTITY_BUILDS: frozenset[str] = frozenset({BUILD_GRCH37, BUILD_GRCH38})

# Variants per resolver round trip when streaming.
_BATCH_SIZE = 1000


def resolve_batch(
    variants: list[Variant],
    *,
    clinvar: ClinVarAnnotator | None,
    gnomad: GnomadAnnotator | None,
    identity: RsidIdentityAnnotator | None = None,
    resolved_coords: ResolvedCoords | None = None,
) -> None:
    """Populate ``ref`` and stamp ``rsid`` in place on one batch (see module docstring).

    ``variants`` must already carry their effective ``build``. Sources
    that are None or not ready are skipped. Resolved coordinates are
    added to ``resolved_coords`` when given.
    """
    if not variants:
        return
    coords: ResolvedCoords = resolved_coords if resolved_coords is not None else {}
    identity_ready = identity is not None and identity.is_ready()
    gnomad_ready = gnomad is not None and gnomad.is_ready()
    if identity_ready:
        assert identity is not None
        _populate_ref_from_identity(identity, variants)
    if clinvar is not None:
        # bulk_resolve_rsids skips any build whose cache is missing or not
        # ready for this version.
        rsidless = [v for v in variants if not v.rsid.startswith("rs")]
        if rsidless:
            for (chrom, pos, ref, alt), rsid in clinvar.bulk_resolve_rsids(rsidless).items():
                coords[rsid] = (chrom, pos, ref, alt)
    if identity_ready:
        assert identity is not None
        _resolve_from_identity(identity, variants, coords)
    if gnomad_ready:
        assert gnomad is not None
        _resolve_from_gnomad(gnomad, variants, coords)
        _populate_ref_from_gnomad(gnomad, variants)


def resolve_rsidless(
    variants: Iterable[Variant],
    *,
    clinvar: ClinVarAnnotator | None,
    gnomad: GnomadAnnotator | None,
    identity: RsidIdentityAnnotator | None = None,
    resolved_coords: ResolvedCoords | None = None,
    batch_size: int = _BATCH_SIZE,
) -> Iterator[Variant]:
    """Stream ``variants`` through ``resolve_batch`` in batches, preserving order.

    For command paths that consume a parser iterator (``stats``, ``extract``,
    ``export plink``). Each variant must carry its effective ``build``.
    """
    batch: list[Variant] = []
    for variant in variants:
        batch.append(variant)
        if len(batch) >= batch_size:
            resolve_batch(
                batch,
                clinvar=clinvar,
                gnomad=gnomad,
                identity=identity,
                resolved_coords=resolved_coords,
            )
            yield from batch
            batch = []
    resolve_batch(
        batch, clinvar=clinvar, gnomad=gnomad, identity=identity, resolved_coords=resolved_coords
    )
    yield from batch


def _populate_ref_from_identity(identity: RsidIdentityAnnotator, variants: list[Variant]) -> None:
    """Set ``Variant.ref`` from the identity table's REF for the call's own rsID.

    Only for calls that carry an rsID, only from the identity rows of that
    rsID at the variant's locus, and only when those rows give a single REF
    of the call's shape: single-base for single-base calls (an indel row's
    anchor REF can share the position), any length otherwise.

    A call whose rsID the identity table does not hold (a withdrawn or
    merged ID) gets no REF, even when another rsID's rows sit at the same
    position: nothing confirms that the other rsID is the same variant.
    Such a call keeps ``ref`` None, and ClinVar matches it by locus only
    when ClinVar holds its rsID.

    A call without an rsID (array vendor probes such as 23andMe
    ``i5000123``) never gets a REF from position alone. The REF at a locus
    says nothing about which strand the vendor reported the probe's alleles
    on, so with an inferred REF the strand-aware ClinVar match would accept
    the complement of a homozygous-reference call as a homozygous
    pathogenic genotype. Keeping ``ref`` None leaves such probes under
    ClinVar's no-REF, no-rsID guard and out of every position resolver.
    """
    pending = [
        v
        for v in variants
        if v.ref is None
        and v.rsid.startswith("rs")
        and v.chromosome
        and v.position > 0
        and v.build in _IDENTITY_BUILDS
    ]
    for build in {v.build for v in pending}:
        group = [v for v in pending if v.build == build]
        rows_at = identity.alleles_at_loci(build, {variant_locus(v) for v in group})
        for v in group:
            refs = {ref for ref, _alt, rsid in rows_at.get(variant_locus(v), []) if rsid == v.rsid}
            if len(v.allele1) <= 1 and len(v.allele2) <= 1:
                refs = {r for r in refs if len(r) == 1}
            if len(refs) == 1:
                v.ref = next(iter(refs))


def _resolve_from_identity(
    identity: RsidIdentityAnnotator, variants: list[Variant], coords: ResolvedCoords
) -> None:
    """Stamp rsIDs on rsID-less calls from the identity table (ADR-0040)."""
    pending = [
        v
        for v in variants
        if not v.rsid.startswith("rs")
        and v.ref is not None
        and (not v.is_no_call or v.is_half_call)
        and v.chromosome
        and v.position > 0
        and v.build in _IDENTITY_BUILDS
    ]
    for build in {v.build for v in pending}:
        group = [v for v in pending if v.build == build]
        rows_at = identity.alleles_at_loci(build, {variant_locus(v) for v in group})
        for v in group:
            rows = rows_at.get(variant_locus(v))
            if not rows:
                continue
            hit = resolve_rsid_at_locus(v.ref, *v.carried_pair, rows)
            if hit is None:
                continue
            rsid, pair = hit
            # Same audit stash as the ClinVar resolver. Debug-only.
            v.original_rsid = v.rsid
            v.rsid = rsid
            if pair is not None:
                coords[rsid] = (v.chromosome, v.position, *pair)


def _resolve_from_gnomad(
    gnomad: GnomadAnnotator, variants: list[Variant], coords: ResolvedCoords
) -> None:
    """Stamp rsIDs from gnomAD's position index on GRCh38 calls still rsID-less.

    A call with no REF is skipped: nothing would confirm that the gnomAD
    row at the position is the same variant (array vendor probes).
    """
    pending = [
        v
        for v in variants
        if v.build == BUILD_GRCH38 and not v.rsid.startswith("rs") and v.ref is not None
    ]
    positions = {
        (v.chromosome, v.position)
        for v in pending
        if (not v.is_no_call or v.is_half_call) and v.chromosome and v.position > 0
    }
    if not positions:
        return
    gnomad_rows = gnomad.bulk_resolve_rsids_from_positions(positions)
    for v in pending:
        rows = gnomad_rows.get((v.chromosome, v.position))
        if not rows:
            continue
        hit = resolve_rsid_at_locus(v.ref, *v.carried_pair, rows)
        if hit is None:
            continue
        rsid, pair = hit
        # Same audit stash as the ClinVar resolver
        # (clinvar.py:bulk_resolve_rsids). Debug-only.
        v.original_rsid = v.rsid
        v.rsid = rsid
        if pair is not None:
            coords[rsid] = (v.chromosome, v.position, *pair)


def _populate_ref_from_gnomad(gnomad: GnomadAnnotator, variants: list[Variant]) -> None:
    """GRCh38 ``Variant.ref`` fallback from gnomAD's rsID coordinates.

    Accepted only when the coordinate is the variant's own locus and its
    REF is unique there. Variants no source covers keep ``ref=None`` and
    degrade gracefully (direct match, no strand flip).
    """
    need_ref = {
        v.rsid
        for v in variants
        if v.build == BUILD_GRCH38 and v.ref is None and v.rsid.startswith("rs")
    }
    if not need_ref:
        return
    coord_map = gnomad.bulk_resolve_coordinates(need_ref)
    for v in variants:
        if v.build != BUILD_GRCH38 or v.ref is not None:
            continue
        refs = {
            ref
            for chrom, pos, ref, _alt in coord_map.get(v.rsid, [])
            if (chrom, pos) == variant_locus(v)
        }
        if len(refs) == 1:
            v.ref = refs.pop()

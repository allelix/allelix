# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""PLINK1 binary format (.bed/.bim/.fam) exporter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.utils.allele import complement, is_strand_ambiguous

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping
    from pathlib import Path

    from allelix.models import Variant

_BED_MAGIC = bytes([0x6C, 0x1B, 0x01])

_CHROM_CODES = {
    "X": "23",
    "Y": "24",
    "MT": "26",
}


def _orient_genotype(
    allele1: str,
    allele2: str,
    ref: str,
    alt: str,
) -> tuple[str, str] | None:
    """Map user alleles to {ref, alt} in a consistent orientation.

    Returns None for palindromic sites, indels, or alleles that don't fit.
    Both alleles are tested in the same orientation — no mixed-strand.
    """
    if len(allele1) != 1 or len(allele2) != 1:
        return None
    if is_strand_ambiguous(ref, alt):
        return None

    pair = {allele1, allele2}
    if pair <= {ref, alt}:
        return (allele1, allele2)

    c1, c2 = complement(allele1), complement(allele2)
    if {c1, c2} <= {ref, alt}:
        return (c1, c2)

    return None


def export_plink(
    variants: Iterator[Variant],
    prefix: Path,
    build: str,
    ref_alt_map: dict[str, tuple[str, str]] | None = None,
) -> tuple[int, int, int, int]:
    """Write .bed/.bim/.fam from parsed variants.

    Args:
        variants: Parsed variant iterator (consumed once).
        prefix: Base path for output files.
        build: Genome build label (informational, not used for liftover).
        ref_alt_map: ``{rsid: (ref, alt)}`` from gnomAD coordinate resolution.
            When provided, uses ref/alt to assign A1/A2 for proper allele coding.
            When None or rsid missing, falls back to ``A2="0"`` for homozygotes.

    Returns:
        ``(variants_written, no_calls_skipped, indels_skipped, monomorphic_count)``

    Note:
        No-call variants and indels (multi-character alleles) are dropped.
        PLINK1 BIM is SNV-only (single-character A1/A2). Indels would
        produce non-standard BIM rows that downstream tools may reject.
    """
    fam_path = prefix.with_suffix(".fam")
    bim_path = prefix.with_suffix(".bim")
    bed_path = prefix.with_suffix(".bed")

    fam_path.write_text("0\tSAMPLE\t0\t0\t0\t-9\n")

    written = 0
    skipped = 0
    indels = 0
    monomorphic = 0

    with bim_path.open("w") as bim_f, bed_path.open("wb") as bed_f:
        bed_f.write(_BED_MAGIC)

        for v in variants:
            if v.is_no_call:
                skipped += 1
                continue

            if len(v.allele1) != 1 or len(v.allele2) != 1:
                indels += 1
                continue

            chrom_code = _CHROM_CODES.get(v.chromosome, v.chromosome)
            a1: str
            a2: str
            bed_code: int

            if ref_alt_map and v.rsid in ref_alt_map:
                ref, alt = ref_alt_map[v.rsid]
                resolved = _orient_genotype(v.allele1, v.allele2, ref, alt)
                if resolved is not None:
                    r1, r2 = resolved
                    a1 = ref
                    a2 = alt
                    a2_count = sum(1 for a in (r1, r2) if a == alt)
                    if a2_count == 0:
                        bed_code = 0b00
                    elif a2_count == 1:
                        bed_code = 0b10
                    else:
                        bed_code = 0b11
                else:
                    a1, a2, bed_code, is_mono = _fallback_coding(v)
                    if is_mono:
                        monomorphic += 1
            else:
                a1, a2, bed_code, is_mono = _fallback_coding(v)
                if is_mono:
                    monomorphic += 1

            bim_f.write(f"{chrom_code}\t{v.rsid}\t0\t{v.position}\t{a1}\t{a2}\n")
            bed_f.write(bytes([bed_code]))
            written += 1

    return written, skipped, indels, monomorphic


def resolve_ref_alt_via_gnomad(
    variant_by_rsid: Mapping[str, Variant], data_dir: Path
) -> dict[str, tuple[str, str]]:
    """Return per-rsid (ref, alt) for PLINK A1/A2 coding via gnomAD.

    For each rsid in ``variant_by_rsid``, looks up gnomAD's coordinate
    record(s). When gnomAD returns a single record, that ref/alt pair is
    used directly. When gnomAD returns multiple records (multi-allelic
    positions), the user's allele pair is matched against each candidate
    via ``_orient_genotype`` and the first match that's also a subset of
    ``{ref, alt}`` wins; if none is a subset, the first orientation-valid
    record wins as a fallback.

    Returns an empty dict if gnomAD isn't ready. Raises on any underlying
    lookup error — the CLI catches and prints the "failed" warning, then
    proceeds with empty (monomorphic) coding.
    """
    from allelix.annotators.gnomad import GnomadAnnotator

    gnomad = GnomadAnnotator(data_dir)
    try:
        if not gnomad.is_ready():
            return {}
        coord_map = gnomad.bulk_resolve_coordinates(set(variant_by_rsid))
        result: dict[str, tuple[str, str]] = {}
        for rsid, coords in coord_map.items():
            if len(coords) == 1:
                _, _, ref, alt = coords[0]
                result[rsid] = (ref, alt)
                continue
            v = variant_by_rsid[rsid]
            pair = {v.allele1, v.allele2}
            for _, _, ref, alt in coords:
                if _orient_genotype(v.allele1, v.allele2, ref, alt) is not None and pair <= {
                    ref,
                    alt,
                }:
                    result[rsid] = (ref, alt)
                    break
            else:
                for _, _, ref, alt in coords:
                    if _orient_genotype(v.allele1, v.allele2, ref, alt) is not None:
                        result[rsid] = (ref, alt)
                        break
        return result
    finally:
        gnomad.close()


def _fallback_coding(v: Variant) -> tuple[str, str, int, bool]:
    """Fallback allele coding when ref/alt is unknown.

    Returns ``(a1, a2, bed_code, is_monomorphic)``.
    """
    if v.is_heterozygous:
        alleles = sorted([v.allele1, v.allele2])
        return alleles[0], alleles[1], 0b10, False

    return v.allele1, "0", 0b00, True

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""PLINK1 binary format (.bed/.bim/.fam) exporter."""

from __future__ import annotations

import contextlib
import os
from typing import TYPE_CHECKING, NamedTuple

from allelix.utils.allele import is_strand_ambiguous

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping
    from pathlib import Path

    from allelix.models import Variant

_BED_MAGIC = bytes([0x6C, 0x1B, 0x01])


class PlinkCounts(NamedTuple):
    """What ``export_plink`` wrote and skipped.

    ``allele_mismatches`` counts calls written as missing because their
    alleles fit neither gnomAD's REF nor ALT at the site.
    """

    written: int
    no_calls: int
    indels: int
    monomorphic: int
    allele_mismatches: int


class EmptyPlinkSetError(ValueError):
    """No variant could be written, so no PLINK set is produced.

    ``no_calls`` and ``indels`` count the variants that were skipped.
    """

    def __init__(self, no_calls: int, indels: int) -> None:
        """Record what was skipped."""
        super().__init__(f"no variant to write ({no_calls} no-calls, {indels} indels skipped)")
        self.no_calls = no_calls
        self.indels = indels


# PLINK1 .bim alleles: single nucleotides only.
_NUCLEOTIDES = frozenset("ACGT")

_CHROM_CODES = {
    "X": "23",
    "Y": "24",
    "MT": "26",
}


# PLINK1 .bed genotype code for a missing call.
_BED_MISSING = 0b01


def _fits_site(allele1: str, allele2: str, ref: str, alt: str) -> bool:
    """True when both single-base alleles are the site's REF or ALT as written.

    The input formats state their strand, so a call is never read on the
    other strand: one that fits neither allele is a mismatch.
    """
    if len(allele1) != 1 or len(allele2) != 1:
        return False
    return {allele1, allele2} <= {ref, alt}


def export_plink(
    variants: Iterator[Variant],
    prefix: Path,
    build: str,
    ref_alt_map: dict[str, tuple[str, str]] | None = None,
) -> PlinkCounts:
    """Write .bed/.bim/.fam from parsed variants.

    Args:
        variants: Parsed variant iterator (consumed once).
        prefix: Base path for output files.
        build: Genome build label (informational, not used for liftover).
        ref_alt_map: ``{rsid: (ref, alt)}`` from gnomAD coordinate resolution.
            When provided, uses ref/alt to assign A1/A2 for proper allele coding.
            When None or rsid missing, falls back to ``A2="0"`` for homozygotes.
            A call whose alleles fit neither REF nor ALT is written as missing.

    Returns:
        A ``PlinkCounts``.

    Note:
        No-call variants and indels are dropped. PLINK1 BIM is SNV-only
        (A1/A2 in ``ACGT``). Multi-character alleles and single-character
        non-nucleotide codes (array ``I``/``D`` indel codes, a VCF ``*``)
        would produce non-standard BIM rows that downstream tools may
        reject; they are counted in ``indels_skipped``.

        The three files are written to ``.tmp`` siblings and renamed into
        place only after ``variants`` is exhausted. An exception raised while
        it is consumed (for example a late build stop, ADR-0044) removes the
        ``.tmp`` files and leaves no output. A killed process can leave
        ``*.tmp`` files behind, but never a final ``.bed``, ``.bim`` or
        ``.fam``.

    Raises:
        EmptyPlinkSetError: no variant could be written (every one was a
            no-call or an indel, or there were none); no files are left.
    """
    final = {ext: prefix.with_suffix(ext) for ext in (".fam", ".bim", ".bed")}
    tmp = {ext: path.with_name(path.name + ".tmp") for ext, path in final.items()}
    try:
        counts = _write_plink(variants, tmp[".fam"], tmp[".bim"], tmp[".bed"], ref_alt_map)
        if counts[0] == 0:
            raise EmptyPlinkSetError(counts[1], counts[2])
        for ext, path in tmp.items():
            os.replace(path, final[ext])
    except BaseException:
        for path in tmp.values():
            with contextlib.suppress(OSError):
                path.unlink()
        raise
    return counts


def _write_plink(
    variants: Iterator[Variant],
    fam_path: Path,
    bim_path: Path,
    bed_path: Path,
    ref_alt_map: dict[str, tuple[str, str]] | None,
) -> PlinkCounts:
    """The body of ``export_plink``: write the three files at the given paths."""
    fam_path.write_text("0\tSAMPLE\t0\t0\t0\t-9\n")

    written = 0
    skipped = 0
    indels = 0
    monomorphic = 0
    mismatches = 0

    with bim_path.open("w") as bim_f, bed_path.open("wb") as bed_f:
        bed_f.write(_BED_MAGIC)

        for v in variants:
            if v.is_no_call:
                skipped += 1
                continue

            if v.allele1 not in _NUCLEOTIDES or v.allele2 not in _NUCLEOTIDES:
                indels += 1
                continue

            chrom_code = _CHROM_CODES.get(v.chromosome, v.chromosome)
            a1: str
            a2: str
            bed_code: int

            if ref_alt_map and v.rsid in ref_alt_map:
                ref, alt = ref_alt_map[v.rsid]
                if not _fits_site(v.allele1, v.allele2, ref, alt):
                    a1, a2, bed_code = ref, alt, _BED_MISSING
                    mismatches += 1
                elif not is_strand_ambiguous(ref, alt):
                    a1 = ref
                    a2 = alt
                    a2_count = sum(1 for a in (v.allele1, v.allele2) if a == alt)
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

    return PlinkCounts(written, skipped, indels, monomorphic, mismatches)


def resolve_ref_alt_via_gnomad(
    variant_by_rsid: Mapping[str, Variant], data_dir: Path
) -> dict[str, tuple[str, str]]:
    """Return per-rsid (ref, alt) for PLINK A1/A2 coding via gnomAD.

    For each rsid in ``variant_by_rsid``, looks up gnomAD's coordinate
    record(s). When gnomAD returns a single record, that ref/alt pair is
    used directly. When gnomAD returns multiple records (multi-allelic
    positions), the first record the user's alleles fit as written wins,
    preferring one that is not palindromic. When none fits, the first
    record is used, and the exporter writes the call as missing.

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
            fitting = [
                (ref, alt)
                for _, _, ref, alt in coords
                if _fits_site(v.allele1, v.allele2, ref, alt)
            ]
            plain = [pair for pair in fitting if not is_strand_ambiguous(*pair)]
            result[rsid] = (plain or fitting or [coords[0][2:4]])[0]
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

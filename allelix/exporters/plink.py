# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""PLINK1 binary format (.bed/.bim/.fam) exporter."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, NamedTuple

from allelix.utils.allele import is_strand_ambiguous

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping
    from pathlib import Path
    from typing import BinaryIO, TextIO

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

# Genotype-file endings taken off the input name for the default prefix,
# longest first so ``x.vcf.gz`` loses ``.vcf.gz``, not just ``.gz``.
_INPUT_ENDINGS = (".vcf.gz", ".vcf.bgz", ".txt", ".csv", ".tsv", ".vcf", ".zip", ".gz")


def _fits_site(allele1: str, allele2: str, ref: str, alt: str) -> bool:
    """True when both single-base alleles are the site's REF or ALT as written.

    The input formats state their strand, so a call is never read on the
    other strand: one that fits neither allele is a mismatch.
    """
    if len(allele1) != 1 or len(allele2) != 1:
        return False
    return {allele1, allele2} <= {ref, alt}


def plink_paths(prefix: Path) -> dict[str, Path]:
    """Return the actual PLINK destinations; shared with input-safety preflight.

    The ending is added to the whole prefix: ``release.v1`` gives
    ``release.v1.bed``, never ``release.bed``.
    """
    return {ext: prefix.with_name(prefix.name + ext) for ext in (".fam", ".bim", ".bed")}


def default_plink_prefix(input_path: Path) -> Path:
    """The input path minus one known genotype-file ending (``my.genome.txt`` → ``my.genome``)."""
    name = input_path.name
    for ending in _INPUT_ENDINGS:
        if name.lower().endswith(ending) and len(name) > len(ending):
            return input_path.with_name(name[: -len(ending)])
    return input_path


def coded_against(v: Variant, pair: tuple[str, str]) -> bool:
    """True when the writer codes ``v`` against ``pair``, False when it uses the fallback.

    A call that does not fit the pair as written is written missing against
    it; a call that fits a palindromic pair (A/T, C/G) gets the fallback.
    """
    return not _fits_site(v.allele1, v.allele2, *pair) or not is_strand_ambiguous(*pair)


def plink_writable(v: Variant) -> bool:
    """True for a call the PLINK set holds: not a no-call, both alleles one of ACGT."""
    return not v.is_no_call and v.allele1 in _NUCLEOTIDES and v.allele2 in _NUCLEOTIDES


RecordKey = tuple[str, str, int, str, str]


def plink_record_key(v: Variant) -> RecordKey:
    """A call's own identity for allele coding: rsID, locus and the two called alleles."""
    return (v.rsid, v.chromosome, v.position, v.allele1, v.allele2)


def _generated_id(v: Variant, emitted: set[str]) -> str:
    """``chr:pos:ref:alt`` for a variant without an ID; ``:2``, ``:3``… if already written."""
    base = f"{v.chromosome}:{v.position}:{v.ref or '.'}:{v.input_alt or '.'}"
    candidate, n = base, 1
    while candidate in emitted:
        n += 1
        candidate = f"{base}:{n}"
    return candidate


def export_plink(
    variants: Iterator[Variant],
    prefix: Path,
    build: str,
    ref_alt_map: dict[str, tuple[str, str]] | None = None,
    *,
    record_alleles: Mapping[RecordKey, tuple[str, str]] | None = None,
    supplied_ids: Iterable[str] = (),
) -> PlinkCounts:
    """Write .bed/.bim/.fam from parsed variants.

    Args:
        variants: Parsed variant iterator (consumed once).
        prefix: Base path for output files.
        build: Genome build label (informational, not used for liftover).
        ref_alt_map: ``{rsid: (ref, alt)}`` from gnomAD (GRCh38 exports).
        record_alleles: ``{plink_record_key(v): (ref, alt)}``, chosen for each
            call (GRCh37 exports, from the rsID identity table). When given,
            it is used instead of ``ref_alt_map``.
        supplied_ids: IDs the input supplies anywhere in the stream; a
            generated ID never takes one of them.
            When provided, uses ref/alt to assign A1/A2 for proper allele coding.
            When None or rsid missing, falls back to ``A2="0"`` for homozygotes.
            A call whose alleles fit neither REF nor ALT is written as missing.

    Returns:
        A ``PlinkCounts``.

    Note:
        A variant without an rsID is written with the ID ``chr:pos:ref:alt``
        (``.`` for an unknown part), with ``:2``, ``:3``… when an earlier
        line already has that ID. No-call variants and indels are dropped. PLINK1 BIM is SNV-only
        (A1/A2 in ``ACGT``). Multi-character alleles and single-character
        non-nucleotide codes (array ``I``/``D`` indel codes, a VCF ``*``)
        would produce non-standard BIM rows that downstream tools may
        reject; they are counted in ``indels_skipped``.

        All three files are staged in owned unique temps and closed before
        publication. A prepublication failure preserves every prior output.
        Publication is atomic per file: a later replacement failure can leave
        already-published complete files alongside unchanged prior files.
        Cleanup removes only owned temps; an uncatchable process kill may
        leave those temps. A bare iterator carries no input path to validate.

    Raises:
        EmptyPlinkSetError: no variant could be written (every one was a
            no-call or an indel, or there were none); no files are left.
    """
    from allelix.reports import staged_output

    final = plink_paths(prefix)
    with contextlib.ExitStack() as owned:
        stages = {ext: owned.enter_context(staged_output(path)) for ext, path in final.items()}
        with contextlib.ExitStack() as streams:
            fam = streams.enter_context(stages[".fam"].open_text(encoding=None))
            bim = streams.enter_context(stages[".bim"].open_text(encoding=None))
            bed = streams.enter_context(stages[".bed"].open_binary())
            counts = _write_plink(
                variants, fam, bim, bed, ref_alt_map, record_alleles, supplied_ids
            )
        if counts[0] == 0:
            raise EmptyPlinkSetError(counts[1], counts[2])
        for stage in stages.values():
            stage.publish()
    return counts


def _write_plink(
    variants: Iterator[Variant],
    fam_f: TextIO,
    bim_f: TextIO,
    bed_f: BinaryIO,
    ref_alt_map: dict[str, tuple[str, str]] | None,
    record_alleles: Mapping[RecordKey, tuple[str, str]] | None = None,
    supplied_ids: Iterable[str] = (),
) -> PlinkCounts:
    """The body of ``export_plink``: write through the three owned streams."""
    fam_f.write("0\tSAMPLE\t0\t0\t0\t-9\n")

    written = 0
    skipped = 0
    indels = 0
    monomorphic = 0
    mismatches = 0
    # IDs a generated one must not take: the supplied ones, plus every ID
    # with a colon written so far (generated IDs always have one).
    emitted: set[str] = set(supplied_ids)

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

        if record_alleles is not None:
            pair = record_alleles.get(plink_record_key(v))
        else:
            pair = ref_alt_map.get(v.rsid) if ref_alt_map else None
        if pair is not None:
            ref, alt = pair
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

        variant_id = v.rsid or _generated_id(v, emitted)
        if ":" in variant_id:
            emitted.add(variant_id)
        bim_f.write(f"{chrom_code}\t{variant_id}\t0\t{v.position}\t{a1}\t{a2}\n")
        bed_f.write(bytes([bed_code]))
        written += 1

    return PlinkCounts(written, skipped, indels, monomorphic, mismatches)


def resolve_ref_alt_via_gnomad(
    variant_by_rsid: Mapping[str, Variant], data_dir: Path
) -> dict[str, tuple[str, str]]:
    """Return per-rsid (ref, alt) for PLINK A1/A2 coding via gnomAD (GRCh38 exports).

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


def resolve_ref_alt_via_identity(
    variants: Iterable[Variant], data_dir: Path
) -> dict[RecordKey, tuple[str, str]] | None:
    """Return GRCh37 (ref, alt) for each call from the rsID identity table, or None without it.

    Keyed by ``plink_record_key``, so two records sharing an rsID (split
    multi-allelic rows, or the same rsID at two positions) each get their own.

    Only the table's rows for the call's own rsID at the call's chromosome and
    position count. One single-base pair is used as written; with several,
    the one the call fits as written, if exactly one does. Anything else is
    left out, so the exporter uses its fallback coding there. No strand is
    complemented: GRCh37 alleles are on the GRCh37 forward strand.
    """
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator

    identity = RsidIdentityAnnotator(data_dir)
    with identity:
        if not identity.is_ready():
            return None
        named = [v for v in variants if v.rsid and plink_writable(v)]
        rows = identity.alleles_at_loci("GRCh37", {(v.chromosome, v.position) for v in named})
    result: dict[RecordKey, tuple[str, str]] = {}
    for v in named:
        pairs = {
            (ref, alt)
            for ref, alt, rsid in rows.get((v.chromosome, v.position), [])
            if rsid == v.rsid and ref in _NUCLEOTIDES and alt in _NUCLEOTIDES
        }
        if len(pairs) > 1:
            pairs = {pair for pair in pairs if _fits_site(v.allele1, v.allele2, *pair)}
        if len(pairs) == 1:
            result[plink_record_key(v)] = pairs.pop()
    return result


def resolve_ref_alt_via_gnomad_grch37(
    variants: Iterable[Variant],
    data_dir: Path,
    have: Mapping[RecordKey, tuple[str, str]],
) -> dict[RecordKey, tuple[str, str]]:
    """GRCh37 pairs from gnomAD for the calls the identity table did not code.

    gnomAD holds GRCh38 alleles. Its pair for the call's rsID is used only
    when it is not strand-ambiguous (A/T, C/G) and the call fits it as
    written. A non-palindromic pair shares no allele with its complement, so
    a site written on the other strand on GRCh38 cannot fit and keeps the
    fallback coding instead of being written missing; a palindromic pair
    cannot be told apart that way and is not trusted. Calls in ``have`` are
    left alone. Empty when gnomAD is not ready.
    """
    from allelix.annotators.gnomad import GnomadAnnotator

    wanted = [
        v for v in variants if v.rsid and plink_writable(v) and plink_record_key(v) not in have
    ]
    if not wanted:
        return {}
    gnomad = GnomadAnnotator(data_dir)
    with gnomad:
        if not gnomad.is_ready():
            return {}
        coords = gnomad.bulk_resolve_coordinates({v.rsid for v in wanted})
    result: dict[RecordKey, tuple[str, str]] = {}
    for v in wanted:
        for _, _, ref, alt in coords.get(v.rsid, []):
            if (
                ref in _NUCLEOTIDES
                and alt in _NUCLEOTIDES
                and not is_strand_ambiguous(ref, alt)
                and _fits_site(v.allele1, v.allele2, ref, alt)
            ):
                result[plink_record_key(v)] = (ref, alt)
                break
    return result


def _fallback_coding(v: Variant) -> tuple[str, str, int, bool]:
    """Fallback allele coding when ref/alt is unknown.

    Returns ``(a1, a2, bed_code, is_monomorphic)``.
    """
    if v.is_heterozygous:
        alleles = sorted([v.allele1, v.allele2])
        return alleles[0], alleles[1], 0b10, False

    return v.allele1, "0", 0b00, True

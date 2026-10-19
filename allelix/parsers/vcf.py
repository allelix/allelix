# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""VCF / gVCF parser.

Handles plain VCF 4.x and gVCF (GATK reference-confidence) files, both
``.vcf`` and ``.vcf.gz``. Streams via the standard library — no pysam
dependency for the base parser. Tabix-indexed random access is a
separate path (see ``allelix[vcf-index]`` extras and
``cli/utility.py extract``).

gVCF distinction
----------------

Plain VCF: absence at a position means reference. Lines are variants
only.

gVCF: explicit reference blocks for positions called as reference-
confident. Detected by ``##ALT=<ID=NON_REF,...>`` in the header or
by ``END=`` INFO tags. Reference blocks never yield a ``Variant`` —
they match nothing in any annotation database and would push the
annotation pipeline through millions of zero-hit queries. Positions
the caller registers with ``VcfParser.watch`` (high-value SNPs)
are classified instead: a block with GQ below the floor or DP 0 makes
them no-calls, a confident block makes them reference, and on a gVCF
a position no record covers is uncovered (ADR-0042).

Quality gate (ADR-0042)
-----------------------------

Rows whose FILTER is not ``PASS`` / ``.`` are skipped unless
``keep_filtered``. A call with GQ below ``min_gq`` (default 20) or DP 0
is yielded as a no-call. The spanning-deletion allele ``*`` is
symbolic. ``VcfParser.quality_stats`` records what the gate did.

Multi-sample handling
---------------------

``VcfParser(sample=None)`` works on single-sample files. For
multi-sample VCFs, the parser must be instantiated with
``VcfParser(sample="<sample_id>")`` or ``parse()`` raises. The CLI
threads the user's ``--sample`` flag through to the constructor.
"""

from __future__ import annotations

import logging
import re
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import dataclass, field
from functools import lru_cache
from typing import TYPE_CHECKING, ClassVar, TextIO

from allelix.models import NO_CALL_MARKER, Variant
from allelix.parsers._helpers import DamagedFileError, normalize_chromosome, open_genotype_text
from allelix.parsers.base import GenotypeMetadata, GenotypeParser
from allelix.parsers.genesforgood import METHOD_ARRAY_IMPUTED
from allelix.utils.build_detect import is_ucsc_mito_name

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

# How many data lines to scan for END= when looking for gVCF markers.
# 100 covers any realistic header-then-data-block layout without
# requiring a full file scan.
_GVCF_SNIFF_LIMIT = 100

# Symbolic alleles indicate non-SNV/indel records (structural variants
# in plain VCF, reference blocks in gVCF). All skipped at the parser
# level — Allelix v2.0 annotates only SNVs and small indels.
_SYMBOLIC_ALT_PREFIX = "<"
# VCF 4.3 spanning-deletion allele: this haplotype's base is removed by an
# upstream deletion. It names no sequence, so it is treated as symbolic.
_SPANNING_DELETION = "*"


def _is_breakend(allele: str) -> bool:
    """VCF breakend notation (``G]chr2:321681]``, ``[chr17:198982[A``, ``G.``, ``.A``)."""
    return "[" in allele or "]" in allele or (len(allele) > 1 and "." in (allele[0], allele[-1]))


def _is_symbolic(allele: str) -> bool:
    """True for symbolic ALTs (``<NON_REF>``, ``<DEL>``, ``<*>``), ``*`` and breakends."""
    return (
        allele.startswith(_SYMBOLIC_ALT_PREFIX)
        or allele == _SPANNING_DELETION
        or _is_breakend(allele)
    )


# gVCF reference-block ALTs: symbolic, but not structural variants.
_REFERENCE_BLOCK_ALTS = frozenset({"<NON_REF>", "<*>"})


@lru_cache(maxsize=1024)
def _structural_kind(alt: str) -> str | None:
    """The structural-variant kind of an all-symbolic ALT column, or None.

    ``<DEL>`` / ``<DUP:TANDEM>`` -> ``DEL`` / ``DUP:TANDEM``; breakends ->
    ``BND``. None for ALTs with a real allele, reference-only rows and
    gVCF reference blocks. Cached: ALT strings repeat.
    """
    alts = [] if alt in (".", "") else alt.split(",")
    kinds = []
    for a in alts:
        if a in _REFERENCE_BLOCK_ALTS or a == _SPANNING_DELETION:
            continue
        if a.startswith(_SYMBOLIC_ALT_PREFIX):
            kinds.append(a.strip("<>"))
        elif _is_breakend(a):
            kinds.append("BND")
        else:
            return None
    return kinds[0] if kinds else None


@lru_cache(maxsize=1024)
def _cannot_yield_call(alt: str) -> bool:
    """True when an ALT column can never produce a Variant and has no ``*``.

    Reference-only (``.``) and all-symbolic ALTs (gVCF ``<NON_REF>`` blocks,
    structural variants). Cached: a gVCF repeats the same few ALT strings
    millions of times.
    """
    alts = [] if alt in (".", "") else alt.split(",")
    return all(_is_symbolic(a) for a in alts) and _SPANNING_DELETION not in alts


# match a ``##contig=<ID=chrN,...>`` line declaring any standard
# human chromosome with the ``chr`` prefix. Standard names only — alt
# contigs (``GL00*``, ``hs37d5``, ``NC_*``) don't disambiguate the
# build. The terminator ``[,>]`` keeps us from matching prefixes like
# ``ID=chr1_KI270706v1_random`` (an alt contig) when only chr1 is
# present as a standard contig.
_CONTIG_ID_RE = re.compile(r"[<,]\s*ID=([^,>]+)")
_CONTIG_LENGTH_RE = re.compile(r"[<,]\s*length=(\d+)")


# GATK's convention: a genotype quality below 20 is low confidence.
DEFAULT_MIN_GQ = 20

# FILTER values that mean the caller did not reject the row.
_UNFILTERED = frozenset({"PASS", ".", ""})

# what the file says about a position. A variant call is the
# strongest evidence, then a no-call, then a confident reference call; a
# position covered by several records takes the strongest.
CALLED = "called"
NO_CALL = "no_call"
REFERENCE = "reference"
_STATE_RANK = {REFERENCE: 1, NO_CALL: 2, CALLED: 3}

# Why a position is a no-call.
REASON_CALLER_NO_CALL = "no_call"  # the caller emitted a missing genotype (./.)
REASON_FILTERED = "filtered"  # the caller's FILTER rejected the record
REASON_LOW_GQ_CALL = "low_gq_call"  # a variant call below the GQ floor or with DP 0
REASON_LOW_GQ_BLOCK = "low_gq_block"  # a reference block or call below the floor or DP 0
REASON_SPANNING_DELETION = "spanning_deletion"  # a haplotype is the deletion allele *
REASON_SYMBOLIC = "symbolic_allele"  # the call is a structural / symbolic allele
REASON_NO_COVERAGE = "no_coverage"  # gVCF: no record covers the position

Locus = tuple[str, int]
LocusState = tuple[str, str]  # (state, reason); reason is "" unless NO_CALL


@dataclass
class VcfQualityStats:
    """What the quality gate saw and did during one ``parse()`` run.

    ``filtered_rows``: rows whose FILTER was not ``PASS`` or ``.``.
    ``filtered_by_label``: the same rows counted by each FILTER label
    (a row with ``LowQual;VQSRTranche…`` counts once under each).
    ``kept_filtered``: True when those rows were kept (``keep_filtered``).
    ``low_quality_calls``: calls reported as no-calls because their GQ
    was below ``min_gq`` or their DP was 0.
    ``spanning_deletion_calls`` / ``spanning_deletion_sites``: calls with
    the deletion allele ``*`` on one or both haplotypes, which are not
    assessed, and each one as ``chrom:pos REF>ALT GT`` (ADR-0047).
    ``half_calls``: calls with one allele called and the other missing
    (GT ``./1``). ClinVar and the GWAS Catalog assess the called allele;
    ClinPGx and SNPedia, which need both, do not (ADR-0047).
    ``malformed_rows``: rows with fewer than eight columns or a POS that
    is not a number. ``no_genotype_rows``: rows with no GT for the
    sample (a sites-only file has nothing else). ``unusable_genotype_calls``:
    calls whose GT names an allele the record does not have (also a GT
    such as ``0/1`` on a record whose ALT is ``.``).
    ``symbolic_allele_calls``: calls whose GT names a symbolic ALT (such as
    ``<NON_REF>``) on a record that also has a real ALT. A reference block
    (``<NON_REF>`` or ``<*>`` alone) is not checked: a gVCF holds them by
    the million (ADR-0047).
    A call written as a no-call (``./.``) is yielded as a no-call, as
    before; it is not a dropped call and is not counted here.
    ``trimmed_matches``: calls ClinVar matched only after trimming an
    allele the record wrote untrimmed; set by the pipeline (ADR-0047).
    ``clinvar_ref_mismatches``: one-base records at a ClinVar position
    whose REF is none of ClinVar's one-base REFs there, the complement
    included (mis-coordinated input); set by the pipeline.
    ``reference_carried_calls``: hom-ref calls (a ``0/0`` call or a
    reference block) at a site where the build's reference carries a
    ClinVar variant allele, assessed as homozygous for that allele.
    ``reference_carried_not_in_file``: such sites a plain VCF has no
    record at. Both set by the pipeline.
    ``symbolic_rows`` / ``symbolic_kinds``: structural-variant and
    copy-number records (all ALTs symbolic, such as ``<DEL>``, ``<DUP>``,
    ``<CNV>``, or breakends), counted whatever their FILTER, and their
    kinds. gVCF reference blocks are not counted. None of them is
    annotated; a file of nothing else is the wrong deliverable.
    ``is_gvcf``: the file has reference blocks (``<NON_REF>`` in the
    header or an ``END=`` record), so every position it does not cover
    is genuinely uncovered.
    ``locus_states`` / ``rsid_states``: the ``(state, reason)`` of each
    watched position / rsID (``VcfParser.watch``) that a record covers.
    """

    filtered_rows: int = 0
    filtered_by_label: Counter[str] = field(default_factory=Counter)
    kept_filtered: bool = False
    low_quality_calls: int = 0
    min_gq: int = DEFAULT_MIN_GQ
    spanning_deletion_calls: int = 0
    spanning_deletion_sites: list[str] = field(default_factory=list)
    half_calls: int = 0
    malformed_rows: int = 0
    no_genotype_rows: int = 0
    unusable_genotype_calls: int = 0
    symbolic_allele_calls: int = 0
    trimmed_matches: int = 0
    clinvar_ref_mismatches: int = 0
    reference_carried_calls: int = 0
    reference_carried_not_in_file: int = 0
    symbolic_rows: int = 0
    symbolic_kinds: Counter[str] = field(default_factory=Counter)
    is_gvcf: bool = False
    locus_states: dict[Locus, LocusState] = field(default_factory=dict)
    rsid_states: dict[str, LocusState] = field(default_factory=dict)


def strongest_state(current: LocusState | None, new: LocusState | None) -> LocusState | None:
    """The stronger of two position states (called > no-call > reference)."""
    if new is None:
        return current
    if current is None or _STATE_RANK[new[0]] > _STATE_RANK[current[0]]:
        return new
    return current


class MultiSampleError(ValueError):
    """Raised when a multi-sample VCF is parsed without a sample selection."""


class SampleNotFoundError(ValueError):
    """Raised when ``--sample <ID>`` does not match any column in the VCF."""


def format_sample_list(samples: list[str], max_shown: int = 10) -> str:
    """Render a sample-ID list for an error message, truncating beyond max_shown.

    1000 Genomes multi-sample VCFs carry 3,202 sample columns. Dumping
    every ID into a SampleNotFoundError / MultiSampleError message
    floods the terminal with screens of useless output before the user
    can scroll up to see the actual error. Show the first ``max_shown``
    IDs followed by an "...and N more" tail.
    """
    if len(samples) <= max_shown:
        return ", ".join(samples)
    return ", ".join(samples[:max_shown]) + f", ... and {len(samples) - max_shown} more"


# how many leading non-blank lines ``can_parse`` will scan for
# ``##fileformat=VCF`` before giving up. 10 covers every real-world
# preamble shape we've encountered (BOM-on-its-own-line edge cases,
# tool banners, blank-line prefixes) without making the negative-case
# cost meaningful — non-VCF inputs reject after reading ~10 short
# lines (a handful of KB).
_PREAMBLE_SCAN_LINES = 10


class VcfParser(GenotypeParser):
    """Streaming parser for VCF 4.x and gVCF files.

    Stateless except for the optional ``sample`` constructor argument
    (the chosen sample column for multi-sample VCFs). Reusable across
    files; the same instance can ``parse`` and ``get_metadata`` any
    number of files.

    **Registry singleton vs. constructor state.** Every other parser
    in ``PARSERS`` is fully stateless and used as a singleton via
    ``get_parser_by_name`` / ``detect_parser``. ``VcfParser`` breaks
    that pattern: the registered instance has ``sample=None`` and
    will raise ``MultiSampleError`` on multi-sample VCFs. The CLI
    works around this in ``cli/_helpers._resolve_parser`` by
    constructing a new ``VcfParser(sample=...)`` instance when the
    user's ``--sample`` value applies.

    Callers reaching for the registry directly (``get_parser_by_name
    ("vcf")``) get the sample-less singleton — fine for single-sample
    VCFs, fails loudly for multi-sample. Programmatic callers that
    handle multi-sample VCFs should construct ``VcfParser(sample=...)``
    explicitly rather than going through the registry.
    """

    name: ClassVar[str] = "vcf"
    display_name: ClassVar[str] = "VCF / gVCF"
    file_extensions: ClassVar[list[str]] = [".vcf", ".vcf.gz"]
    url: ClassVar[str] = "https://samtools.github.io/hts-specs/"

    def __init__(
        self,
        sample: str | None = None,
        *,
        keep_filtered: bool = False,
        min_gq: int = DEFAULT_MIN_GQ,
    ) -> None:
        """``sample``: which sample column to read in a multi-sample VCF.

        ``None`` is fine for single-sample files (the lone sample is
        used automatically). For multi-sample files, ``parse()`` raises
        ``MultiSampleError`` when ``sample`` is None. Use
        ``list_samples(file_path)`` to discover available IDs before
        constructing a sample-bound parser.

        Quality gate: rows the caller filtered (FILTER not ``PASS`` or
        ``.``) are skipped unless ``keep_filtered``. A call whose GQ is
        below ``min_gq`` (or whose DP is 0) is yielded as a no-call;
        ``min_gq=0`` turns the GQ floor off. ``quality_stats`` records
        what the gate did during the most recent ``parse()``.
        """
        self._sample = sample
        self._keep_filtered = keep_filtered
        self._min_gq = min_gq
        self.quality_stats = VcfQualityStats(kept_filtered=keep_filtered, min_gq=min_gq)
        self._watch_loci: dict[str, list[int]] = {}
        self._watch_rsids: frozenset[str] = frozenset()
        # Set by the build decision (GRCh37): records named as UCSC hg19
        # names the mitochondrion (chrM, M) are left out and counted in
        # ``ucsc_mito_dropped``, unless that name is declared with the rCRS
        # length (``rcrs_mito_names``).
        self.drop_ucsc_mito = False
        self.rcrs_mito_names: frozenset[str] = frozenset()
        self.ucsc_mito_dropped = 0

    def watch(self, loci: set[Locus], rsids: set[str]) -> None:
        """Record what the file says at these positions and rsIDs.

        ``loci`` are ``(chromosome, position)`` in Allelix's bare form; a
        record covers a position from POS through its ``END`` (reference
        blocks) or the end of its REF. Results land in
        ``quality_stats.locus_states`` / ``rsid_states`` after ``parse()``.
        Pass empty sets to stop watching.
        """
        by_chrom: dict[str, list[int]] = {}
        for chrom, pos in loci:
            by_chrom.setdefault(chrom, []).append(pos)
        self._watch_loci = {chrom: sorted(set(p)) for chrom, p in by_chrom.items()}
        self._watch_rsids = frozenset(rsids)

    # ── Public API ──────────────────────────────────────────────

    def can_parse(self, file_path: Path) -> bool:
        """True if a leading non-blank line is ``##fileformat=VCF...``.

        catches ``UnicodeDecodeError`` and ``DamagedFileError`` (corrupt
        gzip / truncated BGZF) on top of ``OSError`` so format sniffing
        never raises through into ``detect_parser``. A real-world
        ``.vcf.gz`` shouldn't reach the decode path's replacement
        substitution (``_open_vcf`` uses ``errors="replace"``), but
        binary garbage or a truncated tabix index masquerading as a
        ``.vcf.gz`` would otherwise still leak a traceback at the
        ``read`` boundary of the gzip reader.

        scan multiple leading non-blank lines (capped at
        :data:`_PREAMBLE_SCAN_LINES`) instead of insisting
        ``##fileformat=VCF`` appears on line 1. Real-world VCF
        producers occasionally prepend a non-standard line — a
        tool-banner comment, a blank line wrapped by a downstream
        post-processor, an editor's "edited by …" annotation — and
        rejecting the file on first-line-mismatch read as broken
        when the actual format was fine. The scan terminates the
        moment a line starts with ``##fileformat=VCF`` (success) or
        we've consumed the cap without finding it (failure). A line
        starting with ``#CHROM`` (the data-header line) also
        terminates the scan as failure — we're past the metadata
        block and ``##fileformat`` clearly never appeared.

        The UTF-8 BOM that prompted the original bug report is
        stripped structurally by ``_open_vcf``'s ``utf-8-sig`` codec
        (also ), so the BOM never reaches this loop. The
        multi-line scan complements that by covering the other class
        of preamble quirks the BOM hypothesis doesn't explain.
        """
        try:
            with _open_vcf(file_path) as handle:
                scanned = 0
                for line in handle:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    if stripped.startswith("##fileformat=VCF"):
                        return True
                    # ``#CHROM`` is the data-header line; if we hit it
                    # before ``##fileformat`` we know the metadata
                    # block doesn't contain the declaration and the
                    # file isn't a VCF we can claim.
                    if stripped.startswith("#CHROM"):
                        return False
                    scanned += 1
                    if scanned >= _PREAMBLE_SCAN_LINES:
                        return False
        except (OSError, UnicodeDecodeError, DamagedFileError):
            return False
        return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Yield ``Variant`` objects from the VCF, one per non-reference site.

        Reference blocks (gVCF), reference-only lines (plain VCF with
        ALT='.'), symbolic-allele records (``<NON_REF>``, ``<DEL>``,
        ``<CNV>``, ``<*>``), and malformed rows are silently skipped.
        Malformed rows are also logged at WARNING level (counted by the
        CLI's parser-logging wire-up). Build is **not** filled in here;
        the pipeline overrides every variant's build with the detected
        effective build (ADR-0021).
        """
        stats = VcfQualityStats(kept_filtered=self._keep_filtered, min_gq=self._min_gq)
        self.quality_stats = stats
        drop_mito = self.drop_ucsc_mito
        keep_mito = self.rcrs_mito_names
        self.ucsc_mito_dropped = 0
        watch_rsids = self._watch_rsids
        watcher = _LocusWatcher(self._watch_loci) if self._watch_loci else None
        # Re-classifying a watched row must not count it twice.
        scratch = VcfQualityStats(kept_filtered=self._keep_filtered, min_gq=self._min_gq)
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
            stats.is_gvcf = header.has_non_ref_alt
            sample_idx = self._resolve_sample_index(header, file_path)

            def record(covering: list[tuple[list[str], list[Locus]]]) -> None:
                for row, loci in covering:
                    state = _parse_record(row, "", sample_idx, scratch, want_state=True)[1]
                    if state is None:
                        continue
                    for locus in loci:
                        merged = strongest_state(stats.locus_states.get(locus), state)
                        stats.locus_states[locus] = merged or state

            for line in handle:
                stripped = line.rstrip("\n").rstrip("\r")
                if not stripped or stripped.startswith("#"):
                    continue
                cols = stripped.split("\t")
                if not stats.is_gvcf and _info_end(cols) is not None:
                    stats.is_gvcf = True
                rsid = ""
                if watch_rsids and len(cols) > 2 and cols[2] != ".":
                    rsid = _canonical_rsid(cols[2])
                    rsid = rsid if rsid in watch_rsids else ""
                variant, state = _parse_record(
                    cols, stripped, sample_idx, stats, want_state=bool(rsid)
                )
                if rsid and state is not None:
                    merged = strongest_state(stats.rsid_states.get(rsid), state)
                    stats.rsid_states[rsid] = merged or state
                if watcher is not None:
                    covering = watcher.feed(cols)
                    if covering:
                        record(covering)
                if variant is None:
                    continue
                if drop_mito and is_ucsc_mito_name(cols[0]) and cols[0] not in keep_mito:
                    self.ucsc_mito_dropped += 1
                    continue
                yield variant
            if watcher is not None:
                record(watcher.finish())

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Return the canonical header-derivable metadata.

        VCF-specific extras (``vcf_type``, ``samples``) are exposed via
        ``vcf_type()`` and ``list_samples()`` respectively, since
        ``GenotypeMetadata`` is a fixed TypedDict shared by all parsers.

        For multi-sample files, ``sample_id`` is the bound sample (if
        any). For single-sample files it is the only sample's name.
        Returns ``""`` if the bound sample isn't in the file or the
        file is multi-sample with no binding — let the caller decide
        how to handle it; this method never raises.
        """
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
            samples = header.samples
            if self._sample is not None and self._sample in samples:
                sample_id = self._sample
            elif len(samples) == 1:
                sample_id = samples[0]
            else:
                sample_id = ""
            imputed = (
                header.genes_for_good
                and bool(sample_id)
                and _all_phased(handle, samples.index(sample_id))
            )
        meta = GenotypeMetadata(
            format=self.name,
            sample_id=sample_id,
            build=header.build,
            contig_lengths=dict(header.contig_lengths),
        )
        if imputed:
            # Genes for Good's imputed VCF: marked as its text file is.
            meta["method"] = METHOD_ARRAY_IMPUTED
        elif header.imputation_fields:
            meta["method"] = METHOD_VCF_IMPUTED
        return meta

    def validate_sample(self, file_path: Path) -> None:
        """Raise on an unusable sample binding without consuming variant data.

        Runs the same :class:`MultiSampleError` / :class:`SampleNotFoundError`
        checks as :meth:`parse` but stops after reading the header. The CLI's
        pre-flight calls this so the error surfaces before the analysis
        pipeline starts streaming.

        Single-sample files (sites-only files included) never raise.
        """
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
        self._resolve_sample_index(header, file_path)

    def list_samples(self, file_path: Path) -> list[str]:
        """Return the sample IDs from the ``#CHROM`` column header line.

        Empty list if the file has no column header line (an unusual
        but possible corner of the VCF spec). Single-element list for
        single-sample files. Multi-element list for multi-sample.
        """
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
        return list(header.samples)

    def vcf_type(self, file_path: Path) -> str:
        """``"gvcf"`` if reference-block markers are present, else ``"plain"``.

        Detection: ``##ALT=<ID=NON_REF,...>`` in the header is a
        strong gVCF signal. As a fallback, the first
        :data:`_GVCF_SNIFF_LIMIT` data lines are scanned for an
        ``END=`` INFO tag (which gVCF uses to denote reference-block
        end positions). Either signal flips the result to ``"gvcf"``.
        """
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
            if header.has_non_ref_alt:
                return "gvcf"
            for i, line in enumerate(handle):
                if i >= _GVCF_SNIFF_LIMIT:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                # INFO is column index 7 (0-based). Cheap presence check.
                if "END=" in line:
                    return "gvcf"
        return "plain"

    # ── Internals ───────────────────────────────────────────────

    def _resolve_sample_index(self, header: _VcfHeader, file_path: Path) -> int:
        """Pick which sample column to read; raise on ambiguity or mismatch."""
        if not header.samples:
            # No samples — sites-only VCF. Yield no genotype calls.
            # Treat as no sample selected; data-line parsing returns None.
            return -1
        if self._sample is not None:
            try:
                return header.samples.index(self._sample)
            except ValueError as exc:
                msg = (
                    f"Sample {self._sample!r} not found in {file_path.name}. "
                    f"Available samples: {format_sample_list(header.samples)}"
                )
                raise SampleNotFoundError(msg) from exc
        if len(header.samples) > 1:
            msg = (
                f"Multi-sample VCF {file_path.name}: pass --sample <ID> to "
                f"select. Available samples: {format_sample_list(header.samples)}"
            )
            raise MultiSampleError(msg)
        return 0


# ── Header parsing ─────────────────────────────────────────────


class _VcfHeader:
    """Parsed VCF header — what the pipeline needs from the ``##`` lines."""

    __slots__ = (
        "build",
        "contig_lengths",
        "genes_for_good",
        "has_non_ref_alt",
        "imputation_fields",
        "samples",
    )

    def __init__(self) -> None:
        self.samples: list[str] = []
        self.build: str | None = None
        self.has_non_ref_alt: bool = False
        # ``##contig`` ID -> declared length. The build decision reads the
        # lengths (ADR-0044); contig naming (``chr1`` vs ``1``) is not
        # evidence of a build.
        self.contig_lengths: dict[str, int] = {}
        # A ``##Genes_for_Good_version`` line: Genes for Good's own VCF.
        self.genes_for_good: bool = False
        # The imputation fields the header declares (``IMPUTATION_INFO_IDS``).
        self.imputation_fields: set[str] = set()


def _read_header(handle: TextIO) -> _VcfHeader:
    """Consume lines from ``handle`` up to and including the ``#CHROM`` row.

    Returns parsed metadata. After this call the handle is positioned
    at the first data line (or EOF).

    tolerates a small preamble of single-``#`` comment lines
    BEFORE the first ``##fileformat`` declaration, symmetric with
    :meth:`VcfParser.can_parse`'s scan. Some upstream tools and
    post-processors prepend a banner (e.g. ``# Generated by …``)
    that's not a VCF metadata line but is also not the data block.
    The scan caps at :data:`_PREAMBLE_SCAN_LINES` non-blank lines —
    if no ``##`` declaration appears within that window we treat the
    first non-comment / non-blank line as the start of data,
    preserving the original behavior on inputs that genuinely lack
    a header.
    """
    header = _VcfHeader()
    preamble_scanned = 0
    saw_meta = False
    for line in handle:
        stripped = line.rstrip("\n").rstrip("\r")
        if not stripped:
            continue
        if stripped.startswith("##"):
            saw_meta = True
            _absorb_meta_line(stripped, header)
            continue
        if stripped.startswith("#CHROM"):
            cols = stripped.split("\t")
            # Standard layout: 9 fixed columns then sample columns.
            # FORMAT (col 8) is only present when samples exist; allow
            # sites-only files with 8 cols and no FORMAT.
            if len(cols) > 9:
                header.samples = cols[9:]
            return header
        # single-``#`` comment lines before any ``##`` metadata
        # are treated as preamble (banner / post-processor annotation)
        # up to the scan cap. Once we've seen at least one ``##`` line
        # the header block is "real" and a stray non-``##`` line is
        # treated as data per the original behavior.
        if not saw_meta and stripped.startswith("#"):
            preamble_scanned += 1
            if preamble_scanned < _PREAMBLE_SCAN_LINES:
                continue
        # Reaching a data line before #CHROM means no column header.
        # Return what we have; samples is empty.
        return header
    return header


# Genes for Good writes this line in the header of its VCF downloads. Its
# imputed VCF has no imputation field, but every call in it is phased; its
# measured calls are not (the text download says which file is which).
GENES_FOR_GOOD_LINE = "##Genes_for_Good_version"
# The fields imputation tools declare in the VCFs they write: the ``##INFO``
# fields of Minimac (IMPUTED, TYPED, R2) and Beagle (DR2, IMP), and the
# ``##FORMAT`` dosage ``DS`` (Minimac, Beagle, GLIMPSE) and GLIMPSE's ``HS``.
# A header that declares one marks the file imputed; the values are never
# read. FORMAT ``GP`` alone does not count (a variant caller can add it), nor
# does GATK's INFO ``DS`` ("downsampled"), which is another field.
IMPUTATION_INFO_IDS = frozenset({"IMPUTED", "TYPED", "R2", "DR2", "IMP"})
IMPUTATION_FORMAT_IDS = frozenset({"DS", "HS"})
# The ``method`` of a VCF whose header declares one of them.
METHOD_VCF_IMPUTED = "vcf-imputed"
# How many records with a called GT decide whether a file's calls are phased.
_PHASE_CHECK_CALLS = 1000


def _all_phased(handle: TextIO, sample_idx: int) -> bool:
    """True when the first called diploid GTs of the sample are all phased (``|``).

    Reads from the first data line on, up to ``_PHASE_CHECK_CALLS`` such
    calls. False when any is unphased (``/``) or none is found. A haploid
    GT and a GT with a missing allele say nothing and are passed over.
    """
    phased = 0
    for line in handle:
        cols = line.rstrip("\n").rstrip("\r").split("\t")
        if len(cols) <= 9 + sample_idx or cols[0].startswith("#"):
            continue
        keys = cols[8].split(":")
        if not keys or keys[0] != "GT":
            continue
        gt = cols[9 + sample_idx].split(":")[0]
        if "." in re.split(r"[/|]", gt) or ("/" not in gt and "|" not in gt):
            continue
        if "/" in gt:
            return False
        phased += 1
        if phased >= _PHASE_CHECK_CALLS:
            break
    return phased > 0


def _absorb_meta_line(line: str, header: _VcfHeader) -> None:
    """Extract contig lengths, the assembly label and the gVCF marker from a ``##`` line."""
    if line.startswith("##ALT=") and "ID=NON_REF" in line:
        header.has_non_ref_alt = True
        return
    if line.startswith(GENES_FOR_GOOD_LINE):
        header.genes_for_good = True
        return
    for prefix, ids in (("##INFO=", IMPUTATION_INFO_IDS), ("##FORMAT=", IMPUTATION_FORMAT_IDS)):
        if line.startswith(prefix):
            field_id = _CONTIG_ID_RE.search(line)
            if field_id and field_id.group(1).strip() in ids:
                header.imputation_fields.add(prefix + field_id.group(1).strip())
            return
    if line.startswith("##contig="):
        contig_id = _CONTIG_ID_RE.search(line)
        length = _CONTIG_LENGTH_RE.search(line)
        if contig_id and length:
            header.contig_lengths[contig_id.group(1).strip()] = int(length.group(1))
        if "assembly=" in line and header.build is None:
            # First explicit assembly wins.
            header.build = _extract_assembly(line)


def _extract_assembly(contig_line: str) -> str | None:
    """Pull assembly identifier from a ``##contig`` line, normalized."""
    marker = "assembly="
    start = contig_line.find(marker)
    if start == -1:
        return None
    start += len(marker)
    # assembly value runs until ',', '>', or end-of-line
    end = start
    while end < len(contig_line) and contig_line[end] not in (",", ">"):
        end += 1
    raw = contig_line[start:end].strip()
    upper = raw.upper().replace("HG19", "GRCH37").replace("HG38", "GRCH38")
    if upper.startswith("GRCH37"):
        return "GRCh37"
    if upper.startswith("GRCH38"):
        return "GRCh38"
    return None


# ── Data-line parsing ─────────────────────────────────────────


def _parse_data_line(
    line: str, sample_idx: int, stats: VcfQualityStats | None = None
) -> Variant | None:
    """Parse one tab-separated data line into a ``Variant`` or ``None``.

    See ``_parse_record``; this wrapper drops the position state.
    """
    if stats is None:
        stats = VcfQualityStats()
    variant, _state = _parse_record(line.split("\t"), line, sample_idx, stats)
    return variant


def _parse_record(
    cols: list[str],
    line: str,
    sample_idx: int,
    stats: VcfQualityStats,
    *,
    want_state: bool = False,
) -> tuple[Variant | None, LocusState | None]:
    """Parse one data row: the ``Variant`` it yields and what it says about its span.

    The Variant is ``None`` for: malformed rows, rows the caller filtered
    (FILTER not ``PASS`` / ``.``, unless ``stats.kept_filtered``),
    reference-only rows (ALT='.'), symbolic-allele rows (``<NON_REF>``,
    ``<DEL>``, ``<*>``, ``*``, etc.), rows without a genotype, and
    homozygous-reference calls (the user carries only REF — no actionable
    variant). A call below the GQ floor, or with DP 0, is returned as a
    no-call. QUAL is not consulted.

    The state is ``(CALLED | REFERENCE | NO_CALL, reason)``, or None when
    the row says nothing about the sample (malformed, no genotype). A
    reference block or hom-ref call below the GQ floor (or DP 0) is a
    no-call, by the same rule as a variant call. Reference
    blocks and other rows that cannot yield a call are only classified
    when ``want_state`` (they cover a watched position): a gVCF has
    millions of them.
    """
    if len(cols) < 8:
        logger.warning("Malformed VCF line (fewer than 8 columns): %s", line[:80])
        stats.malformed_rows += 1
        return None, None
    chrom = cols[0]
    pos_str = cols[1]
    rsid = cols[2]
    ref = cols[3]
    alt = cols[4]
    filt = cols[6].strip()

    kind = _structural_kind(alt)
    if kind is not None:
        stats.symbolic_rows += 1
        stats.symbolic_kinds[kind] += 1

    if filt not in _UNFILTERED:
        stats.filtered_rows += 1
        stats.filtered_by_label.update(label for label in filt.split(";") if label)
        if not stats.kept_filtered:
            return None, (NO_CALL, REASON_FILTERED)

    try:
        position = int(pos_str)
    except ValueError:
        logger.warning("Invalid VCF position %r", pos_str)
        stats.malformed_rows += 1
        return None, None

    if not want_state and _cannot_yield_call(alt):
        # Reference block, reference-only row or SV: no Variant, and no
        # one is asking what it says. A reference-only row (ALT '.')
        # whose GT names an ALT is counted (ADR-0047).
        if alt in (".", ""):
            _count_alt_less_call(_gt_indices(_sample_fields(cols, sample_idx)), stats)
        return None, None
    alts = [] if alt in (".", "") else alt.split(",")

    fields = _sample_fields(cols, sample_idx)
    gt = _gt_indices(fields)
    if gt is None:
        # Sites-only file, no sample columns, or no GT in FORMAT; or a GT
        # that is neither haploid nor diploid.
        if "GT" in fields:
            stats.unusable_genotype_calls += 1
        else:
            stats.no_genotype_rows += 1
        return None, None
    low = _low_confidence(fields, stats.min_gq)

    if all(i == "0" for i in gt):
        # Homozygous reference: a reference block (gVCF), a reference-only
        # record, or a hom-ref call at a variant site. Nothing to annotate.
        return None, ((NO_CALL, REASON_LOW_GQ_BLOCK) if low else (REFERENCE, ""))

    stars = sum(1 for i in gt if _points_at(i, alts, _SPANNING_DELETION))
    if stars:
        # a haplotype is the deletion allele '*'. No Variant shape says
        # "this allele, other haplotype deleted", so it is not assessed;
        # counted and listed so the report can say so (ADR-0047).
        stats.spanning_deletion_calls += 1
        stats.spanning_deletion_sites.append(f"{chrom}:{position} {ref}>{alt} {fields['GT']}")
    if any(_points_at(i, alts) for i in gt):
        if not stars and not all(_is_symbolic(a) for a in alts):
            # ADR-0047: a call on a symbolic ALT (0/2 on A,<NON_REF>) of a
            # record that also has a real ALT. Structural-variant records
            # are counted as symbolic rows above; a reference block
            # (<NON_REF> or <*> alone) is not checked (ADR-0047).
            stats.symbolic_allele_calls += 1
        reason = REASON_SPANNING_DELETION if stars else REASON_SYMBOLIC
        return None, (NO_CALL, reason)
    if not alts or all(_is_symbolic(a) for a in alts):
        # Reference-only or symbolic-only record (gVCF block, SV) with a
        # genotype that is not hom-ref. Nothing to annotate.
        if all(i == "." for i in gt):
            return None, (NO_CALL, REASON_CALLER_NO_CALL)
        if not alts:
            _count_alt_less_call(gt, stats)
        return None, None

    a1, a2 = _resolve_genotype(fields, ref, alts)
    if a1 is None or a2 is None:
        # a GT index that names no allele of this record.
        stats.unusable_genotype_calls += 1
        return None, None
    state: LocusState = (CALLED, "")
    half = (a1 == NO_CALL_MARKER) != (a2 == NO_CALL_MARKER)
    if NO_CALL_MARKER in (a1, a2) and not half:
        state = (NO_CALL, REASON_CALLER_NO_CALL)
    elif low:
        # the caller's own genotype quality says this call is not
        # trustworthy. Report it as a no-call, not as the genotype.
        stats.low_quality_calls += 1
        a1 = a2 = NO_CALL_MARKER
        state = (NO_CALL, REASON_LOW_GQ_CALL)
    elif half:
        # ADR-0047: one allele called, the other missing (./1). Yielded
        # with the called allele; the position stays a no-call.
        stats.half_calls += 1
        state = (NO_CALL, REASON_CALLER_NO_CALL)

    # A one-index GT (``1`` on a male X, Y or MT): the file declares one
    # copy. ``_resolve_genotype`` doubles it for matching (ADR-0050).
    gt_field = fields.get("GT", "")
    haploid = "/" not in gt_field and "|" not in gt_field and NO_CALL_MARKER not in (a1, a2)
    variant = Variant(
        rsid=_canonical_rsid(rsid),
        chromosome=normalize_chromosome(chrom),
        position=position,
        allele1=a1,
        allele2=a2,
        # build is overridden by the pipeline (ADR-0021); leave default.
        ref=ref,
        ref_from_input=True,
        haploid=haploid,
    )
    return variant, state


def _count_alt_less_call(gt: list[str] | None, stats: VcfQualityStats) -> None:
    """Count a call on a reference-only record (ALT '.') whose GT names an ALT."""
    if gt is not None and any(i not in ("0", ".") for i in gt):
        stats.unusable_genotype_calls += 1


def _gt_indices(fields: dict[str, str]) -> list[str] | None:
    """The GT field's allele indices (haploid calls doubled), or None."""
    gt = fields.get("GT")
    if gt is None:
        return None
    parts = gt.split("|") if "|" in gt else gt.split("/")
    if len(parts) == 1:
        parts = [parts[0], parts[0]]
    return parts if len(parts) == 2 else None


def _points_at(gt_index: str, alts: list[str], allele: str | None = None) -> bool:
    """True when a GT index names ``allele`` among the ALTs, or any symbolic ALT."""
    if not gt_index.isdigit() or not 1 <= int(gt_index) <= len(alts):
        return False
    named = alts[int(gt_index) - 1]
    return named == allele if allele is not None else _is_symbolic(named)


def _info_end(cols: list[str]) -> int | None:
    """The INFO ``END`` of a record (gVCF reference blocks), or None.

    Located with ``str.find`` rather than splitting INFO: it runs on every
    reference block of a gVCF while high-value positions are watched.
    """
    if len(cols) < 8:
        return None
    info = cols[7]
    at = info.find("END=")
    while at > 0 and info[at - 1] != ";":
        # A key that merely ends in "END" (e.g. "BLEND="); look further on.
        at = info.find("END=", at + 4)
    if at < 0:
        return None
    stop = info.find(";", at)
    try:
        return int(info[at + 4 : stop if stop >= 0 else len(info)])
    except ValueError:
        return None


def _row_stop(cols: list[str], start: int) -> int:
    """Last position a record covers: its INFO ``END``, else the end of its REF."""
    end = _info_end(cols)
    return end if end is not None else start + max(len(cols[3]), 1) - 1


_NOTHING: list[tuple[list[str], list[Locus]]] = []  # shared; never mutated


class _LocusWatcher:
    """Find the records covering watched positions at low per-row cost.

    In a position-sorted file the records covering a watched position P all
    start at or before P. The watcher holds the rows that can still cover
    the current position: the rows sharing the latest start, plus earlier
    variant records whose REF still reaches it (a deletion overlapping the
    records after it). gVCF reference blocks tile the genome and do not
    overlap later records, so an older block is let go once a new start
    arrives. Only when the stream moves past a watched position are the
    held rows' spans (``END`` or REF length) checked against it, so most
    rows cost one ``int()`` and a few comparisons. A file that is not
    sorted (a position goes backwards, or a chromosome comes back) switches
    to checking every row's span against every watched position.
    """

    def __init__(self, loci: dict[str, list[int]]) -> None:
        self._loci = loci
        self._names: dict[str, str] = {}
        self._seen: set[str] = set()
        self._raw_chrom: str | None = None
        self._chrom = ""
        self._positions: list[int] = []
        self._idx = 0
        self._held: list[tuple[list[str], int]] = []
        self._last_start = 0
        self._sorted = True

    def feed(self, cols: list[str]) -> list[tuple[list[str], list[Locus]]]:
        """Rows (this one or earlier) that cover watched positions, with those positions."""
        if len(cols) < 8:
            return _NOTHING
        out = _NOTHING
        if cols[0] != self._raw_chrom:
            out = self.finish()
            self._enter(cols[0])
        positions = self._positions
        if not positions:
            return out
        try:
            start = int(cols[1])
        except ValueError:
            return out
        if (
            self._sorted
            and start > self._last_start
            and (self._idx >= len(positions) or positions[self._idx] >= start)
        ):
            # The common row: nothing watched lies behind it. Hold it.
            self._hold_new_start(cols, start)
            return out
        out = list(out)
        if self._sorted and start < self._last_start:
            self._sorted = False
            out += self._advance(None)
        if not self._sorted:
            stop = _row_stop(cols, start)
            lo, hi = bisect_left(positions, start), bisect_right(positions, stop)
            if lo < hi:
                out.append((cols, [(self._chrom, pos) for pos in positions[lo:hi]]))
            return out
        if start > self._last_start:
            out += self._advance(start)
            self._hold_new_start(cols, start)
        else:
            self._held.append((cols, start))
        return out

    def finish(self) -> list[tuple[list[str], list[Locus]]]:
        """Resolve the watched positions left on the current chromosome."""
        if not self._sorted or self._idx >= len(self._positions):
            return _NOTHING
        return self._advance(None)

    def _hold_new_start(self, cols: list[str], start: int) -> None:
        """Let go of held rows that cannot reach ``start``; hold ``cols``."""
        held = self._held
        if len(held) == 1 and "END=" in held[0][0][7]:
            # The common case: the previous row was a reference block.
            held = []
        elif held:
            # Older variant records stay while their REF reaches the new
            # start; older reference blocks (END=) never overlap it.
            held = [h for h in held if h[1] + len(h[0][3]) > start and "END=" not in h[0][7]]
        held.append((cols, start))
        self._held = held
        self._last_start = start

    def _enter(self, raw_chrom: str) -> None:
        chrom = self._names.get(raw_chrom)
        if chrom is None:
            chrom = self._names[raw_chrom] = normalize_chromosome(raw_chrom)
        if chrom in self._seen:
            self._sorted = False
        self._seen.add(chrom)
        self._raw_chrom = raw_chrom
        self._chrom = chrom
        self._positions = self._loci.get(chrom, [])
        self._idx = 0
        self._held = []
        self._last_start = 0

    def _advance(self, upto: int | None) -> list[tuple[list[str], list[Locus]]]:
        """Match watched positions before ``upto`` (all, if None) to the held rows."""
        covering: dict[int, tuple[list[str], list[Locus]]] = {}
        positions = self._positions
        while self._idx < len(positions) and (upto is None or positions[self._idx] < upto):
            pos = positions[self._idx]
            self._idx += 1
            for row, row_start in self._held:
                if row_start <= pos <= _row_stop(row, row_start):
                    covering.setdefault(id(row), (row, []))[1].append((self._chrom, pos))
        return list(covering.values())


def _canonical_rsid(raw: str) -> str:
    """Pick a usable rsID from the ID column; first ``rs``-prefixed if present.

    VCF ID column may be ``.`` (no ID), a single ID, or a list separated
    by ``;`` (multiple cross-references). Returns the first ID starting
    with ``rs`` (case-insensitive), normalized to lowercase ``rs`` prefix
    so downstream SQL lookups (which all assume the lowercase convention)
    match correctly. Non-rs IDs (e.g., COSMIC ``COSV12345``) pass through
    with original case preserved — their own case conventions apply.

    Positional synthetic IDs (1000 Genomes pipelines emit
    ``"22:10519265:CA:C"``) are not real identifiers — return ``""``
    so the variant flows through rsID resolution downstream rather than
    carrying a meaningless string in its rsid field.

    Empty rsID is a real thing — variants not catalogued in dbSNP have
    no rsID — and downstream annotators simply find no matches, which is
    correct.
    """
    if raw == "." or not raw:
        return ""
    parts = raw.split(";")
    for part in parts:
        if part.lower().startswith("rs"):
            return part.lower()
    first = parts[0]
    # Colons are reserved as field separators in VCF/BCF spec, so they
    # only appear in non-identifier strings — positional synthetics,
    # internal pipeline tags. External database IDs (COSMIC, RCV, etc.)
    # don't contain them.
    if ":" in first:
        return ""
    return first


def _sample_fields(cols: list[str], sample_idx: int) -> dict[str, str]:
    """The chosen sample's FORMAT values by key (``{"GT": "0/1", "GQ": "99"}``).

    Empty for sites-only files, rows without sample columns, or an
    out-of-range sample index. Keys past the end of a truncated sample
    field are absent.
    """
    if sample_idx < 0 or len(cols) < 10 or sample_idx + 9 >= len(cols):
        return {}
    keys = cols[8].split(":")
    values = cols[9 + sample_idx].split(":")
    return dict(zip(keys, values, strict=False))


def _int_field(fields: dict[str, str], key: str) -> int | None:
    """An integer FORMAT value, or None when absent, ``.`` or not a number."""
    raw = fields.get(key, ".")
    try:
        return int(float(raw))
    except ValueError:
        return None


def _low_confidence(fields: dict[str, str], min_gq: int) -> bool:
    """True when GQ is below ``min_gq`` or DP is 0 (the call has no reads)."""
    gq = _int_field(fields, "GQ")
    return (gq is not None and gq < min_gq) or _int_field(fields, "DP") == 0


def _resolve_genotype(
    fields: dict[str, str], ref: str, alts: list[str]
) -> tuple[str | None, str | None]:
    """Pick the user's two alleles from the GT field of the chosen sample.

    Returns ``(allele1, allele2)`` or ``(None, None)`` to signal "skip
    this record." Hom-ref calls return ``(None, None)`` because the
    user doesn't carry a non-reference allele — there's nothing for
    the annotators to bind to.

    Handles: phased (``0|1``), unphased (``0/1``), missing (``.`` or
    ``./.``), multi-allelic (``1/2``), haploid (``0`` or ``1`` for
    MT/Y).
    """
    gt = fields.get("GT")
    if gt is None:
        # Sites-only file, no sample columns, or no GT in FORMAT.
        return None, None

    # Split on either '/' (unphased) or '|' (phased); semantics identical
    # for our purposes.
    gt_parts = gt.split("|") if "|" in gt else gt.split("/")
    if len(gt_parts) == 1:
        # Haploid call. Duplicate to make a homozygous diploid for the
        # rest of the pipeline, same convention as MT/Y array parsers.
        gt_parts = [gt_parts[0], gt_parts[0]]
    if len(gt_parts) != 2:
        return None, None

    a1 = _lookup_allele(gt_parts[0], ref, alts)
    a2 = _lookup_allele(gt_parts[1], ref, alts)
    if a1 is None or a2 is None:
        return None, None

    # Hom-ref → skip (user carries no ALT).
    if a1 == ref and a2 == ref:
        return None, None

    return a1, a2


def _lookup_allele(gt_index: str, ref: str, alts: list[str]) -> str | None:
    """Map a GT index (``"0"``, ``"1"``, ...) to an allele string.

    ``"0"`` is REF; ``"1+N"`` is the Nth ALT (0-indexed within ``alts``).
    ``"."`` is no-call (returned as NO_CALL_MARKER). Anything else
    returns None.
    """
    if gt_index == ".":
        return NO_CALL_MARKER
    try:
        idx = int(gt_index)
    except ValueError:
        return None
    if idx == 0:
        return ref
    if 1 <= idx <= len(alts):
        allele = alts[idx - 1]
        if _is_symbolic(allele):
            # GT picked a symbolic allele (e.g. <NON_REF>, or the
            # spanning deletion *). It names no sequence, so the call
            # is not actionable: _resolve_genotype drops the record when
            # either haplotype is symbolic. A '*' therefore never reaches
            # a genotype string, an annotator or an exporter.
            return None
        return allele
    return None


# ── File handling ─────────────────────────────────────────────


def _open_vcf(file_path: Path) -> TextIO:
    r"""Open a ``.vcf`` or ``.vcf.gz`` file in text mode.

    Sniffs the gzip magic bytes ``1f 8b`` so the parser works on any
    extension (e.g., files named ``foo.txt`` containing valid VCF
    gzipped content still parse). The caller is responsible for closing
    the returned handle (it's a context manager).

    ``errors="replace"`` on the text decode. Real-world VCFs
    routinely carry non-UTF-8 bytes in tool headers — a ``©`` (latin-1
    ``0xa9``) in a ``##CL=`` line is extremely common (GATK / DeepVariant
    / bcftools all do it). A strict decode would raise ``UnicodeDecodeError``
    mid-iteration, leaking a Python traceback through ``can_parse`` and
    ``parse``. Replacement renders the bad byte as U+FFFD inside the
    header block — VCF data lines are pure ASCII tab-separated values
    and never reach the substitution path, so genotype fidelity is
    preserved. ``can_parse`` also catches ``UnicodeDecodeError`` as a
    belt-and-braces guard for the case where ``errors="replace"`` is
    ever lifted.

    ``encoding="utf-8-sig"`` (was ``utf-8``). The ``-sig``
    variant transparently strips a leading UTF-8 BOM (``\\xef\\xbb\\xbf``)
    when present and behaves identically to plain ``utf-8`` on
    unprefixed input. Real-world Sequencing.com / DRAGEN VCFs are
    reported to ship with a UTF-8 BOM under some download flows, and
    the BOM defeated ``can_parse``'s ``startswith("##fileformat=VCF")``
    check: ``\\ufeff##fileformat=VCFv4.2`` does not start with
    ``##fileformat=VCF``. Switching the codec is zero-risk on
    BOM-less input (the codec falls through to standard utf-8) and
    fixes the BOM case structurally — no per-call lstrip, no
    encoding-detection guess. Same swap on both the gzip and plain
    text paths so the fix applies regardless of how the file was
    delivered.

    The gzip path is the shared ``open_genotype_text`` reader: a cut or
    corrupt stream raises ``DamagedFileError`` and a BGZF file without its
    end-of-file block is warned about.

    Returns ``typing.TextIO`` rather than the concrete ``io.TextIOBase``
    so mypy-strict callers see the ``__enter__`` / ``__exit__`` protocol
    on the returned handle.
    """
    return open_genotype_text(file_path, errors="replace")

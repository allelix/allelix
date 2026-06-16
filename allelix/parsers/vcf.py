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
by ``END=`` INFO tags in the first hundred-ish data lines. Reference
blocks are **skipped entirely** by ``parse()`` — they match nothing
in any annotation database and would push the annotation pipeline
through millions of zero-hit queries. The "tested-and-reference vs
not-tested" distinction is reserved for a future R-2 enhancement
(see ``CLAUDE.md`` roadmap, v2.1+).

Multi-sample handling
---------------------

``VcfParser(sample=None)`` works on single-sample files. For
multi-sample VCFs, the parser must be instantiated with
``VcfParser(sample="<sample_id>")`` or ``parse()`` raises. The CLI
threads the user's ``--sample`` flag through to the constructor.
"""

from __future__ import annotations

import gzip
import logging
import re
from typing import TYPE_CHECKING, ClassVar, TextIO

from allelix.models import NO_CALL_MARKER, Variant
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

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

# GH #38: match a ``##contig=<ID=chrN,...>`` line declaring any standard
# human chromosome with the ``chr`` prefix. Standard names only — alt
# contigs (``GL00*``, ``hs37d5``, ``NC_*``) don't disambiguate the
# build. The terminator ``[,>]`` keeps us from matching prefixes like
# ``ID=chr1_KI270706v1_random`` (an alt contig) when only chr1 is
# present as a standard contig.
_CHR_PREFIX_CONTIG_RE = re.compile(r"ID=chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|M)[,>]")


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

    def __init__(self, sample: str | None = None) -> None:
        """``sample``: which sample column to read in a multi-sample VCF.

        ``None`` is fine for single-sample files (the lone sample is
        used automatically). For multi-sample files, ``parse()`` raises
        ``MultiSampleError`` when ``sample`` is None. Use
        ``list_samples(file_path)`` to discover available IDs before
        constructing a sample-bound parser.
        """
        self._sample = sample

    # ── Public API ──────────────────────────────────────────────

    def can_parse(self, file_path: Path) -> bool:
        """True if the first non-blank line is ``##fileformat=VCF...``."""
        try:
            with _open_vcf(file_path) as handle:
                for line in handle:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    return stripped.startswith("##fileformat=VCF")
        except OSError:
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
        with _open_vcf(file_path) as handle:
            header = _read_header(handle)
            sample_idx = self._resolve_sample_index(header, file_path)
            for line in handle:
                stripped = line.rstrip("\n").rstrip("\r")
                if not stripped or stripped.startswith("#"):
                    continue
                variant = _parse_data_line(stripped, sample_idx)
                if variant is not None:
                    yield variant

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
        return GenotypeMetadata(
            format=self.name,
            sample_id=sample_id,
            build=header.build or "",
            chr_prefix_observed=header.chr_prefix_observed,
        )

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

    __slots__ = ("build", "chr_prefix_observed", "has_non_ref_alt", "samples")

    def __init__(self) -> None:
        self.samples: list[str] = []
        self.build: str | None = None
        self.has_non_ref_alt: bool = False
        # GH #38: ``chr``-prefixed contig names indicate GRCh38 in modern
        # variant callers (DeepVariant, DRAGEN, GATK HaplotypeCaller).
        # Tertiary build-detection signal when rsIDs and ``##assembly``
        # both fail to converge.
        self.chr_prefix_observed: bool = False


def _read_header(handle: TextIO) -> _VcfHeader:
    """Consume lines from ``handle`` up to and including the ``#CHROM`` row.

    Returns parsed metadata. After this call the handle is positioned
    at the first data line (or EOF).
    """
    header = _VcfHeader()
    for line in handle:
        stripped = line.rstrip("\n").rstrip("\r")
        if not stripped:
            continue
        if stripped.startswith("##"):
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
        # Reaching a data line before #CHROM means no column header.
        # Return what we have; samples is empty.
        return header
    return header


def _absorb_meta_line(line: str, header: _VcfHeader) -> None:
    """Extract build (from ``##contig``) and gVCF marker from a ``##`` line."""
    if line.startswith("##ALT=") and "ID=NON_REF" in line:
        header.has_non_ref_alt = True
        return
    if line.startswith("##contig="):
        # GH #38: capture the chr-prefix signal once any contig declares
        # it. Match ``ID=chr`` followed by any standard chromosome name
        # (1-22, X, Y, M, MT) terminated by ``,`` or ``>`` so we don't
        # false-positive on alt contigs and decoy sequences (``GL00*``,
        # ``hs37d5``, ``NC_*`` — none disambiguate the build the same
        # way). Previously only checked ``chr1`` and ``chrX``; this
        # widening (v2.0.2) catches per-chromosome VCFs and slices that
        # omit chr1.
        if not header.chr_prefix_observed and _CHR_PREFIX_CONTIG_RE.search(line):
            header.chr_prefix_observed = True
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


def _parse_data_line(line: str, sample_idx: int) -> Variant | None:
    """Parse one tab-separated data line into a ``Variant`` or ``None``.

    Returns ``None`` for: malformed rows, reference-only rows
    (ALT='.'), symbolic-allele rows (``<NON_REF>``, ``<DEL>``, ``<*>``,
    etc.), no-genotype-column files when a sample index was expected,
    and homozygous-reference calls (the user carries only REF — no
    actionable variant).
    """
    cols = line.split("\t")
    if len(cols) < 8:
        logger.warning("Malformed VCF line (fewer than 8 columns): %s", line[:80])
        return None
    chrom = cols[0]
    pos_str = cols[1]
    rsid = cols[2]
    ref = cols[3]
    alt = cols[4]
    # cols[5] = QUAL, cols[6] = FILTER, cols[7] = INFO — not consulted here.

    if alt == "." or alt == "":
        # Reference-only record. Nothing to annotate.
        return None

    try:
        position = int(pos_str)
    except ValueError:
        logger.warning("Invalid VCF position %r", pos_str)
        return None

    alts = alt.split(",")
    # If every ALT is symbolic (<NON_REF>, <DEL>, <CNV>, <*>), this is a
    # pure reference block or structural variant record — skip entirely.
    # A mixed ALT like 'A,<NON_REF>' (typical gVCF variant line) keeps
    # going; _lookup_allele rejects GT indices pointing at symbolic
    # alleles so the variant only yields when the user actually carries
    # one of the real alleles.
    if all(a.startswith(_SYMBOLIC_ALT_PREFIX) for a in alts):
        return None

    # Resolve sample's genotype call.
    a1, a2 = _resolve_genotype(cols, sample_idx, ref, alts)
    if a1 is None or a2 is None:
        return None

    return Variant(
        rsid=_canonical_rsid(rsid),
        chromosome=normalize_chromosome(chrom),
        position=position,
        allele1=a1,
        allele2=a2,
        # build is overridden by the pipeline (ADR-0021); leave default.
    )


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


def _resolve_genotype(
    cols: list[str], sample_idx: int, ref: str, alts: list[str]
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
    if sample_idx < 0:
        # Sites-only VCF or no sample resolved. Caller decides; we
        # signal no actionable genotype available.
        return None, None
    if len(cols) < 10:
        # No sample columns. Can't read genotypes.
        return None, None
    if sample_idx + 9 >= len(cols):
        # Out-of-range — should have been caught at header parse.
        return None, None
    fmt = cols[8].split(":")
    try:
        gt_field_idx = fmt.index("GT")
    except ValueError:
        # No GT in FORMAT. Can't read genotypes.
        return None, None
    sample_field = cols[9 + sample_idx]
    sample_parts = sample_field.split(":")
    if gt_field_idx >= len(sample_parts):
        return None, None
    gt = sample_parts[gt_field_idx]

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
        if allele.startswith(_SYMBOLIC_ALT_PREFIX):
            # GT picked a symbolic allele (e.g. <NON_REF>). Treat as a
            # non-actionable call — the parser yields nothing for this
            # haplotype, and the variant is dropped if both haplotypes
            # are symbolic.
            return None
        return allele
    return None


# ── File handling ─────────────────────────────────────────────


def _open_vcf(file_path: Path) -> TextIO:
    """Open a ``.vcf`` or ``.vcf.gz`` file in text mode.

    Sniffs the gzip magic bytes ``1f 8b`` so the parser works on any
    extension (e.g., files named ``foo.txt`` containing valid VCF
    gzipped content still parse). The caller is responsible for closing
    the returned handle (it's a context manager).

    Returns ``typing.TextIO`` rather than the concrete ``io.TextIOBase``
    so mypy-strict callers see the ``__enter__`` / ``__exit__`` protocol
    on the returned handle.
    """
    with open(file_path, "rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(file_path, "rt", encoding="utf-8")
    return open(file_path, encoding="utf-8")

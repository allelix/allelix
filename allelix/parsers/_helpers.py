# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared helpers for parsers with CSV or concatenated-genotype formats.

Used by FTDNA, MyHeritage, and Living DNA parsers. Extracted here to avoid
duplicating the genotype-splitting and CSV-line-splitting logic across
structurally similar formats.
"""

from __future__ import annotations

import logging

from allelix.models import NO_CALL_MARKER

logger = logging.getLogger(__name__)


def split_csv_line(line: str, columns: int | None = None) -> list[str]:
    """Split a comma-delimited line and strip surrounding quotes from each field.

    Implementation is ``line.split(",")`` followed by a per-field
    ``strip().strip('"')``. This is NOT a real CSV parser: a quoted field
    containing a literal comma yields the wrong column count and is
    silently dropped by callers' ``len(parts) != EXPECTED_COLUMNS``
    guard.

    Adequate for FTDNA / MyHeritage / Living DNA because every value in
    those exports is either an rsID, chromosome identifier, integer
    position, or concatenated genotype string — none of which contain
    commas. If a future format ever ships embedded commas in quoted
    fields, swap to ``csv.reader`` rather than relying on this helper.

    Strips both surrounding double quotes (``"rs1"``) and the
    double-double-quote variant some MyHeritage exports produce
    (``""rs1""``) — the latter via two iterations of the trailing
    ``strip('"')``.

    ``columns``: the format's column count. Trailing empty fields beyond
    it are dropped: a spreadsheet re-save pads every row to the widest
    column count (``"rs1","1","82154","AA",,,``), which would otherwise
    read as seven columns and be skipped. Fields up to ``columns`` are
    kept even when empty, so an empty genotype stays a no-call.
    """
    fields = [field.strip().strip('"') for field in line.split(",")]
    return fields if columns is None else strip_trailing_empty(fields, keep=columns)


def strip_trailing_empty(fields: list[str], keep: int = 0) -> list[str]:
    """Drop trailing empty fields (spreadsheet padding such as 'a,b,,,' -> ['a','b']).

    Never shortens the list below ``keep`` fields.
    """
    while len(fields) > keep and fields[-1] == "":
        fields.pop()
    return fields


def is_padding_only(line: str) -> bool:
    """True for a row of nothing but delimiters (``,,,``): an emptied spreadsheet row."""
    return not line.replace(",", "").replace("\t", "").replace('"', "").strip()


def split_genotype(genotype: str) -> tuple[str, str]:
    """Split a concatenated genotype field into two alleles.

    ``"AG"`` -> ``("A", "G")``, ``"--"`` -> ``("-", "-")``,
    ``"A"`` -> ``("A", "A")`` (haploid MT/Y).
    """
    if genotype == "--":
        return NO_CALL_MARKER, NO_CALL_MARKER
    if len(genotype) == 2:
        return genotype[0], genotype[1]
    if len(genotype) == 1:
        return genotype, genotype
    logger.warning("Unexpected genotype format %r — treating as no-call", genotype)
    return NO_CALL_MARKER, NO_CALL_MARKER


def normalize_chromosome(chrom: str) -> str:
    """Normalize a chromosome identifier to Allelix's canonical bare form.

    VCF/UCSC sources commonly prefix chromosomes with ``chr`` (``chr1``,
    ``chrX``, ``chrM``). Allelix's internal representation uses the bare
    form (``1``, ``X``, ``MT``). Returns the bare form regardless of input.

    Mappings:
      - ``chr1`` / ``CHR1`` / ``1`` -> ``1``
      - ``chrX`` / ``chrx`` / ``X`` / ``x`` -> ``X`` (same for Y)
      - ``chrM`` / ``chrm`` / ``chrMT`` / ``chrmt`` / ``M`` / ``MT`` /
        ``m`` / ``mt`` -> ``MT``
      - Unknown contigs (``chrUn_*``, ``GL00*``, etc.) pass through unchanged
        after the optional ``chr`` strip — downstream code decides whether
        to skip them. Case is preserved for unknown contigs.

    Mitochondrial nomenclature varies: GRCh37 chrMT, GRCh38 chrM, some
    VCFs ship bare M. All collapse to MT for Allelix.

    Standard chromosome names (1-22, X, Y, M, MT) are normalized to
    upper case so a pipeline that lowercases the CHROM column doesn't
    silently miss every annotation database lookup. Unknown contigs
    preserve case because their identifiers are opaque and the
    canonical-vs-non-canonical distinction is what matters.
    """
    if not chrom:
        return chrom
    # Strip case-insensitive 'chr' prefix
    stripped = chrom[3:] if chrom[:3].lower() == "chr" else chrom
    upper = stripped.upper()
    # Collapse mitochondrial variants to canonical MT
    if upper in ("M", "MT"):
        return "MT"
    # Standard sex chromosomes normalize to uppercase
    if upper in ("X", "Y"):
        return upper
    # Autosomes (digits) are unchanged by upper(); return as-is
    if stripped.isdigit():
        return stripped
    # Unknown contigs (GL00*, chrUn_*, alts) preserve case
    return stripped

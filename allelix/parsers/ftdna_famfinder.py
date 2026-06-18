# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Parser for Family Tree DNA FamFinder export files.

A third FTDNA file shape — distinct from the comma-delimited CSV
variant in ``ftdna.py`` and the tab-delimited concatenated-RESULT
variant in ``ftdna_illumina.py``. The FamFinder export splits the
diploid call across **separate ALLELE1 / ALLELE2 columns** instead
of concatenating into a single RESULT.

Format reference (per CLAUDE.md "Parser Format Specifications" — FTDNA
section, FamFinder paragraph):

    # Family Tree DNA - FamFinder
    # Build 37
    RSID    CHROMOSOME      POSITION        ALLELE1     ALLELE2
    rs4477212       1       82154   A       A
    rs3094315       1       752566  A       G
    rs9001001       1       100000  -       -

Specifics:
    - Tab-delimited (same delimiter as ``ftdna_illumina.py``).
    - Comment lines start with ``#``.
    - Header line: ``RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2``
      (case-insensitive); 5 columns vs FTDNA Illumina raw's 4.
    - Separate single-character allele columns (vs. concatenated
      RESULT in the Illumina raw variant).
    - Haploid calls on MT/Y use the single-character pattern
      ``A``/``-`` in ALLELE1 with empty / missing ALLELE2.
    - No-calls: ``-`` in one or both allele columns.
    - Build 37 (GRCh37) on every export observed.
    - Detection key: ``famfinder`` substring (case-insensitive) in
      the first ``SNIFF_LINE_LIMIT`` lines AND the canonical 5-column
      header. Requiring both signals avoids false positives against
      any future tab-delimited 5-column variant.

Registration order in ``allelix/parsers/__init__.py``: FamFinder
must precede ``FTDNAIlluminaParser`` so the more-specific 5-column
header check fires before a permissive 4-column parser could
misclassify a malformed FamFinder file.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 5
HEADER_COLUMNS = ("RSID", "CHROMOSOME", "POSITION", "ALLELE1", "ALLELE2")
_MARKER = "famfinder"


def _is_header_line(line: str) -> bool:
    """True if *line* is the canonical tab-delimited FamFinder header."""
    parts = line.strip().split("\t")
    if len(parts) != EXPECTED_COLUMNS:
        return False
    return tuple(p.strip().upper() for p in parts) == HEADER_COLUMNS


def _normalize_allele(raw: str) -> str:
    """Map a single FamFinder allele cell to the Variant.allele convention.

    Empty cell or ``-`` is the no-call signal — both map to ``-`` so
    downstream code (``Variant.is_no_call``) handles them uniformly.
    """
    stripped = raw.strip().upper()
    if not stripped or stripped == "-":
        return "-"
    return stripped


class FTDNAFamFinderParser(GenotypeParser):
    """Parser for FTDNA FamFinder tab-delimited genotype files with separate allele columns."""

    name: ClassVar[str] = "ftdna_famfinder"
    display_name: ClassVar[str] = "Family Tree DNA (FamFinder)"
    file_extensions: ClassVar[list[str]] = [".txt", ".csv"]
    url: ClassVar[str] = "https://www.familytreedna.com"

    def can_parse(self, file_path: Path) -> bool:
        """Detect by ``famfinder`` marker + canonical 5-column header.

        Requires BOTH signals so a malformed FTDNA Illumina raw file
        (4-col tab) or a future tab-delimited 5-col variant cannot
        accidentally match.
        """
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                marker_seen = False
                for _ in range(SNIFF_LINE_LIMIT):
                    line = fh.readline()
                    if not line:
                        return False
                    stripped = line.rstrip("\r\n")
                    if not stripped:
                        continue
                    if _MARKER in stripped.lower():
                        marker_seen = True
                        continue
                    if stripped.startswith("#"):
                        continue
                    # First non-comment, non-blank, non-marker line must
                    # be the header — anything else means this isn't a
                    # FamFinder file (or the file is malformed).
                    return marker_seen and _is_header_line(stripped)
        except (OSError, UnicodeDecodeError):
            return False
        return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream Variant objects, skipping comments and malformed lines."""
        with file_path.open("r", encoding="utf-8") as fh:
            header_seen = False
            for lineno, raw in enumerate(fh, start=1):
                line = raw.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                if not header_seen:
                    if _is_header_line(line):
                        header_seen = True
                        continue
                    # Skip non-comment marker lines (e.g. "Family Tree DNA — FamFinder")
                    # that precede the header.
                    continue

                parts = line.split("\t")
                if len(parts) != EXPECTED_COLUMNS:
                    logger.warning(
                        "Line %d: expected %d columns, got %d — skipping",
                        lineno,
                        EXPECTED_COLUMNS,
                        len(parts),
                    )
                    continue

                rsid = parts[0].strip()
                chrom = parts[1].strip()
                pos_str = parts[2].strip()
                allele1 = _normalize_allele(parts[3])
                allele2 = _normalize_allele(parts[4])

                try:
                    position = int(pos_str)
                except ValueError:
                    logger.warning("Line %d: invalid position %r — skipping", lineno, pos_str)
                    continue

                yield Variant(
                    rsid=rsid,
                    chromosome=normalize_chromosome(chrom),
                    position=position,
                    allele1=allele1,
                    allele2=allele2,
                    build=DEFAULT_BUILD,
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """FamFinder exports carry no structured sample-ID metadata."""
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=DEFAULT_BUILD,
        )

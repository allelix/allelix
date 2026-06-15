# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Parser for Family Tree DNA Illumina raw export files.

A second FTDNA file shape distinct from the comma-delimited CSV
variant handled by ``ftdna.py``. FTDNA exports both formats from the
same chip data; this one uses tab delimiters and an unquoted
``RSID/CHROMOSOME/POSITION/RESULT`` header, typically named
``*.ftdna-illumina.txt``.

Format reference (from real sample files):

    # Family Tree DNA - Family Finder
    # Illumina raw data export
    # Build 37
    RSID    CHROMOSOME      POSITION        RESULT
    rs4477212       1       82154   AA
    rs3094315       1       752566  AG
    rs9001001       1       100000  --

Specifics:
    - Tab-delimited (not CSV — distinguishes from ``ftdna.py``).
    - Comment lines start with ``#``.
    - Header line: ``RSID\\tCHROMOSOME\\tPOSITION\\tRESULT``
      (case-insensitive).
    - RESULT column is concatenated genotype ("AG", "AA").
    - Haploid calls on MT/Y appear as single characters.
    - No-calls represented as ``--``.
    - Build 37 (GRCh37) on every export observed.
    - Detection key: tab-delimited header matching the canonical
      column names within the first 50 non-comment lines.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import normalize_chromosome, split_genotype
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4
HEADER_COLUMNS = ("RSID", "CHROMOSOME", "POSITION", "RESULT")


def _is_header_line(line: str) -> bool:
    """True if *line* is the canonical tab-delimited FTDNA Illumina header."""
    parts = line.strip().split("\t")
    if len(parts) != EXPECTED_COLUMNS:
        return False
    return tuple(p.strip().upper() for p in parts) == HEADER_COLUMNS


class FTDNAIlluminaParser(GenotypeParser):
    """Parser for FTDNA Illumina raw tab-delimited genotype files."""

    name: ClassVar[str] = "ftdna_illumina"
    display_name: ClassVar[str] = "Family Tree DNA (Illumina raw)"
    file_extensions: ClassVar[list[str]] = [".txt"]
    url: ClassVar[str] = "https://www.familytreedna.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize by tab-delimited ``RSID/CHROMOSOME/POSITION/RESULT`` header."""
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                for _ in range(SNIFF_LINE_LIMIT):
                    line = fh.readline()
                    if not line:
                        return False
                    stripped = line.rstrip("\r\n")
                    if not stripped or stripped.startswith("#"):
                        continue
                    return _is_header_line(stripped)
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
                    logger.warning(
                        "Line %d: expected FTDNA Illumina header, got %r — skipping",
                        lineno,
                        line,
                    )
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
                genotype = parts[3].strip()

                try:
                    position = int(pos_str)
                except ValueError:
                    logger.warning("Line %d: invalid position %r — skipping", lineno, pos_str)
                    continue

                allele1, allele2 = split_genotype(genotype)

                yield Variant(
                    rsid=rsid,
                    chromosome=normalize_chromosome(chrom),
                    position=position,
                    allele1=allele1,
                    allele2=allele2,
                    build=DEFAULT_BUILD,
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """FTDNA Illumina raw files carry no sample-ID metadata header."""
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=DEFAULT_BUILD,
        )

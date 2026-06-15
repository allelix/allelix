# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for Family Tree DNA (FTDNA) raw genotype export files.

Format reference (from real sample files and snps package):

    # FTDNA raw data download
    RSID,CHROMOSOME,POSITION,RESULT
    "rs4477212","1","82154","AA"
    "rs3094315","1","752566","AG"
    "rs9001001","1","100000","--"

Specifics:
    - CSV format, comma-delimited.
    - Optional comment lines starting with ``#``.
    - Header line: ``RSID,CHROMOSOME,POSITION,RESULT`` (quoted or unquoted).
    - Data fields are double-quoted; parser strips quotes.
    - RESULT column is concatenated genotype (e.g., "AG" not "A","G").
    - Haploid calls on MT/Y appear as single characters (e.g., "A").
    - No-calls represented as ``--``.
    - Build 37 (most files).
    - Detection key: header line matching ``RSID,CHROMOSOME,POSITION,RESULT``
      (case-insensitive, with or without quotes) within the first 50 lines.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import split_csv_line, split_genotype
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4
HEADER_CANONICAL = "RSID,CHROMOSOME,POSITION,RESULT"


def _is_header_line(line: str) -> bool:
    """True if *line* is the FTDNA column header (quoted or unquoted)."""
    stripped = line.replace('"', "").replace("'", "").strip()
    return stripped.upper() == HEADER_CANONICAL


class FTDNAParser(GenotypeParser):
    """Parser for Family Tree DNA consumer DNA genotype files."""

    name: ClassVar[str] = "ftdna"
    display_name: ClassVar[str] = "Family Tree DNA"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://www.familytreedna.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by its ``RSID,CHROMOSOME,POSITION,RESULT`` header."""
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                for _ in range(SNIFF_LINE_LIMIT):
                    line = fh.readline()
                    if not line:
                        return False
                    line = line.rstrip("\r\n")
                    if not line or line.startswith("#"):
                        continue
                    return _is_header_line(line)
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
                        "Line %d: expected FTDNA header, got %r — skipping",
                        lineno,
                        line,
                    )
                    continue

                parts = split_csv_line(line)
                if len(parts) != EXPECTED_COLUMNS:
                    logger.warning(
                        "Line %d: expected %d columns, got %d — skipping",
                        lineno,
                        EXPECTED_COLUMNS,
                        len(parts),
                    )
                    continue

                rsid, chrom, pos_str, genotype = parts
                try:
                    position = int(pos_str)
                except ValueError:
                    logger.warning("Line %d: invalid position %r — skipping", lineno, pos_str)
                    continue

                allele1, allele2 = split_genotype(genotype)

                yield Variant(
                    rsid=rsid,
                    chromosome=chrom,
                    position=position,
                    allele1=allele1,
                    allele2=allele2,
                    build=DEFAULT_BUILD,
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Extract metadata from header. FTDNA files have no sample ID field."""
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=DEFAULT_BUILD,
        )

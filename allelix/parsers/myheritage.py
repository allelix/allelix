# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for MyHeritage DNA raw genotype export files.

Format reference (from real sample files and snps package):

    # MyHeritage, https://www.myheritage.com
    RSID,CHROMOSOME,POSITION,RESULT
    "rs4477212","1","82154","AA"
    "rs3094315","1","752566","AG"
    "rs9001001","1","100000","--"

Specifics:
    - CSV format, comma-delimited. Structurally identical to FTDNA.
    - Detection key: ``MyHeritage`` in the first comment line.
    - Header line: ``RSID,CHROMOSOME,POSITION,RESULT`` (quoted or unquoted).
    - Data fields are double-quoted; some exports double-double-quote
      fields (``""rs1""``). ``split_csv_line`` handles both.
    - RESULT column is concatenated genotype (e.g., "AG" not "A","G").
    - Haploid calls on MT/Y appear as single characters (e.g., "A").
    - No-calls represented as ``--``.
    - Build: not declared in file; position-based detection required.
      Defaults to GRCh37.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import split_csv_line, split_genotype
from allelix.parsers.base import GenotypeMetadata, GenotypeParser
from allelix.parsers.ftdna import HEADER_CANONICAL, _is_header_line

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SIGNATURE = "MyHeritage"
SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4


class MyHeritageParser(GenotypeParser):
    """Parser for MyHeritage DNA consumer genotype files."""

    name: ClassVar[str] = "myheritage"
    display_name: ClassVar[str] = "MyHeritage DNA"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://www.myheritage.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by ``MyHeritage`` in the first comment line."""
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                first_line = fh.readline()
                return SIGNATURE in first_line
        except (OSError, UnicodeDecodeError):
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
                        "Line %d: expected %s header, got %r — skipping",
                        lineno,
                        HEADER_CANONICAL,
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
        """Extract metadata. MyHeritage files have no sample ID or build field."""
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=DEFAULT_BUILD,
        )

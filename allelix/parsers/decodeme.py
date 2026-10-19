# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for deCODEme raw genotype files (deCODE genetics, service closed 2013).

Format reference (SNPedia; layout confirmed on real files in aggregate):

    Name,Variation,Chromosome,Position,Strand,YourCode
    rs4477212,A/G,1,72017,+,AA
    rs3094315,C/T,1,742429,-,GA
    rs9001001,A/G,M,100,+,--

Specifics:
    - No comment block: the first line is the column line above. It is the
      detection key, matched exactly (case and surrounding space ignored).
    - Comma-separated, 6 columns, CRLF.
    - ``YourCode`` is the genotype, two letters; ``--`` is a no-call.
    - ``Strand`` is ``+`` or ``-`` per row. A ``-`` row's ``YourCode`` (and
      its ``Variation``) are on the minus strand, so the call is complemented
      to the plus strand that every other format and source uses. Checked on
      real files against real 23andMe exports: minus rows match only after
      complementing, plus rows only as given.
    - Chromosomes 1-22, X, Y and ``M`` (read as MT). Males get two letters
      on X and Y; the calls are kept as given.
    - Build 36 (dbSNP 128). The file states no build, so the known-SNP
      positions decide it; a GRCh36 file runs with position lookups off
      (ADR-0025).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import (
    DamagedFileError,
    first_nonblank_line,
    is_haploid_call,
    normalize_chromosome,
    open_genotype_text,
    split_genotype,
)
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

COLUMNS = ("name", "variation", "chromosome", "position", "strand", "yourcode")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_SNIFF_LINES = 20


def _is_column_line(line: str) -> bool:
    return tuple(field.strip().lower() for field in line.split(",")) == COLUMNS


class DeCODEmeParser(GenotypeParser):
    """Parser for deCODEme ``deCODEme_scan.csv`` raw data."""

    name: ClassVar[str] = "decodeme"
    display_name: ClassVar[str] = "deCODEme"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://www.decode.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by its column line, the first non-blank line."""
        try:
            with open_genotype_text(file_path) as fh:
                return _is_column_line(first_nonblank_line(fh, _SNIFF_LINES))
        except (OSError, UnicodeDecodeError, DamagedFileError):
            return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream Variant objects on the plus strand, skipping malformed lines."""
        with open_genotype_text(file_path) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or _is_column_line(line):
                    continue
                parts = [field.strip() for field in line.split(",")]
                if len(parts) != len(COLUMNS):
                    logger.warning(
                        "Line %d: expected %d columns, got %d — skipping",
                        lineno,
                        len(COLUMNS),
                        len(parts),
                    )
                    continue
                rsid, _variation, chrom, pos_str, strand, code = parts
                try:
                    position = int(pos_str)
                except ValueError:
                    logger.warning("Line %d: invalid position %r — skipping", lineno, pos_str)
                    continue
                if strand not in ("+", "-"):
                    logger.warning("Line %d: invalid strand %r — skipping", lineno, strand)
                    continue
                if strand == "-":
                    code = code.translate(_COMPLEMENT)
                allele1, allele2 = split_genotype(code)
                yield Variant(
                    rsid=rsid,
                    chromosome=normalize_chromosome(chrom),
                    position=position,
                    allele1=allele1,
                    allele2=allele2,
                    build=DEFAULT_BUILD,
                    haploid=is_haploid_call(code),
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """DeCODEme files state neither a sample ID nor a build."""
        return GenotypeMetadata(format=self.name, sample_id="", build=None)

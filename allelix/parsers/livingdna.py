# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for Living DNA raw genotype export files.

Format reference (from snps package and H600 Project wiki):

    # Living DNA customer genotype data download file version: 1.0.1
    # This file contains raw genotype data ...
    # Human Genome Reference Build 37 (GRCh37.p13).
    # Genotypes are presented on the forward strand.
    #
    # rsid  chromosome  position  genotype
    rs1801133   1   11856378    AG
    AX-12345678 3   15000000    GG
    1:726912    1   726912  AA

Specifics:
    - Tab-delimited despite ``.csv`` file extension.
    - Detection key: ``Living DNA`` in the first line.
    - Comment lines start with ``#``, including the column header line.
    - Four columns: ``rsid``, ``chromosome``, ``position``, ``genotype``.
    - Concatenated genotype in result column (e.g., "AA", "CT").
    - No-calls represented as ``--``.
    - Build 37 (GRCh37.p13), forward strand.
    - SNP ID types: rs-numbers, ``AX-`` prefixed (Affymetrix),
      ``AFFX-`` prefixed (Affymetrix control probes), and positional
      notation (``CHR:POS``, e.g., ``1:726912``).
    - Y and MT chromosomes delivered as separate files; main file
      has chromosomes 1-22 and X.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import split_genotype
from allelix.parsers.base import GenotypeMetadata, GenotypeParser
from allelix.utils.build_detect import normalize_build_label

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SIGNATURE = "Living DNA"
SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4


class LivingDNAParser(GenotypeParser):
    """Parser for Living DNA consumer genotype files."""

    name: ClassVar[str] = "livingdna"
    display_name: ClassVar[str] = "Living DNA"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://livingdna.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by ``Living DNA`` in the first line."""
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                first_line = fh.readline()
                return SIGNATURE in first_line
        except (OSError, UnicodeDecodeError):
            return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream Variant objects, skipping comments and malformed lines."""
        with file_path.open("r", encoding="utf-8") as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.rstrip("\r\n")
                if not line or line.startswith("#"):
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
        """Extract build from header comments. Living DNA has no sample ID field."""
        build = DEFAULT_BUILD
        with file_path.open("r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\r\n")
                if not line.startswith("#"):
                    break
                normalized = normalize_build_label(line)
                if normalized:
                    build = normalized
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=build,
        )

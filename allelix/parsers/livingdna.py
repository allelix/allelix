# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for Living DNA raw genotype export files.

Format reference (the export layout as publicly described; the file's own header
states GRCh37 and the forward strand; the rows here are invented):

    # Living DNA customer genotype data download file version: 1.0.1
    # This file contains raw genotype data ...
    # Human Genome Reference Build 37 (GRCh37.p13).
    # Genotypes are presented on the forward strand.
    #
    # rsid  chromosome  position  genotype
    rs3999990021    1   3999021     AG
    AX-39999901     4   39999901    TT
    1:2999901       1   2999901     CC

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
      notation (``CHR:POS``, e.g., ``1:2999901``).
    - Y and MT chromosomes delivered as separate files; main file
      has chromosomes 1-22 and X.
"""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import (
    DamagedFileError,
    first_nonblank_line,
    is_haploid_call,
    open_genotype_text,
    split_genotype,
)
from allelix.parsers.base import GenotypeMetadata, GenotypeParser, read_header_build

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SIGNATURE = "Living DNA"
SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4

# only inspect comment lines that look like a build marker.
# Without this filter every comment line is fed to
# ``normalize_build_label`` and a stray date / version digit can override
# the real build line.
_BUILD_MARKER_RE = re.compile(r"\b(build|reference|genome)\b", re.IGNORECASE)


class LivingDNAParser(GenotypeParser):
    """Parser for Living DNA consumer genotype files."""

    name: ClassVar[str] = "livingdna"
    display_name: ClassVar[str] = "Living DNA"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://livingdna.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by ``Living DNA`` in the first non-blank line."""
        try:
            with open_genotype_text(file_path) as fh:
                first_line = first_nonblank_line(fh, SNIFF_LINE_LIMIT)
                return SIGNATURE in first_line
        except (OSError, UnicodeDecodeError, DamagedFileError):
            return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream Variant objects, skipping comments and malformed lines."""
        with open_genotype_text(file_path) as fh:
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
                    haploid=is_haploid_call(genotype),
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Extract build from header comments. Living DNA has no sample ID field.

        Only comment lines that look like an explicit build marker are
        inspected, and the first match wins (``read_header_build``), so a
        download-date comment like ``# downloaded 2038-01-01`` cannot
        retag the file. The format spec puts the build line near the top:

            # Human Genome Reference Build 37 (GRCh37.p13).
        """
        return GenotypeMetadata(
            format=self.name,
            sample_id="",
            build=read_header_build(file_path),
        )

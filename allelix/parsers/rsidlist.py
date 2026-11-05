# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Reader for rsID and genotype lists: two tab-separated columns, no positions.

    rs4477212	AG
    rs3094315	C

Some services export only each marker's ID and its call (IYG files, for
one). The file carries no chromosome, no position and no build, so its
rows are placed on GRCh38 from the rsID identity table when it is
analyzed (``allelix.annotators._placement``); ``has_positions`` tells the
pipeline. Other commands stop on it.

Rules (ADR-0069):
    - Registered last, after every other reader: the shape is plain.
    - Detection is strict. Within the first ``SNIFF_DATA_LINES`` data lines,
      every one has exactly two tab-separated fields (double quotes around a
      field are removed): an ID with no spaces and a valid genotype (one or
      two of A, C, G, T, D, I, or a no-call marker). At least
      ``MIN_RSID_LINES`` lines, and at least ``MIN_RSID_SHARE`` of them,
      have an rsID. Only ``#`` comment lines and blank lines may come
      between them; there is no column line.
    - Calls: two letters are a diploid call, one letter a haploid call;
      ``--``, ``-``, ``00`` and ``0`` are no-calls. Any other call is read
      as a no-call and counted (``unreadable_genotypes``).
    - IDs other than rsIDs (such as ``MT-…`` marker names) are read; they
      cannot be placed and are counted as unplaced by the run.
    - The name says what the file is, not who made it: no marker in the
      file names its source.
"""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, ClassVar

from allelix.models import NO_CALL_MARKER, Variant
from allelix.parsers._helpers import DamagedFileError, open_genotype_text
from allelix.parsers.base import ArrayInputStats, GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SNIFF_DATA_LINES = 50
MIN_RSID_LINES = 3
MIN_RSID_SHARE = 0.8
# Comment and blank lines read before the first data line, at most.
_SNIFF_PREAMBLE = 200

# Placeholder build: the run places every row and sets the build it decided.
_PLACEHOLDER_BUILD = "GRCh38"
_ID = re.compile(r"\S+")
_RSID = re.compile(r"rs\d+")
_CALL = re.compile(r"[ACGTDI]{1,2}|--|-|00|0", re.IGNORECASE)
_NO_CALLS = frozenset({"--", "-", "00", "0"})
_BASES = frozenset("ACGTDI")


def _fields(line: str) -> list[str]:
    return [field.strip().strip('"') for field in line.split("\t")]


class RsidGenotypeListParser(GenotypeParser):
    """Two columns, ID and call, tab-separated, with no positions."""

    name: ClassVar[str] = "rsid-genotype"
    display_name: ClassVar[str] = "rsID and genotype list"
    file_extensions: ClassVar[list[str]] = [".txt", ".tsv"]
    url: ClassVar[str] = ""
    # No positions: the rows are placed from the rsID identity table.
    has_positions: ClassVar[bool] = False

    def __init__(self) -> None:
        """No counts until the first ``parse()``."""
        self._stats: ArrayInputStats | None = None

    def can_parse(self, file_path: Path) -> bool:
        """Every sniffed data line is ``ID<TAB>genotype``; at least 3, and 80%, rsIDs."""
        try:
            with open_genotype_text(file_path) as fh:
                lines = 0
                rsids = 0
                preamble = 0
                for raw in fh:
                    line = raw.rstrip("\r\n")
                    if not line.strip() or line.startswith("#"):
                        if lines == 0:
                            preamble += 1
                            if preamble > _SNIFF_PREAMBLE:
                                return False
                        continue
                    fields = _fields(line)
                    if (
                        len(fields) != 2
                        or not _ID.fullmatch(fields[0])
                        or not _CALL.fullmatch(fields[1])
                    ):
                        return False
                    lines += 1
                    rsids += bool(_RSID.fullmatch(fields[0]))
                    if lines >= SNIFF_DATA_LINES:
                        break
        except (OSError, UnicodeDecodeError, DamagedFileError):
            return False
        return rsids >= MIN_RSID_LINES and rsids >= MIN_RSID_SHARE * lines

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream unplaced Variant objects (chromosome "", position 0)."""
        stats = ArrayInputStats()
        self._stats = stats
        with open_genotype_text(file_path) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.rstrip("\r\n")
                if not line.strip() or line.startswith("#"):
                    continue
                fields = _fields(line)
                if len(fields) != 2 or not fields[0]:
                    stats.wrong_column_rows += 1
                    logger.warning(
                        "Line %d: expected 2 columns, got %d — skipping", lineno, len(fields)
                    )
                    continue
                rsid, call = fields
                allele1, allele2, haploid, unreadable = _alleles(call.upper())
                if unreadable:
                    stats.unreadable_genotypes += 1
                yield Variant(
                    rsid=rsid,
                    chromosome="",
                    position=0,
                    allele1=allele1,
                    allele2=allele2,
                    build=_PLACEHOLDER_BUILD,
                    haploid=haploid,
                )

    def input_stats(self) -> ArrayInputStats | None:
        """Counts from the most recent ``parse()``."""
        return self._stats

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """No sample ID and no build: the file states neither."""
        return GenotypeMetadata(format=self.name, sample_id="", build=None)


def _alleles(call: str) -> tuple[str, str, bool, bool]:
    """``(allele1, allele2, haploid, unreadable)`` for one call field."""
    if call in _NO_CALLS:
        return NO_CALL_MARKER, NO_CALL_MARKER, False, False
    if len(call) == 1 and call in _BASES:
        return call, call, True, False
    if len(call) == 2 and set(call) <= _BASES:
        return call[0], call[1], False, False
    return NO_CALL_MARKER, NO_CALL_MARKER, False, True

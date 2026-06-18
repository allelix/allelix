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
    - Haploid calls on MT/Y use the single-character pattern: the
      called base in ALLELE1 with EITHER the ALLELE2 column entirely
      absent (4-column line) OR an empty ALLELE2 cell (5-column
      line with trailing tab). Both shapes have been observed in
      real exports. The parser normalizes haploid input to a
      hemizygous-looking Variant where ``allele1 == allele2 ==
      <the called base>`` (consistent with the haploid convention
      in ``_helpers.split_genotype``). GH #113.

      Chromosome-aware: the haploid shape (4-column line OR
      5-column line with empty ALLELE2) is accepted ONLY when the
      chromosome is MT or Y. FTDNA does not publish haploid calls
      on any other chromosome, so the same shape on an autosome or
      X is almost certainly a truncated diploid that lost an
      allele in transit — promoting it to hemizygous would silently
      double the surviving allele and synthesize a wrong-zygosity
      genotype. The parser warns-and-skips that case instead.
      Tightened against the v2.2.1 cross-PR review on PR #117.
    - No-calls: ``-`` in one or both allele columns (with both
      columns present and ALLELE1 != empty).
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
                # GH #113: accept the 4-column haploid shape (missing
                # ALLELE2 column for hemizygous MT/Y calls) in addition
                # to the canonical 5-column diploid shape. The docstring
                # already claimed haploid support; the impl didn't.
                # Anything other than 4 or 5 columns is genuinely
                # malformed.
                if len(parts) not in (EXPECTED_COLUMNS - 1, EXPECTED_COLUMNS):
                    logger.warning(
                        "Line %d: expected %d columns (or %d for haploid), got %d — skipping",
                        lineno,
                        EXPECTED_COLUMNS,
                        EXPECTED_COLUMNS - 1,
                        len(parts),
                    )
                    continue

                rsid = parts[0].strip()
                chrom = parts[1].strip()
                pos_str = parts[2].strip()
                allele1 = _normalize_allele(parts[3])
                # Normalize early — the chromosome-aware haploid guard
                # below compares against the canonical form {"MT", "Y"}.
                normalized_chrom = normalize_chromosome(chrom)
                # Haploid: either the 4-column shape (no ALLELE2 column)
                # OR the 5-column shape with an empty ALLELE2 cell. Both
                # are real in FTDNA's MT/Y output. A genuine hemizygous
                # call is represented as allele1 == allele2 (the
                # convention shared by 23andMe / _helpers.split_genotype),
                # NOT allele2 = NO_CALL_MARKER (which would render as a
                # half no-call and surface the GH #113 reporter's
                # "lab-spec MT/Y vanishes" behavior).
                #
                # GH #113 cross-PR review (#117 review): only accept the
                # haploid shape on MT and Y — the only chromosomes FTDNA
                # publishes haploid calls on. A 4-column line (or
                # 5-column-empty-ALLELE2) on an autosome or X is almost
                # certainly a truncated diploid that lost an allele in
                # transit; promoting it to hemizygous would silently
                # double the surviving allele. Warn-skip instead. This
                # closes the trade-off the v2.3 follow-up was opened
                # for; no real MT/Y haploid call is lost (the canonical
                # shape on those chromosomes is unchanged).
                if len(parts) == EXPECTED_COLUMNS - 1:
                    if normalized_chrom not in {"MT", "Y"}:
                        logger.warning(
                            "Line %d: %d-column haploid line on non-haploid "
                            "chromosome %r — almost certainly a truncated "
                            "diploid; skipping",
                            lineno,
                            len(parts),
                            chrom,
                        )
                        continue
                    allele2 = allele1
                else:
                    raw_allele2 = parts[4].strip()
                    if raw_allele2 == "":
                        if normalized_chrom not in {"MT", "Y"}:
                            logger.warning(
                                "Line %d: empty ALLELE2 on non-haploid "
                                "chromosome %r — almost certainly a "
                                "truncated diploid; skipping",
                                lineno,
                                chrom,
                            )
                            continue
                        allele2 = allele1
                    else:
                        allele2 = _normalize_allele(parts[4])

                try:
                    position = int(pos_str)
                except ValueError:
                    logger.warning("Line %d: invalid position %r — skipping", lineno, pos_str)
                    continue

                yield Variant(
                    rsid=rsid,
                    chromosome=normalized_chrom,
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

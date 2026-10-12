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
    - Detection key: ``MyHeritage`` in one of the first lines. A file
      re-saved through a spreadsheet can start with an emptied row
      (``,,,``), drop the ``#`` from comment lines and pad every row with
      trailing commas; all three are tolerated. Tab-delimited rows (a
      spreadsheet "save as text") are read as well.
    - Low-pass WGS exports (kits processed from January 2026) start with
      ``##fileformat=MyHeritage`` and carry ``##method=Low-pass Whole
      Genome Sequencing`` and ``##reference=build37``. The genotypes are
      imputed from ~2x coverage. ``get_metadata`` reports both lines.
    - Header line: ``RSID,CHROMOSOME,POSITION,RESULT`` (quoted or unquoted).
    - Data fields are double-quoted; some exports double-double-quote
      fields (``""rs1""``). ``split_csv_line`` handles both.
    - RESULT column is concatenated genotype (e.g., "AG" not "A","G").
    - Haploid calls on MT/Y appear as single characters (e.g., "A").
    - No-calls represented as ``--``.
    - Build: array exports declare none (position-based detection
      required; defaults to GRCh37). ``##reference=buildNN`` is an
      explicit header build.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import DEFAULT_BUILD, Variant
from allelix.parsers._helpers import (
    is_padding_only,
    split_csv_line,
    split_genotype,
    strip_trailing_empty,
)
from allelix.parsers.base import GenotypeMetadata, GenotypeParser, build_statement
from allelix.parsers.ftdna import HEADER_CANONICAL, PREAMBLE_WARNING_LIMIT, _is_header_line
from allelix.utils.build_detect import normalize_build_label

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

SIGNATURE = "MyHeritage"
SNIFF_LINE_LIMIT = 50
EXPECTED_COLUMNS = 4

# ``##method=`` values -> GenotypeMetadata ``method``.
_METHODS: dict[str, str] = {
    "low-pass whole genome sequencing": "low-pass-wgs-imputed",
}


def _split_row(line: str) -> list[str]:
    """Split a data or header row: tab-delimited when it has tabs, else CSV.

    Quotes and trailing padding past the four columns are removed either way.
    """
    if "\t" in line:
        fields = [field.strip().strip('"') for field in line.split("\t")]
        return strip_trailing_empty(fields, keep=EXPECTED_COLUMNS)
    return split_csv_line(line, EXPECTED_COLUMNS)


def _is_column_header(line: str) -> bool:
    """The ``RSID,CHROMOSOME,POSITION,RESULT`` header, comma- or tab-delimited."""
    return _is_header_line(",".join(_split_row(line)))


def _looks_like_data_row(line: str) -> bool:
    """Four non-empty fields with an integer position: a genotype row."""
    fields = _split_row(line)
    return len(fields) == EXPECTED_COLUMNS and all(fields) and fields[2].isdigit()


class MyHeritageParser(GenotypeParser):
    """Parser for MyHeritage DNA consumer genotype files."""

    name: ClassVar[str] = "myheritage"
    display_name: ClassVar[str] = "MyHeritage DNA"
    file_extensions: ClassVar[list[str]] = [".csv"]
    url: ClassVar[str] = "https://www.myheritage.com"

    def can_parse(self, file_path: Path) -> bool:
        """Recognize the file by ``MyHeritage`` in one of its first lines.

        Scans up to ``SNIFF_LINE_LIMIT`` lines and stops, unclaimed, at the
        column header or the first genotype row: the signature belongs to
        the preamble, so a file without it before its data is never claimed
        (FTDNA exports share the column layout).
        """
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                for _ in range(SNIFF_LINE_LIMIT):
                    line = fh.readline()
                    if not line:
                        return False
                    if SIGNATURE in line:
                        return True
                    if _is_column_header(line) or _looks_like_data_row(line):
                        return False
        except (OSError, UnicodeDecodeError):
            return False
        return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Stream Variant objects, skipping comments and malformed lines."""
        with file_path.open("r", encoding="utf-8") as fh:
            header_seen = False
            preamble_warnings = 0
            for lineno, raw in enumerate(fh, start=1):
                line = raw.rstrip("\r\n")
                if not line or line.startswith("#") or is_padding_only(line):
                    continue
                if not header_seen:
                    if _is_column_header(line):
                        header_seen = True
                        continue
                    # A spreadsheet re-save drops the '#' from comment lines;
                    # they land here. Warn about the first few only.
                    preamble_warnings += 1
                    if preamble_warnings <= PREAMBLE_WARNING_LIMIT:
                        logger.warning(
                            "Line %d: expected %s header, got %r — skipping",
                            lineno,
                            HEADER_CANONICAL,
                            line,
                        )
                    continue

                parts = _split_row(line)
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
        """Extract metadata from the ``##key=value`` preamble.

        No sample ID in either export. Low-pass WGS exports declare
        ``##method`` and ``##reference=buildNN``: the method becomes
        ``method`` and the build ``build`` (the file's own label; the
        build decision weighs it against the positions, ADR-0044). Array
        exports state the build in a prose comment line instead, which is
        read the same way; ``build`` is None when the file states none.
        """
        meta = GenotypeMetadata(format=self.name, sample_id="", build=None)
        try:
            with file_path.open("r", encoding="utf-8") as fh:
                for _ in range(SNIFF_LINE_LIMIT):
                    line = fh.readline()
                    if not line:
                        break
                    stripped = line.strip()
                    if not stripped.startswith("##") or "=" not in stripped:
                        if _is_column_header(stripped):
                            break
                        # Array exports state the build in prose: "# ... with
                        # respect to the human reference build 37." (also
                        # without the "#" once a spreadsheet re-saves it).
                        if meta["build"] is None:
                            meta["build"] = build_statement(stripped)
                        continue
                    key, value = stripped[2:].split("=", 1)
                    key, value = key.strip().lower(), value.strip()
                    if key == "method" and value:
                        meta["method"] = _METHODS.get(value.lower(), value)
                    elif key == "reference":
                        build = normalize_build_label(value)
                        if build is not None:
                            meta["build"] = build
        except (OSError, UnicodeDecodeError):
            pass
        return meta

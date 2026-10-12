# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Abstract base class for genotype file parsers."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, ClassVar, NotRequired, TypedDict

from allelix.utils.build_detect import normalize_build_label

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from allelix.models import Variant


class GenotypeMetadata(TypedDict):
    """File-level metadata extracted by `GenotypeParser.get_metadata`.

    Header-derivable fields only. SNP count is intentionally NOT here — the
    only reliable source for it is `parse()`, since lines that look like data
    in a header scan may fail validation. Callers that need a count should
    use `sum(1 for _ in parser.parse(file_path))`.

    Keys:
        format: Parser name (matches `GenotypeParser.name`).
        sample_id: Vendor sample identifier, or "" if not present in the file.
        build: The build the file itself states (a header line such as
            ``# Build 37`` or ``##reference=build37``, a VCF ``assembly=``),
            normalized to ``GRCh36`` / ``GRCh37`` / ``GRCh38``, or None when
            the file states none. Parsers never fill in a default here: the
            format default is applied only by the build decision (ADR-0044).
        contig_lengths: ``{contig ID: length}`` from VCF ``##contig`` lines
            that declare a length. Optional — only VCF / gVCF parsers set it.
        method: How the genotypes were produced, when the file says so
            (e.g. ``"low-pass-wgs-imputed"`` for MyHeritage low-pass whole
            genome sequencing). Optional; absent for files that do not
            declare it.
    """

    format: str
    sample_id: str
    build: str | None
    contig_lengths: NotRequired[dict[str, int]]
    method: NotRequired[str]


# Comment lines that can carry a build statement. Other comment lines are
# skipped so a date (``2038-01-01``) or a count cannot read as a build.
_BUILD_MARKER_RE = re.compile(r"\b(build|reference|assembly|genome|grch\d\d|hg\d\d)\b", re.I)
HEADER_SCAN_LIMIT = 50


def header_build_label(lines: Iterable[str]) -> str | None:
    """The build stated in a file's leading ``#`` comment lines, or None.

    Reads until the first non-comment line (at most ``HEADER_SCAN_LIMIT``
    lines). The first line that names a build marker and normalizes to one
    build wins, e.g. ``# Build 37``, ``# Human Genome Reference Build 37
    (GRCh37.p13).``, ``#... human reference build 37.1 coordinates``.
    """
    for count, raw in enumerate(lines):
        line = raw.rstrip("\r\n").lstrip("\ufeff")
        if count >= HEADER_SCAN_LIMIT or (line and not line.startswith("#")):
            break
        build = build_statement(line)
        if build is not None:
            return build
    return None


def build_statement(line: str) -> str | None:
    """The build one header line states, or None (marker word required)."""
    if not _BUILD_MARKER_RE.search(line):
        return None
    return normalize_build_label(line)


def read_header_build(file_path: Path) -> str | None:
    """``header_build_label`` over ``file_path``; None when it cannot be read."""
    try:
        with file_path.open("r", encoding="utf-8", errors="replace") as fh:
            return header_build_label(fh)
    except OSError:
        return None


class GenotypeParser(ABC):
    """Base class for all genotype file parsers.

    Subclasses define metadata as class attributes and implement the three
    abstract methods. Parsers are stateless — `can_parse` and `parse` may be
    called repeatedly on different files.

    Attributes:
        name: Lowercase identifier used by the registry and CLI (e.g., "myhappygenes").
        display_name: Human-readable name for reports ("MyHappyGenes (Tempus)").
        file_extensions: Common file extensions (e.g., [".txt"]). Informational only;
            auto-detection uses `can_parse`, not extension matching.
        url: Vendor URL.
    """

    name: ClassVar[str]
    display_name: ClassVar[str]
    file_extensions: ClassVar[list[str]]
    url: ClassVar[str]

    @abstractmethod
    def can_parse(self, file_path: Path) -> bool:
        """Sniff the file to determine if this parser handles it.

        Must be fast — examines header/structural lines only, not the full file.
        Used by the auto-detection registry.

        Args:
            file_path: Path to the candidate genotype file.

        Returns:
            True if this parser recognizes the format.
        """
        ...

    @abstractmethod
    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Yield normalized Variant objects from the file.

        Streaming: yields one variant at a time. Never loads the whole file.
        Malformed individual lines log a warning and are skipped — they do not
        abort the whole parse.

        Args:
            file_path: Path to the genotype file.

        Yields:
            One Variant per data row in the file.
        """
        ...

    @abstractmethod
    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Extract header-derivable file metadata. Must be cheap (no full parse).

        Args:
            file_path: Path to the genotype file.

        Returns:
            A `GenotypeMetadata` dict.
        """
        ...

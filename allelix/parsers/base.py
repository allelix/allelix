# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Abstract base class for genotype file parsers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, ClassVar, NotRequired, TypedDict

if TYPE_CHECKING:
    from collections.abc import Iterator
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
        build: Reference genome build (e.g., "GRCh37").
        chr_prefix_observed: True when the file uses ``chr``-prefixed contig
            names (``chr1`` / ``chrX``). Optional — only VCF / gVCF parsers
            populate it; consumer-array exports always use bare names and
            don't set it. GH #38: used as a tertiary build-detection signal
            after rsID matching and ``##assembly`` tag checks both fail.
            ``chr``-prefixed contigs overwhelmingly indicate GRCh38.
    """

    format: str
    sample_id: str
    build: str
    chr_prefix_observed: NotRequired[bool]


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

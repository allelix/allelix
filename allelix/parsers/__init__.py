# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser registry. Auto-detection tries each registered parser; first match wins."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.parsers.ancestrydna import AncestryDNAParser
from allelix.parsers.base import GenotypeParser
from allelix.parsers.ftdna import FTDNAParser
from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser
from allelix.parsers.livingdna import LivingDNAParser
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.myheritage import MyHeritageParser
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser

if TYPE_CHECKING:
    from pathlib import Path

PARSERS: list[GenotypeParser] = [
    MyHappyGenesParser(),
    TwentyThreeAndMeParser(),
    AncestryDNAParser(),
    LivingDNAParser(),
    MyHeritageParser(),
    FTDNAParser(),
    FTDNAIlluminaParser(),
    VcfParser(),
]


class ParserNotFoundError(ValueError):
    """Raised when no parser can handle a file or a named parser does not exist."""


def get_parser_by_name(name: str) -> GenotypeParser:
    """Look up a parser by its `name` attribute.

    Args:
        name: Lowercase parser identifier (e.g., "myhappygenes").

    Raises:
        ParserNotFoundError: If no registered parser has that name.
    """
    for parser in PARSERS:
        if parser.name == name:
            return parser
    available = ", ".join(p.name for p in PARSERS)
    raise ParserNotFoundError(f"Unknown parser {name!r}. Available: {available}")


def detect_parser(file_path: Path) -> GenotypeParser:
    """Auto-detect the parser for a file. First match wins.

    Args:
        file_path: Path to the genotype file.

    Raises:
        ParserNotFoundError: If no parser recognizes the format.
    """
    for parser in PARSERS:
        if parser.can_parse(file_path):
            return parser
    raise ParserNotFoundError(
        f"No parser recognized {file_path.name!r}. Try forcing a format with --format <name>."
    )


__all__ = [
    "PARSERS",
    "GenotypeParser",
    "ParserNotFoundError",
    "detect_parser",
    "get_parser_by_name",
]

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser registry. Auto-detection tries each registered parser; first match wins."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.parsers._helpers import (
    DamagedFileError,
    InputContainerError,
    check_container,
    quiet_input_notices,
)
from allelix.parsers.ancestrydna import AncestryDNAParser
from allelix.parsers.base import GenotypeParser
from allelix.parsers.decodeme import DeCODEmeParser
from allelix.parsers.ftdna import FTDNAParser
from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser
from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser
from allelix.parsers.genesforgood import GenesForGoodParser
from allelix.parsers.livingdna import LivingDNAParser
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.myheritage import MyHeritageParser
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser

if TYPE_CHECKING:
    from pathlib import Path

PARSERS: list[GenotypeParser] = [
    MyHappyGenesParser(),
    GenesForGoodParser(),
    TwentyThreeAndMeParser(),
    AncestryDNAParser(),
    LivingDNAParser(),
    DeCODEmeParser(),
    MyHeritageParser(),
    FTDNAParser(),
    FTDNAFamFinderParser(),  # 5-col tab; must precede FTDNAIlluminaParser
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
        DamagedFileError: A damaged zip download.
        InputContainerError: A zip with several files, a password or an
            unreadable compression method.
    """
    check_container(file_path)
    with quiet_input_notices():
        for parser in PARSERS:
            if parser.can_parse(file_path):
                return parser
    raise ParserNotFoundError(
        f"Allelix doesn't recognise the format of {file_path.name!r}. Use the original "
        "file exactly as your testing company provided it, not edited or re-saved. If "
        "you're sure of its source, --format <name> reads it anyway, at your own risk: "
        "Allelix can't check a file it doesn't recognise, and the results may be wrong. "
        'See "Forcing a format" in docs/formats.md.'
    )


__all__ = [
    "PARSERS",
    "DamagedFileError",
    "GenotypeParser",
    "InputContainerError",
    "ParserNotFoundError",
    "detect_parser",
    "get_parser_by_name",
]

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Array files wrapped in gzip, a byte order mark or a leading blank line read as the plain file.

Each wrapped copy is built in ``tmp_path`` from a committed synthetic
fixture, changed in one way only. Detection, the parsed variants and the
header metadata must all match the plain file's.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from allelix.parsers import PARSERS, detect_parser
from allelix.parsers._helpers import first_nonblank_line, open_genotype_text

if TYPE_CHECKING:
    from collections.abc import Callable

    from allelix.parsers.base import GenotypeParser

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"

# (fixture file, parser name it must be detected as)
ARRAY_FIXTURES: list[tuple[str, str]] = [
    ("mock_23andme.txt", "23andme"),
    ("mock_23andme_v5_preamble.txt", "23andme"),
    ("mock_ancestrydna.txt", "ancestrydna"),
    ("mock_ftdna.csv", "ftdna"),
    ("mock_ftdna_illumina.txt", "ftdna_illumina"),
    ("mock_ftdna_famfinder.txt", "ftdna_famfinder"),
    ("mock_livingdna.csv", "livingdna"),
    ("mock_decodeme.csv", "decodeme"),
    ("mock_myheritage.csv", "myheritage"),
    ("myheritage_lowpass_wgs.csv", "myheritage"),
    ("mock_myhappygenes.txt", "myhappygenes"),
    ("mock_genesforgood.txt", "genesforgood"),
    ("mock_genesforgood_imputed.txt", "genesforgood"),
]

BOM = b"\xef\xbb\xbf"


def _bom(data: bytes) -> bytes:
    return BOM + data


def _gzip(data: bytes) -> bytes:
    return gzip.compress(data)


def _blank_first_line(data: bytes) -> bytes:
    return b"\n" + data


def _bom_inside_gzip(data: bytes) -> bytes:
    return gzip.compress(BOM + data)


WRAPPINGS: dict[str, Callable[[bytes], bytes]] = {
    "bom": _bom,
    "gzip": _gzip,
    "blank_first_line": _blank_first_line,
    "bom_inside_gzip": _bom_inside_gzip,
}


def _wrapped_copy(tmp_path: Path, fixture: str, wrapping: str) -> Path:
    plain = FIXTURES_DIR / fixture
    out = tmp_path / f"{wrapping}_{fixture}"
    out.write_bytes(WRAPPINGS[wrapping](plain.read_bytes()))
    return out


@pytest.mark.parametrize("wrapping", sorted(WRAPPINGS))
@pytest.mark.parametrize(("fixture", "parser_name"), ARRAY_FIXTURES)
class TestWrappedArrayFile:
    def test_plain_fixture_detects_as_expected(
        self, fixture: str, parser_name: str, wrapping: str
    ) -> None:
        # Guards the table itself: the plain file must already detect as listed.
        assert detect_parser(FIXTURES_DIR / fixture).name == parser_name

    def test_detected_as_same_parser(
        self, tmp_path: Path, fixture: str, parser_name: str, wrapping: str
    ) -> None:
        wrapped = _wrapped_copy(tmp_path, fixture, wrapping)
        assert detect_parser(wrapped).name == parser_name

    def test_parse_yields_same_variants(
        self, tmp_path: Path, fixture: str, parser_name: str, wrapping: str
    ) -> None:
        plain = FIXTURES_DIR / fixture
        wrapped = _wrapped_copy(tmp_path, fixture, wrapping)
        parser = detect_parser(plain)
        expected = list(parser.parse(plain))
        assert expected
        assert list(parser.parse(wrapped)) == expected

    def test_metadata_matches(
        self, tmp_path: Path, fixture: str, parser_name: str, wrapping: str
    ) -> None:
        plain = FIXTURES_DIR / fixture
        wrapped = _wrapped_copy(tmp_path, fixture, wrapping)
        parser = detect_parser(plain)
        assert parser.get_metadata(wrapped) == parser.get_metadata(plain)


class TestOpenGenotypeText:
    def test_plain_text(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_bytes(b"rs1\t1\t100\tAA\n")
        with open_genotype_text(path) as fh:
            assert fh.read() == "rs1\t1\t100\tAA\n"

    def test_bom_is_dropped(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_bytes(BOM + b"#header\n")
        with open_genotype_text(path) as fh:
            assert fh.readline() == "#header\n"

    def test_gzip_sniffed_whatever_the_extension(self, tmp_path: Path) -> None:
        path = tmp_path / "f.csv"
        path.write_bytes(gzip.compress(BOM + b"RSID,CHROMOSOME,POSITION,RESULT\n"))
        with open_genotype_text(path) as fh:
            assert fh.readline() == "RSID,CHROMOSOME,POSITION,RESULT\n"

    def test_strict_by_default(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_bytes(b"#\xa9 header\n")
        with open_genotype_text(path) as fh, pytest.raises(UnicodeDecodeError):
            fh.read()

    def test_errors_replace_on_request(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_bytes(b"#\xa9 header\n")
        with open_genotype_text(path, errors="replace") as fh:
            assert fh.read() == "#� header\n"


class TestFirstNonblankLine:
    def test_skips_blank_and_whitespace_lines(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_text("\n  \n\t\n#AncestryDNA\n", encoding="utf-8")
        with open_genotype_text(path) as fh:
            assert first_nonblank_line(fh, 50) == "#AncestryDNA\n"

    def test_empty_when_limit_reached(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_text("\n\n\n#AncestryDNA\n", encoding="utf-8")
        with open_genotype_text(path) as fh:
            assert first_nonblank_line(fh, 3) == ""

    def test_empty_file(self, tmp_path: Path) -> None:
        path = tmp_path / "f.txt"
        path.write_text("", encoding="utf-8")
        with open_genotype_text(path) as fh:
            assert first_nonblank_line(fh, 50) == ""


@pytest.mark.parametrize("parser", [p for p in PARSERS if p.name != "vcf"], ids=lambda p: p.name)
def test_cut_gzip_is_not_claimed(tmp_path: Path, parser: GenotypeParser) -> None:
    """A gzip stream cut before its first line ends: detection says no, it does not raise."""
    data = gzip.compress((FIXTURES_DIR / "mock_ftdna.csv").read_bytes())
    path = tmp_path / "cut.csv.gz"
    path.write_bytes(data[:12])
    assert parser.can_parse(path) is False

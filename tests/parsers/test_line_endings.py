# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Every input format reads LF, CRLF and CR-only line endings, plain, gzipped and zipped.

Each copy is built in ``tmp_path`` from a committed synthetic fixture: its
line endings are first made LF, then rewritten to the ending under test,
then optionally wrapped (gzip, a one-file zip, or a gzipped file in a zip).
Detection, the parsed variants and the header metadata must all match the
LF copy's.
"""

from __future__ import annotations

import gzip
import io
import re
import zipfile
from pathlib import Path

import pytest

from allelix.parsers import detect_parser

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"

# (fixture file, parser name it must be detected as)
FIXTURES: list[tuple[str, str]] = [
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
    ("mock_vcf.vcf", "vcf"),
    ("mock_gvcf.g.vcf", "vcf"),
]

ENDINGS: dict[str, bytes] = {"lf": b"\n", "crlf": b"\r\n", "cr": b"\r"}


def _with_ending(data: bytes, ending: bytes) -> bytes:
    lines = re.split(rb"\r\n|\r|\n", data)
    return ending.join(lines)


WRAPPINGS = ["plain", "gzip", "zip", "zip-gzip"]


def _zipped(name: str, data: bytes) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(name, data)
    return buffer.getvalue()


def _copy(tmp_path: Path, fixture: str, ending: str, wrapping: str = "plain") -> Path:
    data = _with_ending((FIXTURES_DIR / fixture).read_bytes(), ENDINGS[ending])
    name = f"{ending}_{fixture}"
    if wrapping in ("gzip", "zip-gzip"):
        data = gzip.compress(data)
        name += ".gz"
    if wrapping in ("zip", "zip-gzip"):
        data = _zipped(name, data)
        name += ".zip"
    out = tmp_path / name
    out.write_bytes(data)
    return out


@pytest.mark.parametrize("wrapping", WRAPPINGS)
@pytest.mark.parametrize("ending", sorted(ENDINGS))
@pytest.mark.parametrize(("fixture", "parser_name"), FIXTURES)
class TestLineEndings:
    def test_copy_has_only_that_ending(
        self, tmp_path: Path, fixture: str, parser_name: str, ending: str, wrapping: str
    ) -> None:
        # Guards the matrix itself: the copy really carries the ending under test.
        data = _copy(tmp_path, fixture, ending).read_bytes()
        assert data.count(ENDINGS[ending]) >= 3
        rest = data.replace(ENDINGS[ending], b"")
        assert b"\r" not in rest
        assert b"\n" not in rest

    def test_detected_as_expected(
        self, tmp_path: Path, fixture: str, parser_name: str, ending: str, wrapping: str
    ) -> None:
        path = _copy(tmp_path, fixture, ending, wrapping)
        assert detect_parser(path).name == parser_name

    def test_parse_yields_same_variants(
        self, tmp_path: Path, fixture: str, parser_name: str, ending: str, wrapping: str
    ) -> None:
        lf_dir = tmp_path / "lf"
        lf_dir.mkdir()
        baseline = _copy(lf_dir, fixture, "lf")
        path = _copy(tmp_path, fixture, ending, wrapping)
        parser = detect_parser(baseline)
        expected = list(parser.parse(baseline))
        assert expected
        assert list(parser.parse(path)) == expected

    def test_metadata_matches(
        self, tmp_path: Path, fixture: str, parser_name: str, ending: str, wrapping: str
    ) -> None:
        lf_dir = tmp_path / "lf"
        lf_dir.mkdir()
        baseline = _copy(lf_dir, fixture, "lf")
        path = _copy(tmp_path, fixture, ending, wrapping)
        parser = detect_parser(baseline)
        assert parser.get_metadata(path) == parser.get_metadata(baseline)

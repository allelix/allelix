# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the deCODEme parser.

The fixture is synthetic: the real column line, invented rows, and four
known-SNP rows at their GRCh36 positions so the build check has evidence.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from allelix.parsers import detect_parser
from allelix.parsers.decodeme import DeCODEmeParser
from allelix.reports._pipeline import run_analysis
from allelix.utils.build_detect import detect_build

FIXTURE = Path(__file__).parent.parent / "fixtures" / "mock_decodeme.csv"


@pytest.fixture
def parser() -> DeCODEmeParser:
    return DeCODEmeParser()


def _calls(parser: DeCODEmeParser, path: Path = FIXTURE) -> dict[str, tuple[str, str, str]]:
    return {v.rsid: (v.chromosome, v.allele1, v.allele2) for v in parser.parse(path)}


class TestDetection:
    def test_detected(self) -> None:
        assert detect_parser(FIXTURE).name == "decodeme"

    def test_column_line_case_and_spaces_ignored(
        self, parser: DeCODEmeParser, tmp_path: Path
    ) -> None:
        path = tmp_path / "scan.csv"
        path.write_text(
            "name, variation ,CHROMOSOME,Position,Strand,yourcode\nrs1,A/G,1,100,+,AG\n"
        )
        assert parser.can_parse(path)

    @pytest.mark.parametrize(
        "fixture",
        ["mock_23andme.txt", "mock_ftdna.csv", "mock_myheritage.csv", "mock_livingdna.csv"],
    )
    def test_other_formats_not_claimed(self, parser: DeCODEmeParser, fixture: str) -> None:
        assert not parser.can_parse(FIXTURE.parent / fixture)


class TestStrand:
    def test_plus_rows_are_kept(self, parser: DeCODEmeParser) -> None:
        calls = _calls(parser)
        assert calls["rs9200001"] == ("1", "A", "G")
        assert calls["rs1801133"] == ("1", "C", "T")

    def test_minus_rows_are_complemented_to_plus(self, parser: DeCODEmeParser) -> None:
        calls = _calls(parser)
        assert calls["rs1801131"] == ("1", "A", "C")  # file: TG on the minus strand
        assert calls["rs9200002"] == ("2", "C", "T")  # file: GA
        assert calls["rs9200006"] == ("Y", "A", "A")  # file: TT
        assert calls["rs9200004"] == ("7", "T", "A")  # palindromic: file AT

    def test_no_call_stays_a_no_call(self, parser: DeCODEmeParser) -> None:
        assert _calls(parser)["rs9200003"] == ("3", "-", "-")

    def test_bad_strand_is_skipped(self, parser: DeCODEmeParser, tmp_path: Path) -> None:
        path = tmp_path / "scan.csv"
        path.write_text(
            "Name,Variation,Chromosome,Position,Strand,YourCode\n"
            "rs1,A/G,1,100,?,AG\nrs2,A/G,1,200,+,AG\nrs3,A/G,1,x,+,AG\nrs4,A/G,1\n"
        )
        assert list(_calls(parser, path)) == ["rs2"]


class TestChromosomesAndBuild:
    def test_m_is_mt(self, parser: DeCODEmeParser) -> None:
        assert _calls(parser)["rs9200007"][0] == "MT"

    def test_no_build_stated(self, parser: DeCODEmeParser) -> None:
        assert parser.get_metadata(FIXTURE) == {
            "format": "decodeme",
            "sample_id": "",
            "build": None,
        }

    def test_seeds_decide_grch36(self, parser: DeCODEmeParser) -> None:
        seeds = detect_build(parser.parse(FIXTURE))
        assert (seeds.build, seeds.is_confident) == ("GRCh36", True)

    def test_run_takes_the_grch36_path(self, parser: DeCODEmeParser) -> None:
        result = run_analysis(FIXTURE, parser, [])
        assert result.build == "GRCh36"
        assert result.build_diagnostics is not None
        assert result.build_diagnostics.source == "seeds"

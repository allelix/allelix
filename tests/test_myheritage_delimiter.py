# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""MyHeritage exports re-saved through a spreadsheet.

A spreadsheet round-trip pads every row with trailing commas, can push the
``MyHeritage`` signature off line 1 (the first row becomes ``,,,``) and can
drop the ``#`` from comment lines. A "save as text" round-trip writes
tab-delimited rows. All fixtures are header plus synthetic rows.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from allelix.parsers import PARSERS, detect_parser
from allelix.parsers._helpers import split_csv_line
from allelix.parsers.ftdna import FTDNAParser
from allelix.parsers.myheritage import MyHeritageParser

FIXTURES = Path(__file__).resolve().parent / "fixtures"

_EXPECTED = [
    ("rs4477212", "1", 82154, "A", "A"),
    ("rs3094315", "1", 752566, "A", "G"),
    ("rs9001001", "1", 100000, "-", "-"),
]


def _calls(path: Path) -> list[tuple]:
    return [
        (v.rsid, v.chromosome, v.position, v.allele1, v.allele2)
        for v in MyHeritageParser().parse(path)
    ]


@pytest.mark.parametrize(
    "fixture",
    [
        "myheritage_spreadsheet_resaved.csv",  # signature on line 2, '#' dropped
        "myheritage_spreadsheet_resaved_hashed.csv",  # '#' kept, every row padded
        "myheritage_tab_delimited.txt",  # tab-delimited rows
        "mock_myheritage.csv",  # the native export, unchanged
    ],
)
class TestVariants:
    def test_detected_as_myheritage_only(self, fixture: str) -> None:
        """No other parser claims it, so routing never depends on registry order."""
        path = FIXTURES / fixture
        assert MyHeritageParser().can_parse(path)
        assert not FTDNAParser().can_parse(path)
        assert [p.name for p in PARSERS if p.can_parse(path)] == ["myheritage"]
        assert detect_parser(path).name == "myheritage"

    def test_parses(self, fixture: str) -> None:
        calls = _calls(FIXTURES / fixture)
        if fixture == "mock_myheritage.csv":
            assert calls  # native fixture content is pinned elsewhere
        else:
            assert calls == _EXPECTED


class TestBoundaries:
    def test_native_ftdna_still_routes_to_ftdna(self) -> None:
        path = FIXTURES / "mock_ftdna.csv"
        assert not MyHeritageParser().can_parse(path)
        assert detect_parser(path).name == "ftdna"

    def test_native_ftdna_illumina_still_routes_to_it(self) -> None:
        path = FIXTURES / "mock_ftdna_illumina.txt"
        assert not MyHeritageParser().can_parse(path)
        assert detect_parser(path).name == "ftdna_illumina"

    def test_signature_after_the_data_is_not_claimed(self, tmp_path: Path) -> None:
        path = tmp_path / "x.csv"
        path.write_text('RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n# MyHeritage\n')
        assert not MyHeritageParser().can_parse(path)

    def test_padding_beyond_the_columns_is_stripped_but_an_empty_genotype_is_kept(
        self,
    ) -> None:
        assert split_csv_line('"rs1","1","100","AA",,,', 4) == ["rs1", "1", "100", "AA"]
        assert split_csv_line('"rs1","1","100",""', 4) == ["rs1", "1", "100", ""]
        assert split_csv_line("rs1,1,100,", 4) == ["rs1", "1", "100", ""]

    def test_preamble_warnings_are_capped(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        preamble = "".join(f"Comment line {i} without a hash.,,,\n" for i in range(20))
        path = tmp_path / "long.csv"
        path.write_text(
            "MyHeritage DNA raw data.,,,\n"
            + preamble
            + 'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n'
        )
        with caplog.at_level("WARNING"):
            calls = list(MyHeritageParser().parse(path))
        assert len(calls) == 1
        assert sum("expected" in r.getMessage() for r in caplog.records) == 5

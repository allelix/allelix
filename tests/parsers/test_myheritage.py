# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the MyHeritage DNA parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.parsers.myheritage import MyHeritageParser

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def parser() -> MyHeritageParser:
    return MyHeritageParser()


def _write(tmp_path: Path, contents: str, name: str = "f.csv") -> Path:
    f = tmp_path / name
    f.write_text(contents, encoding="utf-8")
    return f


class TestParserAttributes:
    def test_required_metadata(self, parser: MyHeritageParser) -> None:
        assert parser.name == "myheritage"
        assert parser.display_name == "MyHeritage DNA"
        assert ".csv" in parser.file_extensions
        assert parser.url


class TestCanParse:
    def test_recognizes_real_fixture(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        assert parser.can_parse(mock_myheritage_path) is True

    def test_rejects_ftdna_file(self, parser: MyHeritageParser, mock_ftdna_path: Path) -> None:
        assert parser.can_parse(mock_ftdna_path) is False

    def test_rejects_23andme_file(self, parser: MyHeritageParser, mock_23andme_path: Path) -> None:
        assert parser.can_parse(mock_23andme_path) is False

    def test_rejects_ancestrydna_file(
        self, parser: MyHeritageParser, mock_ancestrydna_path: Path
    ) -> None:
        assert parser.can_parse(mock_ancestrydna_path) is False

    def test_rejects_empty_file(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "")
        assert parser.can_parse(f) is False

    def test_rejects_missing_file(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        assert parser.can_parse(tmp_path / "does_not_exist.csv") is False

    def test_rejects_binary_file(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = tmp_path / "binary.csv"
        f.write_bytes(b"\x00\x01\x02\xff\xfe")
        assert parser.can_parse(f) is False


class TestParse:
    def test_yields_variants_from_fixture(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        assert len(variants) == 11

    def test_known_mthfr_variant(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        mthfr = next(v for v in variants if v.rsid == "rs1801133")
        assert mthfr.chromosome == "1"
        assert mthfr.position == 11856378
        assert mthfr.allele1 == "A"
        assert mthfr.allele2 == "G"
        assert mthfr.is_heterozygous

    def test_concatenated_genotype_split(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        comt = next(v for v in variants if v.rsid == "rs4680")
        assert comt.allele1 == "A"
        assert comt.allele2 == "G"
        assert comt.is_heterozygous

    def test_homozygous_variant(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        v = next(v for v in variants if v.rsid == "rs1065852")
        assert v.allele1 == "C"
        assert v.allele2 == "C"
        assert not v.is_heterozygous

    def test_no_calls_preserved(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        no_call = next(v for v in variants if v.rsid == "rs9001001")
        assert no_call.is_no_call
        assert no_call.allele1 == "-"
        assert no_call.allele2 == "-"

    def test_mt_chromosome_haploid(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        mt = next(v for v in variants if v.rsid == "rs9001002")
        assert mt.chromosome == "MT"
        assert mt.allele1 == "A"
        assert mt.allele2 == "A"

    def test_x_chromosome(self, parser: MyHeritageParser, mock_myheritage_path: Path) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        x = next(v for v in variants if v.rsid == "rs9001003")
        assert x.chromosome == "X"
        assert x.allele1 == "A"
        assert x.allele2 == "G"

    def test_y_chromosome_haploid(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        variants = list(parser.parse(mock_myheritage_path))
        y = next(v for v in variants if v.rsid == "rs9001004")
        assert y.chromosome == "Y"
        assert y.allele1 == "A"
        assert y.allele2 == "A"

    def test_strips_double_quotes(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"
        assert variants[0].position == 100

    def test_double_double_quoted_fields(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        """MyHeritage's "extra quotes" variant wraps fields in double-double-quotes."""
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n""rs1"",""1"",""100"",""AG""\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "G"

    def test_unquoted_data_rows(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,AG\nrs2,2,200,CT\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 2

    def test_skips_comment_lines(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# MyHeritage\n# extra\nRSID,CHROMOSOME,POSITION,RESULT\n"
            '# data comment\n"rs1","1","100","AG"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1

    def test_skips_blank_lines(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n\n"rs1","1","100","AG"\n\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1

    def test_skips_malformed_column_count(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n'
            "malformed\n"
            '"rs2","2","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs1", "rs2"]

    def test_skips_invalid_position(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","BAD","AG"\n'
            '"rs2","2","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs2"]

    def test_handles_crlf_line_endings(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\r\nRSID,CHROMOSOME,POSITION,RESULT\r\n"rs1","1","100","AG"\r\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].allele1 == "A"

    def test_streaming_returns_iterator(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        result = parser.parse(mock_myheritage_path)
        assert iter(result) is result or hasattr(result, "__next__")


class TestMetadata:
    def test_default_build_37(self, parser: MyHeritageParser, mock_myheritage_path: Path) -> None:
        meta = parser.get_metadata(mock_myheritage_path)
        assert meta["build"] == "GRCh37"
        assert meta["format"] == "myheritage"
        assert meta["sample_id"] == ""

    def test_metadata_has_no_snp_count_field(
        self, parser: MyHeritageParser, mock_myheritage_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_myheritage_path)
        assert "snp_count" not in meta


class TestEdgeCases:
    def test_header_only_file_yields_zero_variants(
        self, parser: MyHeritageParser, tmp_path: Path
    ) -> None:
        f = _write(tmp_path, "# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n")
        assert list(parser.parse(f)) == []

    def test_empty_genotype_is_no_call(self, parser: MyHeritageParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100",""\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call

    def test_three_char_genotype_is_no_call(
        self, parser: MyHeritageParser, tmp_path: Path
    ) -> None:
        f = _write(
            tmp_path,
            '# MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","ABC"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call


class TestAutoDetection:
    def test_registry_detects_myheritage(self, mock_myheritage_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_myheritage_path)
        assert parser.name == "myheritage"

    def test_registry_lookup_by_name(self) -> None:
        from allelix.parsers import get_parser_by_name

        parser = get_parser_by_name("myheritage")
        assert isinstance(parser, MyHeritageParser)

    def test_ftdna_not_confused_with_myheritage(self, mock_ftdna_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_ftdna_path)
        assert parser.name == "ftdna"

    def test_myheritage_detected_before_ftdna(self, mock_myheritage_path: Path) -> None:
        """MyHeritage shares FTDNA's header; registry order must prefer MyHeritage."""
        from allelix.parsers.ftdna import FTDNAParser

        ftdna = FTDNAParser()
        assert ftdna.can_parse(mock_myheritage_path) is True

        from allelix.parsers import detect_parser

        parser = detect_parser(mock_myheritage_path)
        assert parser.name == "myheritage"

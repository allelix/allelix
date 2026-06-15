# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the FTDNA parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.parsers.ftdna import FTDNAParser

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def parser() -> FTDNAParser:
    return FTDNAParser()


def _write(tmp_path: Path, contents: str, name: str = "f.csv") -> Path:
    f = tmp_path / name
    f.write_text(contents, encoding="utf-8")
    return f


class TestParserAttributes:
    def test_required_metadata(self, parser: FTDNAParser) -> None:
        assert parser.name == "ftdna"
        assert parser.display_name == "Family Tree DNA"
        assert ".csv" in parser.file_extensions
        assert parser.url


class TestCanParse:
    def test_recognizes_real_fixture(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        assert parser.can_parse(mock_ftdna_path) is True

    def test_recognizes_quoted_header(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, '"RSID","CHROMOSOME","POSITION","RESULT"\n"rs1","1","100","AG"\n')
        assert parser.can_parse(f) is True

    def test_recognizes_unquoted_header(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\n")
        assert parser.can_parse(f) is True

    def test_recognizes_header_after_comments(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# comment line\n# another\nRSID,CHROMOSOME,POSITION,RESULT\n",
        )
        assert parser.can_parse(f) is True

    def test_rejects_myhappygenes_file(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# MyHappyGenes [TEMPUS]\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tA\tG\n",
        )
        assert parser.can_parse(f) is False

    def test_rejects_23andme_file(self, parser: FTDNAParser, mock_23andme_path: Path) -> None:
        assert parser.can_parse(mock_23andme_path) is False

    def test_rejects_ancestrydna_file(
        self, parser: FTDNAParser, mock_ancestrydna_path: Path
    ) -> None:
        assert parser.can_parse(mock_ancestrydna_path) is False

    def test_rejects_empty_file(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "")
        assert parser.can_parse(f) is False

    def test_rejects_missing_file(self, parser: FTDNAParser, tmp_path: Path) -> None:
        assert parser.can_parse(tmp_path / "does_not_exist.csv") is False

    def test_rejects_binary_file(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = tmp_path / "binary.csv"
        f.write_bytes(b"\x00\x01\x02\xff\xfe")
        assert parser.can_parse(f) is False

    def test_rejects_wrong_header(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "NAME,CHROM,POS,GENO\nrs1,1,100,AG\n")
        assert parser.can_parse(f) is False

    def test_case_insensitive_header(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "rsid,chromosome,position,result\nrs1,1,100,AG\n")
        assert parser.can_parse(f) is True


class TestParse:
    def test_yields_variants_from_fixture(
        self, parser: FTDNAParser, mock_ftdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        assert len(variants) == 11

    def test_known_mthfr_variant(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        mthfr = next(v for v in variants if v.rsid == "rs1801133")
        assert mthfr.chromosome == "1"
        assert mthfr.position == 11856378
        assert mthfr.allele1 == "A"
        assert mthfr.allele2 == "G"
        assert mthfr.is_heterozygous

    def test_concatenated_genotype_split(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        comt = next(v for v in variants if v.rsid == "rs4680")
        assert comt.allele1 == "A"
        assert comt.allele2 == "G"
        assert comt.is_heterozygous

    def test_homozygous_variant(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        v = next(v for v in variants if v.rsid == "rs1065852")
        assert v.allele1 == "C"
        assert v.allele2 == "C"
        assert not v.is_heterozygous

    def test_no_calls_preserved(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        no_call = next(v for v in variants if v.rsid == "rs9001001")
        assert no_call.is_no_call
        assert no_call.allele1 == "-"
        assert no_call.allele2 == "-"

    def test_mt_chromosome_haploid(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        mt = next(v for v in variants if v.rsid == "rs9001002")
        assert mt.chromosome == "MT"
        assert mt.allele1 == "A"
        assert mt.allele2 == "A"

    def test_x_chromosome(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        x = next(v for v in variants if v.rsid == "rs9001003")
        assert x.chromosome == "X"
        assert x.allele1 == "A"
        assert x.allele2 == "G"

    def test_y_chromosome_haploid(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        variants = list(parser.parse(mock_ftdna_path))
        y = next(v for v in variants if v.rsid == "rs9001004")
        assert y.chromosome == "Y"
        assert y.allele1 == "A"
        assert y.allele2 == "A"

    def test_strips_double_quotes(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"
        assert variants[0].position == 100

    def test_unquoted_data_rows(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "RSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,AG\nrs2,2,200,CT\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 2
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "G"

    def test_skips_comment_lines(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# comment\nRSID,CHROMOSOME,POSITION,RESULT\n"
            '# another comment\n"rs1","1","100","AG"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_skips_blank_lines(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n\n"rs1","1","100","AG"\n\n"rs2","2","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 2

    def test_skips_malformed_column_count(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n'
            "this is malformed\n"
            '"rs2","2","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs1", "rs2"]

    def test_skips_invalid_position(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","NOT_A_NUMBER","AG"\n'
            '"rs2","2","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs2"]

    def test_handles_crlf_line_endings(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\r\n"rs1","1","100","AG"\r\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "G"

    def test_streaming_returns_iterator(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        result = parser.parse(mock_ftdna_path)
        assert iter(result) is result or hasattr(result, "__next__")


class TestMetadata:
    def test_default_build_37(self, parser: FTDNAParser, mock_ftdna_path: Path) -> None:
        meta = parser.get_metadata(mock_ftdna_path)
        assert meta["build"] == "GRCh37"
        assert meta["format"] == "ftdna"
        assert meta["sample_id"] == ""

    def test_metadata_has_no_snp_count_field(
        self, parser: FTDNAParser, mock_ftdna_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_ftdna_path)
        assert "snp_count" not in meta


class TestEdgeCases:
    def test_header_only_file_yields_zero_variants(
        self, parser: FTDNAParser, tmp_path: Path
    ) -> None:
        f = _write(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\n")
        assert parser.can_parse(f) is True
        assert list(parser.parse(f)) == []

    def test_empty_genotype_is_no_call(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100",""\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call

    def test_three_char_genotype_is_no_call(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","ABC"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call

    def test_double_double_quoted_fields(self, parser: FTDNAParser, tmp_path: Path) -> None:
        """Some exports wrap fields in extra quotes."""
        f = _write(
            tmp_path,
            'RSID,CHROMOSOME,POSITION,RESULT\n""rs1"","1","100","AG"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_no_header_in_file_rejects(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "# just a comment\n# no header line\n")
        assert parser.can_parse(f) is False

    def test_data_before_header_warns_and_skips(self, parser: FTDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            '"rs1","1","100","AG"\nRSID,CHROMOSOME,POSITION,RESULT\n"rs2","1","200","CT"\n',
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs2"


class TestAutoDetection:
    def test_registry_detects_ftdna(self, mock_ftdna_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_ftdna_path)
        assert parser.name == "ftdna"

    def test_registry_lookup_by_name(self) -> None:
        from allelix.parsers import get_parser_by_name

        parser = get_parser_by_name("ftdna")
        assert isinstance(parser, FTDNAParser)

    def test_mhg_not_confused_with_ftdna(self, mock_mhg_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_mhg_path)
        assert parser.name == "myhappygenes"

    def test_23andme_not_confused_with_ftdna(self, mock_23andme_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_23andme_path)
        assert parser.name == "23andme"

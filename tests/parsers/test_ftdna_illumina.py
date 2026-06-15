# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the FTDNA Illumina raw (tab-delimited) parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def parser() -> FTDNAIlluminaParser:
    return FTDNAIlluminaParser()


def _write(tmp_path: Path, contents: str, name: str = "f.txt") -> Path:
    f = tmp_path / name
    f.write_text(contents, encoding="utf-8")
    return f


class TestParserAttributes:
    def test_required_metadata(self, parser: FTDNAIlluminaParser) -> None:
        assert parser.name == "ftdna_illumina"
        assert parser.display_name == "Family Tree DNA (Illumina raw)"
        assert ".txt" in parser.file_extensions
        assert parser.url


class TestCanParse:
    def test_recognizes_real_fixture(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        assert parser.can_parse(mock_ftdna_illumina_path) is True

    def test_recognizes_header_after_comments(
        self, parser: FTDNAIlluminaParser, tmp_path: Path
    ) -> None:
        f = _write(
            tmp_path,
            "# Family Tree DNA\n# Build 37\nRSID\tCHROMOSOME\tPOSITION\tRESULT\nrs1\t1\t100\tAG\n",
        )
        assert parser.can_parse(f) is True

    def test_case_insensitive_header(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "rsid\tchromosome\tposition\tresult\nrs1\t1\t100\tAG\n")
        assert parser.can_parse(f) is True

    def test_rejects_csv_form(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        """Comma-delimited FTDNA (the other variant) must not match."""
        f = _write(tmp_path, "RSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,AG\n", "f.csv")
        assert parser.can_parse(f) is False

    def test_rejects_myhappygenes(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        """MHG also tab-delimited but with different columns. Must not match."""
        f = _write(
            tmp_path,
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\nrs1\t1\t100\tA\tG\n",
        )
        assert parser.can_parse(f) is False

    def test_rejects_random_text(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "Some random text\nNo header here\n")
        assert parser.can_parse(f) is False

    def test_rejects_empty_file(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "")
        assert parser.can_parse(f) is False


class TestParse:
    def test_yields_expected_variants(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        variants = list(parser.parse(mock_ftdna_illumina_path))
        rsids = [v.rsid for v in variants]
        assert "rs1801133" in rsids
        assert "rs4680" in rsids

    def test_het_genotype_split(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_illumina_path)}
        # rs1801133 is AG (heterozygous)
        v = by_rsid["rs1801133"]
        assert v.allele1 == "A"
        assert v.allele2 == "G"
        assert v.is_heterozygous

    def test_hom_genotype_split(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_illumina_path)}
        # rs429358 is TT (homozygous)
        v = by_rsid["rs429358"]
        assert v.allele1 == "T"
        assert v.allele2 == "T"
        assert not v.is_heterozygous

    def test_no_call_handled(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_illumina_path)}
        # rs9001001 is -- (no-call)
        v = by_rsid["rs9001001"]
        assert v.is_no_call

    def test_haploid_call_duplicated(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        """A single-char RESULT becomes a homozygous diploid (MT/Y convention)."""
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_illumina_path)}
        v = by_rsid["rs9001003"]
        assert v.allele1 == "A"
        assert v.allele2 == "A"
        assert v.chromosome == "MT"

    def test_default_build_is_grch37(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        variants = list(parser.parse(mock_ftdna_illumina_path))
        assert all(v.build == "GRCh37" for v in variants)

    def test_skips_invalid_position(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "RSID\tCHROMOSOME\tPOSITION\tRESULT\n"
            "rs_bad\t1\tNOT_A_NUMBER\tAG\n"
            "rs_good\t1\t100\tAG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs_good"

    def test_skips_wrong_column_count(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "RSID\tCHROMOSOME\tPOSITION\tRESULT\n"
            "rs_short\t1\t100\n"  # 3 cols
            "rs_long\t1\t100\tAG\tEXTRA\n"  # 5 cols
            "rs_good\t1\t100\tAG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs_good"

    def test_comments_skipped(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# A comment\n"
            "# Another\n"
            "RSID\tCHROMOSOME\tPOSITION\tRESULT\n"
            "# Inline comment after header\n"
            "rs1\t1\t100\tAG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_chromosome_normalized(self, parser: FTDNAIlluminaParser, tmp_path: Path) -> None:
        """Defensive normalize_chromosome consistency with the VCF parser.

        Real FTDNA Illumina files use plain numeric chromosomes today,
        but a non-canonical export (chr-prefix, lowercase 'chrm', etc.)
        would otherwise silently miss rsID lookups downstream. The
        parser passes every chromosome through normalize_chromosome
        defensively.
        """
        f = _write(
            tmp_path,
            "RSID\tCHROMOSOME\tPOSITION\tRESULT\n"
            "rs_chr_prefix\tchr1\t100\tAG\n"
            "rs_chr_x\tchrX\t200\tAG\n"
            "rs_chr_m\tchrm\t300\tA\n"
            "rs_canonical\t22\t400\tCT\n",
        )
        by_rsid = {v.rsid: v for v in parser.parse(f)}
        assert by_rsid["rs_chr_prefix"].chromosome == "1"
        assert by_rsid["rs_chr_x"].chromosome == "X"
        assert by_rsid["rs_chr_m"].chromosome == "MT"
        assert by_rsid["rs_canonical"].chromosome == "22"


class TestGetMetadata:
    def test_metadata_basics(
        self, parser: FTDNAIlluminaParser, mock_ftdna_illumina_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_ftdna_illumina_path)
        assert meta["format"] == "ftdna_illumina"
        assert meta["build"] == "GRCh37"
        assert meta["sample_id"] == ""


class TestRegistryIntegration:
    def test_registered_in_parsers_list(self) -> None:
        from allelix.parsers import PARSERS, get_parser_by_name

        names = [p.name for p in PARSERS]
        assert "ftdna_illumina" in names
        assert isinstance(get_parser_by_name("ftdna_illumina"), FTDNAIlluminaParser)

    def test_auto_detect_picks_ftdna_illumina(self, mock_ftdna_illumina_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_ftdna_illumina_path)
        assert isinstance(parser, FTDNAIlluminaParser)

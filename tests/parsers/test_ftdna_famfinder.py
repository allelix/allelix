# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the FTDNA FamFinder (tab-delimited, separate allele columns) parser."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def parser() -> FTDNAFamFinderParser:
    return FTDNAFamFinderParser()


def _write(tmp_path: Path, contents: str, name: str = "f.txt") -> Path:
    f = tmp_path / name
    f.write_text(contents, encoding="utf-8")
    return f


_FAMFINDER_HEADER = (
    "# Family Tree DNA - FamFinder\n# Build 37\nRSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n"
)


class TestParserAttributes:
    def test_required_metadata(self, parser: FTDNAFamFinderParser) -> None:
        assert parser.name == "ftdna_famfinder"
        assert parser.display_name == "Family Tree DNA (FamFinder)"
        assert ".txt" in parser.file_extensions
        assert parser.url


class TestCanParse:
    def test_recognizes_real_fixture(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        assert parser.can_parse(mock_ftdna_famfinder_path) is True

    def test_marker_case_insensitive(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# FAMFINDER export\n"
            "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n"
            "rs1\t1\t100\tA\tG\n",
        )
        assert parser.can_parse(f) is True

    def test_case_insensitive_header(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# famfinder\nrsid\tchromosome\tposition\tallele1\tallele2\nrs1\t1\t100\tA\tG\n",
        )
        assert parser.can_parse(f) is True

    def test_rejects_without_marker(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        """Header alone (no famfinder marker) must NOT match — protects against
        any future tab-delimited 5-column format that isn't a FamFinder file."""
        f = _write(
            tmp_path,
            "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\nrs1\t1\t100\tA\tG\n",
        )
        assert parser.can_parse(f) is False

    def test_rejects_ftdna_illumina_raw(
        self, parser: FTDNAFamFinderParser, mock_ftdna_illumina_path: Path
    ) -> None:
        """The 4-col concatenated-RESULT sibling must not match (mutually exclusive)."""
        assert parser.can_parse(mock_ftdna_illumina_path) is False

    def test_rejects_csv_ftdna(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        """Comma-delimited FTDNA standard variant must not match."""
        f = _write(
            tmp_path,
            '# Family Tree DNA - FamFinder\nRSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,"AG"\n',
            "f.csv",
        )
        # Even with the marker, the CSV body doesn't have the 5-col tab header.
        assert parser.can_parse(f) is False

    def test_rejects_random_text(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "Some random text\nNo header here\n")
        assert parser.can_parse(f) is False

    def test_rejects_empty_file(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "")
        assert parser.can_parse(f) is False


class TestParse:
    def test_yields_expected_variants(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        variants = list(parser.parse(mock_ftdna_famfinder_path))
        rsids = [v.rsid for v in variants]
        assert "rs1801133" in rsids
        assert "rs4680" in rsids

    def test_het_genotype_split(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_famfinder_path)}
        v = by_rsid["rs1801133"]
        assert v.allele1 == "A"
        assert v.allele2 == "G"
        assert v.is_heterozygous

    def test_hom_genotype_split(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_famfinder_path)}
        v = by_rsid["rs429358"]
        assert v.allele1 == "T"
        assert v.allele2 == "T"
        assert not v.is_heterozygous

    def test_full_no_call_handled(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_famfinder_path)}
        # rs9001001 is - / - (full no-call)
        v = by_rsid["rs9001001"]
        assert v.is_no_call

    def test_partial_no_call_handled(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        """A / - is a partial no-call (one allele missing)."""
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_famfinder_path)}
        v = by_rsid["rs9001002"]
        assert v.allele1 == "A"
        assert v.allele2 == "-"

    def test_haploid_on_mt(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        """Haploid MT call (A/A in the FamFinder convention) parses as homozygous."""
        by_rsid = {v.rsid: v for v in parser.parse(mock_ftdna_famfinder_path)}
        v = by_rsid["rs9001003"]
        assert v.allele1 == "A"
        assert v.allele2 == "A"
        assert v.chromosome == "MT"

    def test_default_build_is_grch37(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        variants = list(parser.parse(mock_ftdna_famfinder_path))
        assert all(v.build == "GRCh37" for v in variants)

    def test_skips_invalid_position(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER + "rs_bad\t1\tNOT_A_NUMBER\tA\tG\n" + "rs_good\t1\t100\tA\tG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs_good"

    def test_skips_wrong_column_count(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER
            + "rs_short\t1\t100\tA\n"  # 4 cols
            + "rs_long\t1\t100\tA\tG\tEXTRA\n"  # 6 cols
            + "rs_good\t1\t100\tA\tG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs_good"

    def test_comments_skipped(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# A comment with famfinder marker\n"
            "# Another\n"
            "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n"
            "# Inline comment after header\n"
            "rs1\t1\t100\tA\tG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_chromosome_normalized(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        """Defensive normalize_chromosome consistency with the sibling FTDNA parsers."""
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER
            + "rs_chr_prefix\tchr1\t100\tA\tG\n"
            + "rs_chr_x\tchrX\t200\tA\tG\n"
            + "rs_chr_m\tchrm\t300\tA\tA\n"
            + "rs_canonical\t22\t400\tC\tT\n",
        )
        by_rsid = {v.rsid: v for v in parser.parse(f)}
        assert by_rsid["rs_chr_prefix"].chromosome == "1"
        assert by_rsid["rs_chr_x"].chromosome == "X"
        assert by_rsid["rs_chr_m"].chromosome == "MT"
        assert by_rsid["rs_canonical"].chromosome == "22"

    def test_empty_allele_cell_treated_as_no_call(
        self, parser: FTDNAFamFinderParser, tmp_path: Path
    ) -> None:
        """An empty cell (rather than ``-``) in an allele column maps to ``-``."""
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER + "rs_empty\t1\t100\tA\t\n",
        )
        variants = list(parser.parse(f))
        # 5-col split with trailing tab gives "" — that's allele2 here.
        # The parser should normalize to "-".
        assert len(variants) == 1
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "-"


class TestGetMetadata:
    def test_metadata_basics(
        self, parser: FTDNAFamFinderParser, mock_ftdna_famfinder_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_ftdna_famfinder_path)
        assert meta["format"] == "ftdna_famfinder"
        assert meta["build"] == "GRCh37"
        assert meta["sample_id"] == ""


class TestRegistryIntegration:
    def test_registered_in_parsers_list(self) -> None:
        from allelix.parsers import PARSERS, get_parser_by_name

        names = [p.name for p in PARSERS]
        assert "ftdna_famfinder" in names
        assert isinstance(get_parser_by_name("ftdna_famfinder"), FTDNAFamFinderParser)

    def test_auto_detect_picks_famfinder(self, mock_ftdna_famfinder_path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_ftdna_famfinder_path)
        assert isinstance(parser, FTDNAFamFinderParser)

    def test_auto_detect_still_picks_illumina_for_illumina(
        self,
        mock_ftdna_illumina_path,
    ) -> None:
        """FamFinder registration must not steal Illumina-raw detection."""
        from allelix.parsers import detect_parser
        from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser

        parser = detect_parser(mock_ftdna_illumina_path)
        assert isinstance(parser, FTDNAIlluminaParser)

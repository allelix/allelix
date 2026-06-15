# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the Living DNA parser.

Reviewer test instructions
--------------------------
The synthetic fixture covers all documented edge cases (AX-, AFFX-,
CHR:POS SNP IDs, no-calls, X chromosome, haploid). For real-data
validation:

1. **openSNP archive** (20 GB, December 2017 snapshot):
   ``https://archive.org/details/opensnp_data_dumps``
   Living DNA had very few users pre-2017 — files may be absent or
   scarce. Search the extracted archive for files containing
   "Living DNA" in the first line::

       find opensnp/ -name "*.csv" -exec head -1 {} \\; | grep -l "Living DNA"

   If found, run::

       allelix stats <file>           # should detect as "Living DNA"
       allelix analyze <file>         # should parse and annotate

2. **snps package** (format reference fixtures):
   ``pip install snps`` then inspect ``snps/tests/input/`` for Living
   DNA sample files. These are BSD-licensed format references, not
   real genetic data.

3. **Manual format validation**: Living DNA files are tab-delimited
   despite ``.csv`` extension. Verify with::

       head -20 <file>               # comments start with #
       awk -F'\\t' 'NF!=4 && !/^#/' <file> | head  # should be empty

4. **Non-rs SNP IDs**: Real Living DNA files contain AX- and AFFX-
   prefixed probe IDs and CHR:POS positional notation. These should
   pass through as-is in the rsid field. Verify::

       allelix extract --snps AX-12345678 <file>  # should not crash
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.parsers.livingdna import LivingDNAParser

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def parser() -> LivingDNAParser:
    return LivingDNAParser()


def _write(tmp_path: Path, contents: str, name: str = "f.csv") -> Path:
    f = tmp_path / name
    f.write_text(contents, encoding="utf-8")
    return f


class TestParserAttributes:
    def test_required_metadata(self, parser: LivingDNAParser) -> None:
        assert parser.name == "livingdna"
        assert parser.display_name == "Living DNA"
        assert ".csv" in parser.file_extensions
        assert parser.url


class TestCanParse:
    def test_recognizes_real_fixture(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        assert parser.can_parse(mock_livingdna_path) is True

    def test_rejects_23andme_file(self, parser: LivingDNAParser, mock_23andme_path: Path) -> None:
        assert parser.can_parse(mock_23andme_path) is False

    def test_rejects_ftdna_file(self, parser: LivingDNAParser, mock_ftdna_path: Path) -> None:
        assert parser.can_parse(mock_ftdna_path) is False

    def test_rejects_ancestrydna_file(
        self, parser: LivingDNAParser, mock_ancestrydna_path: Path
    ) -> None:
        assert parser.can_parse(mock_ancestrydna_path) is False

    def test_rejects_myheritage_file(
        self, parser: LivingDNAParser, mock_myheritage_path: Path
    ) -> None:
        assert parser.can_parse(mock_myheritage_path) is False

    def test_rejects_empty_file(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "")
        assert parser.can_parse(f) is False

    def test_rejects_missing_file(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        assert parser.can_parse(tmp_path / "does_not_exist.csv") is False

    def test_rejects_binary_file(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = tmp_path / "binary.csv"
        f.write_bytes(b"\x00\x01\x02\xff\xfe")
        assert parser.can_parse(f) is False


class TestParse:
    def test_yields_variants_from_fixture(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        assert len(variants) == 12

    def test_known_mthfr_variant(self, parser: LivingDNAParser, mock_livingdna_path: Path) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        mthfr = next(v for v in variants if v.rsid == "rs1801133")
        assert mthfr.chromosome == "1"
        assert mthfr.position == 11856378
        assert mthfr.allele1 == "A"
        assert mthfr.allele2 == "G"
        assert mthfr.is_heterozygous

    def test_concatenated_genotype_split(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        comt = next(v for v in variants if v.rsid == "rs4680")
        assert comt.allele1 == "A"
        assert comt.allele2 == "G"
        assert comt.is_heterozygous

    def test_homozygous_variant(self, parser: LivingDNAParser, mock_livingdna_path: Path) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        v = next(v for v in variants if v.rsid == "rs1065852")
        assert v.allele1 == "C"
        assert v.allele2 == "C"
        assert not v.is_heterozygous

    def test_no_calls_preserved(self, parser: LivingDNAParser, mock_livingdna_path: Path) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        no_call = next(v for v in variants if v.rsid == "rs9001001")
        assert no_call.is_no_call
        assert no_call.allele1 == "-"
        assert no_call.allele2 == "-"

    def test_x_chromosome(self, parser: LivingDNAParser, mock_livingdna_path: Path) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        x = next(v for v in variants if v.rsid == "rs9001003")
        assert x.chromosome == "X"
        assert x.allele1 == "A"
        assert x.allele2 == "G"

    def test_ax_prefixed_probe_id(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        ax = next(v for v in variants if v.rsid == "AX-12345678")
        assert ax.chromosome == "3"
        assert ax.allele1 == "G"
        assert ax.allele2 == "G"

    def test_affx_prefixed_probe_id(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        affx = next(v for v in variants if v.rsid == "AFFX-SP-000001")
        assert affx.chromosome == "7"
        assert affx.allele1 == "C"
        assert affx.allele2 == "T"
        assert affx.is_heterozygous

    def test_positional_notation_snp_id(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        variants = list(parser.parse(mock_livingdna_path))
        pos = next(v for v in variants if v.rsid == "1:726912")
        assert pos.chromosome == "1"
        assert pos.position == 726912
        assert pos.allele1 == "A"
        assert pos.allele2 == "A"

    def test_skips_comment_lines(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\n# extra\n# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_skips_blank_lines(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\n# rsid\tchromosome\tposition\tgenotype\n"
            "\nrs1\t1\t100\tAG\n\nrs2\t2\t200\tCT\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 2

    def test_skips_malformed_column_count(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\nrs1\t1\t100\tAG\nthis is malformed\nrs2\t2\t200\tCT\n",
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs1", "rs2"]

    def test_skips_invalid_position(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\nrs1\t1\tNOT_A_NUMBER\tAG\nrs2\t2\t200\tCT\n",
        )
        variants = list(parser.parse(f))
        assert [v.rsid for v in variants] == ["rs2"]

    def test_handles_crlf_line_endings(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\r\nrs1\t1\t100\tAG\r\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "G"

    def test_streaming_returns_iterator(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        result = parser.parse(mock_livingdna_path)
        assert iter(result) is result or hasattr(result, "__next__")


class TestMetadata:
    def test_extracts_build_37_from_header(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_livingdna_path)
        assert meta["build"] == "GRCh37"
        assert meta["format"] == "livingdna"
        assert meta["sample_id"] == ""

    def test_detects_build_38(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\n# Human Genome Reference Build 38 (GRCh38).\nrs1\t1\t100\tAG\n",
        )
        meta = parser.get_metadata(f)
        assert meta["build"] == "GRCh38"

    def test_defaults_to_grch37(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "# Living DNA\nrs1\t1\t100\tAG\n")
        meta = parser.get_metadata(f)
        assert meta["build"] == "GRCh37"

    def test_date_comment_does_not_override_build(
        self, parser: LivingDNAParser, tmp_path: Path
    ) -> None:
        """GH #16: a `2038-…` download-date comment must not silently
        retag the file as GRCh38. Only build-marker lines (containing
        ``build``, ``reference``, or ``genome``) are inspected, and the
        normalizer rejects digits inside year strings via word
        boundaries."""
        f = _write(
            tmp_path,
            "# Living DNA\n"
            "# downloaded 2038-01-01\n"
            "# Human Genome Reference Build 37 (GRCh37.p13).\n"
            "rs1\t1\t100\tAG\n",
        )
        meta = parser.get_metadata(f)
        assert meta["build"] == "GRCh37"

    def test_first_build_marker_line_wins(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        """GH #16: previously the LAST normalize-match won (no break),
        so a stray later comment could overwrite the real header. Now
        the first valid build-marker line wins."""
        f = _write(
            tmp_path,
            "# Living DNA\n"
            "# Human Genome Reference Build 37 (GRCh37.p13).\n"
            "# Bogus trailer: build 38\n"  # would have overwritten before
            "rs1\t1\t100\tAG\n",
        )
        meta = parser.get_metadata(f)
        assert meta["build"] == "GRCh37"

    def test_metadata_has_no_snp_count_field(
        self, parser: LivingDNAParser, mock_livingdna_path: Path
    ) -> None:
        meta = parser.get_metadata(mock_livingdna_path)
        assert "snp_count" not in meta


class TestEdgeCases:
    def test_header_only_file_yields_zero_variants(
        self, parser: LivingDNAParser, tmp_path: Path
    ) -> None:
        f = _write(
            tmp_path,
            "# Living DNA\n# rsid\tchromosome\tposition\tgenotype\n",
        )
        assert parser.can_parse(f) is True
        assert list(parser.parse(f)) == []

    def test_empty_genotype_is_no_call(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "# Living DNA\nrs1\t1\t100\t\n")
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call

    def test_three_char_genotype_is_no_call(self, parser: LivingDNAParser, tmp_path: Path) -> None:
        f = _write(tmp_path, "# Living DNA\nrs1\t1\t100\tABC\n")
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].is_no_call


class TestAutoDetection:
    def test_registry_detects_livingdna(self, mock_livingdna_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_livingdna_path)
        assert parser.name == "livingdna"

    def test_registry_lookup_by_name(self) -> None:
        from allelix.parsers import get_parser_by_name

        parser = get_parser_by_name("livingdna")
        assert isinstance(parser, LivingDNAParser)

    def test_23andme_not_confused_with_livingdna(self, mock_23andme_path: Path) -> None:
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_23andme_path)
        assert parser.name == "23andme"

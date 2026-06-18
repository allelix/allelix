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
        """GH #113: 4 columns is the haploid shape and is now accepted
        (see test_haploid_4_column_line). Genuinely wrong column counts
        — 3 columns (missing both alleles) or 6+ (extra data) — are
        still skipped."""
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER
            + "rs_too_short\t1\t100\n"  # 3 cols — both alleles missing
            + "rs_long\tMT\t100\tA\tG\tEXTRA\n"  # 6 cols
            + "rs_good\t1\t100\tA\tG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs_good"

    def test_haploid_4_column_line(self, parser: FTDNAFamFinderParser, tmp_path: Path) -> None:
        """GH #113: a 4-column line (missing ALLELE2 column entirely)
        is a haploid MT/Y call in the FamFinder convention. Parser
        emits a hemizygous Variant where allele1 == allele2 == the
        called base — NOT a half no-call.

        Pre-#113 the parser warned and skipped the line, dropping
        every haploid MT/Y call.
        """
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER
            + "rs_haploid_mt\tMT\t100\tA\n"  # 4 cols, real call
            + "rs_haploid_y\tY\t200\tC\n"  # 4 cols, real call
            + "rs_haploid_nc\tMT\t300\t-\n"  # 4 cols, no-call
            + "rs_diploid\t1\t400\tA\tG\n",
        )
        by_rsid = {v.rsid: v for v in parser.parse(f)}
        assert len(by_rsid) == 4
        assert by_rsid["rs_haploid_mt"].allele1 == "A"
        assert by_rsid["rs_haploid_mt"].allele2 == "A"
        assert by_rsid["rs_haploid_mt"].chromosome == "MT"
        assert by_rsid["rs_haploid_y"].allele1 == "C"
        assert by_rsid["rs_haploid_y"].allele2 == "C"
        # No-call shape: both alleles end up as the no-call marker via
        # haploid-doubling of "-". Variant.is_no_call must still fire.
        assert by_rsid["rs_haploid_nc"].is_no_call
        # Diploid lines on the same file still work normally.
        assert by_rsid["rs_diploid"].allele1 == "A"
        assert by_rsid["rs_diploid"].allele2 == "G"

    def test_haploid_shape_on_autosome_skipped_as_truncated_diploid(
        self, parser: FTDNAFamFinderParser, tmp_path: Path
    ) -> None:
        """GH #113 cross-PR review (PR #117): the haploid shape (4-col
        line OR 5-col with empty ALLELE2) is accepted ONLY on MT and Y.
        FTDNA does not publish haploid calls on autosomes or X — the
        same shape on chr1 is almost certainly a truncated diploid
        that lost an allele in transit. Promoting it to hemizygous
        would silently double the surviving allele and synthesize a
        wrong-zygosity genotype. The parser warns-and-skips instead.

        Pre-review the parser accepted the autosomal 4-col line as
        haploid → `(A, A)`, silently mis-reporting a truncated `A?`
        as a homozygous AA call. Post-review the row is dropped with
        a warning, matching the canonical published FamFinder shape.
        """
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER
            + "rs_chr1_truncated\t1\t100\tA\n"  # 4-col on autosome
            + "rs_chr22_truncated\t22\t200\tG\n"  # 4-col on autosome
            + "rs_chrx_truncated\tX\t300\tC\n"  # 4-col on X
            + "rs_chr1_empty\t1\t400\tA\t\n"  # 5-col empty-ALLELE2 on autosome
            + "rs_haploid_mt\tMT\t500\tA\n"  # 4-col on MT — kept
            + "rs_diploid\t1\t600\tA\tG\n",  # canonical diploid — kept
        )
        by_rsid = {v.rsid: v for v in parser.parse(f)}
        # Only MT haploid + canonical diploid survive.
        assert set(by_rsid.keys()) == {"rs_haploid_mt", "rs_diploid"}
        assert by_rsid["rs_haploid_mt"].allele1 == "A"
        assert by_rsid["rs_haploid_mt"].allele2 == "A"
        assert by_rsid["rs_diploid"].allele1 == "A"
        assert by_rsid["rs_diploid"].allele2 == "G"

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

    def test_empty_allele2_cell_treated_as_haploid(
        self, parser: FTDNAFamFinderParser, tmp_path: Path
    ) -> None:
        """GH #113: a 5-column line with an empty trailing ALLELE2 cell
        is also the haploid shape — same biology as the 4-column line,
        just from a writer that emitted the trailing tab. Both must
        parse the same way: hemizygous, not half no-call.

        Pre-#113 the empty ALLELE2 was normalized to ``-``, producing
        a synthetic half no-call (e.g. ``A/-``) that the analyze
        pipeline correctly abstained on — losing every haploid MT/Y
        call to a no-call sink.
        """
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER + "rs_empty\tMT\t100\tA\t\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        # Hemizygous, not half no-call.
        assert variants[0].allele1 == "A"
        assert variants[0].allele2 == "A"
        assert not variants[0].is_no_call

    def test_empty_allele1_cell_still_no_call(
        self, parser: FTDNAFamFinderParser, tmp_path: Path
    ) -> None:
        """Symmetry guard: an empty ALLELE1 (the lead allele) doesn't
        get the haploid treatment — there's no called base to
        double. That row genuinely is a no-call, and normalize
        keeps it as ``-/-``."""
        f = _write(
            tmp_path,
            _FAMFINDER_HEADER + "rs_lead_empty\t1\t100\t\tG\n",
        )
        variants = list(parser.parse(f))
        assert len(variants) == 1
        assert variants[0].allele1 == "-"
        assert variants[0].allele2 == "G"
        assert variants[0].is_no_call


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

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Unit tests for shared parser helpers."""

from __future__ import annotations

import pytest

from allelix.parsers._helpers import normalize_chromosome


class TestNormalizeChromosome:
    """`chr` prefix stripping and mitochondrial collapse to Allelix canonical form."""

    @pytest.mark.parametrize("chrom", ["1", "2", "10", "22", "X", "Y", "MT"])
    def test_already_canonical_passes_through(self, chrom: str):
        assert normalize_chromosome(chrom) == chrom

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("chr1", "1"),
            ("chr2", "2"),
            ("chr10", "10"),
            ("chr22", "22"),
            ("chrX", "X"),
            ("chrY", "Y"),
        ],
    )
    def test_chr_prefix_stripped(self, raw: str, expected: str):
        assert normalize_chromosome(raw) == expected

    @pytest.mark.parametrize("raw", ["CHR1", "Chr1"])
    def test_chr_prefix_case_insensitive(self, raw: str):
        assert normalize_chromosome(raw) == "1"

    @pytest.mark.parametrize(
        "raw",
        ["MT", "M", "chrM", "chrMT", "ChrM", "chrm", "chrmt", "mt", "m", "CHRMT"],
    )
    def test_mitochondrial_collapses_to_mt(self, raw: str):
        """All mitochondrial variants (M, MT, chrM, chrMT, case variants) become MT.

        GRCh37 uses chrMT, GRCh38 uses chrM, bare M and bare MT both
        appear in real VCFs. Some pipelines lowercase the entire CHROM
        column. Allelix collapses every variant to MT regardless of case.
        """
        assert normalize_chromosome(raw) == "MT"

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("x", "X"),
            ("y", "Y"),
            ("chrx", "X"),
            ("chry", "Y"),
        ],
    )
    def test_sex_chromosomes_normalize_case(self, raw: str, expected: str):
        """Lowercase X/Y normalize to uppercase canonical form.

        Same failure mode as lowercased mito: a pipeline that lowercases
        the entire CHROM column would otherwise silently miss every
        rsID lookup for sex-chromosome variants.
        """
        assert normalize_chromosome(raw) == expected

    def test_unknown_contig_passes_through_after_chr_strip(self):
        """Unmapped contigs (chrUn_*, GL00*, alts) pass through.

        Downstream code decides whether to skip — the helper doesn't
        whitelist.
        """
        assert normalize_chromosome("chrUn_KI270742v1") == "Un_KI270742v1"
        assert normalize_chromosome("GL000209.1") == "GL000209.1"

    def test_empty_string_passes_through(self):
        assert normalize_chromosome("") == ""

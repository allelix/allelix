# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the Genes for Good parser.

Fixtures are synthetic: a placeholder banner and invented rows.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from allelix.parsers import detect_parser
from allelix.parsers.genesforgood import GenesForGoodParser
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.reports import input_method_label, is_imputed
from allelix.utils.build_detect import detect_build

FIXTURES = Path(__file__).parent.parent / "fixtures"
MEASURED = FIXTURES / "mock_genesforgood.txt"
IMPUTED = FIXTURES / "mock_genesforgood_imputed.txt"


@pytest.fixture
def parser() -> GenesForGoodParser:
    return GenesForGoodParser()


class TestDetection:
    @pytest.mark.parametrize("path", [MEASURED, IMPUTED])
    def test_detected_as_genes_for_good(self, path: Path) -> None:
        assert detect_parser(path).name == "genesforgood"

    def test_23andme_does_not_claim_it(self) -> None:
        assert not TwentyThreeAndMeParser().can_parse(MEASURED)

    def test_does_not_claim_a_23andme_export(self, parser: GenesForGoodParser) -> None:
        assert not parser.can_parse(FIXTURES / "mock_23andme.txt")

    def test_banner_must_be_the_first_line(
        self, parser: GenesForGoodParser, tmp_path: Path
    ) -> None:
        moved = tmp_path / "moved.txt"
        moved.write_text("# a comment\n" + MEASURED.read_text())
        assert not parser.can_parse(moved)


class TestParse:
    def test_rows(self, parser: GenesForGoodParser) -> None:
        variants = list(parser.parse(MEASURED))
        assert len(variants) == 11
        by_id = {v.rsid: v for v in variants}
        assert (by_id["rs9100002"].allele1, by_id["rs9100002"].allele2) == ("D", "I")
        assert (by_id["rs9100004"].allele1, by_id["rs9100004"].allele2) == ("A", "A")
        assert by_id["rs9100006"].chromosome == "MT"
        assert "gfg_placeholder_0001" in by_id

    def test_build_from_seeds_agrees_with_the_label(self, parser: GenesForGoodParser) -> None:
        seeds = detect_build(parser.parse(MEASURED))
        assert (seeds.build, seeds.is_confident) == ("GRCh37", True)


class TestMetadata:
    def test_measured(self, parser: GenesForGoodParser) -> None:
        meta = parser.get_metadata(MEASURED)
        assert meta == {"format": "genesforgood", "sample_id": "", "build": "GRCh37"}

    def test_imputed_is_marked(self, parser: GenesForGoodParser) -> None:
        meta = parser.get_metadata(IMPUTED)
        assert meta["method"] == "array-imputed"
        assert is_imputed(meta["method"])
        assert "most genotypes are imputed" in input_method_label(meta["method"])

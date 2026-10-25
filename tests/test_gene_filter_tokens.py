# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Exact constituent gene matching, with source-specific list separators."""

from pathlib import Path

import pytest

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult


def _ann(source="clinvar", gene="GENE1", **kw):
    return Annotation(
        source=source,
        gene=gene,
        rsid=kw.pop("rsid", "rs990001"),
        magnitude=kw.pop("magnitude", 5),
        category=kw.pop("category", "clinical"),
        significance="synthetic",
        description="synthetic",
        attribution=source,
        genotype_match="A/G",
        **kw,
    )


def _result(*annotations):
    return AnalysisResult(
        Path("synthetic"), "vcf", "VCF", "S", "GRCh37", 0, 0, [], annotations=list(annotations)
    )


@pytest.mark.parametrize("source", ["clinvar", "pharmgkb", "snpedia", "gwas"])
@pytest.mark.parametrize("field", ["GENE1;GENE2", "GENE1, GENE2", " Gene1 ; gene2 , ; "])
@pytest.mark.parametrize("query", ["gene1", " GENE2 ", "OTHER;gene2", "OTHER, gene1"])
def test_constituents_in_source_fields_and_queries(source, field, query) -> None:
    ann = _ann(source, field)
    assert _result(ann).filter(genes={query}) == [ann]
    assert ann.gene == field


@pytest.mark.parametrize("sep", [" - ", " x "])
@pytest.mark.parametrize("source", ["clinvar", "pharmgkb", "snpedia", "gwas"])
def test_spaced_flanking_and_interaction_separators_are_gwas_specific(source, sep) -> None:
    ann = _ann(source, f"GENE1{sep}GENE2")
    assert bool(_result(ann).filter(genes={"GENE2"})) == (source == "gwas")


@pytest.mark.parametrize("separator", [" - ", " x "])
@pytest.mark.parametrize("source", ["clinvar", "pharmgkb", "snpedia", "gwas"])
def test_full_field_query_and_single_gene_query_are_source_symmetric(source, separator) -> None:
    composite = f"GENE1{separator}GENE2"
    ann = _ann(source, composite)
    assert _result(ann).filter(genes={composite}) == [ann]
    constituent = _ann(source, "GENE2")
    assert bool(_result(constituent).filter(genes=iter([composite]))) == (source == "gwas")


@pytest.mark.parametrize("field", ["GENE1-AS1", "HLA-DQA1", "GENE10", "", None])
def test_no_substrings_or_arbitrary_hyphen_splitting(field) -> None:
    result = _result(_ann("gwas", field))
    assert result.filter(genes={"GENE1", "HLA", "DQA1"}) == []
    if field:
        assert len(result.filter(genes={field.lower()})) == 1
    assert result.filter(genes=set()) == []


def test_or_floors_category_and_sorting_still_apply() -> None:
    low = _ann("snpedia", "GENE1;GENE2", magnitude=2, rsid="rs2")
    high = _ann("gwas", "GENE1 - GENE2", magnitude=8, rsid="rs1")
    unrelated = _ann(gene="OTHER", magnitude=6, rsid="rs3")
    result = _result(low, unrelated, high)
    assert result.filter(genes={"GENE2"}, rsids={"RS3"}) == [high, unrelated, low]
    assert result.filter(genes={"GENE2"}, min_magnitude=3) == [high]
    assert result.filter(genes={"GENE2"}, category="trait") == []
    assert result.filter(genes={"GENE2"}, source_min_magnitudes={"gwas": 9}) == [low]
    assert result.filter(
        genes={"GENE2"}, min_magnitude=3, source_min_magnitudes={"snpedia": 1}
    ) == [high, low]


def test_filter_file_normalizes_gene_lists_and_preserves_hyphenated_symbols(
    tmp_path: Path,
) -> None:
    from allelix.cli._helpers import _parse_filter_file

    path = tmp_path / "panel.txt"
    path.write_text(
        "\ufeffgene1; gene2 , ; # genes\nHLA-DQA1\nRS990001 # id\nRSPO1\n", encoding="utf-8"
    )
    genes, rsids = _parse_filter_file(path)
    assert genes == frozenset({"GENE1", "GENE2", "HLA-DQA1", "RSPO1"})
    assert rsids == frozenset({"rs990001"})

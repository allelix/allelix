# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for GWAS MTAG + PheCode rollup (ADR-0024)."""

from __future__ import annotations

from allelix.models import Annotation
from allelix.reports._pipeline import rollup_gwas_duplicates


def _mk(rsid: str, desc: str, mag: float = 9.0, must: bool = False) -> Annotation:
    return Annotation(
        source="gwas",
        rsid=rsid,
        magnitude=mag,
        significance="gwas_association",
        category="trait",
        description=desc,
        attribution="GWAS Catalog",
        genotype_match="AC",
        is_must_include=must,
    )


def test_mtag_twin_collapsed():
    rows = [
        _mk("rs10455872", "GWAS Catalog: Aortic stenosis (p=4.0e-130, gene: LPA)"),
        _mk(
            "rs10455872",
            "GWAS Catalog: Aortic stenosis (MTAG) (p=4.0e-140, gene: LPA)",
        ),
    ]
    out = rollup_gwas_duplicates(rows)
    assert len(out) == 1
    assert "(MTAG)" not in out[0].description


def test_mtag_solo_kept_when_no_plain_twin():
    rows = [_mk("rs99999", "GWAS Catalog: Some trait (MTAG) (p=1.0e-50, gene: X)")]
    assert len(rollup_gwas_duplicates(rows)) == 1


def test_phecode_parent_child_collapsed_strongest_p_wins():
    rows = [
        _mk(
            "rs10455872",
            "GWAS Catalog: Ischemic heart disease (PheCode 411) (p=2.0e-204, gene: LPA)",
        ),
        _mk(
            "rs10455872",
            "GWAS Catalog: Coronary atherosclerosis (PheCode 411.4) (p=1.0e-234, gene: LPA)",
        ),
        _mk(
            "rs10455872",
            "GWAS Catalog: Other chronic IHD (PheCode 411.8) (p=3.0e-160, gene: LPA)",
        ),
    ]
    out = rollup_gwas_duplicates(rows)
    assert len(out) == 1
    assert "411.4" in out[0].description


def test_phecode_distinct_parents_kept_separate():
    rows = [
        _mk(
            "rs10455872",
            "GWAS Catalog: Hyperlipidemia (PheCode 272.1) (p=2.0e-100, gene: LPA)",
        ),
        _mk(
            "rs10455872",
            "GWAS Catalog: Ischemic heart disease (PheCode 411) (p=2.0e-204, gene: LPA)",
        ),
    ]
    assert len(rollup_gwas_duplicates(rows)) == 2


def test_must_include_never_collapsed():
    rows = [
        _mk(
            "rs9271366",
            "GWAS Catalog: MS (PheCode 335) (p=7.0e-184, gene: HLA-DRB1)",
            must=True,
        ),
        _mk(
            "rs9271366",
            "GWAS Catalog: MS (PheCode 335.1) (p=1.0e-50, gene: HLA-DRB1)",
        ),
    ]
    out = rollup_gwas_duplicates(rows)
    must_rsids = [a.rsid for a in out if a.is_must_include]
    assert "rs9271366" in must_rsids


def test_non_gwas_pass_through_untouched():
    rows = [
        Annotation(
            source="clinvar",
            rsid="rs1",
            magnitude=8.0,
            description="X",
            significance="pathogenic",
            attribution="ClinVar",
            category="clinical",
            genotype_match="AA",
        ),
        Annotation(
            source="snpedia",
            rsid="rs2",
            magnitude=3.0,
            description="Y",
            significance="snpedia_bad",
            attribution="SNPedia",
            category="clinical",
            genotype_match="AG",
        ),
    ]
    assert len(rollup_gwas_duplicates(rows)) == 2


def test_real_data_rs10455872_collapses_8_to_5():
    """Reviewer-flagged case: 8 LPA rows collapse to 5 distinct findings."""
    descriptions = [
        "GWAS Catalog: Aortic stenosis (p=4.0e-130, gene: LPA)",
        "GWAS Catalog: Aortic stenosis (MTAG) (p=4.0e-140, gene: LPA)",
        "GWAS Catalog: Hyperlipidemia (PheCode 272.1) (p=2.0e-100, gene: LPA)",
        "GWAS Catalog: Takes medication for coronary artery disease (p=3.0e-121, gene: LPA)",
        "GWAS Catalog: Coronary artery / coronary heart disease (p=5.0e-200, gene: LPA)",
        "GWAS Catalog: Ischemic heart disease (PheCode 411) (p=2.0e-204, gene: LPA)",
        "GWAS Catalog: Other chronic IHD (PheCode 411.8) (p=3.0e-160, gene: LPA)",
        "GWAS Catalog: Coronary atherosclerosis (PheCode 411.4) (p=1.0e-234, gene: LPA)",
    ]
    rows = [_mk("rs10455872", d) for d in descriptions]
    assert len(rollup_gwas_duplicates(rows)) == 5


def test_empty_list_returns_empty():
    assert rollup_gwas_duplicates([]) == []


def test_sort_order_preserved():
    """Output maintains magnitude DESC, rsid ASC sort."""
    rows = [
        _mk("rs222", "GWAS Catalog: Trait A (p=1.0e-50, gene: X)", mag=7.0),
        _mk("rs111", "GWAS Catalog: Trait B (p=1.0e-100, gene: Y)", mag=9.0),
    ]
    out = rollup_gwas_duplicates(rows)
    assert out[0].rsid == "rs111"
    assert out[1].rsid == "rs222"

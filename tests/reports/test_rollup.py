# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for GWAS MTAG + PheCode rollup (ADR-0024).

ADR-0035 PR 3: rollup now reads ``trait`` / ``phecode`` / ``p_value`` from
structured Annotation fields. Tests construct GWAS rows with the structured
fields populated and let the helper render a matching ``description`` so
MTAG-via-description detection still works.
"""

from __future__ import annotations

from allelix.models import Annotation
from allelix.reports._pipeline import rollup_gwas_duplicates


def _mk(
    rsid: str,
    trait: str,
    p_value: float,
    *,
    phecode: str = "",
    mtag: bool = False,
    gene: str = "LPA",
    mag: float = 9.0,
    must: bool = False,
) -> Annotation:
    """Build a GWAS Annotation with structured fields + matching description."""
    phecode_suffix = f" (PheCode {phecode})" if phecode else ""
    mtag_suffix = " (MTAG)" if mtag else ""
    description = (
        f"GWAS Catalog: {trait}{mtag_suffix}{phecode_suffix} (p={p_value:.1e}, gene: {gene})"
    )
    return Annotation(
        source="gwas",
        rsid=rsid,
        magnitude=mag,
        significance="gwas_association",
        category="trait",
        description=description,
        attribution="GWAS Catalog",
        genotype_match="AC",
        is_must_include=must,
        trait=trait,
        p_value=p_value,
        phecode=phecode,
    )


def test_mtag_twin_collapsed():
    rows = [
        _mk("rs10455872", "Aortic stenosis", 4.0e-130),
        _mk("rs10455872", "Aortic stenosis", 4.0e-140, mtag=True),
    ]
    out = rollup_gwas_duplicates(rows)
    assert len(out) == 1
    assert "(MTAG)" not in out[0].description


def test_mtag_solo_kept_when_no_plain_twin():
    rows = [_mk("rs99999", "Some trait", 1.0e-50, mtag=True, gene="X")]
    assert len(rollup_gwas_duplicates(rows)) == 1


def test_phecode_parent_child_collapsed_strongest_p_wins():
    rows = [
        _mk("rs10455872", "Ischemic heart disease", 2.0e-204, phecode="411"),
        _mk("rs10455872", "Coronary atherosclerosis", 1.0e-234, phecode="411.4"),
        _mk("rs10455872", "Other chronic IHD", 3.0e-160, phecode="411.8"),
    ]
    out = rollup_gwas_duplicates(rows)
    assert len(out) == 1
    assert "411.4" in out[0].description


def test_phecode_distinct_parents_kept_separate():
    rows = [
        _mk("rs10455872", "Hyperlipidemia", 2.0e-100, phecode="272.1"),
        _mk("rs10455872", "Ischemic heart disease", 2.0e-204, phecode="411"),
    ]
    assert len(rollup_gwas_duplicates(rows)) == 2


def test_must_include_never_collapsed():
    rows = [
        _mk(
            "rs9271366",
            "MS",
            7.0e-184,
            phecode="335",
            gene="HLA-DRB1",
            must=True,
        ),
        _mk("rs9271366", "MS", 1.0e-50, phecode="335.1", gene="HLA-DRB1"),
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
    rows = [
        _mk("rs10455872", "Aortic stenosis", 4.0e-130),
        _mk("rs10455872", "Aortic stenosis", 4.0e-140, mtag=True),
        _mk("rs10455872", "Hyperlipidemia", 2.0e-100, phecode="272.1"),
        _mk("rs10455872", "Takes medication for coronary artery disease", 3.0e-121),
        _mk("rs10455872", "Coronary artery / coronary heart disease", 5.0e-200),
        _mk("rs10455872", "Ischemic heart disease", 2.0e-204, phecode="411"),
        _mk("rs10455872", "Other chronic IHD", 3.0e-160, phecode="411.8"),
        _mk("rs10455872", "Coronary atherosclerosis", 1.0e-234, phecode="411.4"),
    ]
    assert len(rollup_gwas_duplicates(rows)) == 5


def test_empty_list_returns_empty():
    assert rollup_gwas_duplicates([]) == []


def test_sort_order_preserved():
    """Output maintains magnitude DESC, rsid ASC sort."""
    rows = [
        _mk("rs222", "Trait A", 1.0e-50, gene="X", mag=7.0),
        _mk("rs111", "Trait B", 1.0e-100, gene="Y", mag=9.0),
    ]
    out = rollup_gwas_duplicates(rows)
    assert out[0].rsid == "rs111"
    assert out[1].rsid == "rs222"

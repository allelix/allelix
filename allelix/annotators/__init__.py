# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Annotator registry. Unlike parsers, ALL annotators run on every variant."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.base import Annotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.clinvar import CLINVAR_SUPPORTED_BUILDS, ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator

if TYPE_CHECKING:
    from pathlib import Path


def get_annotators(
    data_dir: Path,
    clinvar_builds: tuple[str, ...] = CLINVAR_SUPPORTED_BUILDS,
    *,
    include_benign: bool = False,
    gwas_filter_traits: bool = True,
    cadd_full: bool = False,
) -> list[Annotator]:
    """Construct all registered annotators bound to the given data directory.

    `clinvar_builds` selects which ClinVar builds are managed by this
    process. Default is both GRCh37 and GRCh38 (per ADR-0021). The CLI
    narrows it via `db update --build grch37|grch38`.

    `include_benign` passes through to ClinVarAnnotator. Default False
    suppresses Benign/Likely_benign annotations (ADR-0008 amendment).

    `gwas_filter_traits` passes through to GWASCatalogAnnotator. Default
    True excludes common-trait noise categories (ADR-0024 amendment).

    `cadd_full` enables CADD full mode (tabix queries against the
    complete 81 GB CADD file). Requires ``pysam`` and a local copy.

    ADR-0023: ClinVar's `reference_for(rsid, build)` is wired into
    ClinPGx and SNPedia as the primary hom-ref suppression filter — the
    REF allele lookup universally determines whether the user is
    homozygous reference (and thus a non-finding for that variant).
    """
    clinvar = ClinVarAnnotator(data_dir, builds=clinvar_builds, include_benign=include_benign)
    pharmgkb = PharmGKBAnnotator(data_dir, clinvar_ref_provider=clinvar.reference_for)
    gwas = GWASCatalogAnnotator(data_dir, filter_traits=gwas_filter_traits)
    snpedia = SNPediaAnnotator(data_dir, clinvar_ref_provider=clinvar.reference_for)
    gnomad = GnomadAnnotator(data_dir)
    alphamissense = AlphaMissenseAnnotator(data_dir)
    cadd = CaddAnnotator(data_dir, full_mode=cadd_full)
    return [clinvar, pharmgkb, gwas, snpedia, gnomad, alphamissense, cadd]


_ANNOTATOR_CLASSES: dict[str, type[Annotator]] = {
    cls.name: cls
    for cls in [
        ClinVarAnnotator,
        PharmGKBAnnotator,
        GWASCatalogAnnotator,
        SNPediaAnnotator,
        GnomadAnnotator,
        AlphaMissenseAnnotator,
        CaddAnnotator,
    ]
}


def get_annotator_class(name: str) -> type[Annotator] | None:
    """Return the annotator class for a given source name, or None."""
    return _ANNOTATOR_CLASSES.get(name)


__all__ = [
    "AlphaMissenseAnnotator",
    "Annotator",
    "CaddAnnotator",
    "ClinVarAnnotator",
    "GWASCatalogAnnotator",
    "GnomadAnnotator",
    "PharmGKBAnnotator",
    "SNPediaAnnotator",
    "get_annotator_class",
    "get_annotators",
]

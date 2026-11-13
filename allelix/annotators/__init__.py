# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Annotator registry. Unlike parsers, ALL annotators run on every variant."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.base import Annotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.clinvar import CLINVAR_SUPPORTED_BUILDS, ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.reference_fallback import ReferenceCarriedFallback
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator

if TYPE_CHECKING:
    from collections.abc import Collection

# The cache directory of a switched-off source's stand-in (ADR-0067). It sits
# under the null device, which is a file, so no cache can ever exist there:
# the stand-in reads exactly as a source whose cache is absent.
SWITCHED_OFF_DIR = Path(os.devnull) / "allelix-switched-off"


def _absent_connection(cls: type[Annotator]) -> object:
    """What a source's ``_connection`` gives without a cache, without opening one."""
    if cls is ClinVarAnnotator:
        return lambda *_args, **_kwargs: None

    def refuse(*_args: object, **_kwargs: object) -> None:
        msg = f"{cls.name} is switched off"
        raise FileNotFoundError(msg)

    return refuse


def switched_off_stand_in[A: Annotator](cls: type[A], **kwargs: object) -> A:
    """A stand-in for a switched-off source: it reads as an absent cache (ADR-0067).

    Bound to ``SWITCHED_OFF_DIR``, where no cache can exist, and given the
    absent cache's connection answer, so it never opens anything.
    """
    stand_in = cls(SWITCHED_OFF_DIR, **kwargs)
    setattr(stand_in, "_connection", _absent_connection(cls))  # noqa: B010
    return stand_in


def configured_source[A: Annotator](
    name: str, cls: type[A], data_dir: Path, **kwargs: object
) -> A:
    """``cls`` bound to ``data_dir``, or its stand-in when the config turns ``name`` off.

    For the commands that read a source as a resolver or exporter (``stats``,
    ``extract``, ``export plink``): a source switched off in the data
    directory's config is not consulted there either (ADR-0067). The source
    is named by the caller, not read from ``cls``, so ``cls`` may be any
    factory that builds it.
    """
    from allelix.config import load_config

    if not load_config(data_dir, create=False).is_enabled(name, _ANNOTATOR_CLASSES):
        return switched_off_stand_in(cls, **kwargs)
    return cls(data_dir, **kwargs)


def get_annotators(
    data_dir: Path,
    clinvar_builds: tuple[str, ...] = CLINVAR_SUPPORTED_BUILDS,
    *,
    include_benign: bool = False,
    gwas_filter_traits: bool = True,
    cadd_full: bool = False,
    switched_off: Collection[str] = (),
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

    ADR-0023: ClinVar's `reference_at(rsid, build, locus)` is wired into
    ClinPGx and SNPedia as the primary hom-ref suppression filter — the
    REF allele lookup universally determines whether the user is
    homozygous reference (and thus a non-finding for that variant).

    `switched_off` names the sources the user turned off. ADR-0067: such a
    source is never consulted as a dependency: every other source is wired
    to a stand-in bound to `SWITCHED_OFF_DIR`, which reads as an absent
    cache, so its dependents take their unavailable path.
    """
    rsid_identity = RsidIdentityAnnotator(data_dir)
    identity = (
        switched_off_stand_in(RsidIdentityAnnotator)
        if "rsid_identity" in switched_off
        else rsid_identity
    )
    clinvar = ClinVarAnnotator(
        data_dir,
        builds=clinvar_builds,
        include_benign=include_benign,
        other_build_refs=identity.normal_alleles_on_other_build,
    )
    clinvar_ref = (
        switched_off_stand_in(
            ClinVarAnnotator,
            builds=clinvar_builds,
            include_benign=include_benign,
            other_build_refs=identity.normal_alleles_on_other_build,
        )
        if "clinvar" in switched_off
        else clinvar
    )
    gnomad = GnomadAnnotator(data_dir)
    gnomad_ref = switched_off_stand_in(GnomadAnnotator) if "gnomad" in switched_off else gnomad
    # ClinVar's REF where its cache for the call's build is readable; at a
    # reference-carried site without it, the other build's ClinVar or
    # gnomAD's major allele, else nothing (2.5.2).
    references = ReferenceCarriedFallback(clinvar_ref, identity, gnomad_ref)
    pharmgkb = PharmGKBAnnotator(
        data_dir,
        clinvar_ref_provider=references.reference_at,
        palindromic_provider=identity.palindromic_rsids,
        normal_allele_provider=references.normal_allele_for,
        reference_prefetch=references.prefetch,
        forward_site_provider=identity.forward_sites,
    )
    gwas = GWASCatalogAnnotator(
        data_dir,
        filter_traits=gwas_filter_traits,
        forward_site_provider=identity.forward_sites,
    )
    snpedia = SNPediaAnnotator(
        data_dir,
        clinvar_ref_provider=references.reference_at,
        normal_allele_provider=references.normal_allele_for,
        forward_site_provider=identity.forward_sites,
        reference_prefetch=references.prefetch,
    )
    alphamissense = AlphaMissenseAnnotator(data_dir)
    cadd = CaddAnnotator(data_dir, full_mode=cadd_full)
    return [clinvar, pharmgkb, gwas, snpedia, gnomad, alphamissense, cadd, rsid_identity]


_ANNOTATOR_CLASS_LIST: list[type[Annotator]] = [
    ClinVarAnnotator,
    PharmGKBAnnotator,
    GWASCatalogAnnotator,
    SNPediaAnnotator,
    GnomadAnnotator,
    AlphaMissenseAnnotator,
    CaddAnnotator,
    RsidIdentityAnnotator,
]
_ANNOTATOR_CLASSES: dict[str, type[Annotator]] = {cls.name: cls for cls in _ANNOTATOR_CLASS_LIST}


def get_annotator_class(name: str) -> type[Annotator] | None:
    """Return the annotator class for a given source name, or None."""
    return _ANNOTATOR_CLASSES.get(name)


__all__ = [
    "SWITCHED_OFF_DIR",
    "AlphaMissenseAnnotator",
    "Annotator",
    "CaddAnnotator",
    "ClinVarAnnotator",
    "GWASCatalogAnnotator",
    "GnomadAnnotator",
    "PharmGKBAnnotator",
    "RsidIdentityAnnotator",
    "SNPediaAnnotator",
    "configured_source",
    "get_annotator_class",
    "get_annotators",
    "switched_off_stand_in",
]

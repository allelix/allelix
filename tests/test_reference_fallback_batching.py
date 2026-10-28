# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The reference-carried fallback reads its caches in bulk, not once per call."""

from __future__ import annotations

import pytest

from allelix.annotators.reference_fallback import ReferenceCarriedFallback
from allelix.models import Variant

GRCH37, GRCH38 = "GRCh37", "GRCh38"
SITES = 300  # More than one gnomAD and one ClinVar query chunk.


class _ClinVar:
    def __init__(self, builds: set[str]) -> None:
        self.builds = builds
        self.reads = 0
        self.bulk = 0

    def reads_build(self, build: str) -> bool:
        self.reads += 1
        return build in self.builds

    def snv_alleles_many(self, rsids, build):
        self.bulk += 1
        # Odd sites: the other build's ClinVar names this build's REF as the variant.
        return {r: {("C", "T")} for r in rsids if int(r[2:]) % 2}


class _Identity:
    def __init__(self) -> None:
        self.ready = 0

    def is_ready(self) -> bool:
        self.ready += 1
        return True

    def reference_differing_sites(self):
        return {f"rs{i}": ("1", 1000 + i, "T", 2000 + i, "C") for i in range(SITES)}


class _Gnomad:
    def __init__(self) -> None:
        self.ready = 0
        self.bulk = 0

    def is_ready(self) -> bool:
        self.ready += 1
        return True

    def sites_frequencies(self, sites):
        self.bulk += 1
        # Every third site: the GRCh38 REF is major; the rest abstain (not listed).
        return {s: {"T": 0.1} for s in sites if (s[1] - 2000) % 3 == 0}


def _variants(build: str) -> list[Variant]:
    pos = 1000 if build == GRCH37 else 2000
    calls = [Variant(f"rs{i}", "1", pos + i, "T", "T", build=build) for i in range(SITES)]
    return [
        *calls,
        *(Variant(f"rs{9000 + i}", "1", 5000 + i, "A", "A", build) for i in range(SITES)),
    ]


def _answers(fallback: ReferenceCarriedFallback, variants: list[Variant]):
    return [
        (
            fallback.reference_at(v.rsid, v.build, ("1", v.position)),
            fallback.normal_allele_for(v.rsid, v.build, ("1", v.position)),
        )
        for v in variants
    ]


@pytest.mark.parametrize("clinvar_builds", [set(), {GRCH38}], ids=["noclinvar", "no37"])
def test_prefetch_decides_a_batch_with_one_query_per_source(clinvar_builds):
    clinvar, identity, gnomad = _ClinVar(clinvar_builds), _Identity(), _Gnomad()
    fallback = ReferenceCarriedFallback(clinvar, identity, gnomad)  # type: ignore[arg-type]
    variants = _variants(GRCH37)
    fallback.prefetch(variants)
    batched = _answers(fallback, variants)
    assert gnomad.bulk == 1
    assert clinvar.bulk == (1 if clinvar_builds else 0)
    # Readiness is read once per run, not per call (each check opens the cache).
    assert identity.ready == 1 and gnomad.ready == 1 and clinvar.reads <= 2

    single = ReferenceCarriedFallback(_ClinVar(clinvar_builds), _Identity(), _Gnomad())  # type: ignore[arg-type]
    assert _answers(single, variants) == batched
    assert fallback.abstentions() == single.abstentions()


def test_per_call_answers_without_prefetch_are_unchanged():
    clinvar, identity, gnomad = _ClinVar({GRCH38}), _Identity(), _Gnomad()
    fallback = ReferenceCarriedFallback(clinvar, identity, gnomad)  # type: ignore[arg-type]
    got = dict(
        zip(
            [v.rsid for v in _variants(GRCH37)], _answers(fallback, _variants(GRCH37)), strict=True
        )
    )
    # Other-build ClinVar first, then gnomAD's major allele, else abstain.
    assert got["rs1"] == ("C", "C")
    assert got["rs0"] == ("C", "C")
    assert got["rs2"] == (None, None)
    assert got["rs9000"] == (None, None)
    assert identity.ready == 1 and gnomad.ready == 1

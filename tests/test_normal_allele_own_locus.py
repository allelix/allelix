# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Both normal-allele consumers use the supplied locus, never another placement.

Synthetic VariationIDs deliberately have different normal alleles at two
placements of one rsID. Cross-build matching is by VariationID, as before;
its result must stay attached to the input placement. Legacy callers that
omit the locus are not assigned new behavior by this file.
"""

import pytest
from test_homref_ref_by_locus import (
    _analyze,
    _cv,
    _vcf,
    _write_clinvar,
    _write_empty_gwas,
    _write_identity,
    _write_pharm,
    _write_snpedia,
)

from allelix.annotators.clinvar import ClinVarAnnotator

RSID = "rs900910009"
BUILDS = ("GRCh37", "GRCh38")


def _data(tmp_path, build, reverse, *, second_carried):
    data = tmp_path / "data"
    data.mkdir()
    other = "GRCh38" if build == "GRCh37" else "GRCh37"
    own = [
        _cv(RSID, "1", 5000, "T", "T", allele_id=91),
        _cv(RSID, "1", 1000, "T", "T", allele_id=92)
        if second_carried
        else _cv(RSID, "1", 1000, "A", "G", allele_id=92),
    ]
    remote = [_cv(RSID, "1", 6000, "C", "T", allele_id=91)]
    if second_carried:
        remote.append(_cv(RSID, "1", 2000, "G", "T", allele_id=92))
    _write_clinvar(data, build, list(reversed(own)) if reverse else own)
    _write_clinvar(data, other, list(reversed(remote)) if reverse else remote)
    return data


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
def test_normal_allele_does_not_leak_into_an_ordinary_locus(tmp_path, build, reverse):
    data = _data(tmp_path, build, reverse, second_carried=False)
    with ClinVarAnnotator(data) as clinvar:
        assert clinvar.normal_allele_for(RSID, build, ("1", 1000)) is None
        assert clinvar.reference_at(RSID, build, ("1", 1000)) == "A"
        assert clinvar.normal_allele_for(RSID, build, ("1", 5000)) == "C"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
def test_report_alt_uses_the_ordinary_calls_ref(tmp_path, build, reverse):
    data = _data(tmp_path, build, reverse, second_carried=False)
    _write_identity(data, [(RSID, "1", 1000, "A", "G", 1000, "A", "G", 0)])
    _write_pharm(data, [(RSID, "AG", "synthetic")])
    _write_snpedia(data, [(RSID, "A", "G", 4.0, "synthetic", "plus")])
    _write_empty_gwas(data)
    source = _vcf(
        tmp_path / "input.vcf",
        [("1", 1000, RSID, "A", "G", 50, "PASS", ".", "GT:GQ", "0/1:99")],
        build=build,
    )
    report = _analyze(source, data)
    rows = [a for a in report["annotations"] if a["source"] in ("snpedia", "pharmgkb")]
    assert len(rows) == 2, "both genotype sources must actually match"
    assert {a["source"] for a in rows} == {"snpedia", "pharmgkb"}
    assert {a["alt"] for a in rows} == {"G"}


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("method", ["reference_at", "normal_allele_for"])
def test_reference_carried_answers_stay_at_their_own_loci(tmp_path, build, reverse, method):
    data = _data(tmp_path, build, reverse, second_carried=True)
    with ClinVarAnnotator(data) as clinvar:
        lookup = getattr(clinvar, method)
        assert lookup(RSID, build, ("1", 1000)) == "G"
        assert lookup(RSID, build, ("1", 5000)) == "C"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("method", ["reference_at", "normal_allele_for"])
@pytest.mark.parametrize("locus", [("1", 9999), ("2", 5000)])
def test_absent_placement_supplies_no_normal_allele(tmp_path, build, method, locus):
    data = _data(tmp_path, build, False, second_carried=False)
    with ClinVarAnnotator(data) as clinvar:
        assert getattr(clinvar, method)(RSID, build, locus) is None


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("method", ["reference_at", "normal_allele_for"])
def test_single_reference_carried_site_keeps_its_answer(tmp_path, build, method):
    other = "GRCh38" if build == "GRCh37" else "GRCh37"
    _write_clinvar(tmp_path, build, [_cv(RSID, "1", 5000, "T", "T", allele_id=91)])
    _write_clinvar(tmp_path, other, [_cv(RSID, "1", 6000, "C", "T", allele_id=91)])
    with ClinVarAnnotator(tmp_path) as clinvar:
        assert getattr(clinvar, method)(RSID, build, ("1", 5000)) == "C"

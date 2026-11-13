"""A stored REF equality cannot prefilter the normalized reference decision."""

from contextlib import ExitStack

import pytest
from test_annotation_parity_and_named_disclosures import _gnomad_site, _identity

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.reference_fallback import ReferenceCarriedFallback
from allelix.annotators.rsid_identity import RsidIdentityAnnotator, snv_form
from allelix.models import Variant

RSID = "rs900503999"
PLAIN = (RSID, "1", 1000, "A", "C", 2000, "G", "T", 0)
PADDED = (RSID, "1", 1000, "AG", "CG", 1999, "AG", "AT", 0)


def _caches(tmp_path):
    plain, padded = tmp_path / "plain", tmp_path / "padded"
    for directory, row in ((plain, PLAIN), (padded, PADDED)):
        directory.mkdir()
        _identity(directory, [row])
        _gnomad_site(directory, RSID, "1", 2000, "G", "T", 0.8)
    assert PADDED[3] == PADDED[6]
    assert snv_form(*PADDED[2:5]) == PLAIN[2:5]
    assert snv_form(*PADDED[5:8]) == PLAIN[5:8]
    return plain, padded


def test_equal_stored_refs_can_have_different_trimmed_reference_bases(tmp_path):
    plain, padded = _caches(tmp_path)
    with RsidIdentityAnnotator(plain) as identity:
        expected = identity.reference_differing_sites()
    assert expected == {RSID: ("1", 1000, "A", 2000, "G", False)}
    with RsidIdentityAnnotator(padded) as identity:
        assert identity.reference_differing_sites() == expected, (
            "Candidate selection must not discard a normalized REF difference"
        )


@pytest.mark.parametrize("build,pos,normal", (("GRCh37", 1000, "A"), ("GRCh38", 2000, "T")))
@pytest.mark.parametrize("prefetch", (False, True))
def test_equal_stored_refs_keep_the_gnomad_normal_decision(tmp_path, build, pos, normal, prefetch):
    plain, padded = _caches(tmp_path)
    answers = []
    for directory in (plain, padded):
        with ExitStack() as stack:
            cv = stack.enter_context(ClinVarAnnotator(directory))
            identity = stack.enter_context(RsidIdentityAnnotator(directory))
            gnomad = stack.enter_context(GnomadAnnotator(directory))
            fallback = ReferenceCarriedFallback(cv, identity, gnomad)
            if prefetch:
                fallback.prefetch([Variant(RSID, "1", pos, normal, normal, build=build)])
            answers.append(
                (
                    fallback.reference_at(RSID, build, ("1", pos)),
                    fallback.normal_allele_for(RSID, build, ("1", pos)),
                )
            )
    assert answers[0] == (normal, normal), "The unpadded control decides the normal allele"
    assert answers[1] == answers[0], "Padding cannot erase a confident normal allele"


@pytest.mark.parametrize("build,pos", (("GRCh37", 1000), ("GRCh38", 2000)))
def test_coordinate_and_rsid_lookups_already_agree_on_this_normalization(tmp_path, build, pos):
    plain, padded = _caches(tmp_path)
    answers = []
    for directory in (plain, padded):
        with RsidIdentityAnnotator(directory) as identity:
            answers.append(
                (
                    identity.forward_sites({RSID}, build),
                    identity.alleles_at_loci(build, {("1", pos)}),
                    identity.loci_of_rsids({RSID}),
                )
            )
    assert all(answers[0])
    assert answers[1] == answers[0]

"""Reverse-only single-base site evidence must not conflate deletion alleles."""

import pytest
from test_annotation_parity_and_named_disclosures import _identity, _pharm
from test_clinpgx_reverse_homozygotes import _call, _data, _find

from allelix.models import Variant
from allelix.utils.allele import trim_allele_pair


@pytest.mark.parametrize("build,pos", (("GRCh37", 1000), ("GRCh38", 2000)))
@pytest.mark.parametrize("mode", ("single", "batch"))
def test_trimming_a_deletion_does_not_turn_its_anchor_into_another_source_allele(
    tmp_path, build, pos, mode
):
    rsid = "rs900502999"
    _pharm(tmp_path, rsid, "SYNGENE", ["CC"], level="3")
    _identity(
        tmp_path,
        [(rsid, "1", 1000, "GAC", alt, 2000, "GAC", alt, 0) for alt in ("G", "GC", "GACC")],
    )
    # The input deletes one base. The source's complemented G is a different
    # identity allele (GAC>G, deleting two bases), not this call's anchor.
    assert trim_allele_pair(pos, "GAC", "GC") == (pos, "GA", "G")
    assert trim_allele_pair(pos, "GAC", "G") != (pos, "GA", "G")
    padded = Variant(rsid, "1", pos, "GC", "GC", build=build, ref="GAC", ref_from_input=True)
    trimmed = Variant(rsid, "1", pos, "G", "G", build=build, ref="GA", ref_from_input=True)
    assert not _find(tmp_path, padded, mode), "The sequence allele is not the CC source genotype"
    assert not _find(tmp_path, trimmed, mode), (
        "Removing shared sequence padding cannot make the source's different deletion a match"
    )


@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize("mode", ("single", "batch"))
def test_reverse_homozygous_snv_still_uses_the_site_evidence(tmp_path, build, mode):
    _data(tmp_path)
    rows = _find(tmp_path, _call(build), mode)
    assert len(rows) == 1

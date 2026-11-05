"""Matching insertion evidence keeps real source findings; unrelated records do not."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect, sql


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("gt", ("0/0", "0|0", "0"))
@pytest.mark.parametrize("evidence", ("matching", "anchor_snv", "other_insertion"))
def test_insertion_target_findings_need_matching_record_evidence(
    tmp_path, build, named, gt, evidence
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    column = "alt37" if build == "GRCh37" else "alt38"
    sql(data, "rsid_identity.sqlite", f"UPDATE rsid_identity SET {column}='GA'")
    alt = {"matching": "GA", "anchor_snv": "T", "other_insertion": "GT"}[evidence]
    path = _vcf(
        tmp_path / "input.vcf",
        [(RSID if named else ".", "1", pos, "G", alt, gt, ".")],
    )
    sources = ("gwas", "pharmgkb", "snpedia")
    result = analyze(data, path, build, sources)
    if evidence == "matching":
        assert {a.source for a in result.annotations} == set(sources)
        assert all(a.rsid == RSID for a in result.annotations)
        assert all(a.genotype_match in ("G", "GG", "G/G") for a in result.annotations)
    else:
        assert not result.annotations

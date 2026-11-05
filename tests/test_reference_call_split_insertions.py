"""A carried insertion is not a reference allele merely because its anchor survives."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect, sql


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("padded", (False, True))
@pytest.mark.parametrize("target_alt", ("A", "GA"))
@pytest.mark.parametrize("gt", ("1/1", "0/0"))
def test_split_insertion_does_not_invent_a_reference_genotype(
    tmp_path, build, reverse, padded, target_alt, gt
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    sql(
        data,
        "rsid_identity.sqlite",
        "UPDATE rsid_identity SET alt37=?, alt38=?",
        (target_alt, target_alt),
    )
    at, ref, alt = (pos - 1, "CG", "CGA") if padded else (pos, "G", "GA")
    # The SNV control has its own SNV evidence; the insertion control
    # relies on the matching second record when both calls are reference.
    first_alt = "A" if target_alt == "A" else "GT"
    rows = [(".", "1", pos, "G", first_alt, "0/0", "."), (".", "1", at, ref, alt, gt, ".")]
    if reverse:
        rows.reverse()
    sources = ("gwas", "pharmgkb", "snpedia")
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, sources)
    if gt == "1/1":
        assert not [a for a in result.annotations if a.genotype_match in ("G", "GG", "G/G")]
    else:
        assert {a.source for a in result.annotations} == set(sources)

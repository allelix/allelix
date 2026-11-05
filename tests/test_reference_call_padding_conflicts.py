"""Only the changed site of a padded SNV can displace a separate reference call."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect, signatures

SOURCES = ("gwas", "pharmgkb", "snpedia")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("shape", ("prefix", "suffix"))
@pytest.mark.parametrize("gt", ("0/1", "1/1"))
@pytest.mark.parametrize("changes_target", (False, True))
def test_reference_eligibility_is_the_same_beside_padded_or_trimmed_calls(
    tmp_path, build, reverse, shape, gt, changes_target
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    reference = (".", "1", pos, "G", "A", "0/0", ".")
    if changes_target:
        padded = (pos - 1, "CG", "CA") if shape == "prefix" else (pos, "GC", "AC")
        trimmed = (pos, "G", "A")
    elif shape == "prefix":
        padded = (pos - 1, "CGT", "CGA")
        trimmed = (pos + 1, "T", "A")
    else:
        padded = (pos - 1, "TG", "AG")
        trimmed = (pos - 1, "T", "A")

    results = []
    for label, (at, ref, alt) in (("padded", padded), ("trimmed", trimmed)):
        rows = [reference, (".", "1", at, ref, alt, gt, ".")]
        if reverse:
            rows.reverse()
        results.append(analyze(data, _vcf(tmp_path / f"{label}.vcf", rows), build, SOURCES))

    if not changes_target:
        assert {a.source for a in results[1].annotations} == set(SOURCES)
    assert signatures(results[0]) == signatures(results[1])
    if changes_target:
        assert not [
            a
            for result in results
            for a in result.annotations
            if a.genotype_match in ("G", "GG", "G/G")
        ]

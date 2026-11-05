"""The ploidy warning counts records in scope, rather than loci or every reference row."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect, warning_count


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("source", ("gwas", "snpedia", "clinvar"))
def test_target_references_and_nonreference_records_share_one_count(
    tmp_path, build, reverse, source
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    rows = [
        (".", "1", pos, "G", "A", "0", "."),
        (".", "1", pos, "G", "A", "1/1", "."),
        # Non-target reference does not make this diploid call contradictory.
        (".", "1", pos + 100, "G", "A", "0", "."),
        (".", "1", pos + 100, "G", "A", "1/1", "."),
        # The two non-reference records count even away from all targets.
        (".", "1", pos + 200, "G", "A", "1", "."),
        (".", "1", pos + 200, "G", "A", "0/1", "."),
    ]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, (source,))
    assert warning_count(result, "contradictory_calls") == 4

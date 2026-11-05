"""A rebuilt reference call keeps its record when the next record is on another chromosome."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0", "0/0"])
def test_record_number_survives_a_chromosome_change(tmp_path, build, gt):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    rows = [
        ("rs900260801", "1", 100, "C", "T", "0/1", "."),
        (RSID, "1", pos, "G", "A", gt, "."),
        ("rs900260802", "2", 500, "C", "T", "0/1", "."),
    ]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    found = [a for a in result.annotations if a.source == "clinvar" and a.rsid == RSID]
    assert found
    assert all(a.record_index == 1 for a in found)

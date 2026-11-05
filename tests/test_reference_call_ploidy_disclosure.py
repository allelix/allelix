"""Reference ploidy counts only at targets; non-reference records count everywhere."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect, warning_count


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("at_target", (False, True))
@pytest.mark.parametrize("sources", ((), ("pharmgkb",)), ids=("no_sources", "pharmgkb"))
@pytest.mark.parametrize(
    "gts,expected",
    [
        (("0", "0/0"), 2),
        (("0", "1/1"), 2),
        (("0/0", "1"), 2),
        (("0/0", "0/0"), 0),
        (("1", "0/1"), 2),
        (("0", "0/0", "0/0"), 3),
    ],
)
def test_reference_record_ploidy_counts_only_at_targets(
    tmp_path, build, reverse, at_target, sources, gts, expected
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    if not at_target:
        pos += 100
    rows = [(".", "1", pos, "G", "A", gt, ".") for gt in gts]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, sources)
    has_reference = any(set(gt.replace("|", "/").split("/")) == {"0"} for gt in gts)
    if has_reference and not (at_target and sources):
        expected = 0
    assert warning_count(result, "contradictory_calls") == expected

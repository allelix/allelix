"""Padding a base heterozygote must not bypass the shared identity guard."""

import pytest
from validated_call_helpers import BUILDS, RSID, analyze, call, setup, sql, warning_count


@pytest.mark.parametrize("rsid", (RSID, "."), ids=("named", "recovered"))
@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ("plain", "prefix", "suffix", "both"))
@pytest.mark.parametrize("gt", ("0/1", "1/0", "0|1", "1|0"))
@pytest.mark.parametrize("misplaced", (False, True), ids=("placed", "conflict"))
def test_padding_does_not_hide_a_base_call_identity_conflict(
    tmp_path, build, shape, gt, misplaced, rsid
):
    data = tmp_path / "db"
    pos = setup(data, build)
    if misplaced and rsid == ".":
        # Another source can recover the rsID at a contradictory locus.
        sql(data, f"clinvar.{build}.sqlite", "UPDATE clinvar_variants SET position=position+100")
    path = call(
        tmp_path / "input.vcf",
        pos + (100 if misplaced else 0),
        rsid=rsid,
        shape=shape,
        gt=gt,
    )
    result = analyze(data, path, build)
    assert warning_count(result, "identity_locus_conflict") == int(misplaced)
    if misplaced:
        assert not result.annotations
    else:
        assert any(a.source == "pharmgkb" for a in result.annotations)

"""Resolved two-rsID split pairs count the two calls actually assessed."""

import pytest
from split_site_contract_helpers import BUILDS, OTHER, RSID, run, signature, write_vcf
from test_split_pair_two_position_rsids import _counts, _pair, _two_rsid_data


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("auto", [False, True], ids=["explicit", "auto"])
def test_named_unnamed_and_roundtrip_have_the_same_call_total(tmp_path, build, auto):
    data = _two_rsid_data(tmp_path / "data")
    annotated = tmp_path / "annotated.vcf"
    unnamed, _ = run(
        data,
        write_vcf(tmp_path / "unnamed.vcf", build, _pair((".", "."))),
        tmp_path / "unnamed.json",
        build,
        auto=auto,
        extra=("--vcf-out", str(annotated)),
    )
    written, _ = run(
        data,
        write_vcf(tmp_path / "written.vcf", build, _pair((RSID, OTHER))),
        tmp_path / "written.json",
        build,
        auto=auto,
    )
    reread, _ = run(data, annotated, tmp_path / "reread.json", build, auto=auto)
    assert signature(unnamed) == signature(written) == signature(reread)
    assert _counts(unnamed) == _counts(written) == _counts(reread) == (2, 2)
    assert written["input"]["total_variants"] == reread["input"]["total_variants"] == 2
    assert unnamed["input"]["total_variants"] == 2, (
        "Route S assesses two half calls, just like written rsIDs and the annotated VCF; "
        "do not retain the parser's provisional one-pair total"
    )

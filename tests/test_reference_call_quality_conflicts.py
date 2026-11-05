"""A rejected reference genotype cannot overrule the parser's no-call decision."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect

from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("kind", ("low_gq", "zero_dp", "site_filter", "sample_filter"))
@pytest.mark.parametrize("override", (False, True))
def test_reference_blocking_uses_the_quality_decision(tmp_path, build, reverse, kind, override):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    path = _vcf(
        tmp_path / "input.vcf",
        [(".", "1", pos, "G", "A", "0/0", "."), (".", "1", pos, "G", "T", "0/0", ".")],
    )
    lines = path.read_text().splitlines()
    header = [line for line in lines if line.startswith("#")]
    rows = [line for line in lines if not line.startswith("#")]
    cols = rows[1].split("\t")
    if kind == "low_gq":
        cols[8:] = ["GT:GQ", "0/0:5"]
    elif kind == "zero_dp":
        cols[8:] = ["GT:DP", "0/0:0"]
    elif kind == "site_filter":
        cols[6] = "LowQual"
    else:
        cols[8:] = ["GT:FT", "0/0:LowAD"]
    rows[1] = "\t".join(cols)
    if reverse:
        rows.reverse()
    path.write_text("\n".join(header + rows) + "\n")
    parser = VcfParser(keep_filtered=override, min_gq=0 if override else 20)
    result = analyze(data, path, build, ("gwas", "pharmgkb", "snpedia"), parser=parser)
    allowed = override and kind != "zero_dp"
    expected = {"gwas", "pharmgkb", "snpedia"} if allowed else set()
    assert {a.source for a in result.annotations} == expected

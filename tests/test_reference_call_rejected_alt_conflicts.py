"""A rejected ALT call makes all declared alleles uncertain for reference evidence."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect

from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ("0/1", "1/1"))
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("kind", ("low_gq", "zero_dp", "site_filter", "sample_filter"))
@pytest.mark.parametrize("override", (False, True))
def test_rejected_alt_cannot_supply_certain_padding(tmp_path, build, gt, reverse, kind, override):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    path = _vcf(
        tmp_path / "input.vcf",
        [(".", "1", pos, "G", "A", "0/0", "."), (".", "1", pos - 1, "CG", "TG,CA", gt, ".")],
    )
    lines = path.read_text().splitlines()
    cols = lines[-1].split("\t")
    if kind == "low_gq":
        cols[8:] = ["GT:GQ", f"{gt}:5"]
    elif kind == "zero_dp":
        cols[8:] = ["GT:DP", f"{gt}:0"]
    elif kind == "site_filter":
        cols[6] = "LowQual"
    else:
        cols[8:] = ["GT:FT", f"{gt}:LowAD"]
    lines[-1] = "\t".join(cols)
    if reverse:
        lines[-2:] = reversed(lines[-2:])
    path.write_text("\n".join(lines) + "\n")
    result = analyze(
        data,
        path,
        build,
        ("gwas", "pharmgkb", "snpedia"),
        parser=VcfParser(keep_filtered=override, min_gq=0 if override else 20),
    )
    expected = {"gwas", "pharmgkb", "snpedia"} if override and kind != "zero_dp" else set()
    assert {a.source for a in result.annotations} == expected

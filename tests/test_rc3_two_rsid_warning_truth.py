"""Two known ALT copies routed by rsID do not make the genotype unknowable."""

import pytest
from rc3_verification_helpers import OTHER, RS, analyze, caches, warning, write_vcf


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("written", [False, True])
@pytest.mark.parametrize("gts", [("1/0", "0/1"), ("1/.", "./1"), ("1|0", "0|1")])
def test_two_rsid_routing_warning_does_not_claim_the_genotype_is_unknown(
    tmp_path, build, written, gts
):
    pairs = [(RS, "A"), (OTHER, "T")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, pairs)
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [
            dict(id=RS if written else ".", alt="A", gt=gts[0], ps=7),
            dict(id=OTHER if written else ".", alt="T", gt=gts[1], ps=7),
        ],
    )
    report, console = analyze(data, source, tmp_path / "report.json", build)
    notices = warning(report, "split_site_calls")
    assert len(notices) == 1 and notices[0]["count"] == 2
    assert report["input"]["vcf_quality"]["split_site_half_calls"] == 2
    message = " ".join(notices[0]["message"].lower().split())
    assert "genotype there cannot be told" not in message, (
        "The two ALT copies are known; only their different rsIDs prevent pairing"
    )
    assert "genotype there cannot be told" not in console.lower()
    assert "clinvar" in message and "gwas" in message
    assert "half" in message, "Keep the disclosure that these are assessed as half calls"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_contradictory_site_remains_counted_and_disclosed(tmp_path, build):
    pairs = [(RS, "A"), (RS, "C")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, [])
    source = write_vcf(
        tmp_path / "input.vcf", build, [dict(alt="A", gt="1/1"), dict(alt="C", gt="0/1")]
    )
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    assert warning(report, "split_site_calls")[0]["count"] == 2
    assert not report["annotations"]
    assert report["input"]["vcf_quality"]["split_site_half_calls"] == 0

"""The split warning must not claim a disabled source assessed a call."""

import pytest
from split_site_contract_helpers import BUILDS, data_dir, record, run, write_vcf


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("disabled", ("clinvar", "gwas", "both"))
def test_split_warning_and_half_note_do_not_claim_disabled_sources_assessed(
    tmp_path, build, disabled
):
    data = data_dir(tmp_path / "data")
    off = ("clinvar", "gwas") if disabled == "both" else (disabled,)
    config = data / "config.toml"
    text = config.read_text()
    for source in off:
        text = text.replace(source + "=true", source + "=false")
    config.write_text(text)
    source = write_vcf(tmp_path / "input.vcf", build, [record("A", "0/1"), record("C", "0/1")])
    report, text = run(data, source, tmp_path / "out.json", build)
    assert all(s["state"] == "skipped" for s in report["sources"] if s["name"] in off)
    warning = next(
        w["message"] for w in report["report_warnings"] if w["code"] == "split_site_calls"
    )
    assert "ClinVar and the GWAS Catalog assessed" not in warning
    assert "ClinVar and the GWAS Catalog assessed" not in " ".join(text.split())
    assert report["input"]["vcf_quality"]["split_site_half_calls"] == 2

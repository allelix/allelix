"""An unnamed contradicting split site is named by ClinVar's P/LP row at its position."""

import pytest
from split_site_contract_helpers import BUILDS, RSID, by_source, data_dir, record, run, write_vcf


@pytest.mark.parametrize("build", BUILDS)
def test_unnamed_contradicting_site_is_named_by_the_clinvar_row_at_its_position(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=False, alts="AC")
    rows = [record("A", "1/1", rsid="."), record("C", "0/1", rsid=".")]
    report, text = run(
        data,
        write_vcf(tmp_path / "input.vcf", build, rows),
        tmp_path / "out.json",
        build,
        extra=("--show-unchecked",),
    )
    assert all(not found for found in by_source(report).values())
    named = [r for r in report.get("high_value_unchecked", []) if r.get("reason") == "split_site"]
    assert [r.get("rsid") for r in named] == [RSID]
    assert RSID in text and "SYNGENE" in text


@pytest.mark.parametrize("build", BUILDS)
def test_unnamed_contradicting_site_without_a_plp_row_is_only_counted(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=False, alts="AC", pathogenic=False)
    rows = [record("A", "1/1", rsid="."), record("C", "0/1", rsid=".")]
    report, _text = run(
        data, write_vcf(tmp_path / "input.vcf", build, rows), tmp_path / "out.json", build
    )
    named = report.get("high_value_unchecked", [])
    assert not [r for r in named if r.get("reason") == "split_site"]
    warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warnings) == 1 and warnings[0]["count"] == 2

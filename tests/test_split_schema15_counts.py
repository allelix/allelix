"""Schema 15 keeps assessed standalone and split-site half calls disjoint."""

import json

import pytest
from split_site_contract_helpers import BUILDS, OTHER, POS, data_dir, record, run, write_vcf

from allelix.reports.diff import load_previous_report

SHAPES = {
    "standalone": ([record("A", "1/.")], (1, 0, 0)),
    "exact": ([record("A", "1/0"), record("C", "0/1")], (0, 0, 0)),
    "vt_exact": ([record("A", "1/."), record("C", "./1")], (0, 0, 0)),
    "phased_vt_exact": ([record("A", "1|.", ps="3"), record("C", ".|1", ps="3")], (0, 0, 0)),
    "agreeing": ([record("A", "0/1"), record("C", "0/1")], (0, 2, 2)),
    "half_and_het": ([record("A", "1/."), record("C", "0/1")], (0, 2, 2)),
    "rejected_sibling": ([record("A", "0/1"), record("C", "0/1", filt="q10")], (0, 1, 2)),
    "all_rejected": ([record("A", "./."), record("C", "./.")], (0, 0, 2)),
    "contradicting": ([record("A", "1/."), record("C", "1/1")], (0, 0, 2)),
}


def _counts(report):
    quality = report["input"]["vcf_quality"]
    assert "split_site_half_calls" in quality, "Schema 15 must publish the separate assessed count"
    warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warnings) <= 1
    return (
        quality["half_calls"],
        quality["split_site_half_calls"],
        warnings[0]["count"] if warnings else 0,
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("kind", SHAPES)
def test_half_counts_are_after_site_decision_and_never_overlap(tmp_path, build, kind):
    data = data_dir(tmp_path / "data")
    rows, expected = SHAPES[kind]
    report, _text = run(
        data, write_vcf(tmp_path / "in.vcf", build, rows), tmp_path / "out.json", build
    )
    assert _counts(report) == expected


@pytest.mark.parametrize("build", BUILDS)
def test_standalone_and_split_halves_in_one_file_have_separate_counts(tmp_path, build):
    data = data_dir(tmp_path / "data")
    rows = [
        record("A", "0/1"),
        record("C", "0/1"),
        record("A", "1/.", rsid=OTHER, pos=POS[build] + 100),
    ]
    report, text = run(
        data, write_vcf(tmp_path / "in.vcf", build, rows), tmp_path / "out.json", build
    )
    assert _counts(report) == (1, 2, 2)
    assert "split" in text.lower() and "ClinVar" in text and "GWAS" in text


@pytest.mark.parametrize("build", BUILDS)
def test_new_reports_write_schema15(tmp_path, build):
    data = data_dir(tmp_path / "data")
    report, _ = run(
        data,
        write_vcf(tmp_path / "in.vcf", build, [record("A", "0/1")]),
        tmp_path / "out.json",
        build,
    )
    assert report["schema_version"] == "15"


@pytest.mark.parametrize("schema", [str(n) for n in range(1, 16)])
def test_diff_accepts_current_and_historical_schema_versions(tmp_path, schema):
    report = tmp_path / "baseline.json"
    report.write_text(json.dumps({"schema_version": schema, "annotations": []}))
    assert load_previous_report(report)["annotations"] == []


def test_diff_refuses_schema16(tmp_path):
    report = tmp_path / "baseline.json"
    report.write_text(json.dumps({"schema_version": "16", "annotations": []}))
    with pytest.raises(ValueError, match=r"schema|version"):
        load_previous_report(report)

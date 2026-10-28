# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Explicit writes reject unwritable destinations with actionable CLI errors."""

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("report_format", ["json", "html"])
@pytest.mark.parametrize("blocked", ["file", "directory"])
def test_report_permission_failure_precedes_analysis(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    rc_no_analysis,
    command,
    report_format,
    blocked,
):
    case = rc_case
    output = case.out / f"report.{report_format}"
    if blocked == "file":
        output.write_text("existing report must survive\n")
        output.chmod(0o444)
    else:
        case.out.chmod(0o555)
    before = rc_snapshot(case.root)
    result = case.invoke([*case.report_args(command), "--output", output])
    assert not rc_no_analysis, result.output
    rc_clean_error(result, output if blocked == "file" else case.out)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_vcf_permission_failure_precedes_report_write(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    rc_no_analysis,
    suffix,
):
    case = rc_case
    blocked = case.root / "readonly"
    blocked.mkdir()
    blocked.chmod(0o555)
    output = blocked / ("annotated" + suffix)
    before = rc_snapshot(case.root)
    result = case.invoke(
        [*case.report_args("analyze"), "--output", case.out / "report.json", "--vcf-out", output]
    )
    assert not rc_no_analysis, result.output
    rc_clean_error(result, blocked)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("blocked", ["directory", ".bed", ".bim", ".fam"])
def test_plink_permission_failure_precedes_reading_calls(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    rc_no_analysis,
    blocked,
):
    case = rc_case
    prefix = case.out / "export"
    path = case.out if blocked == "directory" else prefix.with_suffix(blocked)
    if blocked != "directory":
        path.write_bytes(b"previous output\n")
    path.chmod(0o555 if blocked == "directory" else 0o444)
    before = rc_snapshot(case.root)
    result = case.invoke(
        [
            "export",
            "plink",
            case.source,
            "--data-dir",
            case.data,
            "--build",
            "grch37",
            "--output-prefix",
            prefix,
        ]
    )
    assert not rc_no_analysis, result.output
    rc_clean_error(result, path)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("blocked", ["file", "directory"])
def test_config_set_permission_error_preserves_config(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    blocked,
):
    case = rc_case
    path = case.data / "config.toml" if blocked == "file" else case.data
    path.chmod(0o444 if blocked == "file" else 0o555)
    before = rc_snapshot(case.root)
    result = case.invoke(["config", "set", "--data-dir", case.data, "sources.clinvar", "false"])
    rc_clean_error(result, path)
    assert rc_snapshot(case.root) == before


def test_db_clean_unwritable_directory_is_a_clean_error(rc_case, rc_snapshot, rc_clean_error):
    case = rc_case
    case.data.chmod(0o555)
    before = rc_snapshot(case.root)
    result = case.invoke(["db", "clean", "--data-dir", case.data, "--yes"])
    rc_clean_error(result, case.data)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize(
    "blocked",
    [
        "directory",
        "clinvar.GRCh37.sqlite",
        "clinvar.GRCh38.sqlite",
        "clinvar-build.log",
    ],
)
def test_db_update_unwritable_destination_precedes_download(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    monkeypatch,
    blocked,
):
    case = rc_case
    calls = []

    def forbidden(*args, **kwargs):
        calls.append("source setup")
        raise AssertionError("unwritable data directory reached source setup")

    monkeypatch.setattr(ClinVarAnnotator, "setup", forbidden)
    path = case.data if blocked == "directory" else case.data / blocked
    if blocked == "clinvar-build.log":
        path.write_text("previous build log\n")
    path.chmod(0o555 if blocked == "directory" else 0o444)
    before = rc_snapshot(case.root)
    result = case.invoke(
        [
            "db",
            "update",
            "--data-dir",
            case.data,
            "--only",
            "clinvar",
            "--force",
        ]
    )
    assert not calls, result.output
    rc_clean_error(result, path)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_existing_readonly_vcf_is_refused_without_other_writes(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    rc_no_analysis,
    suffix,
):
    case = rc_case
    output = case.out / ("annotated" + suffix)
    output.write_bytes(b"existing VCF bytes\n")
    output.chmod(0o444)
    before = rc_snapshot(case.root)
    result = case.invoke(
        [*case.report_args("analyze"), "--output", case.out / "report.json", "--vcf-out", output]
    )
    assert not rc_no_analysis, result.output
    rc_clean_error(result, output)
    assert rc_snapshot(case.root) == before

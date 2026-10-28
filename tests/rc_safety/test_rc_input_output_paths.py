# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Every explicit file input is protected from every applicable output."""

import os
from pathlib import Path

import pytest

REPORTS = ("analyze", "methylation", "pharmacogenomics")
INPUT_OUTPUTS = [
    (command, input_kind, output_kind)
    for command in REPORTS
    for input_kind in (("main", "diff", "filter") if command == "analyze" else ("main", "diff"))
    for output_kind in (("report", "vcf") if command == "analyze" else ("report",))
]


@pytest.mark.parametrize("command,input_kind,output_kind", INPUT_OUTPUTS)
@pytest.mark.parametrize("alias", ["same", "relative", "symlink", "hardlink", "parent-symlink"])
@pytest.mark.parametrize("report_format", ["json", "html"])
def test_file_input_output_aliases_are_refused_before_analysis(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    monkeypatch,
    command,
    input_kind,
    output_kind,
    alias,
    report_format,
):
    case = rc_case
    args = case.report_args(command)
    target = case.source
    if input_kind == "diff":
        target = case.root / "baseline.vcf"
        result = case.invoke(
            [*case.report_args(command), "--output", target, "--report-format", "json"]
        )
        assert result.exit_code == 0, (result.output, repr(result.exception))
        args += ["--diff", target]
    elif input_kind == "filter":
        target = case.root / "selection.vcf"
        target.write_text(case.rsid + "\n")
        args += ["--filter-file", target]
    if alias == "same":
        destination = target
    elif alias == "relative":
        destination = Path(os.path.relpath(target, case.root))
    elif alias == "parent-symlink":
        parent = case.root / "linked-parent"
        parent.symlink_to(target.parent, target_is_directory=True)
        destination = parent / target.name
    else:
        destination = case.root / "alias.vcf"
        if alias == "symlink":
            destination.symlink_to(target)
        else:
            os.link(target, destination)
    if output_kind == "report":
        args += ["--output", destination, "--report-format", report_format]
    else:
        args += [
            "--output",
            case.out / f"report.{report_format}",
            "--report-format",
            report_format,
            "--vcf-out",
            destination,
        ]
    calls = []

    def forbidden(*args, **kwargs):
        calls.append("analysis")
        raise AssertionError("input/output alias reached analysis")

    from allelix.cli import _helpers

    monkeypatch.setattr(_helpers, "run_analysis", forbidden)
    before = rc_snapshot(case.root)
    result = case.invoke(args)
    assert calls == [], result.output
    rc_clean_error(result, destination)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("command", REPORTS)
def test_distinct_baseline_is_read_without_modification(rc_case, command):
    case = rc_case
    baseline = case.root / "baseline.json"
    initial = case.invoke([*case.report_args(command), "--output", baseline])
    assert initial.exit_code == 0, initial.output
    before = baseline.read_bytes()
    output = case.out / "current.json"
    result = case.invoke([*case.report_args(command), "--diff", baseline, "--output", output])
    assert result.exit_code == 0, (result.output, repr(result.exception))
    assert baseline.read_bytes() == before
    assert output.is_file()

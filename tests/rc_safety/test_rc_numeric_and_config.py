# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Thresholds and configuration fail cleanly without corrupting saved inputs."""

import json

import click
import pytest

from allelix.cli import main
from allelix.config import load_config

REPORTS = ("analyze", "methylation", "pharmacogenomics")
FLOORS = ("--min-magnitude", "--gwas-min-magnitude", "--snpedia-min-magnitude")
CONFIG_CONSUMERS = ("config get", "config set", "config show", "db status", "db update", *REPORTS)


def test_every_float_option_is_in_the_finite_threshold_matrix():
    found = set()

    def walk(command, prefix):
        for param in command.params:
            if isinstance(param, click.Option) and isinstance(
                param.type, click.types.FloatParamType
            ):
                found.add((prefix, param.opts[0]))
        if isinstance(command, click.Group):
            for name, subcommand in command.commands.items():
                walk(subcommand, f"{prefix} {name}".strip())

    walk(main, "")
    assert found <= {(command, option) for command in REPORTS for option in FLOORS}, found


@pytest.mark.parametrize("command", REPORTS)
@pytest.mark.parametrize("option", FLOORS)
@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "NaN", "Infinity", "-Infinity"])
def test_nonfinite_threshold_is_a_usage_error_before_analysis(
    rc_case,
    rc_snapshot,
    rc_no_analysis,
    command,
    option,
    value,
):
    case = rc_case
    before = rc_snapshot(case.root)
    result = case.invoke(
        [*case.report_args(command), option, value, "--output", case.out / "r.json"]
    )
    assert not rc_no_analysis, result.output
    assert result.exit_code == 2, (result.output, repr(result.exception))
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "error" in result.output.lower() and option in result.output
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("command", REPORTS)
@pytest.mark.parametrize("option", FLOORS)
@pytest.mark.parametrize("value", ["-1", "0", "2.5"])
def test_finite_thresholds_still_produce_standard_json(rc_case, command, option, value):
    case = rc_case
    output = case.out / "report.json"
    result = case.invoke([*case.report_args(command), option, value, "--output", output])
    assert result.exit_code == 0, (result.output, repr(result.exception))

    def reject(token):
        raise AssertionError(f"nonstandard JSON number {token}")

    report = json.loads(output.read_text(), parse_constant=reject)
    assert report["schema_version"] == "11"


@pytest.mark.parametrize("command", CONFIG_CONSUMERS)
@pytest.mark.parametrize("contents", [b"[broken", b"\xff\x00"])
def test_config_consumers_name_invalid_file_without_traceback_or_writes(
    rc_case,
    rc_snapshot,
    rc_clean_error,
    command,
    contents,
):
    case = rc_case
    path = case.data / "config.toml"
    path.write_bytes(contents)
    if command in REPORTS:
        args = [*case.report_args(command), "--output", case.out / "report.json"]
    else:
        args = [*command.split(), "--data-dir", case.data]
        if command == "config get":
            args += ["sources.clinvar"]
        elif command == "config set":
            args += ["sources.clinvar", "false"]
        elif command == "db update":
            args += ["--only", "clinvar"]
    before = rc_snapshot(case.root)
    result = case.invoke(args)
    rc_clean_error(result, path)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize(
    "key",
    [
        "sources.clinvra",
        "sources.未知",
        "sources.",
        "license.clinvra",
        "license.未知",
        "license.",
    ],
)
@pytest.mark.parametrize("operation,value", [("get", None), ("set", "true"), ("set", "false")])
def test_unknown_source_key_is_rejected_without_changing_config(
    rc_case,
    rc_snapshot,
    key,
    operation,
    value,
):
    case = rc_case
    before = rc_snapshot(case.root)
    args = ["config", operation, "--data-dir", case.data, key]
    if value is not None:
        args.append(value)
    result = case.invoke(args)
    assert result.exit_code == 1, (result.output, repr(result.exception))
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "unknown" in result.output.lower(), result.output
    assert rc_snapshot(case.root) == before
    load_config(case.data)


@pytest.mark.parametrize("key", ["sources.clinvar", "license.cadd", "options.cadd_full"])
@pytest.mark.parametrize("value", ["true", "false"])
def test_known_config_keys_remain_writable(rc_case, key, value):
    case = rc_case
    result = case.invoke(["config", "set", "--data-dir", case.data, key, value])
    assert result.exit_code == 0, (result.output, repr(result.exception))
    queried = case.invoke(["config", "get", "--data-dir", case.data, key])
    assert queried.exit_code == 0, queried.output
    assert value in queried.output.lower(), queried.output
    load_config(case.data)

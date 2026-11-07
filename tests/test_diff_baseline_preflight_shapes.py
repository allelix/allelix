# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Unusable baseline shapes stop cleanly before update or analysis.

Moving the loader earlier is insufficient if it accepts values the diff
engine cannot use. Preserve legacy baselines; reject malformed input here.
All inputs are synthetic. The analysis sentinel avoids doing unnecessary
analysis on the known-broken path, after the update spies have fired.
"""

import json

import pytest
from click.testing import CliRunner
from test_early_output_checks import _offline, analyze_inputs  # noqa: F401

from allelix.cli import _helpers, main


def _report(**fields):
    return {
        "schema_version": "12",
        "annotations": [{"source": "clinvar", "rsid": "rs1", **fields}],
    }


BAD = [
    pytest.param([], id="array-root"),
    pytest.param(None, id="null-root"),
    pytest.param("report", id="string-root"),
    pytest.param(42, id="number-root"),
    pytest.param({"schema_version": [], "annotations": []}, id="array-schema"),
    pytest.param(_report(source=[]), id="array-source"),
    pytest.param(_report(rsid=["rs1"]), id="array-rsid"),
    pytest.param(_report(condition=["condition"]), id="array-condition"),
    pytest.param(_report(alt=["A"]), id="array-alt"),
    pytest.param(_report(magnitude="not a number"), id="nonnumeric-magnitude"),
]


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("payload", BAD)
def test_unusable_baseline_precedes_update_and_analysis(
    tmp_path,
    analyze_inputs,  # noqa: F811
    monkeypatch,
    command,
    payload,
):
    data, genome, calls = analyze_inputs
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(payload))
    output = tmp_path / "report.json"

    def analysis_must_not_run(*args, **kwargs):
        calls.append("analysis")
        raise AssertionError("unusable baseline reached analysis")

    monkeypatch.setattr(_helpers, "run_analysis", analysis_must_not_run)
    result = CliRunner().invoke(
        main,
        [
            command,
            str(genome),
            "--data-dir",
            str(data),
            "--output",
            str(output),
            "--diff",
            str(baseline),
        ],
    )
    assert calls == [], (calls, result.output, result.exception)
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit), (result.output, result.exception)
    assert "Error:" in result.output
    assert not output.exists()

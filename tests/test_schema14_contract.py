# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The release writes schema 14 and reads only the supported baseline versions."""

import json
from pathlib import Path

import pytest

from allelix.reports._pipeline import AnalysisResult
from allelix.reports.diff import load_previous_report
from allelix.reports.json_report import render_json


def test_written_report_is_schema14_and_can_be_read_back(tmp_path: Path) -> None:
    result = AnalysisResult(
        file_path=Path("synthetic.txt"),
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes (Tempus)",
        sample_id="synthetic",
        build="GRCh37",
        total_variants=1,
        skipped_count=0,
        annotators_used=[],
        annotations=[],
    )
    output = tmp_path / "report.json"
    render_json(result, output_path=output)
    payload = json.loads(output.read_text())
    assert payload["schema_version"] == "14"
    assert load_previous_report(output) == payload


@pytest.mark.parametrize("version", [str(n) for n in range(1, 15)])
def test_diff_accepts_each_released_or_current_schema(tmp_path: Path, version: str) -> None:
    payload = {"schema_version": version, "annotations": []}
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(payload))
    assert load_previous_report(baseline) == payload


@pytest.mark.parametrize("version", ["15", "99"])
def test_diff_refuses_future_schemas(tmp_path: Path, version: str) -> None:
    baseline = tmp_path / "future.json"
    baseline.write_text(json.dumps({"schema_version": version, "annotations": []}))
    with pytest.raises(ValueError, match="schema version"):
        load_previous_report(baseline)

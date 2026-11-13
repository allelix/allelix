# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A 2.5.2 report (schema 11) is a ``--diff`` baseline for a 2.6 report (schema 12).

``tests/data/report-2.5.2-schema11.json`` was written by the v2.5.2 release's own
``analyze``, unedited, on the synthetic VCF and ClinVar cache below (the ``data``
fixture and ``_row()`` of ``test_reference_carried_safety``).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.reports.diff import load_previous_report
from allelix.reports.json_report import SCHEMA_VERSION
from tests import test_reference_carried_safety as reference
from tests.test_reference_carried_safety import _row, _write

data = reference.data  # the synthetic ClinVar cache the baseline was written on

BASELINE = Path(__file__).parent / "data" / "report-2.5.2-schema11.json"


def test_baseline_is_a_2_5_2_schema_11_report() -> None:
    payload = load_previous_report(BASELINE)
    assert payload["schema_version"] == "11"
    assert payload["allelix_version"] == "2.5.2"
    assert [(a["rsid"], a["genotype_match"]) for a in payload["annotations"]] == [
        ("rs900230001", "TT")
    ]


def _diff(data: Path, tmp_path: Path, row: str) -> dict:
    path = _write(tmp_path / "reference.vcf", [row])
    output = tmp_path / "diff.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--data-dir",
            str(data),
            "--no-update",
            "--diff",
            str(BASELINE),
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(output.read_text())
    assert payload["schema_version"] == SCHEMA_VERSION == "14"
    return payload["diff"]


def test_same_input_diffs_clean_against_the_2_5_2_report(data: Path, tmp_path: Path) -> None:
    diff = _diff(data, tmp_path, _row())
    assert (diff["new"], diff["removed"], diff["changed"]) == ([], [], [])


@pytest.mark.parametrize("ref", ["A", "C"])
def test_a_finding_gone_since_2_5_2_is_removed(data: Path, tmp_path: Path, ref: str) -> None:
    diff = _diff(data, tmp_path, _row(ref=ref))
    assert diff["new"] == [] and diff["changed"] == []
    assert [a["rsid"] for a in diff["removed"]] == ["rs900230001"]

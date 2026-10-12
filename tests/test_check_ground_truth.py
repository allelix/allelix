# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""test_data/check_ground_truth.py reads every key of HG002_GROUND_TRUTH.yaml."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml

TEST_DATA = Path(__file__).resolve().parents[1] / "test_data"
sys.path.insert(0, str(TEST_DATA))

import check_ground_truth  # noqa: E402


def _truth() -> dict:
    return yaml.safe_load((TEST_DATA / "HG002_GROUND_TRUTH.yaml").read_text())


def test_every_key_in_the_yaml_is_read() -> None:
    assert check_ground_truth._schema_errors(_truth()) == []


@pytest.mark.parametrize(
    ("block", "key"),
    [
        (None, "extract"),
        ("floors.giab_grch38_benchmark", "total_annotations_mag_0_min"),
        ("floors.giab_grch38_benchmark", "expected_build_detection"),
        ("analyze.0", "expected_genotype"),
    ],
)
def test_an_unread_key_is_an_error(block: str | None, key: str) -> None:
    truth = _truth()
    target = truth
    for part in (block or "").split(".") if block else []:
        target = target[int(part)] if part.isdigit() else target[part]
    target[key] = 1
    errors = check_ground_truth._schema_errors(truth)
    assert any(repr(key) in e for e in errors), errors


def test_a_floor_without_build_source_is_an_error() -> None:
    truth = _truth()
    del truth["floors"]["giab_grch37_benchmark"]["build_source"]
    assert check_ground_truth._schema_errors(truth) == [
        "floors giab_grch37_benchmark: missing key 'build_source'"
    ]


def _report(**input_fields: object) -> dict:
    fields = {"file": "HG002_GRCh38_benchmark.vcf.gz", "build": "GRCh38"}
    fields.update(input_fields)
    return {"input": fields, "annotations": []}


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({}, None),
        ({"build_source": "contig_lengths"}, "input.build_source is 'contig_lengths'"),
        ({"build": "GRCh37"}, "input.build is 'GRCh37'"),
        ({"file": "other.vcf.gz"}, "input.file is 'other.vcf.gz'"),
    ],
)
def test_identity(change: dict, message: str | None) -> None:
    floors = _truth()["floors"]["giab_grch38_benchmark"]
    failures: list[str] = []
    check_ground_truth._check_identity(
        _report(**{"build_source": "override", **change}), floors, failures
    )
    if message is None:
        assert failures == []
    else:
        assert len(failures) == 1
        assert failures[0].startswith(f"identity: {message}, expected ")


def test_main_exits_two_on_an_unread_key(tmp_path: Path, monkeypatch, capsys) -> None:
    truth = _truth()
    truth["floors"]["giab_grch38_benchmark"]["default_filter_written_min"] = 240
    copy = tmp_path / "check_ground_truth.py"
    copy.write_text((TEST_DATA / "check_ground_truth.py").read_text())
    (tmp_path / "HG002_GROUND_TRUTH.yaml").write_text(yaml.safe_dump(truth))
    report = tmp_path / "r.json"
    report.write_text(json.dumps(_report(build_source="override")))
    monkeypatch.setattr(sys, "argv", [str(copy), str(report), "giab_grch38_benchmark"])
    monkeypatch.setattr(check_ground_truth, "__file__", str(copy))
    assert check_ground_truth.main() == 2
    assert "unknown key 'default_filter_written_min'" in capsys.readouterr().err

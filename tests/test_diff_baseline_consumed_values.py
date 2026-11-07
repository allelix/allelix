# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A baseline must be usable by both the diff engine and its output surfaces.

Synthetic removed findings expose null numeric/significance values and
rendered fields not read by compute_diff itself. Missing legacy fields and
null locus fields remain valid. No product changes or real data.
"""

import json

import pytest
from click.testing import CliRunner
from test_early_output_checks import _offline, analyze_inputs  # noqa: F401

from allelix.cli import _helpers, main

BAD = [
    pytest.param("magnitude", None, id="null-magnitude"),
    pytest.param("significance", None, id="null-significance"),
    pytest.param("magnitude", 10**400, id="unrepresentable-magnitude"),
    pytest.param("gene", ["GENE"], id="array-gene"),
    pytest.param("attribution", ["ClinVar"], id="array-attribution"),
    pytest.param("generated_at", {"date": "2026-11-05"}, id="object-date"),
]


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("field,value", BAD)
def test_unusable_consumed_value_stops_before_work(
    tmp_path,
    analyze_inputs,  # noqa: F811
    monkeypatch,
    command,
    field,
    value,
):
    data, genome, calls = analyze_inputs
    entry = {"source": "clinvar", "rsid": "rs900001", "magnitude": 1.0}
    payload = {"schema_version": "12", "annotations": [entry]}
    (payload if field == "generated_at" else entry)[field] = value
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(payload))

    def must_not_analyze(*args, **kwargs):
        calls.append("analysis")
        raise AssertionError("unusable baseline reached analysis")

    monkeypatch.setattr(_helpers, "run_analysis", must_not_analyze)
    result = CliRunner().invoke(
        main,
        [command, str(genome), "--data-dir", str(data), "--diff", str(baseline)],
    )
    assert calls == [], (calls, result.output, result.exception)
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit), (result.output, result.exception)
    assert "Error:" in result.output


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("kind", ["legacy-omissions", "null-locus", "integer-magnitude"])
def test_supported_baselines_still_run(tmp_path, analyze_inputs, command, kind):  # noqa: F811
    data, genome, calls = analyze_inputs
    entry = {"source": "clinvar", "rsid": "rs900001"}
    if kind == "null-locus":
        entry.update(chromosome=None, position=None, ref=None, alt=None)
    elif kind == "integer-magnitude":
        entry["magnitude"] = 4
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps({"schema_version": "12", "annotations": [entry]}))
    result = CliRunner().invoke(
        main,
        [command, str(genome), "--data-dir", str(data), "--diff", str(baseline)],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert "prepare" in calls

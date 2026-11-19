"""Consumed baseline build metadata must not crash after analysis."""

import json

import pytest
from click.testing import CliRunner
from test_diff_cross_build_contract import _baseline, _setup, _source

from allelix.cli import main


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("invalid", [[], {}, ["GRCh37"], {"name": "GRCh37"}])
def test_malformed_baseline_build_is_a_diagnostic_not_a_typeerror(tmp_path, build, invalid):
    data = _setup(tmp_path)
    baseline = _baseline(tmp_path, data, build)
    payload = json.loads(baseline.read_text())
    payload["input"]["build"] = invalid
    baseline.write_text(json.dumps(payload))
    source = _source(tmp_path, "new", build)
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            build,
            "--diff",
            str(baseline),
        ],
    )
    assert not isinstance(result.exception, (TypeError, AttributeError)), (
        result.output,
        result.exception,
    )
    # Either reject malformed metadata clearly, or consistently treat it as unknown.
    if result.exit_code:
        assert result.exit_code == 1 and "Error" in result.output and "build" in result.output
    else:
        assert "builds differ" not in result.output

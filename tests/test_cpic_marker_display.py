# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An opaque entity tag must not supply a cache-health warning."""

import pytest
from click.testing import CliRunner

from allelix.cli import main
from tests.test_http_freshness_fields import _update
from tests.test_http_freshness_fields import remote as remote


@pytest.mark.parametrize("command", ["update", "status"])
def test_opaque_tag_is_not_a_cpic_health_marker(tmp_path, remote, command):
    remote.headers["ETag"] = '"release|cpic:unavailable"'
    data = tmp_path / "data"
    output = _update(data, remote)
    if command == "status":
        result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])
        assert result.exit_code == 0, result.output
        output = result.output.lower()
    assert "non-finding filter degraded" not in output, output
    assert "no cpic" not in output, output

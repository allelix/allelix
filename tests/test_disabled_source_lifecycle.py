# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report startup must not refresh or reinterpret a switched-off cache."""

import os

import pytest
from click.testing import CliRunner
from test_disabled_source_reports import CHOICES, FLOORS, _configure, _prepare
from validated_call_helpers import sql

from allelix.annotators import get_annotators
from allelix.cli import main


@pytest.mark.parametrize("source,method", CHOICES)
@pytest.mark.parametrize("operation", ["refresh", "prepare"])
def test_report_startup_does_not_consult_disabled_source(
    tmp_path, monkeypatch, source, method, operation
):
    data, path = _prepare(tmp_path, "GRCh38")
    options = _configure(data, source, method)
    registry = get_annotators(data)
    calls = []
    for annotator in registry:
        cls = type(annotator)

        def remote(self):
            calls.append((self.name, "refresh"))
            return None

        def prepare(self):
            calls.append((self.name, "prepare"))
            return False

        # No test makes a remote request or rebuilds a cache.
        monkeypatch.setattr(cls, "fetch_remote_signal", remote)
        monkeypatch.setattr(cls, "prepare_cache", prepare)
    pattern = "clinvar.*.sqlite" if source == "clinvar" else source + ".sqlite"
    for cache in data.glob(pattern):
        if operation == "prepare":
            sql(data, cache.name, "UPDATE database_versions SET local_version_tag='obsolete'")
        os.utime(cache, (1, 1))
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--build",
            "grch38",
            "--data-dir",
            str(data),
            "--output",
            str(tmp_path / "report.json"),
            *FLOORS,
            *options,
        ],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert all(name != source for name, _operation in calls), calls

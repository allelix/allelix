# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Complete synthetic reports retain their pre-selection-change bytes."""

import hashlib

import pytest
from test_disabled_source_reports import FLOORS, _bytes, _configure, _prepare

from allelix.reports import json_report

# Captured with every cache installed before position-free input support.
# The placement field was added as null; every other report byte is unchanged.
# The schema-12 recapture changes only schema_version from 11 to 12, and the
# schema-13 recapture only from 12 to 13 (2.6.1).
# The version label is fixed so a release number alone does not invalidate the control.
# All report fields, including provenance and warnings, participate in the digest.
EXPECTED = {
    "GRCh37/ordinary/default": "63fe8d0b3e78ca5ef82c271dc98742a61fad1d7dc673cbf6f7802cd3311487dc",
    "GRCh37/ordinary/all-on": "f71dd455af0541ee71258453e651d958168fc09d3700866e8a6af76fd713319a",
    "GRCh37/ref/default": "754130074fa50d3494bdf625efb7363e82f5ac8f999575dcc0e3ed49481858c3",
    "GRCh37/ref/all-on": "4636df8b5c26b1de28e6708064770a065d758a093b6e23ef6c457df3873ca4ec",
    "GRCh38/ordinary/default": "2f8a2ab1e5688cc82dcd845d6d6d3d80858b4004cc9bb5d21fd1ea59e083923a",
    "GRCh38/ordinary/all-on": "5398aaba4ab1d2b89bd7f6e4a4f1560f3bb32ae3b37fe9b2a013149d8e4fd457",
    "GRCh38/ref/default": "b7c2e1f6b62b5fa0e40eddd7b833b59ecd594c6bb73c278fef6ebf5a07037833",
    "GRCh38/ref/all-on": "8de89cd92c9d9ee944a0a7cd9ade01d68c6a9ea4272baaf6170ad2b92b13274f",
}


@pytest.mark.parametrize("case", EXPECTED)
def test_default_and_all_on_reports_keep_their_complete_bytes(tmp_path, monkeypatch, case):
    monkeypatch.setattr(json_report, "__version__", "synthetic")
    build, scenario, mode = case.split("/")
    data, path = _prepare(tmp_path, build, "reference-all" if scenario == "ref" else scenario)
    if mode == "all-on":
        _configure(data)
    _, raw = _bytes(tmp_path, data, path, build, "analyze", FLOORS, "report")
    # _bytes removes generated_at only. No fields or rows may be dropped here.
    assert hashlib.sha256(raw).hexdigest() == EXPECTED[case]

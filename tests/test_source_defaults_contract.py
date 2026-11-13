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
# schema-13 recapture only from 12 to 13 (2.6.1); the schema-14 recapture
# changes only 13 to 14 (2.7.0), verified against all eight prior digests.
# The version label is fixed so a release number alone does not invalidate the control.
# All report fields, including provenance and warnings, participate in the digest.
EXPECTED = {
    "GRCh37/ordinary/default": "45afba9234aa5cc7d64a43a0d082a7e65c4d5de60f255e765aa861d030963afb",
    "GRCh37/ordinary/all-on": "4f3956bb9cc857dc08d1d969e7f83b6daf87df3b7ad67050ea229ee02ff22661",
    "GRCh37/ref/default": "bbfc7e5cfdf55b3e92134128e72a86ebfc1e1cac4a48c73bb94e350dc1075ccc",
    "GRCh37/ref/all-on": "ea1a14f80aced387ec5dfd2fc9f6c3193d76e80abd43d1afa8571df27ec8fc50",
    "GRCh38/ordinary/default": "f007383b22b30320a129a309d071454ee006a91838b4987a211d2411a032f39e",
    "GRCh38/ordinary/all-on": "a5fb576d8c85e1f8b9e2fc0a75017ee3510f68e2573896ee6329efe402376523",
    "GRCh38/ref/default": "dc3a067ae4ecc9e69b866f972e6914310649faaf448156fb630cd3c63fbadbcc",
    "GRCh38/ref/all-on": "00e3b22045ff6f18ccecadc2693c917c3bb74ce83b2a2c70140a0d5e4ff19588",
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

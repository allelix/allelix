# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Complete synthetic reports retain their pre-selection-change bytes."""

import hashlib

import pytest
from test_disabled_source_reports import FLOORS, _bytes, _configure, _prepare

from allelix.reports import json_report

# Captured with every cache installed before position-free input support.
# The placement field was added as null; every other report byte is unchanged.
# The schema-12 recapture changes only schema_version from 11 to 12.
# The version label is fixed so a release number alone does not invalidate the control.
# All report fields, including provenance and warnings, participate in the digest.
EXPECTED = {
    "GRCh37/ordinary/default": "88fc9e3a7a8bbcc7f339d44319fc909e6f15550f6d363a3f119836fad123be5d",
    "GRCh37/ordinary/all-on": "2c9d6b10758bd8c1547333c9ea3a6cb2853de7bff28e4fd8774a1cd4b0ad1b6e",
    "GRCh37/ref/default": "ff81276cc2e0e89da5c56582da317988201bd8c4fc172e6813a096f43c8fcaba",
    "GRCh37/ref/all-on": "6a33f3aaf5b5744e14f52d377580a198885e580f9cf7798017577355ec92ee09",
    "GRCh38/ordinary/default": "91b917a988768af611c65094bc14d71d59131b5b12a1bfe27ed6abab02f861d3",
    "GRCh38/ordinary/all-on": "389cab744015fc62a04a1abfffe49af7a6ab294a0e66031c279a6deac2d448b3",
    "GRCh38/ref/default": "1ebe78a64762d100bfe841dc9ea05ff987a9ff42c9087dcc92dbe7ee5a99a080",
    "GRCh38/ref/all-on": "fd156ac32944cb1a7630cd0a7968e3b1b1ada5201ef63ebde4c59349854b2c85",
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

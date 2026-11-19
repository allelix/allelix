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
# changes only 13 to 14 (2.7.0), verified against all eight prior digests;
# the schema-15 recapture (2.7.1) changes only 14 to 15 and adds
# input.vcf_quality.split_site_half_calls = 0, verified the same way.
# The version label is fixed so a release number alone does not invalidate the control.
# All report fields, including provenance and warnings, participate in the digest.
EXPECTED = {
    "GRCh37/ordinary/default": "71fcb5348801649b6c793edb25b7ca674bb7c92181d86528b1c52155e5dd1a62",
    "GRCh37/ordinary/all-on": "25ce0859f9e5ee40b66d40192ccb2a31a1d43fe09eb59c42a60cbff3ec10101c",
    "GRCh37/ref/default": "f626329ae8a6bcf27fdfa42a9703312ac82ee43031f55b6bc6436ba8750b0e5e",
    "GRCh37/ref/all-on": "f6642bec6148432621e14d68cdbc642fa1923b1f277c028612671ad761c85304",
    "GRCh38/ordinary/default": "ccf4ea9d283e42ce58fdfde4d019a09aa9bf6dd039123c03d88055c75fc041c8",
    "GRCh38/ordinary/all-on": "b84455b4e596ba31f9a5ad8d5fb4480a99fc24d738419dd8506e24f3786bc634",
    "GRCh38/ref/default": "765df1298ceb06348611e4485e1c0d46600a48c868ba1f24f08b696e36c4552d",
    "GRCh38/ref/all-on": "5db6d516339098080e3c9c23ea947b7b1116c51e633ae047edc10cd6f53c6cf8",
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

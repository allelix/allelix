# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A rebuilt finding keeps the number of the eligible source record."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0", "0/0"])
@pytest.mark.parametrize("tail", ["none", "same-chrom", "next-chrom"])
def test_reference_origin_when_watcher_releases(tmp_path, build, gt, tail):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    rows = [(RSID, "1", pos, "G", "A", gt, ".")]
    if tail != "none":
        rows.append(
            ("rs900260802", "1" if tail == "same-chrom" else "2", pos + 20, "C", "T", "0/1", ".")
        )
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    found = [a for a in result.annotations if a.source == "clinvar"]
    assert found
    assert all(a.record_index == 0 for a in found)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("mismatch_first", [False, True])
def test_reference_origin_rejects_mismatch(tmp_path, build, mismatch_first):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    matching = (RSID, "1", pos, "G", "A", "0/0", ".")
    mismatching = (RSID, "1", pos, "C", "T", "0/0", ".")
    rows = [mismatching, matching] if mismatch_first else [matching, mismatching]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    found = [a for a in result.annotations if a.source == "clinvar"]
    assert found
    assert all(a.record_index == int(mismatch_first) for a in found)

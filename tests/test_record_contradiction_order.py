# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Contradiction counts describe records, independent of their input order."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, setup, warning_count

from allelix.reports import _pipeline


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch_size", [1, 1000])
@pytest.mark.parametrize("interleaved", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_a_different_locus_does_not_hide_a_ploidy_contradiction(
    tmp_path, monkeypatch, build, batch_size, interleaved, reverse
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    haploid = (RSID, "1", pos, "G", "A", "1", ".")
    diploid = (RSID, "1", pos, "G", "A", "0/1", ".")
    spacer = ("rs900260702", "1", pos + 20, "G", "A", "0/1", ".")
    rows = [haploid, spacer, diploid] if interleaved else [haploid, diploid, spacer]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    assert warning_count(result, "contradictory_calls") == 2

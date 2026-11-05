# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Padding does not hide a contradiction at the same biological site."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, setup, warning_count

from allelix.reports import _pipeline


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch_size", [1, 1000])
@pytest.mark.parametrize("padded", [False, True])
def test_padding_preserves_the_same_site_contradiction_count(
    tmp_path, monkeypatch, build, batch_size, padded
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    first = (
        (RSID, "1", pos - 1, "CG", "CA", "1", ".")
        if padded
        else (RSID, "1", pos, "G", "A", "1", ".")
    )
    second = (RSID, "1", pos, "G", "A", "0/1", ".")
    result = analyze(data, _vcf(tmp_path / "input.vcf", [first, second]), build)
    assert warning_count(result, "contradictory_calls") == 2

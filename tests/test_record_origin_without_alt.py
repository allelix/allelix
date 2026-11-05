# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A padded finding with no reported effect allele still has its own origin."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, setup, sql

from allelix.reports import _pipeline


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["plain", "prefix", "suffix", "both"])
@pytest.mark.parametrize("gt", ["1", "1/1"])
@pytest.mark.parametrize("batch_size", [1, 1000])
def test_unknown_effect_allele_keeps_the_padded_calls_origin(
    tmp_path, monkeypatch, build, shape, gt, batch_size
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    sql(data, "gwas.sqlite", "UPDATE gwas_associations SET risk_allele=NULL")
    ref, alt, shift = {
        "plain": ("G", "A", 0),
        "prefix": ("CG", "CA", 1),
        "suffix": ("GC", "AC", 0),
        "both": ("CGC", "CAC", 1),
    }[shape]
    rows = [
        (RSID, "1", pos - shift, ref, alt, gt, "."),
        ("rs900260702", "1", pos + 20, "G", "A", "0/1", "."),
    ]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    found = [a for a in result.annotations if a.source == "gwas"]
    assert found
    for annotation in found:
        assert annotation.record_index == 0
        assert annotation.haploid == (gt == "1")
        diploid = alt + alt if len(alt) == 1 else f"{alt}/{alt}"
        assert annotation.genotype_match == (alt if gt == "1" else diploid)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Record provenance for recovered alleles, reference blocks and rejected half calls."""

from __future__ import annotations

import contextlib
import sqlite3

import pytest

from allelix.annotators.clinvar import clinvar_db_filename
from tests import test_vcf_reference_findings as b


@pytest.fixture
def data(tmp_path):
    return b.data.__wrapped__(tmp_path)


@pytest.mark.parametrize("rsid", [".", "rs900252005"])
@pytest.mark.parametrize("reverse", [False, True])
def test_each_resolved_record_keeps_only_its_own_allele(data, tmp_path, rsid, reverse):
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh38"))) as conn:
        b._insert(
            conn, "rs900252005", 5000, "G", "T", gene="SECOND", sig="Likely_pathogenic", allele=90
        )
        conn.commit()
    genotypes = ["0/2", "0/1"] if reverse else ["0/1", "0/2"]
    rows = [b._row(5000, rsid, "G", "A,T", gt) for gt in genotypes]
    out = b._run(data, b._input(tmp_path, rows), tmp_path, "GRCh38")
    values = [
        b._value("rs900252005", "ALTGENE"),
        b._value("rs900252005", "SECOND", "clinvar_likely_pathogenic"),
    ]
    assert len(out[0]["annotations"]) == 2
    for gt, info in zip(genotypes, out[2], strict=True):
        expected = (values[0], ".") if gt == "0/1" else (".", values[1])
        assert info[b.ALT_FIELD] == expected


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("alt", [".", "<NON_REF>"])
def test_end_block_cannot_borrow_reference_stamp(data, tmp_path, reverse, alt):
    rows = [
        b._row(4000, "rs900252004", "T", alt, "0/0"),
        b._row(4000, "rs900252004", "T", alt, "0/0", info="END=4001;KEEP=x"),
    ]
    if reverse:
        rows.reverse()
    out = b._run(data, b._input(tmp_path, rows), tmp_path, "GRCh37")
    assert out[2][int(reverse)].get(b.REF_FIELD) == b._value("rs900252004")
    assert out[2][int(not reverse)] == {"KEEP": "x"}


@pytest.mark.parametrize("alt", [".", "C"])
@pytest.mark.parametrize("gt", ["0/.", ".|0"])
@pytest.mark.parametrize("gq,dp", [(3, 20), (99, 0)])
def test_reference_half_call_quality_failure_is_reported(data, tmp_path, alt, gt, gq, dp):
    rows = [
        b._row(5000, "rs900252005", "G", "A", "0/1"),
        b._row(169519049, "rs6025", "T", alt, gt, gq=gq, dp=dp),
    ]
    out = b._run(data, b._input(tmp_path, rows), tmp_path, "GRCh37")
    assert not any(a["rsid"] == "rs6025" for a in out[0]["annotations"])
    assert out[0]["input"]["vcf_quality"]["low_quality_calls"] == 1
    assert any(
        v["rsid"] == "rs6025" and v["reason"] == "low_gq_call"
        for v in out[0]["high_value_no_calls"]
    )

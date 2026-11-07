# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""SNPedia reads 23andMe's own i-probe pages for a file the 23andMe parser read.

A 23andMe ``i`` ID names 23andMe's own probe, and SNPedia's page for it is
written for that probe, so the page answers which indel the D/I call names
(the question ADR-0070 is about). A D/I call at an ``i`` ID in a file read
by the 23andMe parser is matched to its page again, as in 2.5.2. Everywhere
else ADR-0070 stands: a D/I call at an rsID, any D/I call from another
parser (Genes for Good writes 23andMe-style files, but its ``i`` IDs are not
23andMe's probes), and every other source. The D/I count is unchanged.

Fixtures are the synthetic ones of ``test_di_calls_not_interpreted``.
"""

from __future__ import annotations

import pytest
from test_di_calls_not_interpreted import (
    _CHIP_CALLS,
    _DI_RSIDS,
    DI_PLUS,
    IPROBE,
    _analyze_json,
    _array,
    _both_paths,
    _call,
    _di_rows,
    _snpedia,
    data,  # noqa: F401  (pytest fixture)
    identity,  # noqa: F401  (pytest fixture)
)


@pytest.mark.parametrize("genotype", ("DI", "DD"))
def test_the_switch_lets_a_d_i_call_reach_its_i_probe_page(
    data,  # noqa: F811
    identity,  # noqa: F811
    genotype,
):
    annotator = _snpedia(data, identity)
    annotator.read_vendor_probe_pages(True)
    single, batch = _both_paths(annotator, _call(IPROBE, genotype))
    expected = [f"SNPedia: i-ID page {genotype}"]
    assert [a.description for a in single] == expected
    assert [a.description for a in batch] == expected


@pytest.mark.parametrize("genotype", ("DI", "DD", "II"))
def test_the_switch_leaves_a_d_i_call_at_an_rsid_unread(
    data,  # noqa: F811
    identity,  # noqa: F811
    genotype,
):
    annotator = _snpedia(data, identity)
    annotator.read_vendor_probe_pages(True)
    assert _both_paths(annotator, _call(DI_PLUS, genotype)) == ([], [])


def test_the_switch_is_off_by_default(
    data,  # noqa: F811
    identity,  # noqa: F811
):
    assert _both_paths(_snpedia(data, identity), _call(IPROBE, "DI")) == ([], [])


def test_a_23andme_file_gets_its_i_probe_page_and_nothing_else_for_d_i_calls(
    tmp_path,
    data,  # noqa: F811
):
    report = _analyze_json(
        data, _array(tmp_path / "a.txt", _CHIP_CALLS), tmp_path / "r.json", "--build", "grch37"
    )
    assert _di_rows(report, _DI_RSIDS) == [("snpedia", IPROBE)]
    assert report["input"]["array_quality"]["di_calls"] == sum(
        1 for _rsid, gt in _CHIP_CALLS if set(gt) & {"D", "I"}
    )


def test_another_parser_does_not_read_i_probe_pages(
    tmp_path,
    data,  # noqa: F811
):
    report = _analyze_json(
        data,
        _array(tmp_path / "a.txt", _CHIP_CALLS),
        tmp_path / "r.json",
        "--build",
        "grch37",
        "--format",
        "genesforgood",
    )
    assert _di_rows(report, _DI_RSIDS) == []

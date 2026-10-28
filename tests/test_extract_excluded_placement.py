# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An excluded mitochondrial placement cannot settle an autosomal extract request."""

from __future__ import annotations

import pysam
import pytest
from test_extract_reference_calls import _cache, _extract, _row


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("length", [16569, 16571])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
def test_excluded_best_placement_cannot_stop_before_a_better_eligible_call(
    tmp_path, monkeypatch, build, length, reverse, indexed
):
    places = [("MT", 100, 100, "G", "A"), ("1", 200, 1000, "G", "A"), ("2", 300, 1100, "G", "A")]
    data = _cache(tmp_path, places, reverse=reverse)
    p1, p2 = (200, 300) if build == "GRCh37" else (1000, 1100)
    lengths = (249250621, 243199373) if build == "GRCh37" else (248956422, 242193529)
    source = tmp_path / "mixed.vcf"
    rows = [_row(100, "0/0", chrom="MT"), _row(p2, "1/1", chrom="2")]
    rows += [_row(p, "0/1", chrom="2", name="rs900252299") for p in range(p2 + 1, p2 + 1101)]
    rows.append(_row(p1, "0/1", chrom="1"))
    source.write_text(
        "##fileformat=VCFv4.2\n"
        + f"##contig=<ID=MT,length={length}>\n"
        + f"##contig=<ID=1,length={lengths[0]}>\n"
        + f"##contig=<ID=2,length={lengths[1]}>\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(rows)
    )
    compressed = source.with_suffix(".vcf.gz")
    pysam.tabix_compress(str(source), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    expected = ("MT", 100, "G", "G") if length == 16569 else ("1", p1, "G", "A")
    call, _ = _extract(monkeypatch, compressed if indexed else source, data, build)
    assert (call.chromosome, call.position, call.allele1, call.allele2) == expected

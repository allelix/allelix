# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extract keeps malformed-row, build-exclusion and header-order contracts."""

from __future__ import annotations

import pysam
import pytest
from test_extract_reference_calls import RSID, _cache, _extract, _input, _row


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("columns", [1, 2, 4, 7])
def test_extract_skips_short_rows_before_a_valid_call(
    tmp_path, monkeypatch, build, pos, compressed, columns
):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    malformed = "\t".join(["1", "20", ".", "G", "A", ".", "PASS"][:columns]) + "\n"
    source = _input(tmp_path, build, [malformed, _row(pos, "0/1")], indexed=False)
    if not compressed:
        source = source.with_suffix("")
    call, text = _extract(monkeypatch, source, data, build)
    assert (call.position, call.allele1, call.allele2) == (pos, "G", "A")
    assert "malformed" in text.lower()


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("cache_reverse", [False, True])
@pytest.mark.parametrize("header_reverse", [False, True])
@pytest.mark.parametrize("body_reverse", [False, True])
def test_extract_obeys_header_order_when_body_contig_order_differs(
    tmp_path, monkeypatch, build, cache_reverse, header_reverse, body_reverse
):
    placements = [("1", 200, 1000, "G", "A"), ("2", 300, 1100, "G", "A")]
    data = _cache(tmp_path, placements, reverse=cache_reverse)
    body = list(reversed(placements)) if body_reverse else placements
    header = ("2", "1") if header_reverse else ("1", "2")
    rows = []
    for number, (chrom, pos37, pos38, _ref, _alt) in enumerate(body):
        pos = pos37 if build == "GRCh37" else pos38
        rows.append(_row(pos, "0/1" if chrom == "1" else "1/1", chrom=chrom))
        if number == 0:
            # Exceed the resolver's batch so its read-ahead cannot mask early stopping.
            rows.extend(
                _row(p, "0/1", chrom=chrom, name="rs900252299") for p in range(pos + 1, pos + 1101)
            )
    # Body is sorted within each contig and contigs stay together: a valid real index.
    source = _input(tmp_path, build, rows, indexed=True, contigs=header)
    indexed, _ = _extract(monkeypatch, source, data, build)
    source.with_suffix(source.suffix + ".tbi").unlink()
    sequential, _ = _extract(monkeypatch, source, data, build)
    expected_chrom = header[0]
    expected_pos = next(p[1 if build == "GRCh37" else 2] for p in placements if p[0] == header[0])
    for call in (indexed, sequential):
        assert (call.chromosome, call.position) == (expected_chrom, expected_pos)
        assert (call.allele1, call.allele2) == (
            ("G", "A") if expected_chrom == "1" else ("A", "A")
        )


@pytest.mark.parametrize(
    "build,chrom,length,excluded",
    [
        ("GRCh37", "chrM", None, True),
        ("GRCh37", "M", None, True),
        ("GRCh37", "chrM", 16571, True),
        ("GRCh38", "chrM", 16571, True),
        ("GRCh37", "MT", 16571, True),
        ("GRCh37", "chrM", 16569, False),
        ("GRCh37", "MT", 16569, False),
        ("GRCh38", "MT", 16569, False),
    ],
)
@pytest.mark.parametrize("shape", ["reference", "unnamed-reference", "block", "non-reference"])
@pytest.mark.parametrize("indexed", [False, True])
def test_extract_reference_candidates_obey_mitochondrial_build_exclusions(
    tmp_path, monkeypatch, build, chrom, length, excluded, shape, indexed
):
    data = _cache(tmp_path, [("MT", 200, 1000, "G", "A")])
    pos = 200 if build == "GRCh37" else 1000
    source = tmp_path / "mito.vcf"
    declaration = "" if length is None else f"##contig=<ID={chrom},length={length}>\n"
    row = _row(
        pos - 10 if shape == "block" else pos,
        "0/1" if shape == "non-reference" else "0/0",
        chrom=chrom,
        name="." if shape in ("block", "unnamed-reference") else RSID,
        ref="A" if shape == "block" else "G",
        alt="<NON_REF>" if shape == "block" else "A",
        info=f"END={pos + 10}" if shape == "block" else ".",
    )
    source.write_text(
        "##fileformat=VCFv4.2\n"
        + declaration
        + '##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + row
    )
    if indexed:
        compressed = source.with_suffix(source.suffix + ".gz")
        pysam.tabix_compress(str(source), str(compressed))
        pysam.tabix_index(str(compressed), preset="vcf")
        source = compressed
    call, _text = _extract(monkeypatch, source, data, build, missing=excluded)
    if not excluded:
        assert (call.chromosome, call.position) == ("MT", pos)
        assert (call.allele1, call.allele2) == (
            ("G", "A") if shape == "non-reference" else ("G", "G")
        )

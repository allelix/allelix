# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Several gnomAD-mapped placements: indexed and sequential extract agree (ADR-0054).

When gnomAD maps one rsID to more than one placement and the file holds a
record at each, both paths prefer usable mapped calls, then header contig order, position
and input line order. Cache insertion order never decides. A mapped
no-call remains preferred over an unmapped call. All rows are synthetic.
"""

from __future__ import annotations

import contextlib
import sqlite3

import pysam
import pytest
from click.testing import CliRunner

from allelix.cli import main, utility
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA

RSID = "rs900252201"
_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=248956422>\n"
    "##contig=<ID=2,length=242193529>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def _cache(tmp_path, placements):
    data = tmp_path / "data"
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
            "VALUES (?,?,'G','A',?,0.01)",
            [(chrom, pos, RSID) for chrom, pos in placements],
        )
        conn.execute(
            "INSERT INTO database_versions VALUES "
            "('gnomad','test://synthetic','1','2026-10-25',?,NULL,?)",
            (len(placements), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()
    return data


def _input(tmp_path, rows):
    path = tmp_path / "input.vcf"
    path.write_text(
        _HEADER
        + "".join(
            f"{chrom}\t{pos}\t{RSID}\tG\tA\t99\tPASS\t.\tGT:GQ\t{gt}:99\n"
            for chrom, pos, gt in rows
        )
    )
    gz = tmp_path / "input.vcf.gz"
    pysam.tabix_compress(str(path), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)
    return gz


def _both(monkeypatch, data, source):
    """The extracted record's locus and no-call flag, with and without the index."""
    captures = []
    monkeypatch.setattr(
        utility,
        "_render_extract_table",
        lambda _path, _wanted, found: captures.append(
            {k: (v.chromosome, v.position, v.is_no_call) for k, v in found.items()}
        ),
    )
    args = ["extract", str(source), "--snps", RSID, "--data-dir", str(data), "--build", "grch38"]
    assert utility._try_tabix_extract(source, {RSID}, None, data, "GRCh38") is not None
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    index = source.with_suffix(source.suffix + ".tbi")
    hidden = index.with_suffix(".off")
    index.rename(hidden)
    try:
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output
    finally:
        hidden.rename(index)
    assert len(captures) == 2
    return captures


@pytest.mark.parametrize("cache_order", ["ascending", "descending"])
@pytest.mark.parametrize(
    "rows,expected",
    [
        ([("1", 100, "0/1"), ("1", 200, "1/1")], ("1", 100, False)),
        ([("1", 100, "./."), ("1", 200, "1/1")], ("1", 200, False)),
        ([("1", 100, "0/1"), ("1", 200, "./.")], ("1", 100, False)),
        ([("1", 200, "1/1"), ("2", 50, "0/1")], ("1", 200, False)),
    ],
    ids=["two-calls", "mapped-no-call-first", "mapped-no-call-last", "two-chromosomes"],
)
def test_usable_call_then_input_order_breaks_a_tie_between_mapped_placements(
    tmp_path, monkeypatch, cache_order, rows, expected
):
    placements = sorted({(chrom, pos) for chrom, pos, _gt in rows})
    if cache_order == "descending":
        placements.reverse()
    data = _cache(tmp_path, placements)
    indexed, sequential = _both(monkeypatch, data, _input(tmp_path, rows))
    assert indexed == sequential
    assert indexed[RSID] == expected

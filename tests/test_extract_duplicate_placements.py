# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Indexed and sequential extract choose the same duplicate placement.

Among records yielded by the parser, prefer a gnomAD-mapped placement;
otherwise keep the first matching record in file order. A mapped no-call
is still a candidate. Confident 0/0 calls participate with their REF genotype.
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


def _cache(tmp_path, mapped):
    data = tmp_path / "data"
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            ("1", mapped, "G", "A", RSID, 0.01),
        )
        conn.execute(
            "INSERT INTO database_versions VALUES ('gnomad','test://synthetic','1','2026-10-25',1,NULL,?)",
            (f"sv:{GNOMAD_SCHEMA_VERSION}",),
        )
        conn.commit()
    return data


def _input(tmp_path, rows, samples="S"):
    path = tmp_path / "input.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n##contig=<ID=chr1,length=248956422>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        + samples
        + "\n"
        + "".join(rows)
    )
    gz = tmp_path / "input.vcf.gz"
    pysam.tabix_compress(str(path), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)
    return gz


def _row(pos, gt, *, rsid=RSID, chrom="1", second=None):
    row = f"{chrom}\t{pos}\t{rsid}\tG\tA\t99\tPASS\t.\tGT:GQ\t{gt}:99"
    return row + (f"\t{second}:99" if second else "") + "\n"


def _pack(found):
    return {
        k: (v.chromosome, v.position, v.allele1, v.allele2, v.ref, v.input_alt, v.is_no_call)
        for k, v in found.items()
    }


def _both(tmp_path, monkeypatch, data, source, *, sample=None, wanted=RSID):
    captures = []
    monkeypatch.setattr(
        utility,
        "_render_extract_table",
        lambda _path, _wanted, found: captures.append(_pack(found)),
    )
    args = ["extract", str(source), "--snps", wanted, "--data-dir", str(data), "--build", "grch38"]
    if sample:
        args.extend(["--sample", sample])
    # Prove the real .tbi is usable; no fake tabix class or mocked fetch.
    fast = utility._try_tabix_extract(source, set(wanted.split(",")), sample, data, "GRCh38")
    assert fast is not None
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    index = source.with_suffix(source.suffix + ".tbi")
    hidden = index.with_suffix(".off")
    index.rename(hidden)
    try:
        assert utility._try_tabix_extract(source, {RSID}, sample, data, "GRCh38") is None
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output
    finally:
        hidden.rename(index)
    assert len(captures) == 2
    return captures


@pytest.mark.parametrize("chrom", ["1", "chr1"])
@pytest.mark.parametrize("gt", ["0/1", "1/1", "1|0", "1|1", "0/0"])
def test_later_mapped_placement_has_the_same_result_with_and_without_index(
    tmp_path, monkeypatch, chrom, gt
):
    data = _cache(tmp_path, 200)
    src = _input(tmp_path, [_row(100, "./.", chrom=chrom), _row(200, gt, chrom=chrom)])
    indexed, sequential = _both(tmp_path, monkeypatch, data, src)
    assert indexed == sequential
    assert indexed[RSID][1] == 200
    assert indexed[RSID][-1] is False
    if gt == "0/0":
        assert indexed[RSID][2:4] == ("G", "G")


@pytest.mark.parametrize(
    "mapped,first,last,expected,nocall",
    [
        (100, "./.", "1/1", 100, True),
        (100, "0/1", "1/1", 100, False),
        (200, "1/1", "./.", 200, True),
        (300, "./.", "1/1", 100, True),
    ],
)
def test_duplicate_policy_keeps_mapped_no_calls_and_file_order_fallback(
    tmp_path, monkeypatch, mapped, first, last, expected, nocall
):
    data = _cache(tmp_path, mapped)
    src = _input(tmp_path, [_row(100, first), _row(200, last)])
    indexed, sequential = _both(tmp_path, monkeypatch, data, src)
    assert indexed == sequential
    assert indexed[RSID][1] == expected
    assert indexed[RSID][-1] is nocall


@pytest.mark.parametrize("sample,expected,nocall", [("S", 200, False), ("SECOND", 200, False)])
def test_duplicate_selection_uses_the_same_sample_in_both_paths(
    tmp_path, monkeypatch, sample, expected, nocall
):
    data = _cache(tmp_path, 200)
    src = _input(
        tmp_path,
        [_row(100, "./.", second="1/1"), _row(200, "0/1", second="0/0")],
        samples="S\tSECOND",
    )
    indexed, sequential = _both(tmp_path, monkeypatch, data, src, sample=sample)
    assert indexed == sequential
    assert indexed[RSID][1] == expected
    assert indexed[RSID][-1] is nocall
    if sample == "SECOND":
        assert indexed[RSID][2:4] == ("G", "G")


def test_pending_unmapped_request_does_not_change_the_mapped_choice(tmp_path, monkeypatch):
    data = _cache(tmp_path, 200)
    other = "rs900252202"
    src = _input(tmp_path, [_row(100, "./."), _row(200, "1/1"), _row(300, "0/1", rsid=other)])
    indexed, sequential = _both(tmp_path, monkeypatch, data, src, wanted=RSID + "," + other)
    assert indexed == sequential
    assert indexed[RSID][1] == 200
    assert indexed[other][1] == 300

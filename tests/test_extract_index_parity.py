# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Mixed indexed/fallback lookup parity with synthetic VCFs and caches."""

from __future__ import annotations

import contextlib
import gzip
import sqlite3
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from allelix.cli import main, utility
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.parsers.vcf import VcfParser, _parse_data_line


def _row(pos, rsid, *, chrom="1", ref="A", alt="G", gt="0/1:99:10", filt="PASS", info="."):
    return f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t{filt}\t{info}\tGT:GQ:DP\t{gt}\n"


def _write(path, rows, *, samples=("S",), contigs=None):
    if contigs is None:
        contigs = list(dict.fromkeys(row.split("\t")[0] for row in rows))
    text = "##fileformat=VCFv4.2\n"
    for chrom in contigs:
        length = 16569 if chrom in ("M", "MT", "chrM", "chrMT") else 248956422
        text += f"##contig=<ID={chrom},length={length}>\n"
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    text += "".join(rows)
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        path.write_text(text)
    return path


def _cache(path, coords):
    path.mkdir(exist_ok=True)
    with contextlib.closing(sqlite3.connect(path / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,0.1)",
            coords,
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
            "VALUES ('gnomad','test://','synthetic','2001-01-01',?,'sv:1')",
            (len(coords),),
        )
        conn.commit()
    return path


class _Index:
    def __init__(self, path):
        opener = gzip.open if Path(path).suffix == ".gz" else open
        with opener(path, "rt") as handle:
            self.rows = [line.rstrip() for line in handle if not line.startswith("#")]
        self.contigs = tuple(dict.fromkeys(row.split("\t")[0] for row in self.rows))
        self.queries = []
        self.closed = False

    def fetch(self, chrom, start, stop):
        self.queries.append((chrom, start, stop))
        if chrom not in self.contigs:
            raise ValueError("contig absent from index")
        return [
            row
            for row in self.rows
            if row.split("\t")[0] == chrom
            and int(row.split("\t")[1]) - 1 < stop
            and int(row.split("\t")[1]) + len(row.split("\t")[3]) - 1 > start
        ]

    def close(self):
        self.closed = True


def _pack(found):
    return {
        rsid: (v.chromosome, v.position, v.allele1, v.allele2, v.ref, v.build, v.haploid)
        for rsid, v in found.items()
    }


@pytest.mark.parametrize("mapped", ["hit", "miss", "empty"])
def test_all_requested_but_unfound_ids_share_one_fallback_pass(tmp_path, monkeypatch, mapped):
    path = _write(tmp_path / "mixed.vcf", [_row(100, "rs900500001"), _row(200, "rs900500002")])
    data = _cache(
        tmp_path / "data", [("1", 100 if mapped == "hit" else 110, "A", "G", "rs900500001")]
    )
    if mapped == "empty":
        from allelix.annotators.gnomad import GnomadAnnotator

        monkeypatch.setattr(
            GnomadAnnotator, "bulk_resolve_coordinates", lambda *a: {"rs900500001": []}
        )
    wanted = {"rs900500001", "rs900500002", "rs900500099"}
    original = utility._sequential_extract
    calls = []

    def sequential(parser, source, pending, *args):
        calls.append(set(pending))
        return original(parser, source, pending, *args)

    monkeypatch.setattr(utility, "_sequential_extract", sequential)
    indexed = utility._execute_tabix_extract(
        path, wanted, None, data, SimpleNamespace(TabixFile=_Index), "GRCh38"
    )
    assert indexed is not None
    assert _pack(indexed) == _pack(original(VcfParser(), path, wanted, data, "GRCh38"))
    assert calls == [wanted - {"rs900500001"} if mapped == "hit" else wanted]


def test_pending_fallback_preserves_an_indexed_no_call(tmp_path, monkeypatch):
    path = _write(
        tmp_path / "mixed.vcf",
        [
            _row(100, "rs900500001", gt="./.:99:10"),
            _row(110, "rs900500001", gt="1/1:99:10"),
            _row(200, "rs900500002"),
        ],
    )
    data = _cache(tmp_path / "data", [("1", 100, "A", "G", "rs900500001")])
    original = utility._sequential_extract
    pending = []

    def sequential(parser, source, wanted, *args):
        pending.append(set(wanted))
        return original(parser, source, wanted, *args)

    monkeypatch.setattr(utility, "_sequential_extract", sequential)
    found = utility._execute_tabix_extract(
        path,
        {"rs900500001", "rs900500002"},
        None,
        data,
        SimpleNamespace(TabixFile=_Index),
        "GRCh38",
    )
    assert found is not None and set(found) == {"rs900500001", "rs900500002"}
    assert found["rs900500001"].is_no_call and found["rs900500001"].position == 100
    assert pending == [{"rs900500002"}]


@pytest.mark.parametrize(
    "actual,query",
    [("1", "chr1"), ("chr1", "1"), ("M", "MT"), ("MT", "chrM"), ("chrM", "MT"), ("chrMT", "M")],
)
def test_fetch_uses_actual_index_contigs_in_both_alias_directions(tmp_path, actual, query):
    path = _write(tmp_path / "alias.vcf", [_row(100, "rs900500001", chrom=actual)])
    index = _Index(path)
    found = utility._fetch_rows(index, query, 100, 0, _parse_data_line)
    assert len(found) == 1 and found[0].rsid == "rs900500001"
    assert index.queries == [(actual, 99, 100)]


def test_missing_index_contig_does_not_guess_a_name(tmp_path):
    path = _write(tmp_path / "one.vcf", [_row(100, "rs900500001")], contigs=["1", "chr2"])
    index = _Index(path)
    assert utility._fetch_rows(index, "2", 100, 0, _parse_data_line) == []
    assert index.queries == []


def test_alias_traversal_follows_index_order_and_filters_overlaps(tmp_path):
    path = _write(
        tmp_path / "aliases.vcf",
        [
            _row(99, "rs900500099", chrom="chrM", ref="AA", alt="A"),
            _row(100, "rs900500001", chrom="chrM"),
            _row(100, "rs900500002", chrom="MT"),
        ],
    )
    index = _Index(path)
    found = utility._fetch_rows(index, "M", 100, 0, _parse_data_line)
    assert [v.rsid for v in found] == ["rs900500001", "rs900500002"]
    assert index.queries == [("chrM", 99, 100), ("MT", 99, 100)]


@pytest.mark.parametrize(
    "gt,alt,filt",
    [
        ("0/1:99:10", "G", "PASS"),
        ("0/1:99:10", "G", "q10"),
        ("0/1:3:10", "G", "PASS"),
        ("0/1:99:0", "G", "PASS"),
        ("./1:99:10", "G", "PASS"),
        ("./.:99:10", "G", "PASS"),
        ("1:99:10", "G", "PASS"),
        ("1/2:99:10", "G,C", "PASS"),
        ("0/0:99:10", "G", "PASS"),
        ("0/1:99:10", "<DEL>", "PASS"),
        ("0/0:99:10", "<NON_REF>", "PASS"),
    ],
)
def test_mixed_indexed_fallback_preserves_quality_rules(tmp_path, gt, alt, filt):
    path = _write(
        tmp_path / "quality.vcf",
        [
            _row(100, "rs900500001", gt=gt, alt=alt, filt=filt, info="END=105"),
            _row(200, "rs900500002", gt=gt, alt=alt, filt=filt, info="END=205"),
            _row(300, "rs900500090"),
        ],
    )
    data = _cache(tmp_path / "data", [("1", 100, "A", "G", "rs900500001")])
    wanted = {"rs900500001", "rs900500002"}
    indexed = utility._execute_tabix_extract(
        path, wanted, None, data, SimpleNamespace(TabixFile=_Index), "GRCh38"
    )
    assert indexed is not None
    assert _pack(indexed) == _pack(
        utility._sequential_extract(VcfParser(), path, wanted, data, "GRCh38")
    )


def test_repeated_interval_queries_are_fetched_once(tmp_path):
    path = _write(
        tmp_path / "same-site.vcf", [_row(100, "rs900500001"), _row(100, "rs900500002", alt="C")]
    )
    data = _cache(
        tmp_path / "data",
        [("1", 100, "A", alt, r) for alt, r in (("G", "rs900500001"), ("C", "rs900500002"))],
    )
    opened = []

    def index(source):
        obj = _Index(source)
        opened.append(obj)
        return obj

    found = utility._execute_tabix_extract(
        path,
        {"rs900500001", "rs900500002"},
        None,
        data,
        SimpleNamespace(TabixFile=index),
        "GRCh38",
    )
    assert found is not None and len(found) == 2
    assert opened[0].queries.count(("1", 99, 100)) == 1
    assert opened[0].closed


@pytest.mark.parametrize("placed", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
def test_cli_plain_and_indexed_mixed_lookup_have_identical_tables(
    tmp_path, monkeypatch, placed, indexed
):
    rows = [_row(100, "rs900500001"), _row(200, "rs900500002")]
    plain = _write(tmp_path / "same.vcf", rows)
    compressed = _write(tmp_path / "same.vcf.gz", rows)
    if indexed:
        compressed.with_suffix(".gz.tbi").write_bytes(b"synthetic index double")
    data = _cache(
        tmp_path / "data", [("1", 100, "A", "G", "rs900500001" if placed else "rs900500088")]
    )
    monkeypatch.setitem(
        sys.modules, "pysam", SimpleNamespace(TabixFile=_Index, set_verbosity=lambda _level: 0)
    )
    outputs = []
    for path in (plain, compressed):
        result = CliRunner(env={"COLUMNS": "140"}).invoke(
            main,
            [
                "extract",
                str(path),
                "--snps",
                "rs900500001,rs900500002,rs900500099",
                "--data-dir",
                str(data),
                "--build",
                "grch38",
            ],
        )
        assert result.exit_code == 0, result.output
        outputs.append(" ".join(result.output.replace(path.name, "source").split()))
    assert outputs[0] == outputs[1]
    assert outputs[1].count("rs900500001") == outputs[1].count("rs900500002") == 1


@pytest.mark.parametrize("sample", [None, "absent", "S2"])
def test_mixed_indexed_sample_binding_matches_stream(tmp_path, sample):
    rows = [
        _row(pos, rsid).rstrip() + "\t1/1:99:10\n"
        for pos, rsid in ((100, "rs900500001"), (200, "rs900500002"))
    ]
    path = _write(tmp_path / "samples.vcf", rows, samples=("S1", "S2"))
    data = _cache(tmp_path / "data", [("1", 100, "A", "G", "rs900500001")])
    wanted = {"rs900500001", "rs900500002"}
    found = utility._execute_tabix_extract(
        path, wanted, sample, data, SimpleNamespace(TabixFile=_Index), "GRCh38"
    )
    if sample != "S2":
        assert found is None
    else:
        assert found is not None
        assert _pack(found) == _pack(
            utility._sequential_extract(VcfParser(sample="S2"), path, wanted, data, "GRCh38")
        )
        assert all(v.allele1 == v.allele2 == "G" for v in found.values())


def test_indexed_rsidless_hit_and_unplaced_named_hit_match_stream(tmp_path):
    path = _write(tmp_path / "rsidless.vcf", [_row(100, "."), _row(200, "rs900500002")])
    data = _cache(tmp_path / "data", [("1", 100, "A", "G", "rs900500001")])
    wanted = {"rs900500001", "rs900500002"}
    found = utility._execute_tabix_extract(
        path, wanted, None, data, SimpleNamespace(TabixFile=_Index), "GRCh38"
    )
    assert found is not None and set(found) == wanted
    assert _pack(found) == _pack(
        utility._sequential_extract(VcfParser(), path, wanted, data, "GRCh38")
    )


@pytest.mark.parametrize("actual", ["1", "chr1", "M", "MT", "chrM", "chrMT"])
def test_real_tabix_mixed_lookup_matches_stream(tmp_path, actual):
    pysam = pytest.importorskip("pysam")
    plain = _write(
        tmp_path / "real.vcf",
        [
            _row(100, "rs900500001", chrom=actual),
            _row(200, "rs900500002", chrom=actual),
        ],
    )
    compressed = tmp_path / "real.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    canonical = "MT" if actual in ("M", "MT", "chrM", "chrMT") else "1"
    data = _cache(tmp_path / "data", [(canonical, 100, "A", "G", "rs900500001")])
    wanted = {"rs900500001", "rs900500002"}
    indexed = utility._try_tabix_extract(compressed, wanted, None, data, "GRCh38")
    assert indexed is not None
    assert _pack(indexed) == _pack(
        utility._sequential_extract(VcfParser(), compressed, wanted, data, "GRCh38")
    )

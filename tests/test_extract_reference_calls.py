# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extract explicit/reference-block calls with one placement and quality policy.

Mapped beats unmapped; usable mapped beats rejected/no-call mapped; ties use
header contig order, numeric POS, then line order. No cache insertion tie-break.
GRCh37 uses its identity placements, never GRCh38 gnomAD coordinates. Presence
of a real index must not change the result; GRCh37 may use sequential fallback.
All variants and database rows are synthetic.
"""

from __future__ import annotations

import contextlib
import sqlite3

import pysam
import pytest
from click.testing import CliRunner

from allelix.cli import main, utility
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GNOMAD_SCHEMA

RSID = "rs900252291"


def _cache(root, placements, *, reverse=False):
    """Each placement is chromosome, GRCh37 POS, GRCh38 POS, REF, ALT."""
    data = root / "data"
    data.mkdir()
    ordered = list(reversed(placements)) if reverse else list(placements)
    write_identity_db(
        [(RSID, chrom, p37, ref, alt, p38, ref, alt, 0) for chrom, p37, p38, ref, alt in ordered],
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://synthetic",
        downloaded_at="2001-01-01",
    )
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,0.01)",
            [(chrom, p38, ref, alt, RSID) for chrom, _p37, p38, ref, alt in ordered],
        )
        conn.execute(
            "INSERT INTO database_versions VALUES "
            "('gnomad','test://synthetic','1','2026-10-25',?,NULL,?)",
            (len(ordered), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()
    return data


def _row(
    pos,
    gt="0/0",
    *,
    chrom="1",
    name=RSID,
    ref="G",
    alt="A",
    filt="PASS",
    gq=50,
    dp=20,
    info=".",
    second=None,
):
    sample = f"{gt}:{gq}:{dp}"
    if second is not None:
        sample += f"\t{second}:50:20"
    return f"{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t{filt}\t{info}\tGT:GQ:DP\t{sample}\n"


def _input(root, build, rows, *, indexed, contigs=("1", "2"), samples="S"):
    lengths = {
        "GRCh37": {"1": 249250621, "2": 243199373},
        "GRCh38": {"1": 248956422, "2": 242193529},
    }
    source = root / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        + "".join(f"##contig=<ID={c},length={lengths[build][c]}>\n" for c in contigs)
        + '##ALT=<ID=NON_REF,Description="Any non-reference allele">\n'
        + '##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
        + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + samples
        + "\n"
        + "".join(rows)
    )
    compressed = root / "input.vcf.gz"
    pysam.tabix_compress(str(source), str(compressed))
    if indexed:
        pysam.tabix_index(str(compressed), preset="vcf")
        with pysam.TabixFile(str(compressed)) as index:
            assert list(index.fetch(contigs[0]))
    return compressed


def _extract(monkeypatch, source, data, build, *, sample=None, missing=False):
    captured = []
    original = utility._render_extract_table

    def render(path, wanted, found):
        captured.append(dict(found))
        original(path, wanted, found)

    monkeypatch.setattr(utility, "_render_extract_table", render)
    args = [
        "extract",
        str(source),
        "--snps",
        RSID,
        "--data-dir",
        str(data),
        "--build",
        build.lower(),
    ]
    if sample:
        args += ["--sample", sample]
    result = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(captured) == 1
    if missing:
        assert RSID not in captured[0], result.output
        return None, result.output
    assert RSID in captured[0], result.output
    return captured[0][RSID], result.output


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "index-present"])
@pytest.mark.parametrize(
    "shape",
    [
        "reference",
        "phased-reference",
        "alt-dot",
        "bogus-before-reference",
        "filtered-reference",
        "low-gq-reference",
        "zero-dp-reference",
        "block-interior",
        "block-at-end",
        "low-gq-block",
        "filtered-block",
    ],
)
def test_extract_includes_reference_calls(tmp_path, monkeypatch, build, pos, indexed, shape):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    no_call = shape in {
        "filtered-reference",
        "low-gq-reference",
        "zero-dp-reference",
        "low-gq-block",
        "filtered-block",
    }
    if "block" in shape:
        # The block's literal anchor is A; the requested interior site's reference is G.
        # Copying the block's anchor base would silently report the wrong genotype.
        end = pos if shape == "block-at-end" else pos + 10
        rows = [
            _row(
                pos - 10,
                name=".",
                ref="A",
                alt="<NON_REF>",
                info=f"END={end}",
                gq=5 if shape == "low-gq-block" else 50,
                filt="LowQual" if shape == "filtered-block" else "PASS",
            )
        ]
    else:
        rows = [
            _row(
                pos,
                gt="0|0" if shape == "phased-reference" else "0/0",
                alt="." if shape == "alt-dot" else "A",
                gq=5 if shape == "low-gq-reference" else 50,
                dp=0 if shape == "zero-dp-reference" else 20,
                filt="LowQual" if shape == "filtered-reference" else "PASS",
            )
        ]
        if shape == "bogus-before-reference":
            rows.insert(0, _row(100, "./."))
    source = _input(tmp_path, build, rows, indexed=indexed)
    call, text = _extract(monkeypatch, source, data, build)
    assert (call.chromosome, call.position) == ("1", pos)
    assert call.is_no_call is no_call
    assert call.is_heterozygous is False
    if no_call:
        assert call.allele1 not in {"A", "C", "G", "T"}
        assert "yes" in next(line for line in text.splitlines() if RSID in line)
    else:
        assert (call.allele1, call.allele2) == ("G", "G")
        assert "G/G" in text
        row = next(line for line in text.splitlines() if RSID in line)
        assert row.count("no") >= 2  # Het? no, No-call? no in the actual terminal table.


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("cache_reverse", [False, True])
@pytest.mark.parametrize("file_reverse", [False, True])
@pytest.mark.parametrize(
    "calls",
    [
        "both-usable",
        "first-no-call",
        "first-filtered",
        "first-low-gq",
        "both-usable-alt",
        "first-no-call-alt",
    ],
)
def test_mapped_quality_then_header_order_not_cache_order(
    tmp_path, monkeypatch, build, indexed, cache_reverse, file_reverse, calls
):
    placements = [("1", 200, 1000, "G", "A"), ("2", 300, 1100, "G", "A")]
    data = _cache(tmp_path, placements, reverse=cache_reverse)
    order = list(reversed(placements)) if file_reverse else placements
    rows = []
    for i, (chrom, p37, p38, _ref, _alt) in enumerate(order):
        pos = p37 if build == "GRCh37" else p38
        rows.append(
            _row(
                pos,
                (
                    "./."
                    if i == 0 and calls in ("first-no-call", "first-no-call-alt")
                    else ("0/1" if i == 0 else "1/1")
                    if calls.endswith("-alt")
                    else "0/0"
                ),
                chrom=chrom,
                filt="LowQual" if i == 0 and calls == "first-filtered" else "PASS",
                gq=5 if i == 0 and calls == "first-low-gq" else 50,
            )
        )
    source = _input(
        tmp_path, build, rows, indexed=indexed, contigs=tuple(item[0] for item in order)
    )
    call, _text = _extract(monkeypatch, source, data, build)
    selected = order[0 if calls in ("both-usable", "both-usable-alt") else 1]
    assert (call.chromosome, call.position) == (
        selected[0],
        selected[1 if build == "GRCh37" else 2],
    )
    assert not call.is_no_call
    if calls == "both-usable-alt":
        assert {call.allele1, call.allele2} == {"G", "A"}
    elif calls == "first-no-call-alt":
        assert (call.allele1, call.allele2) == ("A", "A")
    else:
        assert (call.allele1, call.allele2) == ("G", "G")


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
def test_line_order_breaks_a_tie_at_the_same_mapped_position(
    tmp_path, monkeypatch, build, pos, indexed, reverse
):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    rows = [_row(pos, "0/0"), _row(pos, "1/1")]
    if reverse:
        rows.reverse()
    source = _input(tmp_path, build, rows, indexed=indexed)
    call, _text = _extract(monkeypatch, source, data, build)
    allele = "A" if reverse else "G"
    assert (call.allele1, call.allele2) == (allele, allele)
    assert not call.is_no_call


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True])
def test_a_mapped_no_call_still_beats_an_unmapped_usable_call(
    tmp_path, monkeypatch, build, pos, indexed
):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    source = _input(tmp_path, build, [_row(100, "1/1"), _row(pos, "./.")], indexed=indexed)
    call, _text = _extract(monkeypatch, source, data, build)
    assert call.position == pos
    assert call.is_no_call


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True])
def test_reference_call_uses_the_selected_sample(tmp_path, monkeypatch, build, pos, indexed):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    source = _input(
        tmp_path, build, [_row(pos, "1/1", second="0/0")], indexed=indexed, samples="S\tSECOND"
    )
    call, text = _extract(monkeypatch, source, data, build, sample="SECOND")
    assert (call.allele1, call.allele2) == ("G", "G")
    assert "G/G" in text


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("side", ["before", "after"])
def test_a_reference_block_does_not_cover_a_position_outside_its_span(
    tmp_path, monkeypatch, build, pos, indexed, side
):
    data = _cache(tmp_path, [("1", 200, 1000, "G", "A")])
    start, end = (pos - 20, pos - 1) if side == "before" else (pos + 1, pos + 20)
    source = _input(
        tmp_path,
        build,
        [_row(start, name=".", ref="A", alt="<NON_REF>", info=f"END={end}")],
        indexed=indexed,
    )
    _extract(monkeypatch, source, data, build, missing=True)


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("indexed", [False, True])
def test_explicit_reference_call_needs_no_reference_cache(
    tmp_path, monkeypatch, build, pos, indexed
):
    data = tmp_path / "data"
    data.mkdir()
    source = _input(tmp_path, build, [_row(pos)], indexed=indexed)
    call, _text = _extract(monkeypatch, source, data, build)
    assert (call.allele1, call.allele2) == ("G", "G")
    assert not call.is_no_call

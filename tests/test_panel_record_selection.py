# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Panel coverage uses a valid allele placement and the best usable observation."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest

from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import RSID_IDENTITY_SCHEMA_VERSION
from allelix.databases.schema import RSID_IDENTITY_SCHEMA
from allelix.models import Variant
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from allelix.reports.high_value import HighValueSNP, scan_no_calls
from allelix.reports.json_report import render_json

IDENTITY = [
    ("rs6025", "1", 169519049, "T", "C", 169549811, "C", "T", 0),
    ("rs900252101", "1", 1000, "A", "G", 2000, "A", "G", 0),
    ("rs45578832", "1", 11850974, "G", "A", 11790917, "G", "A", 0),
    ("rs751611452", "1", 11850974, "GCA", "G", 11790917, "GCA", "G", 0),
]


@pytest.fixture
def data(tmp_path):
    root = tmp_path / "data"
    root.mkdir()
    with contextlib.closing(sqlite3.connect(root / "rsid_identity.sqlite")) as conn:
        conn.executescript(RSID_IDENTITY_SCHEMA)
        conn.executemany("INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)", IDENTITY)
        conn.execute(
            "INSERT INTO database_versions VALUES ('rsid_identity','test://synthetic','1','2026-10-25',?,NULL,?)",
            (len(IDENTITY), f"sv:{RSID_IDENTITY_SCHEMA_VERSION}"),
        )
        conn.commit()
    return root


def _row(pos, rsid, ref, alt, gt, *, chrom="1", filt="PASS", gq=99, dp=20, second=None):
    text = f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t50\t{filt}\t.\tGT:GQ:DP\t{gt}:{gq}:{dp}"
    return text + (f"\t{second}" if second else "") + "\n"


def _run(tmp_path, data, rows, build, panel, *, sample=None, return_result=False):
    src = tmp_path / "input.vcf"
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
    if sample:
        header += "\tSECOND"
    src.write_text(header + "\n" + "".join(rows))
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(
            src,
            VcfParser(sample=sample),
            [],
            build_override=build,
            identity=identity,
            panel_rsids=frozenset(panel),
            high_value_rsids=set(panel) if return_result else None,
        )
    coverage = result.panel_coverage([])
    assert coverage["requested"] == len(panel)
    assert sum(coverage["counts"].values()) == len(panel)
    out = tmp_path / "out.json"
    render_json(result, output_path=out)
    assert json.loads(out.read_text())["panel_coverage"] == coverage
    if return_result:
        return result
    return {r["rsid"]: r for r in coverage["variants"]}


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("good_gt", ["0/1", "1/1"])
@pytest.mark.parametrize(
    "bad_gt,options",
    [
        ("0/1", {"filt": "LowQual"}),
        ("0/1", {"gq": 3}),
        ("0/1", {"dp": 0}),
        ("./.", {}),
    ],
)
def test_usable_panel_call_wins_over_a_rejected_duplicate(
    tmp_path, data, build, reverse, good_gt, bad_gt, options
):
    pos, ref, alt = (169519049, "T", "C") if build == "GRCh37" else (169549811, "C", "T")
    good = _row(pos, "rs6025", ref, alt, good_gt)
    bad = _row(pos, "rs6025", ref, alt, bad_gt, **options)
    rows = [good, bad] if reverse else [bad, good]
    site = _run(tmp_path, data, rows, build, {"rs6025"})["rs6025"]
    expected = f"{ref}/{alt}" if good_gt == "0/1" else f"{alt}/{alt}"
    assert (site["state"], site["genotype"], site["reason"]) == (
        "typed_no_annotation",
        expected,
        None,
    )


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "gt,options,reason",
    [
        ("0/1", {"filt": "LowQual"}, "filtered"),
        ("0/1", {"gq": 3}, "low_gq_call"),
        ("0/1", {"dp": 0}, "low_gq_call"),
        ("./.", {}, "no_call"),
    ],
)
def test_a_rejected_record_alone_remains_a_no_call(tmp_path, data, build, gt, options, reason):
    pos = 1000 if build == "GRCh37" else 2000
    site = _run(
        tmp_path, data, [_row(pos, "rs900252101", "A", "G", gt, **options)], build, {"rs900252101"}
    )["rs900252101"]
    assert (site["state"], site["reason"]) == ("no_call", reason)


def _assert_attribution(result, consumer, rsid, genotype):
    """Observe coverage or the scan's response to an independent no-call sentinel."""
    if consumer == "panel":
        site = next(r for r in result.panel_coverage([])["variants"] if r["rsid"] == rsid)
        expected = (
            ("typed_no_annotation", genotype, None) if genotype else ("not_on_chip", None, None)
        )
        assert (site["state"], site["genotype"], site["reason"]) == expected
    else:
        # The sentinel is not an input record. A call attributed to this rsID
        # suppresses its warning; a rejected record must leave it intact.
        sentinel = Variant(rsid, "1", 1, "-", "-")
        warnings = scan_no_calls(
            [*result.hv_variants, sentinel],
            {rsid: HighValueSNP(rsid, "SYNTH", "", "Synthetic site")},
            result.high_value_vcf_states,
        )
        assert [w.snp.rsid for w in warnings] == ([] if genotype else [rsid])


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt", ["0/0", "1/1"])
@pytest.mark.parametrize("mismatch", ["id", "ref", "alt", "position", "chromosome"])
def test_placement_and_ref_define_record_attribution(
    tmp_path, data, build, gt, mismatch, consumer
):
    pos = 1000 if build == "GRCh37" else 2000
    rsid, ref, alt, chrom = "rs900252101", "A", "G", "1"
    if mismatch == "id":
        rsid = "rs900252999"
    elif mismatch == "ref":
        ref = "C"
    elif mismatch == "alt":
        alt = "T"
    elif mismatch == "position":
        pos += 50
    else:
        chrom = "2"
    result = _run(
        tmp_path,
        data,
        [_row(pos, rsid, ref, alt, gt, chrom=chrom)],
        build,
        {"rs900252101"},
        return_result=True,
    )
    genotype = ("A/A" if gt == "0/0" else f"{alt}/{alt}") if mismatch in {"id", "alt"} else None
    _assert_attribution(result, consumer, "rs900252101", genotype)


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt", ["0/0", "1/1"])
@pytest.mark.parametrize("same_id", [True, False])
@pytest.mark.parametrize("cache", ["missing_rsid", "unavailable"])
def test_no_identity_placement_falls_back_to_id(
    tmp_path, data, build, gt, same_id, cache, consumer
):
    target = "rs900252999"
    if cache == "unavailable":
        (data / "rsid_identity.sqlite").unlink()
    record_id = target if same_id else "rs900252998"
    result = _run(
        tmp_path,
        data,
        [_row(9000, record_id, "C", "T", gt, chrom="2")],
        build,
        {target},
        return_result=True,
    )
    genotype = ("C/C" if gt == "0/0" else "T/T") if same_id else None
    _assert_attribution(result, consumer, target, genotype)


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "case",
    ["matching", "other_id", "i_id", "position", "chromosome", "no_placement", "outside_pair"],
)
def test_text_chip_attribution_keeps_rsid_matching(tmp_path, data, build, case, consumer):
    target = "rs900252101"
    pos = 1000 if build == "GRCh37" else 2000
    rsid, chrom = target, "1"
    if case == "other_id":
        rsid = "rs900252999"
    elif case == "i_id":
        rsid = "i900252999"
    elif case == "position":
        pos += 50
    elif case == "chromosome":
        chrom = "2"
    elif case == "no_placement":
        target = rsid = "rs900252999"
    genotype = "TT" if case == "outside_pair" else "GG"
    src = tmp_path / "chip.txt"
    src.write_text(
        "# This data file generated by 23andMe at: synthetic\n"
        "# rsid\tchromosome\tposition\tgenotype\n"
        f"{rsid}\t{chrom}\t{pos}\t{genotype}\n"
    )
    parser = TwentyThreeAndMeParser()
    assert parser.can_parse(src)
    parsed = list(parser.parse(src))
    assert len(parsed) == 1 and parsed[0].ref is None
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(
            src,
            parser,
            [],
            build_override=build,
            identity=identity,
            panel_rsids=frozenset({target}),
            high_value_rsids={target},
        )
    # Text-chip annotation can match by rsID despite mismatching coordinates
    # or an observed genotype outside the identity pair. Position alone does
    # not turn an unrelated rsID or an i-number into the target rsID.
    expected = None if case in {"other_id", "i_id"} else "/".join(genotype)
    _assert_attribution(result, consumer, target, expected)


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt,expected", [("0/0", "A/A"), ("1/2", "G/T"), ("2/2", "T/T")])
@pytest.mark.parametrize("reverse", [False, True])
def test_multiallelic_identity_rows_attribute_one_record(
    tmp_path, data, build, gt, expected, reverse, consumer
):
    rows = [IDENTITY[1], ("rs900252101", "1", 1000, "A", "T", 2000, "A", "T", 0)]
    with contextlib.closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn:
        conn.execute("DELETE FROM rsid_identity WHERE rsid='rs900252101'")
        conn.executemany(
            "INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)", rows[::-1] if reverse else rows
        )
        conn.commit()
    pos = 1000 if build == "GRCh37" else 2000
    result = _run(
        tmp_path,
        data,
        [_row(pos, "rs900252999", "A", "G,T", gt)],
        build,
        {"rs900252101"},
        return_result=True,
    )
    _assert_attribution(result, consumer, "rs900252101", expected)


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("lowercase", ["input", "identity"])
def test_reference_comparison_is_case_insensitive(tmp_path, data, build, lowercase, consumer):
    if lowercase == "identity":
        with contextlib.closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn:
            conn.execute("UPDATE rsid_identity SET ref37='a', ref38='a' WHERE rsid='rs900252101'")
            conn.commit()
    pos = 1000 if build == "GRCh37" else 2000
    result = _run(
        tmp_path,
        data,
        [_row(pos, "rs900252999", "a" if lowercase == "input" else "A", "G", "1/1")],
        build,
        {"rs900252101"},
        return_result=True,
    )
    _assert_attribution(result, consumer, "rs900252101", "G/G")


@pytest.mark.parametrize("consumer", ["panel", "high_value"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("second", [False, True])
def test_each_of_two_identity_placements_can_supply_the_call(
    tmp_path, data, build, second, consumer
):
    with contextlib.closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn:
        conn.execute(
            "INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)",
            ("rs900252101", "2", 3000, "C", "T", 4000, "C", "T", 0),
        )
        conn.commit()
    pos = (
        (3000 if build == "GRCh37" else 4000) if second else (1000 if build == "GRCh37" else 2000)
    )
    ref, alt, chrom = ("C", "T", "2") if second else ("A", "G", "1")
    result = _run(
        tmp_path,
        data,
        [_row(pos, "rs900252999", ref, alt, "1/1", chrom=chrom)],
        build,
        {"rs900252101"},
        return_result=True,
    )
    _assert_attribution(result, consumer, "rs900252101", f"{alt}/{alt}")


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("rsid", ["rs900252101", "."])
@pytest.mark.parametrize("gt,expected", [("0/0", "A/A"), ("0/1", "A/G"), ("1/1", "G/G")])
def test_the_matching_panel_placement_is_counted(tmp_path, data, build, rsid, gt, expected):
    pos = 1000 if build == "GRCh37" else 2000
    site = _run(tmp_path, data, [_row(pos, rsid, "A", "G", gt)], build, {"rs900252101"})[
        "rs900252101"
    ]
    assert (site["state"], site["genotype"], site["reason"]) == (
        "typed_no_annotation",
        expected,
        None,
    )


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt", ["0/0", "0/1", "1/1", "./."])
@pytest.mark.parametrize("rsid", ["rs751611452", "."])
def test_a_deletion_at_the_same_start_is_not_coverage_of_the_snp(tmp_path, data, build, gt, rsid):
    pos = 11850974 if build == "GRCh37" else 11790917
    sites = _run(
        tmp_path, data, [_row(pos, rsid, "GCA", "G", gt)], build, {"rs45578832", "rs751611452"}
    )
    snp = sites["rs45578832"]
    assert (snp["state"], snp["genotype"], snp["reason"]) == ("not_on_chip", None, None)
    assert sites["rs751611452"]["state"] == ("no_call" if gt == "./." else "typed_no_annotation")


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
def test_an_unrelated_deletion_cannot_override_the_real_snp_call(tmp_path, data, build, reverse):
    pos = 11850974 if build == "GRCh37" else 11790917
    rows = [_row(pos, "rs751611452", "GCA", "G", "./."), _row(pos, "rs45578832", "G", "A", "1/1")]
    if reverse:
        rows.reverse()
    site = _run(tmp_path, data, rows, build, {"rs45578832"})["rs45578832"]
    assert (site["state"], site["genotype"], site["reason"]) == (
        "typed_no_annotation",
        "A/A",
        None,
    )


@pytest.mark.parametrize("sample,expected", [("S", "G/G"), ("SECOND", "A/G")])
def test_duplicate_selection_uses_the_selected_sample(tmp_path, data, sample, expected):
    rows = [
        _row(2000, "rs900252101", "A", "G", "0/1", gq=3, second="0/1:99:20"),
        _row(2000, "rs900252101", "A", "G", "1/1", second="1/1:3:20"),
    ]
    site = _run(tmp_path, data, rows, "GRCh38", {"rs900252101"}, sample=sample)["rs900252101"]
    assert (site["state"], site["genotype"], site["reason"]) == (
        "typed_no_annotation",
        expected,
        None,
    )


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
def test_filter_file_cli_reports_the_usable_f5_call(tmp_path, data, build, reverse):
    from click.testing import CliRunner

    from allelix.cli import main

    pos, ref, alt = (169519049, "T", "C") if build == "GRCh37" else (169549811, "C", "T")
    rows = [
        _row(pos, "rs6025", ref, alt, "0/1", filt="LowQual"),
        _row(pos, "rs6025", ref, alt, "1/1"),
    ]
    if reverse:
        rows.reverse()
    _run(tmp_path, data, rows, build, {"rs6025"})
    panel = tmp_path / "panel.txt"
    panel.write_text("rs6025\n")
    output = tmp_path / "cli.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(tmp_path / "input.vcf"),
            "--build",
            build.lower(),
            "--filter-file",
            str(panel),
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    site = json.loads(output.read_text())["panel_coverage"]["variants"][0]
    assert (site["state"], site["genotype"], site["reason"]) == (
        "typed_no_annotation",
        f"{alt}/{alt}",
        None,
    )

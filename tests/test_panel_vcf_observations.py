# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Offline coverage observations must not become annotation inputs."""

from __future__ import annotations

import contextlib
import io
import json
import sqlite3
from pathlib import Path

import pytest
from rich.console import Console

from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import RSID_IDENTITY_SCHEMA_VERSION
from allelix.databases.schema import RSID_IDENTITY_SCHEMA
from allelix.models import Annotation
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from allelix.reports.diff import compute_diff, load_previous_report
from allelix.reports.json_report import render_json

PANEL = frozenset(f"rs99000{i}" for i in range(1, 6))
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def _write(tmp_path: Path, rows: str, header: str = HEADER) -> Path:
    path = tmp_path / "synthetic.vcf"
    path.write_text(header + rows)
    return path


def _row(i: int, gt: str, *, gq: int = 99, dp: int = 20, filt: str = "PASS") -> str:
    return f"1\t{i * 100}\trs99000{i}\tA\tG\t50\t{filt}\t.\tGT:GQ:DP\t{gt}:{gq}:{dp}\n"


def _run(path: Path, parser: VcfParser | None = None, panel: frozenset[str] | None = PANEL):
    return run_analysis(
        path, parser or VcfParser(), [], build_override="GRCh37", panel_rsids=panel
    )


def _sites(result):
    return {v["rsid"]: v for v in result.panel_coverage([])["variants"]}


def test_five_explicit_records_are_distinct_from_absence_and_findings(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        _row(1, "0/0") + _row(2, "0/1", gq=3) + _row(3, "0/1", filt="q10") + _row(4, "./."),
    )
    parser = VcfParser()
    result = _run(path, parser)
    assert result.annotations == result.hv_variants == []
    assert result.high_value_vcf_states is None
    assert parser.quality_stats.rsid_states == {}
    assert result.total_variants == 2  # Existing annotation stream remains unchanged.
    assert result.panel_coverage([])["counts"] == {
        "annotated": 0,
        "typed_no_annotation": 1,
        "no_call": 3,
        "not_on_chip": 1,
    }
    sites = _sites(result)
    assert result.genotyped_panel_rsids == PANEL - {"rs990005"}
    assert sites["rs990001"]["genotype"] == "A/A"
    assert sites["rs990001"]["reason"] is None
    assert [sites[f"rs99000{i}"]["reason"] for i in (2, 3, 4)] == [
        "low_gq_call",
        "filtered",
        "no_call",
    ]
    assert all(sites[f"rs99000{i}"]["genotype"] != "A/G" for i in (2, 3, 4))
    assert sites["rs990005"]["genotype"] is sites["rs990005"]["reason"] is None


@pytest.mark.parametrize(
    "options,expected",
    [
        ({}, ("no_call", "no_call", "no_call", "no_call")),
        ({"min_gq": 0}, ("typed_no_annotation", "no_call", "no_call", "no_call")),
        ({"keep_filtered": True}, ("no_call", "no_call", "typed_no_annotation", "no_call")),
    ],
)
def test_quality_options_reuse_parser_policy(tmp_path: Path, options, expected) -> None:
    path = _write(
        tmp_path,
        _row(1, "0/0", gq=3) + _row(2, "0/1", dp=0) + _row(3, "0/1", filt="q10") + _row(4, "./1"),
    )
    sites = _sites(_run(path, VcfParser(**options)))
    assert tuple(sites[f"rs99000{i}"]["state"] for i in range(1, 5)) == expected
    assert sites["rs990004"]["genotype"] == "G/."


@pytest.mark.parametrize(
    "first,second,gt,state,reason",
    [
        (_row(1, "0/0"), _row(1, "0/1", gq=3), "A/A", "typed_no_annotation", None),
        # 2.5.2: a usable record replaces an earlier filtered or no-call one.
        (_row(1, "0/1", filt="q10"), _row(1, "0/1"), "A/G", "typed_no_annotation", None),
        (_row(1, "./1"), _row(1, "0/0"), "A/A", "typed_no_annotation", None),
    ],
)
def test_duplicate_observation_is_one_consistent_record(
    tmp_path: Path,
    first,
    second,
    gt,
    state,
    reason,
) -> None:
    site = _sites(_run(_write(tmp_path, first + second)))["rs990001"]
    assert (site["genotype"], site["state"], site["reason"]) == (gt, state, reason)


def test_reuse_and_no_panel_reset_bounded_observations(tmp_path: Path) -> None:
    parser = VcfParser()
    first = _run(_write(tmp_path, _row(1, "0/0")), parser)
    assert _sites(first)["rs990001"]["state"] == "typed_no_annotation"
    result = _run(_write(tmp_path, "".join(_row(i, "0/0") for i in range(6, 200))), parser)
    assert result.genotyped_panel_rsids == frozenset()
    assert parser.panel_observations == {}
    assert len(parser.panel_observations) <= len(PANEL)
    assert _run(_write(tmp_path, _row(1, "0/0")), parser, None).panel_coverage([]) is None
    assert parser.panel_observations == {}


def test_selected_sample_is_used_for_reference_observations(tmp_path: Path) -> None:
    header = HEADER.rstrip() + "\tOTHER\n"
    path = _write(tmp_path, _row(1, "0/0").rstrip() + "\t0/1:2:20\n", header)
    assert _sites(_run(path, VcfParser(sample="S")))["rs990001"]["genotype"] == "A/A"
    assert _sites(_run(path, VcfParser(sample="OTHER")))["rs990001"]["reason"] == "low_gq_call"


@pytest.mark.parametrize("gt", ["0/0", "0/1", "./."])
def test_excluded_mt_is_present_but_never_typed(tmp_path: Path, gt: str) -> None:
    header = HEADER.replace("#CHROM", "##contig=<ID=chrM,length=16571>\n#CHROM")
    mt = f"chrM\t100\trs990001\tA\tG\t50\tPASS\t.\tGT\t{gt}\n"
    result = _run(_write(tmp_path, mt + _row(2, "0/1"), header))
    site = _sites(result)["rs990001"]
    assert (site["state"], site["genotype"], site["reason"]) == ("no_call", None, "excluded_mt")
    assert "rs990001" in result.genotyped_panel_rsids
    assert result.annotations == []


def test_plain_reference_only_row_and_haploid_reference(tmp_path: Path) -> None:
    rows = _row(1, "0/0").replace("\tG\t50", "\t.\t50") + _row(2, "0")
    sites = _sites(_run(_write(tmp_path, rows)))
    assert sites["rs990001"]["genotype"] == sites["rs990002"]["genotype"] == "A/A"


def test_block_does_not_invent_interior_panel_presence(tmp_path: Path) -> None:
    block = "1\t90\t.\tA\t<NON_REF>\t50\tPASS\tEND=500\tGT:GQ\t0/0:99\n"
    result = _run(_write(tmp_path, block + _row(4, "0/0")))
    assert result.panel_coverage([])["counts"]["not_on_chip"] == 4
    assert _sites(result)["rs990004"]["genotype"] == "A/A"


@pytest.mark.parametrize(
    "build,pos,expected",
    [
        ("GRCh37", 100, "typed_no_annotation"),
        ("GRCh37", 200, "not_on_chip"),
        ("GRCh38", 200, "typed_no_annotation"),
        ("GRCh38", 100, "not_on_chip"),
    ],
)
def test_identity_placements_use_only_effective_build(
    tmp_path: Path, build, pos, expected
) -> None:
    identity = _identity(tmp_path)
    path = _write(
        tmp_path,
        _row(1, "0/0").replace("\t100\t", f"\t{pos}\t").replace("rs990001", "."),
        HEADER.replace("##contig=<ID=1,length=249250621>\n", ""),
    )
    result = run_analysis(
        path, VcfParser(), [], build_override=build, identity=identity, panel_rsids=PANEL
    )
    site = _sites(result)["rs990001"]
    assert site["state"] == expected
    assert site["genotype"] == ("A/A" if expected == "typed_no_annotation" else None)
    assert result.annotations == result.hv_variants == []
    assert result.high_value_vcf_states is None


def _identity(tmp_path: Path, rows=None) -> RsidIdentityAnnotator:
    rows = rows or [("rs990001", "1", 100, "A", "G", 200, "A", "G", 0)]
    with contextlib.closing(sqlite3.connect(tmp_path / "rsid_identity.sqlite")) as conn, conn:
        conn.executescript(RSID_IDENTITY_SCHEMA)
        conn.executemany("INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)", rows)
        conn.execute(
            "INSERT INTO database_versions VALUES ('rsid_identity','test://','1',"
            "'2026-10-20',?,NULL,?)",
            (len(rows), f"sv:{RSID_IDENTITY_SCHEMA_VERSION}"),
        )
    identity = RsidIdentityAnnotator(tmp_path)
    assert identity.is_ready()
    return identity


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("build,pos", [("GRCh37", 100), ("GRCh38", 200)])
@pytest.mark.parametrize("case", ["ambiguous", "same_locus", "other_build_ambiguous"])
def test_panel_placement_ambiguity_is_order_and_build_independent(
    tmp_path: Path,
    reverse,
    build,
    pos,
    case,
) -> None:
    first = ("rs990001", "1", 100, "A", "G", 200, "A", "G", 0)
    second_pos37 = 300 if case == "ambiguous" else 100
    second_pos38 = 200 if case == "same_locus" else 400
    second = ("rs990001", "1", second_pos37, "A", "T", second_pos38, "A", "T", 0)
    rows = [second, first] if reverse else [first, second]
    identity = _identity(tmp_path, rows)
    usable = case == "same_locus" or (case == "other_build_ambiguous" and build == "GRCh37")
    with identity:
        sites = identity.forward_sites({"rs990001"}, build)
        assert bool(sites) == usable
        if usable:
            assert sites["rs990001"] == ("1", pos, frozenset({"A", "G", "T"}))
    path = _write(
        tmp_path,
        _row(1, "0/0").replace("\t100\t", f"\t{pos}\t").replace("rs990001", "."),
        HEADER.replace("##contig=<ID=1,length=249250621>\n", ""),
    )
    result = run_analysis(
        path, VcfParser(), [], build_override=build, identity=identity, panel_rsids=PANEL
    )
    # Since 2.5.2 panel attribution uses every identity placement on the
    # build (``allelix.attribution``), not ``forward_sites``, so a record at
    # either placement with its REF counts, in any row order.
    site = _sites(result)["rs990001"]
    assert site["state"] == "typed_no_annotation"
    assert site["genotype"] == "A/A"
    assert result.annotations == result.hv_variants == []


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("build,pos", [("GRCh37", 100), ("GRCh38", 200)])
def test_explicit_panel_id_survives_ambiguous_identity(
    tmp_path: Path, reverse, build, pos
) -> None:
    rows = [
        ("rs990001", "1", 100, "A", "G", 200, "A", "G", 0),
        ("rs990001", "1", 300, "A", "G", 400, "A", "G", 0),
    ]
    identity = _identity(tmp_path, list(reversed(rows)) if reverse else rows)
    path = _write(
        tmp_path,
        _row(1, "0/0").replace("\t100\t", f"\t{pos}\t"),
        HEADER.replace("##contig=<ID=1,length=249250621>\n", ""),
    )
    result = run_analysis(
        path, VcfParser(), [], build_override=build, identity=identity, panel_rsids=PANEL
    )
    site = _sites(result)["rs990001"]
    assert (site["state"], site["genotype"]) == ("typed_no_annotation", "A/A")


def test_identity_panel_placement_does_not_infer_block_interior(tmp_path: Path) -> None:
    identity = _identity(tmp_path)
    block = "1\t90\t.\tA\t<NON_REF>\t50\tPASS\tEND=500\tGT:GQ\t0/0:99\n"
    result = run_analysis(
        _write(tmp_path, block + _row(4, "0/1")),
        VcfParser(),
        [],
        build_override="GRCh37",
        identity=identity,
        panel_rsids=PANEL,
    )
    assert _sites(result)["rs990001"]["state"] == "not_on_chip"


def test_array_html_keeps_existing_column_layout(tmp_path: Path) -> None:
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.html import render_html

    result = AnalysisResult(
        Path("synthetic"), "23andme", "23andMe", "S", "GRCh37", 0, 0, [], panel_rsids=PANEL
    )
    path = tmp_path / "array.html"
    render_html(result, output_path=path)
    assert "<th>Genotype</th><th>State</th><th>Annotations</th>" in path.read_text()
    assert "<th>Reason</th>" not in path.read_text()


def test_all_supported_renderers_expose_reference_and_rejection_reasons(tmp_path: Path) -> None:
    from allelix.reports.html import render_html
    from allelix.reports.terminal import render_terminal

    result = _run(
        _write(
            tmp_path,
            _row(1, "0/0") + _row(2, "0/1", gq=3) + _row(3, "0/1", filt="q10") + _row(4, "./."),
        )
    )
    html = tmp_path / "report.html"
    render_html(result, output_path=html)
    text = io.StringIO()
    render_terminal(result, Console(file=text, width=250, color_system=None))
    for output in (html.read_text(), text.getvalue()):
        assert "A/A" in output
        for label in (
            "Low genotype quality or zero depth",
            "Caller FILTER rejected",
            "Caller no-call",
        ):
            assert label in output
    out = tmp_path / "report.json"
    render_json(result, output_path=out)
    assert json.loads(out.read_text())["panel_coverage"] == result.panel_coverage([])


def test_schema10_actual_writer_snapshot_diff_roundtrip(tmp_path: Path) -> None:
    result = _run(_write(tmp_path, _row(1, "0/0")))
    result.annotations.append(
        Annotation(
            source="clinvar",
            rsid="rs990099",
            gene="GENE1;GENE2",
            magnitude=5,
            category="clinical",
            significance="synthetic",
            description="synthetic",
            attribution="ClinVar",
            genotype_match="A/G",
        )
    )
    path = tmp_path / "report.json"
    render_json(result, output_path=path)
    payload = json.loads(path.read_text())
    assert payload["schema_version"] == "11"
    previous = load_previous_report(path)
    assert previous["panel_coverage"]["variants"][0]["genotype"] == "A/A"
    assert not compute_diff(
        result.annotations, previous["annotations"], previous["generated_at"]
    ).has_changes


@pytest.mark.parametrize("length,reason", [(None, "excluded_ucsc_mito"), (16569, None)])
def test_ucsc_mito_eligibility_matches_existing_parser_policy(
    tmp_path: Path, length, reason
) -> None:
    header = HEADER
    if length:
        header = header.replace("#CHROM", f"##contig=<ID=chrM,length={length}>\n#CHROM")
    mt = "chrM\t100\trs990001\tA\tG\t50\tPASS\t.\tGT:GQ\t0/0:99\n"
    site = _sites(_run(_write(tmp_path, mt + _row(2, "0/1"), header)))["rs990001"]
    assert site["reason"] == reason
    assert site["state"] == ("no_call" if reason else "typed_no_annotation")


@pytest.mark.parametrize("version", [str(i) for i in range(1, 12)])
def test_snapshot_legacy_and_current_compatibility(tmp_path: Path, version: str) -> None:
    path = tmp_path / "baseline.json"
    path.write_text(json.dumps({"schema_version": version, "annotations": []}))
    assert load_previous_report(path)["schema_version"] == version


def test_future_schema_rejected(tmp_path: Path) -> None:
    path = tmp_path / "future.json"
    path.write_text('{"schema_version":"12","annotations":[]}')
    with pytest.raises(ValueError, match="schema version"):
        load_previous_report(path)

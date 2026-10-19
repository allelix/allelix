# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ADR-0046: ClinVar variants without an rsID, matched for VCF input only.

The cache is built with the real loader from synthetic TSVs. Two variants
have no rsID (RS# -1) and one has rs103; each test follows a no-rsID
finding through one path that used to key on the rsID alone.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from rich.console import Console

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.databases.manager import (
    _SS_REQUIRED_COLUMNS,
    _VS_REQUIRED_COLUMNS,
    load_clinvar_tsv,
)
from allelix.models import (
    Variant,
    annotation_label,
    annotation_to_public_dict,
    public_dict_label,
)
from allelix.parsers import detect_parser
from allelix.reports._pipeline import AnalysisResult, run_analysis
from allelix.reports.diff import compute_diff
from allelix.reports.html import _group_annotations, render_html
from allelix.reports.json_report import render_json
from allelix.reports.terminal import render_terminal
from allelix.reports.vcf import _format_primary

if TYPE_CHECKING:
    from allelix.models import Annotation

BUILD = "GRCh38"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
MOCK_VCF = FIXTURES / "mock_vcf.vcf"

# (VariationID, RS#, position): two variants without an rsID, one with.
_RECORDS = (("1", "-1", 1_000), ("2", "-1", 2_000), ("3", "103", 3_000))


def _vs_line(vid: str, rs: str, pos: int) -> str:
    cols = [""] * 43
    cols[1] = "single nucleotide variant"
    cols[2] = f"NM_000001.1(GENE{vid}):c.{vid}A>G"
    cols[4] = f"GENE{vid}"
    cols[6] = "Pathogenic"
    cols[9] = rs
    cols[16] = BUILD
    cols[18] = "1"
    cols[24] = "criteria provided, single submitter"
    cols[30] = vid
    cols[31] = str(pos)
    cols[32] = "A"
    cols[33] = "G"
    return "\t".join(cols)


def _ss_line(vid: str) -> str:
    cols = [""] * 16
    cols[0] = vid
    cols[1] = "Pathogenic"
    cols[5] = f"C000000{vid}:Condition {vid}"
    cols[6] = "criteria provided, single submitter"
    cols[10] = f"SCV{int(vid):09d}.1"
    cols[15] = "yes"
    return "\t".join(cols)


def _header(columns: dict[int, str], width: int) -> str:
    return "#" + "\t".join(columns.get(i, f"col{i}") for i in range(width))


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    vs = tmp_path / "vs.tsv"
    ss = tmp_path / "ss.tsv"
    vs.write_text(
        "\n".join([_header(_VS_REQUIRED_COLUMNS, 43)] + [_vs_line(*r) for r in _RECORDS]) + "\n"
    )
    ss.write_text(
        "\n".join([_header(_SS_REQUIRED_COLUMNS, 16)] + [_ss_line(r[0]) for r in _RECORDS]) + "\n"
    )
    db = tmp_path / clinvar_db_filename(BUILD)
    load_clinvar_tsv(vs, ss, db, BUILD, source_url="test://norsid", remote_signal="md5:norsid")
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.execute(
            "UPDATE database_versions SET name = ? WHERE name = 'clinvar'",
            (clinvar_record_name(BUILD),),
        )
        conn.commit()
    return tmp_path


def _vcf_call(pos: int, rsid: str = "") -> Variant:
    """A VCF record with ID "." (rsid ""): REF read from the file."""
    return Variant(rsid, "1", pos, "A", "G", build=BUILD, ref="A", ref_from_input=True)


def _annotate(data_dir: Path, *variants: Variant) -> list[Annotation]:
    with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
        return list(clinvar.batch_annotate(list(variants)))


def _findings(data_dir: Path) -> list[Annotation]:
    return _annotate(data_dir, _vcf_call(1_000), _vcf_call(2_000))


class TestMatching:
    def test_loader_keeps_the_rows_with_an_empty_rsid(self, data_dir):
        db = data_dir / clinvar_db_filename(BUILD)
        with contextlib.closing(sqlite3.connect(db)) as conn:
            rows = conn.execute("SELECT position, rsid FROM clinvar_variants ORDER BY position")
            assert rows.fetchall() == [(1_000, ""), (2_000, ""), (3_000, "rs103")]

    def test_vcf_call_matches_a_row_without_rsid(self, data_dir):
        (a,) = _annotate(data_dir, _vcf_call(1_000))
        assert (a.significance, a.magnitude, a.rsid) == ("clinvar_pathogenic", 9.0, "")
        assert a.locus == ("1", 1_000)
        assert annotation_label(a) == "1:1000 A>G"

    def test_vcf_call_keeps_its_own_non_rs_id(self, data_dir):
        (a,) = _annotate(data_dir, _vcf_call(1_000, rsid="COSV9"))
        assert a.rsid == "COSV9"
        assert annotation_label(a) == "COSV9"

    def test_array_call_with_a_table_ref_does_not_match(self, data_dir):
        """An array call whose REF came from the identity table (not from
        the input) must not pick up a row nothing ties to its rsID."""
        array_call = Variant("rs999", "1", 1_000, "A", "G", build=BUILD, ref="A")
        assert array_call.ref_from_input is False
        assert _annotate(data_dir, array_call) == []

    def test_row_with_rsid_still_matches_an_array_call(self, data_dir):
        array_call = Variant("rs103", "1", 3_000, "A", "G", build=BUILD, ref="A")
        (a,) = _annotate(data_dir, array_call)
        assert a.rsid == "rs103"

    def test_reference_lookup_ignores_rows_without_rsid(self, data_dir):
        with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
            assert clinvar.reference_for("", BUILD) is None
            assert clinvar.reference_for("rs103", BUILD) == "A"

    def test_resolver_does_not_stamp_an_empty_rsid(self, data_dir):
        call = _vcf_call(1_000)
        with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
            assert clinvar.bulk_resolve_rsids([call]) == {}
        assert (call.rsid, call.original_rsid) == ("", None)

    def test_vcf_parser_marks_ref_from_input(self):
        calls = list(detect_parser(MOCK_VCF).parse(MOCK_VCF))
        assert calls and all(v.ref_from_input for v in calls)

    def test_array_parser_leaves_the_marker_off(self):
        path = FIXTURES / "mock_23andme.txt"
        calls = list(detect_parser(path).parse(path))
        assert calls and not any(v.ref_from_input for v in calls)


class TestReports:
    def _result(self, data_dir, tmp_path) -> AnalysisResult:
        return AnalysisResult(
            file_path=tmp_path / "input.vcf",
            parser_name="vcf",
            parser_display_name="VCF",
            sample_id="S1",
            build=BUILD,
            total_variants=2,
            skipped_count=0,
            annotators_used=[("clinvar", "test")],
            annotations=_findings(data_dir),
        )

    def test_json_names_the_locus(self, data_dir, tmp_path):
        out = tmp_path / "r.json"
        render_json(self._result(data_dir, tmp_path), output_path=out)
        rows = json.loads(out.read_text())["annotations"]
        assert sorted(
            (r["rsid"], r["chromosome"], r["position"], r["ref"], r["alt"]) for r in rows
        ) == [
            ("", "1", 1_000, "A", "G"),
            ("", "1", 2_000, "A", "G"),
        ]
        assert sorted(public_dict_label(r) for r in rows) == ["1:1000 A>G", "1:2000 A>G"]

    def test_html_keeps_two_variants_apart(self, data_dir, tmp_path):
        """Same genotype, no rsID: grouped by locus, not merged."""
        groups = _group_annotations(_findings(data_dir))
        assert sorted(annotation_label(g[0]) for g in groups) == ["1:1000 A>G", "1:2000 A>G"]
        out = tmp_path / "r.html"
        render_html(self._result(data_dir, tmp_path), output_path=out)
        page = out.read_text()
        assert "1:1000 A>G" in page and "1:2000 A>G" in page

    def test_terminal_names_the_locus(self, data_dir, tmp_path):
        console = Console(record=True, width=300)
        render_terminal(self._result(data_dir, tmp_path), console)
        text = console.export_text()
        assert "1:1000 A>G" in text and "1:2000 A>G" in text

    def test_diff_keeps_two_variants_apart(self, data_dir):
        """Two findings with the same text at different loci (a report
        often says "not provided") are two entries, not one."""
        first, second = sorted(_findings(data_dir), key=lambda a: a.locus or ("", 0))
        second = replace(second, condition=first.condition)
        diff = compute_diff([first, second], [annotation_to_public_dict(first)], "")
        assert [annotation_label(a) for a in diff.new] == ["1:2000 A>G"]
        assert diff.removed == []

    def test_vcf_info_writes_a_dot_for_no_rsid(self, data_dir):
        a = _findings(data_dir)[0]
        assert _format_primary(a) == ("clinvar", "clinvar_pathogenic|.|GENE1")

    def test_vcf_out_maps_each_record_by_locus(self, data_dir, tmp_path):
        header = [
            line
            for line in MOCK_VCF.read_text(encoding="utf-8").splitlines(keepends=True)
            if line.startswith("#")
        ]
        vcf = tmp_path / "in.vcf"
        vcf.write_text(
            "".join(header)
            + "1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "1\t2000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        )
        with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
            result = run_analysis(
                vcf, detect_parser(vcf), [clinvar], collect_vcf_write_inputs=True
            )
        assert result.vcf_write_inputs is not None
        by_key = result.vcf_write_inputs.annotations_by_key
        assert sorted(by_key) == [("1", 1_000, "A", "G"), ("1", 2_000, "A", "G")]
        assert [a.gene for a in by_key[("1", 2_000, "A", "G")]] == ["GENE2"]

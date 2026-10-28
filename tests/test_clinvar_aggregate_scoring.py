# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ADR-0045: ClinVar rows scored against ClinVar's own view of the variant.

Caches are built with the real loader from synthetic TSVs, so each test
runs the loader, the annotator and the reports on the same rows.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from rich.console import Console

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.databases.manager import (
    _SS_REQUIRED_COLUMNS,
    _VS_REQUIRED_COLUMNS,
    load_clinvar_tsv,
)
from allelix.models import Variant
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.clinvar_text import (
    detail_fields,
    scored_down_count,
    scored_down_line,
    submitter_split,
    terminal_note,
)
from allelix.reports.html import render_html
from allelix.reports.json_report import SCHEMA_VERSION, render_json
from allelix.reports.terminal import render_terminal

if TYPE_CHECKING:
    from pathlib import Path

BUILD = "GRCh38"
CONFLICTING = "Conflicting classifications of pathogenicity"

# (VariationID, rsID, position, aggregate, [(significance, review, contributes)])
_CR = "criteria provided, single submitter"
_EP = "reviewed by expert panel"
_RECORDS = (
    # 226007-shaped: ClinVar calls it Conflicting, 2 of 5 say Pathogenic.
    (
        "1",
        "rs101",
        1_000,
        CONFLICTING,
        [("Benign", _CR, "yes")] * 3 + [("Pathogenic", _CR, "yes")] * 2,
    ),
    # Every submitter agrees: Pathogenic, not capped.
    ("2", "rs102", 2_000, "Pathogenic", [("Pathogenic", _CR, "yes")] * 3),
    # A Benign aggregate holding a P submission: capped at Benign's 1.0.
    ("3", "rs103", 3_000, "Benign", [("Pathogenic", _CR, "yes"), ("Benign", _EP, "yes")]),
    # F5-shaped: an expert panel's drug response outranks agreeing P calls.
    (
        "4",
        "rs104",
        4_000,
        "drug response",
        [("drug response", _EP, "yes"), ("Pathogenic", _CR, "no"), ("Pathogenic", _CR, "no")],
    ),
    # MTHFR-shaped: outranked P calls conflict with outranked Benign calls.
    (
        "5",
        "rs105",
        5_000,
        "drug response",
        [("drug response", _EP, "yes"), ("Pathogenic", _CR, "no"), ("Benign", _CR, "no")],
    ),
    # Vocabulary: a low-penetrance call.
    (
        "6",
        "rs106",
        6_000,
        "Pathogenic, low penetrance",
        [("Pathogenic, low penetrance", _CR, "yes")],
    ),
    # SOD2 rs4880-shaped: "Conflicting ...; risk factor". The aggregate's
    # germline term caps the disease calls, not the risk factor call.
    (
        "8",
        "rs108",
        8_000,
        f"{CONFLICTING}; risk factor",
        [("Pathogenic", _CR, "yes"), ("Benign", _CR, "yes"), ("risk factor", _CR, "yes")],
    ),
    # A classification Allelix does not score.
    ("7", "rs107", 7_000, "Pathogenic", [("Newly coined classification", _CR, "yes")]),
)


def _vs_line(vid: str, rsid: str, pos: int, aggregate: str) -> str:
    cols = [""] * 43
    cols[1] = "single nucleotide variant"
    cols[2] = f"NM_000001.1(GENE{vid}):c.{vid}A>G"
    cols[4] = f"GENE{vid}"
    cols[6] = aggregate
    cols[9] = rsid.removeprefix("rs")
    cols[16] = BUILD
    cols[18] = "1"
    cols[24] = _CR
    cols[30] = vid
    cols[31] = str(pos)
    cols[32] = "A"
    cols[33] = "G"
    return "\t".join(cols)


def _ss_line(vid: str, sig: str, review: str, contributes: str, scv: int) -> str:
    cols = [""] * 16
    cols[0] = vid
    cols[1] = sig
    cols[5] = "C0000001:Some condition"
    cols[6] = review
    cols[10] = f"SCV{scv:09d}.1"
    cols[15] = contributes
    return "\t".join(cols)


def _header(columns: dict[int, str], width: int) -> str:
    return "#" + "\t".join(columns.get(i, f"col{i}") for i in range(width))


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    vs = tmp_path / "vs.tsv"
    ss = tmp_path / "ss.tsv"
    vs.write_text(
        "\n".join(
            [_header(_VS_REQUIRED_COLUMNS, 43)]
            + [_vs_line(vid, rs, pos, agg) for vid, rs, pos, agg, _subs in _RECORDS]
        )
        + "\n"
    )
    scv = iter(range(1, 1000))
    ss.write_text(
        "\n".join(
            [_header(_SS_REQUIRED_COLUMNS, 16)]
            + [
                _ss_line(vid, sig, review, contributes, next(scv))
                for vid, _rs, _pos, _agg, subs in sorted(_RECORDS, key=lambda r: int(r[0]))
                for sig, review, contributes in subs
            ]
        )
        + "\n"
    )
    db = tmp_path / clinvar_db_filename(BUILD)
    load_clinvar_tsv(vs, ss, db, BUILD, source_url="test://agg", remote_signal="md5:agg")
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.execute(
            "UPDATE database_versions SET name = ? WHERE name = 'clinvar'",
            (clinvar_record_name(BUILD),),
        )
        conn.commit()
    return tmp_path


def _het(rsid: str, pos: int) -> Variant:
    return Variant(rsid, "1", pos, "A", "G", ref="A", build=BUILD)


def _annotate(data_dir: Path, *, include_benign: bool = False):
    with ClinVarAnnotator(data_dir, builds=(BUILD,), include_benign=include_benign) as clinvar:
        out = {rs: clinvar.annotate(_het(rs, pos)) for _vid, rs, pos, _agg, _subs in _RECORDS}
        return out, clinvar.run_counts()["unrecognised_rows"]


class TestScoring:
    def test_conflicting_aggregate_caps_pathogenic_row(self, data_dir):
        rows, _ = _annotate(data_dir)
        (row,) = rows["rs101"]
        assert row.significance == "clinvar_pathogenic"
        assert row.magnitude == 4.0
        assert row.aggregate_capped is True
        assert row.outside_aggregate is False
        assert row.clinvar_aggregate == CONFLICTING
        assert (row.submitters_same, row.submitters_total) == (2, 5)

    def test_benign_rows_carry_their_split(self, data_dir):
        rows, _ = _annotate(data_dir, include_benign=True)
        benign = [a for a in rows["rs101"] if a.significance == "clinvar_benign"]
        assert [(a.submitters_same, a.submitters_total, a.magnitude) for a in benign] == [
            (3, 5, 1.0)
        ]
        assert benign[0].aggregate_capped is False

    def test_agreeing_aggregate_is_not_capped(self, data_dir):
        rows, _ = _annotate(data_dir)
        (row,) = rows["rs102"]
        assert (row.magnitude, row.aggregate_capped) == (9.0, False)
        assert (row.submitters_same, row.submitters_total) == (3, 3)

    def test_benign_aggregate_caps_at_benign(self, data_dir):
        rows, _ = _annotate(data_dir)
        (row,) = rows["rs103"]
        assert (row.significance, row.magnitude, row.aggregate_capped) == (
            "clinvar_pathogenic",
            1.0,
            True,
        )

    def test_outranked_agreeing_submissions_are_scored_as_submitted(self, data_dir):
        rows, _ = _annotate(data_dir)
        by_sig = {a.significance: a for a in rows["rs104"]}
        assert set(by_sig) == {"clinvar_drug_response", "clinvar_pathogenic"}
        p = by_sig["clinvar_pathogenic"]
        assert (p.magnitude, p.aggregate_capped, p.outside_aggregate) == (9.0, False, True)
        assert (p.submitters_same, p.submitters_total) == (2, 2)
        assert p.clinvar_aggregate == "drug response"
        assert by_sig["clinvar_drug_response"].outside_aggregate is False

    def test_outranked_conflicting_submissions_are_capped(self, data_dir):
        rows, _ = _annotate(data_dir)
        p = next(a for a in rows["rs105"] if a.significance == "clinvar_pathogenic")
        assert (p.magnitude, p.aggregate_capped, p.outside_aggregate) == (4.0, True, True)
        assert (p.submitters_same, p.submitters_total) == (1, 2)

    def test_compound_aggregate_caps_only_disease_calls(self, data_dir):
        rows, _ = _annotate(data_dir)
        by_sig = {a.significance: a for a in rows["rs108"]}
        assert (
            by_sig["clinvar_pathogenic"].magnitude,
            by_sig["clinvar_pathogenic"].aggregate_capped,
        ) == (4.0, True)
        risk = by_sig["clinvar_risk_factor"]
        assert (risk.magnitude, risk.aggregate_capped) == (6.0, False)
        assert (risk.submitters_same, risk.submitters_total) == (1, 3)

    def test_low_penetrance_vocabulary(self, data_dir):
        rows, _ = _annotate(data_dir)
        (row,) = rows["rs106"]
        assert (row.significance, row.magnitude) == ("clinvar_pathogenic,_low_penetrance", 7.5)

    def test_unscored_classification_is_hidden_and_counted(self, data_dir):
        rows, unrecognised = _annotate(data_dir)
        assert rows["rs107"] == []
        assert unrecognised == 1


class TestDisplay:
    def _result(self, data_dir, tmp_path) -> AnalysisResult:
        rows, unrecognised = _annotate(data_dir)
        annotations = [a for anns in rows.values() for a in anns]
        return AnalysisResult(
            file_path=tmp_path / "input.vcf",
            parser_name="vcf",
            parser_display_name="VCF",
            sample_id="S1",
            build=BUILD,
            total_variants=len(_RECORDS),
            skipped_count=0,
            annotators_used=[("clinvar", "test")],
            annotations=annotations,
            clinvar_scored_down=scored_down_count(annotations),
            clinvar_unrecognised=unrecognised,
        )

    def test_scored_down_count_counts_capped_pathogenic_rows(self, data_dir, tmp_path):
        # rs101 (Conflicting), rs103 (Benign aggregate), rs105 (outranked
        # conflict), rs108 (Conflicting; risk factor).
        assert self._result(data_dir, tmp_path).clinvar_scored_down == 4
        assert scored_down_line(0) == ""
        assert scored_down_line(3).startswith("3 ClinVar Pathogenic or Likely pathogenic rows")

    def test_text(self, data_dir):
        rows, _ = _annotate(data_dir)
        capped = rows["rs101"][0]
        assert submitter_split(capped) == "2 of 5 submitters: Pathogenic"
        assert terminal_note(capped) == " (2 of 5; ClinVar: conflicting)"
        outranked = next(a for a in rows["rs104"] if a.outside_aggregate)
        assert terminal_note(outranked) == " (2 of 2; outranked; ClinVar: drug response)"
        assert terminal_note(rows["rs102"][0]) == ""
        labels = [label for label, _value in detail_fields(capped)]
        assert labels == ["ClinVar classification", "Submitters", "Scored"]

    def test_json_carries_fields_and_counts(self, data_dir, tmp_path):
        out = tmp_path / "r.json"
        render_json(self._result(data_dir, tmp_path), output_path=out)
        payload = json.loads(out.read_text())
        assert payload["schema_version"] == SCHEMA_VERSION == "11"
        assert payload["clinvar"] == {
            "pathogenic_rows_scored_down": 4,
            "unrecognised_rows_hidden": 1,
            "common_allele_rows_scored_down": 0,
            # The result is built with the count alone (no named variants).
            "scored_down_variants": [],
        }
        row = next(a for a in payload["annotations"] if a["rsid"] == "rs101")
        assert row["clinvar_aggregate"] == CONFLICTING
        assert (row["submitters_same"], row["submitters_total"]) == (2, 5)
        assert row["aggregate_capped"] is True
        assert row["outside_aggregate"] is False

    def test_html_shows_banner_and_detail(self, data_dir, tmp_path):
        out = tmp_path / "r.html"
        render_html(self._result(data_dir, tmp_path), output_path=out)
        page = out.read_text()
        assert "<strong>ClinVar.</strong>" in page
        assert "4 ClinVar Pathogenic or Likely pathogenic rows were scored down" in page
        assert "1 ClinVar row(s) not shown" in page
        assert "2 of 5 submitters: Pathogenic" in page

    def test_terminal_shows_note(self, data_dir, tmp_path):
        console = Console(record=True, width=300)
        render_terminal(self._result(data_dir, tmp_path), console)
        text = console.export_text()
        assert "pathogenic (2 of 5; ClinVar: conflicting)" in text
        assert "pathogenic (2 of 2; outranked; ClinVar: drug response)" in text

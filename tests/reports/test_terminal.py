# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for terminal report rendering."""

from __future__ import annotations

import io
from pathlib import Path

from rich.console import Console

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.terminal import render_terminal, render_terminal_diff


def _ann(**overrides) -> Annotation:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 5.0,
        "description": "test",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "GENE1",
        "condition": "Some condition",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("dummy.txt"),
        parser_name="x",
        parser_display_name="X",
        sample_id="S",
        build="GRCh37",
        total_variants=0,
        skipped_count=0,
        annotators_used=[],
        annotations=annotations,
    )


def _render(annotations, **kwargs) -> tuple[str, int]:
    buf = io.StringIO()
    console = Console(file=buf, force_terminal=False, width=200)
    count = render_terminal(_result(annotations), console=console, **kwargs)
    return buf.getvalue(), count


class TestRenderTerminal:
    def test_empty_list_renders_message(self):
        out, count = _render([])
        assert count == 0
        assert "No annotations" in out

    def test_attribution_column_present(self):
        out, _ = _render([_ann()])
        assert "ClinVar" in out
        assert "Source" in out

    def test_sorts_by_magnitude_descending(self):
        annotations = [
            _ann(rsid="rs_low", magnitude=2.0),
            _ann(rsid="rs_high", magnitude=9.0),
            _ann(rsid="rs_mid", magnitude=5.0),
        ]
        out, _ = _render(annotations)
        assert out.index("rs_high") < out.index("rs_mid") < out.index("rs_low")

    def test_min_magnitude_filter(self):
        annotations = [
            _ann(rsid="rs_skip", magnitude=2.0),
            _ann(rsid="rs_keep", magnitude=8.0),
        ]
        out, count = _render(annotations, min_magnitude=5.0)
        assert count == 1
        assert "rs_keep" in out
        assert "rs_skip" not in out

    def test_category_filter(self):
        annotations = [
            _ann(rsid="rs_clinical", category="clinical"),
            _ann(rsid="rs_pharma", category="pharma"),
        ]
        out, _ = _render(annotations, category="clinical")
        assert "rs_clinical" in out
        assert "rs_pharma" not in out

    def test_genes_filter(self):
        annotations = [
            _ann(rsid="m", gene="MTHFR"),
            _ann(rsid="b", gene="BRCA1"),
        ]
        out, count = _render(annotations, genes={"MTHFR"})
        assert count == 1
        assert "m" in out
        assert "BRCA1" not in out

    def test_review_status_column_present(self):
        out, _ = _render([_ann(review_status="criteria_provided,_single_submitter")])
        assert "Review Status" in out
        assert "criteria_provided" in out

    def test_review_status_dash_when_empty(self):
        out, _ = _render([_ann(review_status="")])
        assert "Review Status" in out

    def test_freq_column_present_when_frequency_set(self):
        out, _ = _render([_ann(allele_frequency=0.35)])
        assert "Freq" in out
        assert "35.00%" in out

    def test_freq_column_absent_when_no_frequency(self):
        out, _ = _render([_ann()])
        assert "Freq" not in out

    def test_freq_rare_variant_format(self):
        out, _ = _render([_ann(allele_frequency=0.00005)])
        assert "<0.01%" in out

    def test_freq_none_shows_dash(self):
        annotations = [
            _ann(rsid="rs1", allele_frequency=0.35),
            _ann(rsid="rs2", allele_frequency=None),
        ]
        out, _ = _render(annotations)
        assert "Freq" in out
        assert "35.00%" in out


class TestZygosityColumn:
    def test_zygosity_column_present(self):
        out, _ = _render([_ann(genotype_match="A/G")])
        assert "Zygosity" in out

    def test_heterozygous_label(self):
        out, _ = _render([_ann(genotype_match="A/G")])
        assert "Heterozygous" in out

    def test_homozygous_label(self):
        out, _ = _render([_ann(genotype_match="A/A")])
        assert "Homozygous" in out

    def test_no_call_label(self):
        out, _ = _render([_ann(genotype_match="A/-")])
        assert "No Call" in out


class TestAlphaMissenseColumn:
    def test_am_column_present_when_data_exists(self):
        out, _ = _render([_ann(am_pathogenicity=0.95, am_class="likely_pathogenic")])
        assert "AM" in out
        assert "0.950" in out

    def test_am_column_absent_when_no_data(self):
        out, _ = _render([_ann()])
        assert "AM" not in out

    def test_am_none_shows_dash(self):
        annotations = [
            _ann(rsid="rs1", am_pathogenicity=0.80, am_class="likely_pathogenic"),
            _ann(rsid="rs2"),
        ]
        out, _ = _render(annotations)
        assert "AM" in out
        assert "0.800" in out

    def test_am_pharmgkb_footnote(self):
        annotations = [
            _ann(
                source="pharmgkb",
                attribution="ClinPGx",
                am_pathogenicity=0.95,
                am_class="likely_pathogenic",
            ),
        ]
        out, _ = _render(annotations)
        assert "protein structure impact only" in out

    def test_am_no_footnote_without_pharmgkb(self):
        out, _ = _render([_ann(am_pathogenicity=0.95, am_class="likely_pathogenic")])
        assert "protein structure impact only" not in out


def _ann_dict(**overrides) -> dict:
    defaults = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "clinvar: test",
        "attribution": "ClinVar",
        "genotype_match": "AG",
        "references": [],
        "condition": "MTHFR deficiency",
        "gene": "MTHFR",
    }
    defaults.update(overrides)
    return defaults


class TestRenderTerminalDiff:
    def test_render_terminal_diff_new_only(self, capsys):
        """Verifies the New Annotations table renders with all columns."""
        from allelix.reports.diff import DiffResult

        diff = DiffResult(
            new=[
                _ann(
                    rsid="rs1801133",
                    gene="MTHFR",
                    magnitude=9.0,
                    review_status="criteria_provided,_single_submitter",
                )
            ],
            previous_generated_at="2026-05-01T00:00:00",
        )
        total = render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert total == 1
        assert "New Annotations (1)" in out
        assert "rs1801133" in out
        assert "Review Status" in out
        assert "criteria_provided" in out

    def test_render_terminal_diff_changed_only(self, capsys):
        """Old Sig / New Sig / Old Mag / New Mag / Review Status columns all present."""
        from allelix.reports.diff import ChangedAnnotation, DiffResult

        diff = DiffResult(
            changed=[
                ChangedAnnotation(
                    current=_ann(magnitude=7.0, review_status="reviewed_by_expert_panel"),
                    previous_significance="old_sig",
                    previous_magnitude=9.0,
                )
            ],
            previous_generated_at="2026-05-01T00:00:00",
        )
        render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert "Changed Annotations (1)" in out
        assert "old_sig" in out
        assert "9.0" in out and "7.0" in out
        assert "Review Status" in out
        assert "reviewed_by_expert_panel" in out

    def test_render_terminal_diff_removed_only(self, capsys):
        from allelix.reports.diff import DiffResult

        diff = DiffResult(
            removed=[_ann_dict()],
            previous_generated_at="2026-05-01T00:00:00",
        )
        render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert "Removed Annotations (1)" in out
        assert "rs1801133" in out

    def test_render_terminal_diff_no_changes(self, capsys):
        from allelix.reports.diff import DiffResult

        diff = DiffResult(previous_generated_at="2026-05-01T00:00:00")
        total = render_terminal_diff(diff, Console(force_terminal=True, width=200))
        assert total == 0
        assert "No changes since previous report." in capsys.readouterr().out

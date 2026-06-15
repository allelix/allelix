# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for terminal report rendering.

GH #9: the terminal table is intentionally bare-min — rsID, Gene
(conditional), Source, Significance, Mag, GT, Condition (conditional).
Enrichment columns (Review Status, Zygosity, Freq, AM, CADD) belong
in the HTML and JSON reports, not in a 12-column-wide table that Rich
squashes to hairline-zero-width on typical terminals. These tests
pin the new layout — the data is still on the Annotation model and
still surfaces in the other reports; just not in the terminal.
"""

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


class TestBareMinColumns:
    """GH #9: enrichment columns are intentionally NOT in the terminal."""

    def test_review_status_not_in_terminal(self):
        """ClinVar review status surfaces in HTML/JSON, not terminal."""
        out, _ = _render([_ann(review_status="criteria_provided,_single_submitter")])
        assert "Review Status" not in out
        assert "criteria_provided" not in out

    def test_zygosity_not_in_terminal(self):
        """Zygosity is derivable from the GT column; not duplicated here."""
        out, _ = _render([_ann(genotype_match="A/G")])
        assert "Zygosity" not in out
        assert "Heterozygous" not in out
        assert "Homozygous" not in out

    def test_freq_not_in_terminal_even_when_set(self):
        out, _ = _render([_ann(allele_frequency=0.35)])
        assert "Freq" not in out
        assert "35.00%" not in out

    def test_am_not_in_terminal_even_when_set(self):
        out, _ = _render([_ann(am_pathogenicity=0.95, am_class="likely_pathogenic")])
        assert "0.950" not in out
        assert "protein structure impact only" not in out

    def test_cadd_not_in_terminal_even_when_set(self):
        out, _ = _render([_ann(cadd_phred=24.3)])
        assert "CADD" not in out
        assert "24.3" not in out


class TestCompactDisplay:
    """GH #9: terminal cells use compact strings to fit narrow terminals."""

    def test_significance_strips_source_prefix(self):
        """Source column already shows ClinVar; significance shouldn't repeat it."""
        out, _ = _render([_ann(significance="clinvar_pathogenic", source="clinvar")])
        assert "pathogenic" in out
        # The full prefixed form should NOT appear:
        assert "clinvar_pathogenic" not in out

    def test_gwas_source_compacted(self):
        """``GWAS Catalog`` always truncates to ``GWAS Ca…`` in narrow
        terminals; render as plain ``GWAS`` instead."""
        out, _ = _render(
            [
                _ann(
                    source="gwas",
                    attribution="GWAS Catalog",
                    significance="gwas_association",
                )
            ]
        )
        assert "GWAS" in out
        assert "GWAS Catalog" not in out

    def test_gt_column_header(self):
        """Genotype column is labeled ``GT`` (matching the web report)."""
        out, _ = _render([_ann(genotype_match="A/G")])
        assert "GT" in out
        # The cell still carries the full genotype string:
        assert "A/G" in out


class TestConditionalColumns:
    """Gene and Condition are dropped when no row has them."""

    def test_gene_column_dropped_when_empty(self):
        out, _ = _render([_ann(gene="")])
        assert "Gene" not in out

    def test_gene_column_present_when_any_row_has_it(self):
        out, _ = _render(
            [
                _ann(rsid="rs1", gene="MTHFR"),
                _ann(rsid="rs2", gene=""),
            ]
        )
        assert "Gene" in out
        assert "MTHFR" in out

    def test_condition_column_dropped_when_empty(self):
        out, _ = _render([_ann(condition="")])
        assert "Condition" not in out

    def test_condition_column_present_when_any_row_has_it(self):
        out, _ = _render(
            [
                _ann(rsid="rs1", condition="MTHFR deficiency"),
                _ann(rsid="rs2", condition=""),
            ]
        )
        assert "Condition" in out
        assert "MTHFR deficiency" in out


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
        """New Annotations table renders with the bare-min column set."""
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
        # Review Status was provided but does NOT surface in terminal (bare-min):
        assert "Review Status" not in out
        assert "criteria_provided" not in out

    def test_render_terminal_diff_changed_only(self, capsys):
        """Changed table keeps Old Sig / New Sig / Old Mag / New Mag (those
        are the whole point of the diff). Review Status is still axed."""
        from allelix.reports.diff import ChangedAnnotation, DiffResult

        diff = DiffResult(
            changed=[
                ChangedAnnotation(
                    current=_ann(magnitude=7.0, review_status="reviewed_by_expert_panel"),
                    previous_significance="clinvar_old_sig",
                    previous_magnitude=9.0,
                )
            ],
            previous_generated_at="2026-05-01T00:00:00",
        )
        render_terminal_diff(diff, Console(force_terminal=True, width=200))
        out = capsys.readouterr().out
        assert "Changed Annotations (1)" in out
        # Significance prefixes are stripped (source = clinvar):
        assert "old_sig" in out
        assert "pathogenic" in out
        assert "9.0" in out and "7.0" in out
        assert "Review Status" not in out
        assert "reviewed_by_expert_panel" not in out

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

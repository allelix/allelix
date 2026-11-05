# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A scored-down finding retains the explanation beside its new magnitude."""

from __future__ import annotations

import io

import pytest
from rich.console import Console

from allelix.models import Annotation
from allelix.reports.diff import ChangedAnnotation, DiffResult
from allelix.reports.terminal import render_terminal_diff
from tests.test_narrow_changed_table import _cells


@pytest.mark.parametrize("width", [60, 70, 79, 80, 120])
@pytest.mark.parametrize("source", ["ClinVar", "ClinPGx", "SNPedia", "GWAS Catalog"])
@pytest.mark.parametrize("gene", ["", "LOC105374351"])
def test_below_floor_explanation_keeps_every_word(width, source, gene):
    row = Annotation(
        source="clinvar",
        rsid="rs2200000000",
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=4.0,
        description="Synthetic finding",
        attribution=source,
        genotype_match="A/C",
        gene=gene,
        condition="Mock condition",
    )
    diff = DiffResult(
        changed=[
            ChangedAnnotation(
                current=row,
                previous_significance="clinvar_pathogenic",
                previous_magnitude=9.0,
                previous_genotype_match="A/G",
                below_display_floor=True,
            )
        ]
    )
    stream = io.StringIO()
    render_terminal_diff(
        diff, Console(file=stream, width=width, color_system=None, force_terminal=False)
    )
    cells = _cells(stream.getvalue())
    for value in ("9.0", "4.0", "belowdisplayfloor", "Mockcondition", "A/G->A/C"):
        assert any(value in cell for cell in cells), (value, stream.getvalue())

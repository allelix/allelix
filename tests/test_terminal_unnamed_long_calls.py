# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An allele in a positional label must not hide the rest of a terminal row."""

from __future__ import annotations

import io
import re

import pytest
from rich.console import Console

from allelix.models import Annotation, annotation_label, annotation_to_public_dict
from allelix.reports.diff import ChangedAnnotation, DiffResult
from allelix.reports.terminal import _print_table, render_terminal_diff


def _cells(output: str) -> list[str]:
    rows = [
        re.split("[│┃]", line)[1:-1] for line in output.splitlines() if "│" in line or "┃" in line
    ]
    assert rows, output
    return ["".join("".join(row[i].split()) for row in rows) for i in range(len(rows[0]))]


@pytest.mark.parametrize("width", [80, 250])
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("haploid", [False, True])
@pytest.mark.parametrize("view", ["analyze", "new", "changed", "removed"])
def test_positional_allele_label_keeps_other_columns(width, named, padded, haploid, view):
    ref, alt = "A", "A" + "CGTA" * 75
    pos = 20_010_383
    if padded:
        ref, alt, pos = "T" + ref + "C", "T" + alt + "C", pos - 1
    genotype = alt if haploid else f"{ref}/{alt}"
    annotation = Annotation(
        source="clinvar",
        rsid="rs900260111" if named else "",
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=9.0,
        description="Synthetic finding",
        attribution="ClinVar",
        genotype_match=genotype,
        gene="MOCK",
        condition="Mock trait",
        locus=("22", pos),
        locus_ref=ref,
        alt=alt,
    )
    stream = io.StringIO()
    console = Console(file=stream, width=width, color_system=None, force_terminal=False)
    if view == "analyze":
        _print_table([annotation], console)
    else:
        diff = DiffResult()
        if view == "new":
            diff.new = [annotation]
        elif view == "changed":
            diff.changed = [
                ChangedAnnotation(
                    current=annotation,
                    previous_significance="clinvar_benign",
                    previous_magnitude=1.0,
                    previous_genotype_match="A/G",
                )
            ]
        else:
            diff.removed = [annotation_to_public_dict(annotation)]
        assert render_terminal_diff(diff, console) == 1
    output = stream.getvalue()
    cells = _cells(output)
    # Check each cell separately: a sequence in the label is not a visible GT.
    assert "MOCK" in cells[1], output
    assert "ClinVar" in cells[2], output
    assert "pathogenic" in cells[4 if view == "changed" else 3], output
    assert "9.0" in cells[6 if view == "changed" else 4], output
    assert genotype in cells[7 if view == "changed" else 5], output
    assert "Mocktrait" in cells[-1], output
    assert "".join(annotation_label(annotation).split()) in cells[0], output
    if view == "changed":
        assert "benign" in cells[3] and "1.0" in cells[5], output
        assert "A/G->" in cells[7], output

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Positional labels must fit with the other cells near a wrapping boundary."""

from __future__ import annotations

import io
import re

import pytest
from rich.console import Console

from allelix.models import Annotation, annotation_label, annotation_to_public_dict
from allelix.reports.diff import ChangedAnnotation, DiffResult
from allelix.reports.terminal import _print_table, render_terminal_diff


def _cells(output: str) -> list[str]:
    rows = []
    for line in output.splitlines():
        if not line.startswith(("│", "┃")):
            continue
        row = re.split("[│┃]", line)[1:]
        if row[-1] == "":
            row.pop()
        rows.append(row)
    assert rows, output
    # A clipped right border alone is not a missing answer. Keep its cell text.
    return ["".join("".join(row[i].split()) for row in rows) for i in range(min(map(len, rows)))]


@pytest.mark.parametrize("width", [80, 120, 250])
@pytest.mark.parametrize("alt_length", [9, 10, 11])
@pytest.mark.parametrize("gene", ["MOCKGENE", "GENESYMBOL"])
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("view", ["analyze", "new", "changed", "removed"])
def test_medium_positional_labels_keep_the_condition(width, alt_length, gene, named, view):
    alt = "C" * alt_length
    annotation = Annotation(
        source="clinvar",
        rsid="rs900260121" if named else "",
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=9.0,
        description="Synthetic finding",
        attribution="ClinVar",
        genotype_match=f"A/{alt}",
        gene=gene,
        condition="Mock trait",
        locus=("22", 20_010_383),
        locus_ref="A",
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
    assert len(cells) == (9 if view == "changed" else 7), output
    assert "Mocktrait" in cells[-1], output
    assert gene in cells[1] and "ClinVar" in cells[2], output
    assert "pathogenic" in cells[4 if view == "changed" else 3], output
    assert "9.0" in cells[6 if view == "changed" else 4], output
    assert annotation.genotype_match in cells[7 if view == "changed" else 5], output
    assert "".join(annotation_label(annotation).split()) in cells[0], output
    if view == "changed":
        assert "benign" in cells[3] and "1.0" in cells[5], output
        assert "A/G->" in cells[7], output

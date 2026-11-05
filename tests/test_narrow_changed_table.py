# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A narrow terminal must retain every answer in a changed finding."""

from __future__ import annotations

import io
import re

import pytest
from rich.console import Console

from allelix.models import Annotation, annotation_label, annotation_to_public_dict
from allelix.reports.diff import ChangedAnnotation, DiffResult
from allelix.reports.terminal import _print_table, render_terminal_diff

GENES = ["MOCKGENEABC"[:length] for length in range(4, 12)] + ["LOC105374351"]


def _cells(output: str) -> list[str]:
    """Join wrapping within columns, retaining text if the right border is clipped."""
    rows = []
    for line in output.splitlines():
        if not line.startswith(("│", "┃")):
            continue
        cells = re.split("[│┃]", line)[1:]
        if cells and cells[-1] == "":
            cells.pop()
        rows.append(cells)
    assert rows, output
    # No column count, width, order or wrapping strategy is prescribed.
    return [
        "".join("".join(row[index].split()) for row in rows if index < len(row))
        for index in range(max(map(len, rows)))
    ]


def _render(width: int, named: bool, long_call: bool, gene: str, view: str) -> None:
    ref, alt = ("ACGT" * 20, "A") if long_call else ("A", "C")
    genotype = f"{ref}/{ref}" if long_call else "A/C"
    finding = Annotation(
        source="clinvar",
        rsid="rs2200000000" if named else "",
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=9.0,
        description="Synthetic finding",
        attribution="ClinVar",
        genotype_match=genotype,
        gene=gene,
        condition="Mock trait",
        locus=("22", 20_010_383),
        locus_ref=ref,
        alt=alt,
    )
    stream = io.StringIO()
    console = Console(file=stream, width=width, color_system=None, force_terminal=False)
    expected = [annotation_label(finding), gene, "ClinVar", "pathogenic", "9.0", genotype]
    expected.append("Mock trait")
    if view == "analyze":
        _print_table([finding], console)
    else:
        diff = DiffResult()
        if view == "new":
            diff.new = [finding]
        elif view == "removed":
            diff.removed = [annotation_to_public_dict(finding)]
        else:
            diff.changed = [
                ChangedAnnotation(
                    current=finding,
                    previous_significance="clinvar_benign",
                    previous_magnitude=1.0,
                    previous_genotype_match="G/T",
                )
            ]
            expected.extend(["benign", "1.0", "G/T"])
        assert render_terminal_diff(diff, console) == 1
    output = stream.getvalue()
    cells = _cells(output)
    missing = [value for value in expected if not any("".join(value.split()) in c for c in cells)]
    assert not missing, f"Missing cell values: {missing}\n{output}"


@pytest.mark.parametrize("width", [60, 70, 80])
@pytest.mark.parametrize("named", [True, False])
@pytest.mark.parametrize("long_call", [False, True])
@pytest.mark.parametrize("gene", GENES)
def test_changed_row_preserves_every_cell(
    width: int, named: bool, long_call: bool, gene: str
) -> None:
    """Changed fields survive below 80 columns and at the control width."""
    _render(width, named, long_call, gene, "changed")


@pytest.mark.parametrize("width", [60, 70, 80])
@pytest.mark.parametrize("named", [True, False])
@pytest.mark.parametrize("long_call", [False, True])
@pytest.mark.parametrize("gene", ["MOCK", "LOC105374351"])
@pytest.mark.parametrize("view", ["analyze", "new", "removed"])
def test_other_report_tables_preserve_every_cell(
    width: int, named: bool, long_call: bool, gene: str, view: str
) -> None:
    """The three other finding tables retain the same field values."""
    _render(width, named, long_call, gene, view)

"""Adding a site to an rsID must not hide other finding cells.

Synthetic contig names are legal input names; rsID-based sources can still
report them. Join wrapped lines within each column, including a last cell
whose right border alone was clipped.
"""

import io
import re
from dataclasses import replace

import pytest
from rich.console import Console
from test_html_site_identity import ann, result

from allelix.reports.terminal import render_terminal


def _cells(text):
    rows = []
    for line in text.splitlines():
        if not line.startswith(("│", "┃")):
            continue
        cells = re.split("[│┃]", line)[1:]
        if cells[-1] == "":
            cells.pop()
        rows.append(cells)
    assert rows, text
    return [
        "".join("".join(row[i].split()) for row in rows if i < len(row))
        for i in range(max(map(len, rows)))
    ]


def _render(annotations, width):
    buffer = io.StringIO()
    console = Console(file=buffer, width=width, color_system=None, force_terminal=False)
    render_terminal(result(annotations), console)
    text = buffer.getvalue()
    return _cells(text), text


def _assert_values(cells, text):
    # All are finding-cell text, not a disclaimer or a banner outside the table.
    for value in ("rs900529001", "LOC105374351", "ClinVar", "pathogenic", "8.0", "AG"):
        assert any(value in cell for cell in cells), f"Missing {value!r}\n{text}"
    for value in ("Firstcondition", "Secondcondition"):
        assert any(value in cell for cell in cells), f"Missing {value!r}\n{text}"


@pytest.mark.parametrize("width", [80, 120, 160])
@pytest.mark.parametrize(
    "chromosome", ["synthetic_contig_" + "x" * 20, "synthetic_contig_" + "x" * 64]
)
@pytest.mark.parametrize("origin", ["locus", "record_key"])
def test_distinct_sites_keep_all_cells_and_coordinates(width, chromosome, origin):
    first = ann(locus=(chromosome, 12345), gene="LOC105374351")
    second = replace(
        first, locus=(chromosome, 12346), record_index=1, condition="Second condition"
    )
    if origin == "record_key":
        # A trimmed source locus must display the input record's position.
        first = replace(first, record_key=(chromosome, 12344, "GA", "G"))
        second = replace(second, record_key=(chromosome, 12345, "GA", "G"))
    cells, text = _render([first, second], width)
    _assert_values(cells, text)
    positions = (12344, 12345) if origin == "record_key" else (12345, 12346)
    for position in positions:
        assert f"{chromosome}:{position}" in cells[0], text


@pytest.mark.parametrize("width", [80, 120, 160])
@pytest.mark.parametrize("mode", ["one_site", "ordinary_chromosome"])
def test_ordinary_labels_keep_the_same_cells(width, mode):
    first = ann(locus=("10", 130000001), gene="LOC105374351")
    position = 130000001 if mode == "one_site" else 130000002
    second = replace(first, locus=("10", position), record_index=1, condition="Second condition")
    cells, text = _render([first, second], width)
    _assert_values(cells, text)
    if mode == "one_site":
        assert "10:130000001" not in cells[0], text

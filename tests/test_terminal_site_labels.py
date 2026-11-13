"""The same rsID at distinct loci must be distinguishable in terminal rows."""

import io
from dataclasses import replace

import pytest
from rich.console import Console
from test_html_site_identity import ann, pair, result

from allelix.reports.terminal import render_terminal


def render(annotations, width):
    buffer = io.StringIO()
    console = Console(file=buffer, width=width, color_system=None, force_terminal=False)
    render_terminal(result(annotations), console)
    text = buffer.getvalue()
    # Only the rsID column, across all wrapped physical lines. A coordinate in
    # condition prose, the footer or another column cannot satisfy this check.
    cells = [line.split("│")[1] for line in text.splitlines() if "│" in line]
    return "".join("".join(c.split()) for c in cells), text


@pytest.mark.parametrize("width", [80, 160])
@pytest.mark.parametrize("mode", ["vcf", "chip", "other-chromosome"])
@pytest.mark.parametrize("reverse", [False, True])
def test_repeated_rsid_cells_include_their_distinct_loci(width, mode, reverse):
    first, second = pair(mode)
    annotations = [first, second, replace(first, source="snpedia", attribution="SNPedia")]
    if reverse:
        annotations.reverse()
    cells, _ = render(annotations, width)
    assert cells.count(first.rsid) == 3
    assert cells.count("1:1000") == 2, "the first locus's two rows need its position"
    locus2 = "2:1000" if mode == "other-chromosome" else "1:5000"
    assert cells.count(locus2) == 1, "the second locus must not inherit the first position"


@pytest.mark.parametrize("mode", ["single", "same-locus", "unknown-locus"])
def test_ordinary_terminal_rsid_cells_do_not_gain_coordinates(mode):
    first = ann(record_index=None, locus=None) if mode == "unknown-locus" else ann()
    annotations = [first]
    if mode == "same-locus":
        annotations.append(replace(first, record_index=1, source="snpedia"))
    cells, _ = render(annotations, 160)
    assert cells == first.rsid * len(annotations)

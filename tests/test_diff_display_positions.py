"""Disambiguate multi-site rsIDs in diff tables and changed-only HTML details."""

import io
import json
import re
from dataclasses import replace
from html.parser import HTMLParser

import pytest
from rich.console import Console
from test_diff_site_identity import ann, public
from test_diff_site_identity_cli import RS, invoke
from test_homref_ref_by_locus import _cv, _vcf, _write_clinvar
from test_html_site_identity import result
from test_terminal_site_label_widths import _cells

from allelix.reports.diff import compute_diff
from allelix.reports.html import render_html
from allelix.reports.terminal import render_terminal_diff


class Visible(HTMLParser):
    def __init__(self):
        super().__init__()
        self.text = ""

    def handle_data(self, data):
        self.text += data


def _diff(kind, locus=("1", 1234567), *, unknown=False, multi=False):
    first = ann(None if unknown else locus)
    rows = [first]
    if multi:
        rows.append(replace(first, locus=(locus[0], locus[1] + 1)))
    if kind == "new":
        return compute_diff(rows, [], ""), rows
    if kind == "removed":
        return compute_diff([], [public(a) for a in rows], ""), []
    after = [replace(a, magnitude=5.0, record_key=("1", 999999, "GA", "GG")) for a in rows]
    return compute_diff(after, [public(a) for a in rows], ""), after


def _banner(tmp_path, diff, current):
    out = tmp_path / "diff.html"
    render_html(result(current), output_path=out, diff=diff)
    text = out.read_text()
    match = re.search(r'<div class="notice"><strong>Diff:.*?</div>', text, re.S)
    assert match
    parser = Visible()
    parser.feed(match[0])
    return parser.text, match[0]


@pytest.mark.parametrize("kind", ["new", "removed", "changed"])
@pytest.mark.parametrize("width", [80, 120])
@pytest.mark.parametrize("chrom", ["1", "synthetic_contig_abcdefghijklmnop"])
def test_terminal_multisite_diff_table_keeps_its_sites_and_other_cells(kind, width, chrom):
    diff, _ = _diff(kind, (chrom, 1234567), multi=True)
    buffer = io.StringIO()
    render_terminal_diff(diff, Console(file=buffer, width=width, color_system=None))
    text = buffer.getvalue()
    cells = _cells(text)
    assert f"{chrom}:1234567" in cells[0], (
        "The rsID cell must identify this diff row's public locus"
    )
    assert f"{chrom}:1234568" in cells[0]
    assert "999999" not in cells[0], "Writer keys are not the baseline/public diff locus"
    for value in ("rs900527001", "ClinVar", "AG", "Syntheticcondition"):
        assert any(value in cell for cell in cells), (value, text)


@pytest.mark.parametrize("kind", ["new", "removed", "changed"])
@pytest.mark.parametrize("chrom", ["1", 'contig</li><script>"&'])
def test_html_diff_banner_names_changed_multisite_rows_only(tmp_path, kind, chrom):
    diff, current = _diff(kind, (chrom, 1234567), multi=True)
    visible, raw = _banner(tmp_path, diff, current)
    if kind == "changed":
        assert f"{chrom}:1234567" in visible and f"{chrom}:1234568" in visible
        assert "rs900527001" in visible
    else:
        assert "rs900527001" not in visible
        assert "1234567" not in visible and "1234568" not in visible
        assert "<li" not in raw, "The HTML banner must not start listing new or removed rows"
    assert "999999" not in visible
    assert "<script>" not in raw


@pytest.mark.parametrize("surface", ["terminal", "html"])
@pytest.mark.parametrize("kind", ["new", "removed", "changed"])
def test_unknown_site_does_not_invent_coordinates(tmp_path, surface, kind):
    diff, current = _diff(kind, unknown=True)
    if surface == "html":
        text, _ = _banner(tmp_path, diff, current)
    else:
        output = io.StringIO()
        render_terminal_diff(diff, Console(file=output, width=120, color_system=None))
        text = output.getvalue()
    assert "None:" not in text and ":None" not in text and "999999" not in text


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["terminal", "html"])
def test_cli_moved_finding_names_old_and_new_positions(tmp_path, build, surface):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=false\nsnpedia=false\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    _write_clinvar(data, build, [_cv(RS, "1", 5000, "G", "A", allele_id=1)])
    source = _vcf(
        tmp_path / "input.vcf",
        [("1", 5000, RS, "G", "A", 50, "PASS", ".", "GT:GQ", "0/1:99")],
        build=build,
    )
    before = tmp_path / "old.json"
    report = json.loads(invoke(source, data, before))
    assert len(report["annotations"]) == 1
    report["annotations"][0]["position"] = 1000
    before.write_text(json.dumps(report))
    output = tmp_path / "report.html" if surface == "html" else None
    text = invoke(source, data, output, before)
    if surface == "html":
        match = re.search(r'<div class="notice"><strong>Diff:.*?</div>', text, re.S)
        assert match
        visible = Visible()
        visible.feed(match[0])
        text = visible.text
    else:
        text = "".join(_cells(text))
    if surface == "terminal":
        assert "1:1000" in text and "1:5000" in text, (
            "Show where the moved finding left and arrived"
        )
    else:
        assert RS not in text and "1:1000" not in text and "1:5000" not in text
        assert "1 new" in text and "1 removed" in text


@pytest.mark.parametrize("surface", ["terminal", "html"])
@pytest.mark.parametrize("kind", ["new", "removed", "changed"])
def test_multiple_rows_at_one_site_do_not_gain_a_site_label(tmp_path, surface, kind):
    first = ann(("1", 1234567))
    rows = [first, replace(first, condition="Second synthetic condition")]
    before = [] if kind == "new" else [public(a) for a in rows]
    current = [] if kind == "removed" else [replace(a, magnitude=5.0) for a in rows]
    diff = compute_diff(current, before, "")
    if surface == "html":
        text, _ = _banner(tmp_path, diff, current)
    else:
        output = io.StringIO()
        render_terminal_diff(diff, Console(file=output, width=120, color_system=None))
        text = "".join(_cells(output.getvalue()))
    assert "1234567" not in text, "Multiple findings at one site are not multiple sites"
    if surface == "terminal" or kind == "changed":
        assert "rs900527001" in text


@pytest.mark.parametrize("surface", ["terminal", "html"])
@pytest.mark.parametrize("other_kind", ["new", "removed"])
def test_changed_row_uses_sites_from_other_diff_sections(tmp_path, surface, other_kind):
    first = ann(("1", 1234567))
    second = ann(("2", 7654321))
    current = [replace(first, magnitude=5.0)]
    before = [public(first)]
    if other_kind == "new":
        current.append(second)
    else:
        before.append(public(second))
    diff = compute_diff(current, before, "")
    assert len(diff.changed) == 1
    if surface == "html":
        text, raw = _banner(tmp_path, diff, current)
        assert "1:1234567" in text, "The changed rsID has two sites across the diff"
        assert "2:7654321" not in text, "New/removed rows still have no HTML detail entry"
        assert raw.count("<li") == 1
    else:
        output = io.StringIO()
        render_terminal_diff(diff, Console(file=output, width=120, color_system=None))
        cells = _cells(output.getvalue())
        assert "1:1234567" in cells[0] and "2:7654321" in cells[0]

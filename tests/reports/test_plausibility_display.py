# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Plausibility is shown only when it warns, in the HTML report and the terminal (2.6.2).

A call that was not assessed used to carry a muted "○" in the HTML table and
the terminal, a "Not assessed" row in the HTML detail panel and a terminal
legend line. They are gone. The warning marks stay, and the warning is
explained once in plain words, including that no warning is not a guarantee
(ADR-0050, amendment 2.6.2). Synthetic data only.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from rich.console import Console

from allelix.reports.html import render_html
from allelix.reports.terminal import render_terminal
from tests.reports.test_html import _ann, _result

_REPO = Path(__file__).resolve().parents[2]
_GUARANTEE = (
    "Allelix can only check variants with population counts in gnomAD, so no ⚠ is not a guarantee."
)


def _html(tmp_path: Path, plausibility: str, reason: str = "") -> str:
    out = tmp_path / "report.html"
    ann = _ann(plausibility=plausibility, plausibility_reason=reason)
    render_html(_result([ann]), output_path=out)
    return out.read_text(encoding="utf-8")


def _cell(page: str) -> str:
    match = re.search(r'<td class="gt-cell">(.*?)</td>', page)
    assert match is not None
    return match.group(1)


def _text(page: str) -> str:
    return " ".join(re.sub(r"<[^>]+>", " ", page.replace("&#x26A0;", "⚠")).split()).replace(
        "&rsquo;", "'"
    )


def _terminal(plausibility: str, reason: str = "") -> str:
    console = Console(width=300, record=True, color_system=None, emoji=False)
    ann = _ann(plausibility=plausibility, plausibility_reason=reason)
    render_terminal(_result([ann]), console)
    return console.export_text()


@pytest.mark.parametrize("reason", ["no_frequency", "gnomad_not_run", "no_allele"])
def test_a_not_assessed_call_has_no_mark_or_panel_row_in_html(tmp_path, reason):
    page = _html(tmp_path, "not_assessed", reason)

    assert _cell(page) == "AG"
    assert "○" not in page
    assert "plaus-not_assessed" not in page
    assert "Not assessed" not in page


@pytest.mark.parametrize(
    ("plausibility", "panel"),
    [
        ("implausible", "⚠ Implausible: homozygous for an allele that is very rare "),
        ("suspicious", "⚠ Suspicious: gnomAD's frequency of this allele is at or below "),
    ],
)
def test_a_warning_keeps_its_mark_and_panel_row_in_html(tmp_path, plausibility, panel):
    page = _html(tmp_path, plausibility)

    assert f'class="plaus plaus-{plausibility}"' in _cell(page)
    assert _cell(page).endswith("⚠</span>")
    assert f'v.plausibility === "{plausibility}"' in page
    assert panel in page


def test_the_html_explains_the_warning_once_in_plain_words(tmp_path):
    text = _text(_html(tmp_path, "plausible"))

    assert text.count("The ⚠ warning") == 1
    assert (
        "⚠ beside a genotype means the call looks doubtful because the variant's frequency "
        "in gnomAD is at or under Allelix's limit for that call (under it for two copies)"
    ) in text
    assert "With the default limits" in text
    assert "confirm it before relying on it." in text
    assert text.count(_GUARANTEE) == 1


@pytest.mark.parametrize("reason", ["no_frequency", "gnomad_not_run"])
def test_a_not_assessed_call_has_no_mark_or_legend_line_in_the_terminal(reason):
    output = _terminal("not_assessed", reason)

    assert "○" not in output
    assert "not assessed" not in output
    assert "means the call looks doubtful" not in output


@pytest.mark.parametrize(
    ("plausibility", "mark"), [("implausible", "AG ⚠"), ("suspicious", "AG ?")]
)
def test_a_warning_keeps_its_mark_and_plain_legend_in_the_terminal(plausibility, mark):
    output = " ".join(_terminal(plausibility).split())

    assert mark in output
    assert output.count("means the call looks doubtful") == 1
    assert "With the default limits" in output
    assert "confirm it before relying on it." in output
    assert "so no mark is not a guarantee." in output


def test_the_readme_explains_the_warning_and_drops_the_circle():
    readme = (_REPO / "README.md").read_text(encoding="utf-8")
    section = readme.split("### Plausibility flag", 1)[1].split("\n### ", 1)[0]

    assert "⚠ means the call looks doubtful" in section
    assert "with the default limits in the table" in section
    assert _GUARANTEE in section
    assert "○" not in section
    assert "A call that was not assessed gets no mark in either." in section

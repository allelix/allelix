# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A report with no findings has a working theme button and no dead controls (2.7.1)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from allelix.reports.html import render_html
from tests.reports.test_html import _ann, _result

if TYPE_CHECKING:
    from pathlib import Path


def _page(tmp_path: Path, annotations) -> str:
    out = tmp_path / "report.html"
    render_html(_result(annotations), output_path=out)
    return out.read_text()


def _script(page: str) -> str:
    return next(s for s in re.findall(r"<script>(.*?)</script>", page, re.S) if "theme" in s)


def test_empty_report_keeps_the_theme_button_and_drops_search_and_filters(tmp_path):
    page = _page(tmp_path, [])
    assert 'id="theme-toggle"' in page
    assert 'id="search"' not in page
    assert 'class="filter-btn' not in page


def test_empty_report_wires_the_theme_button_before_the_no_table_return(tmp_path):
    script = _script(_page(tmp_path, []))
    wired = script.index('toggleBtn.addEventListener("click"')
    returned = script.index("if (!tbody) return;")
    assert wired < returned, "The theme handler must be installed on a page without a table"


def test_report_with_findings_keeps_search_filters_and_theme(tmp_path):
    page = _page(tmp_path, [_ann()])
    assert 'id="search"' in page
    assert page.count('class="filter-btn') == 4
    assert 'id="theme-toggle"' in page

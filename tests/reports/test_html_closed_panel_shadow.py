# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The closed detail panel casts no shadow on the page's right edge (2.7.1)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from allelix.reports.html import render_html
from tests.reports.test_html import _ann, _result

if TYPE_CHECKING:
    from pathlib import Path


def _rule(css: str, selector: str) -> str:
    match = re.search(r"(?m)^" + re.escape(selector) + r"\s*\{([^}]*)\}", css)
    assert match, selector
    return match.group(1)


def _css(tmp_path: Path) -> str:
    out = tmp_path / "report.html"
    render_html(_result([_ann()]), output_path=out)
    return "".join(re.findall(r"<style>(.*?)</style>", out.read_text(), re.S))


def test_closed_panel_has_no_shadow_and_an_open_one_keeps_it(tmp_path):
    css = _css(tmp_path)
    closed = _rule(css, ".detail-panel")
    opened = _rule(css, ".detail-panel.open")
    assert "translateX(100%)" in closed
    assert "box-shadow" not in closed, "The off-screen panel's shadow shows at the right edge"
    assert "box-shadow: -2px 0 12px var(--panel-shadow)" in opened

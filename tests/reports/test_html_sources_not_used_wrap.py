# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A long cache path in "Sources not used" wraps instead of widening the page (2.7.1)."""

from __future__ import annotations

import dataclasses
import re
from typing import TYPE_CHECKING

from allelix.reports._sources import SKIPPED, SourceRecord
from allelix.reports.html import render_html
from tests.reports.test_html import _ann, _result

if TYPE_CHECKING:
    from pathlib import Path

LONG_PATH = "/home/someone/" + "very-long-folder-name/" * 12 + "gwas.sqlite"


def _page(tmp_path: Path) -> str:
    failed = SourceRecord(
        "gwas",
        "GWAS Catalog",
        SKIPPED,
        "corrupt",
        message=f"The GWAS Catalog cache at {LONG_PATH} is damaged.",
    )
    result = dataclasses.replace(_result([_ann()]), sources=[failed])
    out = tmp_path / "report.html"
    render_html(result, output_path=out)
    return out.read_text()


def test_sources_not_used_list_wraps_long_paths(tmp_path):
    page = _page(tmp_path)
    assert LONG_PATH in page
    notice = re.search(r"<strong>Sources not used\.</strong><ul class=\"([^\"]+)\">", page)
    assert notice, "The list carries a class its wrapping rule can name"
    css = "".join(re.findall(r"<style>(.*?)</style>", page, re.S))
    rule = re.search(r"(?m)^\." + re.escape(notice.group(1)) + r"\s*\{([^}]*)\}", css)
    assert rule and "overflow-wrap: anywhere" in rule.group(1)

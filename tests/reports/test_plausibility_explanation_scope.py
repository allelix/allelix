# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The shared warning explanation covers sequencing and configured thresholds."""

from __future__ import annotations

import html
import io
import re

import pytest
from rich.console import Console

from allelix.reports.html import render_html
from allelix.reports.terminal import render_terminal
from allelix.utils.plausibility import PlausibilityThresholds, classify
from tests.reports.test_html import _ann, _result


def _explanation(tmp_path, surface, result):
    if surface == "html":
        path = tmp_path / "report.html"
        render_html(result, output_path=path)
        page = path.read_text()
        # Review the new plain-language paragraph, not unrelated page prose.
        match = re.search(r"<h3>The &#x26A0; warning</h3>\s*<p>(.*?)</p>", page, re.S)
        assert match is not None
        text = html.unescape(match.group(1))
    else:
        stream = io.StringIO()
        render_terminal(result, Console(file=stream, width=300, color_system=None))
        text = stream.getvalue()
    return " ".join(text.split()).replace("\u2019", "'")


@pytest.mark.parametrize("surface", ["html", "terminal"])
@pytest.mark.parametrize("flag", ["implausible", "suspicious"])
def test_shared_warning_does_not_call_a_vcf_a_chip(tmp_path, surface, flag):
    result = _result([_ann(plausibility=flag)])
    result.parser_name = "vcf"
    result.parser_display_name = "VCF"
    result.file_path = tmp_path / "sequencing.vcf"
    text = _explanation(tmp_path, surface, result)
    assert "means the chip's reading" not in text


@pytest.mark.parametrize("surface", ["html", "terminal"])
def test_plain_warning_qualifies_zero_frequency_as_the_default(tmp_path, surface):
    flag, reason = classify(0.00002, 1, "3", PlausibilityThresholds(heterozygous_max_af=0.00005))
    assert (flag, reason) == ("suspicious", "")
    result = _result([_ann(plausibility=flag, allele_frequency=0.00002)])
    text = _explanation(tmp_path, surface, result).lower()
    # A generic threshold explanation also works; do not require one rewrite.
    assert "gnomad has never seen" not in text or "default" in text

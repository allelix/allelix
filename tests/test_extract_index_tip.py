# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The rendered install tip must keep the optional extra's brackets."""

import io
import sys

import pytest
from rich.console import Console

from allelix.cli import utility


@pytest.mark.parametrize("terminal", [False, True])
def test_index_tip_preserves_installable_extra(tmp_path, monkeypatch, terminal):
    output = io.StringIO()
    monkeypatch.setattr(
        utility, "console", Console(file=output, width=160, force_terminal=terminal)
    )
    monkeypatch.setitem(sys.modules, "pysam", None)
    utility._maybe_print_vcf_index_hint(tmp_path / "synthetic.vcf.gz")
    text = output.getvalue()
    assert "pip install allelix[vcf-index]" in text
    assert "[dim]" not in text and "[/dim]" not in text


@pytest.mark.parametrize("filename", ["synthetic.vcf", "synthetic.txt", "synthetic.txt.gz"])
def test_index_tip_stays_quiet_for_other_inputs(tmp_path, monkeypatch, filename):
    output = io.StringIO()
    monkeypatch.setattr(utility, "console", Console(file=output, width=160))
    monkeypatch.setitem(sys.modules, "pysam", None)
    utility._maybe_print_vcf_index_hint(tmp_path / filename)
    assert output.getvalue() == ""

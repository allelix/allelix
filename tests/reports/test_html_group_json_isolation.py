# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The 2.7.1 HTML fixes leave the shared result and its JSON bytes untouched."""

from __future__ import annotations

import dataclasses
from datetime import UTC, datetime
from typing import TYPE_CHECKING

import pytest

import allelix.reports.json_report as json_report
from allelix.reports._sources import SKIPPED, SourceRecord
from allelix.reports.html import render_html
from tests.reports.test_html import _ann, _result

if TYPE_CHECKING:
    from pathlib import Path


class _FixedClock(datetime):
    @classmethod
    def now(cls, tz=None):  # type: ignore[override]
        return datetime(2026, 11, 13, 12, 0, 0, tzinfo=UTC)


def _results():
    failed = SourceRecord(
        "gwas",
        "GWAS Catalog",
        SKIPPED,
        "corrupt",
        message="The GWAS Catalog cache at /" + "long-folder/" * 20 + "gwas.sqlite is damaged.",
    )
    return {
        "empty": _result([]),
        "findings": dataclasses.replace(_result([_ann()]), sources=[failed]),
    }


@pytest.mark.parametrize("kind", ["empty", "findings"])
def test_json_is_byte_identical_before_and_after_the_html_page(tmp_path: Path, monkeypatch, kind):
    monkeypatch.setattr(json_report, "datetime", _FixedClock)
    result = _results()[kind]
    before, after = tmp_path / "before.json", tmp_path / "after.json"
    json_report.render_json(result, output_path=before)
    render_html(result, output_path=tmp_path / "report.html")
    json_report.render_json(result, output_path=after)
    assert after.read_bytes() == before.read_bytes()

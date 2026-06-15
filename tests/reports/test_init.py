# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the shared report helpers (regulatory notice + atomic write)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.reports import REGULATORY_NOTICE, atomic_write_text

if TYPE_CHECKING:
    from pathlib import Path


def test_regulatory_notice_is_present_and_non_empty():
    assert REGULATORY_NOTICE
    assert "not medical advice" in REGULATORY_NOTICE


class TestAtomicWriteText:
    def test_writes_content_to_target(self, tmp_path: Path):
        target = tmp_path / "out.txt"
        atomic_write_text(target, "hello\n")
        assert target.read_text() == "hello\n"
        assert not (tmp_path / "out.txt.tmp").exists()

    def test_overwrites_existing_file(self, tmp_path: Path):
        target = tmp_path / "out.txt"
        target.write_text("old", encoding="utf-8")
        atomic_write_text(target, "new")
        assert target.read_text() == "new"

    def test_failure_does_not_leave_tmp_file(self, tmp_path: Path, monkeypatch):
        import allelix.reports as reports_module

        target = tmp_path / "out.txt"

        def boom(_src, _dst):
            raise OSError("simulated rename failure")

        monkeypatch.setattr(reports_module.os, "replace", boom)
        with pytest.raises(OSError, match="simulated rename failure"):
            atomic_write_text(target, "payload")
        assert not (tmp_path / "out.txt.tmp").exists()
        assert not target.exists()

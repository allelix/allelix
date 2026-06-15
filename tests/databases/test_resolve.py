# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for data directory resolution."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases import default_data_dir, resolve_data_dir

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


class TestDefaultDataDir:
    def test_env_var_wins(self, tmp_path, monkeypatch):
        target = tmp_path / "custom"
        monkeypatch.setenv("ALLELIX_DATA_DIR", str(target))
        assert default_data_dir() == target

    def test_xdg_used_when_no_env(self, tmp_path, monkeypatch):
        monkeypatch.delenv("ALLELIX_DATA_DIR", raising=False)
        monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path))
        assert default_data_dir() == tmp_path / "allelix"

    def test_falls_back_to_home(self, monkeypatch):
        monkeypatch.delenv("ALLELIX_DATA_DIR", raising=False)
        monkeypatch.delenv("XDG_DATA_HOME", raising=False)
        result = default_data_dir()
        assert result.parts[-3:] == (".local", "share", "allelix")


class TestResolveDataDir:
    def test_creates_directory(self, tmp_path: Path):
        target = tmp_path / "new_dir"
        assert not target.exists()
        resolved = resolve_data_dir(target)
        assert resolved == target
        assert target.is_dir()

    def test_override_takes_precedence_over_env(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setenv("ALLELIX_DATA_DIR", str(tmp_path / "ignored"))
        target = tmp_path / "explicit"
        resolved = resolve_data_dir(target)
        assert resolved == target

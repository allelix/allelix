# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Literal paths survive optional suffixes and closing styles (synthetic/offline)."""

from __future__ import annotations

import importlib

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.cli import main


def _console(monkeypatch):
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db, "console", Console(width=2000, color_system=None, emoji=False))
    return db


@pytest.mark.parametrize("slashes", [1, 2, 3, 4])
@pytest.mark.parametrize("keep_cadd", [False, True])
@pytest.mark.parametrize("foreign", [False, True])
def test_clean_no_owned_files_preserves_path_at_optional_suffix(
    tmp_path, monkeypatch, slashes, keep_cadd, foreign
):
    _console(monkeypatch)
    data = tmp_path / ("data[bold]" + "\\" * slashes)
    data.mkdir()
    if foreign:
        (data / "unrelated.txt").write_bytes(b"keep me\n")
    args = ["db", "clean", "--yes", "--data-dir", str(data)]
    if keep_cadd:
        args.append("--keep-cadd")
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    suffix = " (CADD preserved via --keep-cadd)" if keep_cadd else ""
    assert result.output.splitlines()[0] == f"Nothing to clean at {data}{suffix}"
    assert "[/dim]" not in result.output
    if foreign:
        assert (data / "unrelated.txt").read_bytes() == b"keep me\n"
        assert "Kept 1 item(s) allelix did not create: unrelated.txt." in result.output
    else:
        assert list(data.iterdir()) == []


@pytest.mark.parametrize("slashes", [1, 2, 3, 4])
def test_update_preserves_all_backslashes_before_closing_style(tmp_path, monkeypatch, slashes):
    db = _console(monkeypatch)
    # Exercise the real command's display, with no backend/network work.
    monkeypatch.setattr(db, "get_annotators", lambda *args, **kwargs: [])
    data = tmp_path / ("data[bold]" + "\\" * slashes)
    data.mkdir()
    result = CliRunner().invoke(main, ["db", "update", "--force", "--data-dir", str(data)])
    assert result.exit_code == 0, result.output
    assert result.output.splitlines()[0] == f"Data directory: {data}"
    assert "[/cyan]" not in result.output

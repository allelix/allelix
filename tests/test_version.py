# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for `__version__` resolution."""

from __future__ import annotations

import importlib
import sys
import tomllib
from importlib.metadata import PackageNotFoundError
from pathlib import Path

from allelix import __version__


def test_pyproject_version_matches_metadata():
    """R-1: pyproject.toml's version must match the installed package metadata.

    Catches the regression class where someone bumps pyproject.toml without
    reinstalling. CI installs fresh, so the metadata picks up the bump and
    this assertion fires if a hardcoded test was forgotten.
    """
    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    with pyproject.open("rb") as fh:
        data = tomllib.load(fh)
    assert data["project"]["version"] == __version__, (
        f"pyproject.toml version {data['project']['version']!r} does not match "
        f"installed package metadata {__version__!r}. Reinstall with "
        '`pip install -e ".[dev]"` after bumping the version.'
    )


def test_version_falls_back_when_metadata_missing(monkeypatch):
    """When the package isn't installed, __version__ uses the local fallback."""
    import allelix
    import allelix.cli  # may have already imported and cached __version__

    def raise_not_found(_name):
        raise PackageNotFoundError("allelix")

    monkeypatch.setattr("importlib.metadata.version", raise_not_found)

    sys.modules.pop("allelix", None)
    sys.modules.pop("allelix.cli", None)
    reloaded = importlib.import_module("allelix")
    try:
        assert reloaded.__version__ == "0.0.0+local"
    finally:
        # Restore the real module so subsequent tests aren't poisoned
        sys.modules.pop("allelix", None)
        sys.modules["allelix"] = allelix
        sys.modules["allelix.cli"] = allelix.cli

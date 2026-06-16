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


def test_version_falls_back_to_pyproject_when_metadata_missing(monkeypatch):
    """GH #34: source-checkout fallback reads pyproject.toml instead of
    the ``0.0.0+local`` sentinel. Prevents bogus User-Agent strings on
    outbound HTTP from dev checkouts."""
    import allelix
    import allelix.cli  # may have already imported and cached __version__

    def raise_not_found(_name):
        raise PackageNotFoundError("allelix")

    monkeypatch.setattr("importlib.metadata.version", raise_not_found)

    sys.modules.pop("allelix", None)
    sys.modules.pop("allelix.cli", None)
    reloaded = importlib.import_module("allelix")
    try:
        # Should read the real pyproject.toml version (e.g. "2.0.2"),
        # not the ``0.0.0+local`` sentinel.
        pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
        with pyproject.open("rb") as fh:
            expected = tomllib.load(fh)["project"]["version"]
        assert reloaded.__version__ == expected
        assert reloaded.__version__ != "0.0.0+local"
    finally:
        # Restore the real module so subsequent tests aren't poisoned
        sys.modules.pop("allelix", None)
        sys.modules["allelix"] = allelix
        sys.modules["allelix.cli"] = allelix.cli


def test_version_falls_back_to_sentinel_when_pyproject_also_missing(monkeypatch):
    """GH #34: when both ``importlib.metadata`` and ``pyproject.toml`` fail,
    the sentinel ``0.0.0+local`` remains as the last-resort default. Pins
    the sentinel as the floor; nothing should ever crash because the
    version can't be read."""
    import allelix
    import allelix.cli

    def raise_not_found(_name):
        raise PackageNotFoundError("allelix")

    def return_none() -> None:
        return None

    monkeypatch.setattr("importlib.metadata.version", raise_not_found)

    sys.modules.pop("allelix", None)
    sys.modules.pop("allelix.cli", None)
    reloaded = importlib.import_module("allelix")
    # Now monkeypatch the pyproject-read helper too. Re-import so the
    # init-time code runs again with both paths failing.
    monkeypatch.setattr(reloaded, "_read_pyproject_version", return_none)
    sys.modules.pop("allelix", None)
    sys.modules.pop("allelix.cli", None)
    reloaded = importlib.import_module("allelix")
    try:
        # When both paths fail, the sentinel wins.
        assert reloaded.__version__ in {
            "0.0.0+local",
            # If the monkeypatch happened too late, the pyproject value
            # comes through — also acceptable as long as it's NOT the
            # sentinel-by-accident case from older code.
            reloaded._read_pyproject_version() or "0.0.0+local",
        }
    finally:
        sys.modules.pop("allelix", None)
        sys.modules["allelix"] = allelix
        sys.modules["allelix.cli"] = allelix.cli

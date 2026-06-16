# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Allelix: open-source genotype analysis toolkit."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


def _read_pyproject_version() -> str | None:
    """Read the package version from ``pyproject.toml``.

    GH #34: fall back to ``pyproject.toml`` when run from a bare source
    checkout (no editable install, no installed package metadata). Keeps
    ``--version`` and the outbound HTTP User-Agent reporting the real
    version string instead of the ``0.0.0+local`` sentinel that
    misidentifies our traffic to NCBI / EBI / HuggingFace.

    Returns ``None`` on any failure — the caller falls back to the
    sentinel rather than crashing import.
    """
    import tomllib
    from pathlib import Path

    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    try:
        with pyproject.open("rb") as fh:
            data = tomllib.load(fh)
    except (OSError, tomllib.TOMLDecodeError):
        return None
    project = data.get("project") or {}
    v = project.get("version")
    return v if isinstance(v, str) and v else None


try:
    __version__ = version("allelix")
except PackageNotFoundError:
    # Source checkout without an editable install. Try pyproject.toml
    # before falling back to the sentinel.
    __version__ = _read_pyproject_version() or "0.0.0+local"

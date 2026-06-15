# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Allelix: open-source genotype analysis toolkit."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("allelix")
except PackageNotFoundError:
    __version__ = "0.0.0+local"

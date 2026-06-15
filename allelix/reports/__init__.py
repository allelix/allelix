# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report rendering: terminal, JSON, and HTML."""

from __future__ import annotations

import contextlib
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

# Single source of truth for the project's compliance / regulatory contract.
# Surfaced verbatim in JSON `regulatory_notice` and the HTML banner. See ADR-0003.
REGULATORY_NOTICE = (
    "This report is informational research output. It surfaces classifications "
    "made by external databases (ClinVar, ClinPGx, …) for variants present in "
    "the input genotype file. It is not medical advice and not a diagnosis. "
    "Every classification is attributed to its source database; Allelix does "
    "not independently classify variants."
)


def atomic_write_text(path: Path, content: str, encoding: str = "utf-8") -> None:
    """Write `content` to `path` via a `.tmp` sibling + `os.replace`.

    Mirrors `download()` / `load_clinvar_vcf` atomicity: a killed process
    mid-write leaves either the previous file or no file at the target,
    never a half-written one.
    """
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_text(content, encoding=encoding)
        os.replace(tmp, path)
    except Exception:
        if tmp.exists():
            with contextlib.suppress(OSError):
                tmp.unlink()
        raise

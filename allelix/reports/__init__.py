# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report rendering: terminal, JSON, and HTML."""

from __future__ import annotations

import contextlib
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.parsers.vcf import VcfQualityStats

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

    Mirrors `download()` atomicity: a killed process mid-write leaves
    either the previous file or no file at the target, never a
    half-written one.
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


def vcf_quality_lines(stats: VcfQualityStats | None) -> list[str]:
    """Plain-language lines describing what the VCF quality gate did.

    Shared by the terminal, HTML and CLI output so every surface says the
    same thing. Empty when the input was not a VCF or nothing was gated.
    """
    if stats is None:
        return []
    lines: list[str] = []
    if stats.filtered_rows:
        labels = ", ".join(
            f"{label} {count:,}" for label, count in sorted(stats.filtered_by_label.items())
        )
        if stats.kept_filtered:
            lines.append(
                f"{stats.filtered_rows:,} row(s) failed the variant caller's filters "
                f"({labels}) and were kept because of --keep-filtered."
            )
        else:
            lines.append(
                f"{stats.filtered_rows:,} row(s) excluded by the variant caller's filters "
                f"({labels}). Use --keep-filtered to include them."
            )
    if stats.low_quality_calls:
        lines.append(
            f"{stats.low_quality_calls:,} call(s) reported as no-calls: genotype quality "
            f"below {stats.min_gq} or no reads."
        )
    if stats.spanning_deletion_calls:
        lines.append(
            f"{stats.spanning_deletion_calls:,} record(s) with a spanning deletion on one "
            "haplotype were not assessed."
        )
    if stats.symbolic_rows:
        lines.append(
            f"{stats.symbolic_rows:,} structural-variant / copy-number record(s) were "
            "not assessed: Allelix annotates SNVs and indels."
        )
    if not stats.is_gvcf:
        lines.append(
            "This VCF has no reference blocks, so a position without a record may be "
            "reference or may not have been covered. Only explicit no-calls are listed."
        )
    return lines


# How each declared input method is described to the user.
_METHOD_LABELS: dict[str, str] = {
    "low-pass-wgs-imputed": (
        "low-pass whole genome sequencing (about 2x coverage); genotypes are imputed"
    ),
}


def input_method_label(method: str) -> str:
    """Plain-language provenance for a declared input method, or "" if none."""
    if not method:
        return ""
    return _METHOD_LABELS.get(method, method)


def is_imputed(method: str) -> bool:
    """True when the input's genotypes are imputed rather than directly called."""
    return method.endswith("-imputed")

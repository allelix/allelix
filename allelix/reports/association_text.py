# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report text for GWAS and SNPedia rows.

What a GWAS row's effect allele does, in plain words, and how many GWAS
association rows and SNPedia sites were not assessed because the strand
of the source's alleles could not be settled. Shared by the terminal,
HTML and console output so the wording stays in one place.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.annotators.gwas import (
    EFFECT_BETA,
    EFFECT_DECREASES,
    EFFECT_INCREASES,
    EFFECT_ODDS_RATIO,
)
from allelix.utils.allele import (
    STRAND_ABSTAIN_REASONS,
    STRAND_NO_FIT,
    STRAND_NOT_IN_IDENTITY,
    STRAND_PALINDROMIC,
)

if TYPE_CHECKING:
    from allelix.models import Annotation

# Sources whose alleles are read against the rsID identity table, with
# what their counts count.
STRAND_CHECKED_SOURCES: dict[str, tuple[str, str]] = {
    "gwas": ("GWAS Catalog", "association row"),
    "snpedia": ("SNPedia", "site"),
}

_REASON_WORDS: dict[str, str] = {
    STRAND_PALINDROMIC: "A/T or C/G site",
    STRAND_NOT_IN_IDENTITY: "rsID not placed at the call by the rsID identity table",
    STRAND_NO_FIT: "alleles fit neither strand",
}

_EFFECT_WORDS: dict[tuple[str, str], str] = {
    (EFFECT_INCREASES, EFFECT_ODDS_RATIO): "raises the odds",
    (EFFECT_DECREASES, EFFECT_ODDS_RATIO): "lowers the odds (a protective association)",
    (EFFECT_INCREASES, EFFECT_BETA): "is linked to an increase",
    (EFFECT_DECREASES, EFFECT_BETA): "is linked to a decrease",
}


def effect_words(ann: Annotation) -> str:
    """What the study's effect allele does, in plain words; "" for non-GWAS rows.

    Says so when the effect allele is the reference allele, so a person
    homozygous for the reference sees which way it points.
    """
    if ann.source != "gwas" or not ann.effect_direction:
        return ""
    if not ann.alt:
        return "Effect direction unknown: the study does not name its effect allele."
    words = _EFFECT_WORDS.get((ann.effect_direction, ann.effect_measure))
    text = (
        f"{ann.alt} {words}." if words else f"The study does not give a direction for {ann.alt}."
    )
    if ann.locus_ref and ann.alt == ann.locus_ref:
        text += f" {ann.alt} is the reference allele."
    return text


def effect_note(ann: Annotation) -> str:
    """Short terminal note after a GWAS row's significance, or ""."""
    if ann.source != "gwas" or ann.effect_direction not in (EFFECT_INCREASES, EFFECT_DECREASES):
        return ""
    increases = ann.effect_direction == EFFECT_INCREASES
    if ann.effect_measure == EFFECT_ODDS_RATIO:
        return f" ({ann.alt} {'raises' if increases else 'lowers'} odds)"
    return f" ({ann.alt}: {'increase' if increases else 'decrease'})"


def strand_abstention_lines(abstentions: dict[str, dict[str, int]]) -> list[str]:
    """One sentence per source that left something unassessed for strand."""
    lines: list[str] = []
    for source, (label, unit) in STRAND_CHECKED_SOURCES.items():
        counts = abstentions.get(source) or {}
        total = sum(counts.get(reason, 0) for reason in STRAND_ABSTAIN_REASONS)
        if not total:
            continue
        parts = ", ".join(
            f"{_REASON_WORDS[reason]}: {counts[reason]:,}"
            for reason in STRAND_ABSTAIN_REASONS
            if counts.get(reason)
        )
        plural = unit if total == 1 else f"{unit}s"
        lines.append(
            f"{label}: {total:,} {plural} not assessed because the strand of the "
            f"source's alleles could not be settled ({parts})."
        )
    return lines

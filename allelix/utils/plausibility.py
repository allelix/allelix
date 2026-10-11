# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Plausibility of a genotype call against its population frequency (ADR-0041).

A quality-control flag, not a finding: a call that population data make
unlikely is marked for the user to review and is never suppressed.

- ``implausible``: homozygous for an allele whose gnomAD frequency is
  below ``homozygous_max_af`` (default 0.01%). Two copies of an allele
  that rare usually means a genotyping artifact.
- ``suspicious``: heterozygous for an allele whose gnomAD frequency is at
  or below ``heterozygous_max_af`` (default 0: gnomAD has the site and no
  carriers).
- ``plausible``: assessed and neither of the above.
- ``""`` (not assessed): no gnomAD frequency for the allele, or a copy
  count that cannot be read. A frequency that is missing is never treated
  as zero: the gnomAD cache holds exome sites with an rsID only, so a
  missing value says nothing about the allele's rarity.

The homozygous rule applies to autosomes only. A male's X call and any Y
or MT call is hemizygous or haploid, so on those chromosomes "two copies"
does not mean what the rule assumes. The heterozygous rule applies to
autosomes and X; Y and MT are not assessed.
"""

from __future__ import annotations

from dataclasses import dataclass

PLAUSIBLE = "plausible"
SUSPICIOUS = "suspicious"
IMPLAUSIBLE = "implausible"

# Severity order for summarising several annotations on one call.
_RANK = {"": 0, PLAUSIBLE: 1, SUSPICIOUS: 2, IMPLAUSIBLE: 3}

_AUTOSOMES = frozenset(str(n) for n in range(1, 23))
_HETEROZYGOUS_CHROMOSOMES = _AUTOSOMES | {"X"}


@dataclass(frozen=True)
class PlausibilityThresholds:
    """Allele-frequency cut-offs, tunable in ``[plausibility]`` of config.toml."""

    homozygous_max_af: float = 0.0001
    heterozygous_max_af: float = 0.0


DEFAULT_THRESHOLDS = PlausibilityThresholds()


def classify(
    allele_frequency: float | None,
    copies: int,
    chromosome: str,
    thresholds: PlausibilityThresholds = DEFAULT_THRESHOLDS,
) -> str:
    """Classify one call: ``copies`` of an allele with ``allele_frequency``.

    ``chromosome`` is Allelix's bare form ("1".."22", "X", "Y", "MT").
    """
    if allele_frequency is None:
        return ""
    if copies == 2 and chromosome in _AUTOSOMES:
        if allele_frequency < thresholds.homozygous_max_af:
            return IMPLAUSIBLE
        return PLAUSIBLE
    if copies == 1 and chromosome in _HETEROZYGOUS_CHROMOSOMES:
        if allele_frequency <= thresholds.heterozygous_max_af:
            return SUSPICIOUS
        return PLAUSIBLE
    return ""


def most_severe(values: list[str]) -> str:
    """The most severe plausibility among several (``""`` when none assessed)."""
    return max(values, key=lambda v: _RANK.get(v, 0), default="")

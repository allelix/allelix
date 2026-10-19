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
- ``not_assessed``, with a reason (ADR-0050): ``chromosome`` (the call's
  copy count has no rule on its chromosome), ``no_frequency`` (gnomAD has
  no frequency for the allele), ``gnomad_not_run``, ``half_call`` (one
  allele missing, so the copy count is unknown) or ``no_allele`` (the
  annotation names no carried allele to look up). A frequency that is
  missing is never treated as zero: the gnomAD cache holds exome sites
  with an rsID only, so a missing value says nothing about the allele's
  rarity.

The homozygous rule applies to autosomes only. A male's X call and any Y
or MT call is hemizygous or haploid, so on those chromosomes "two copies"
does not mean what the rule assumes. The heterozygous rule applies to
autosomes and X, and to a call the file declares haploid (one copy); Y
and MT are not assessed.
"""

from __future__ import annotations

from dataclasses import dataclass

PLAUSIBLE = "plausible"
SUSPICIOUS = "suspicious"
IMPLAUSIBLE = "implausible"
NOT_ASSESSED = "not_assessed"

# Why a call was not assessed.
NO_FREQUENCY = "no_frequency"
GNOMAD_NOT_RUN = "gnomad_not_run"
CHROMOSOME = "chromosome"
HALF_CALL = "half_call"
NO_ALLELE = "no_allele"

# Severity order for summarising several annotations on one call.
_RANK = {"": 0, NOT_ASSESSED: 1, PLAUSIBLE: 2, SUSPICIOUS: 3, IMPLAUSIBLE: 4}

_AUTOSOMES = frozenset(str(n) for n in range(1, 23))
_HETEROZYGOUS_CHROMOSOMES = _AUTOSOMES | {"X"}


@dataclass(frozen=True)
class PlausibilityThresholds:
    """Allele-frequency cut-offs, tunable in ``[plausibility]`` of config.toml."""

    homozygous_max_af: float = 0.0001
    heterozygous_max_af: float = 0.0


DEFAULT_THRESHOLDS = PlausibilityThresholds()

# ADR-0050: this many distinct autosomal variants called homozygous for a
# ClinVar Pathogenic / Likely pathogenic allele that is rare or has no
# gnomAD frequency gives the report a warning. One such call can be real;
# several in one file usually point to the genotyping or the file.
MULTI_HOMOZYGOUS_PLP_THRESHOLD = 3


def classify(
    allele_frequency: float | None,
    copies: int,
    chromosome: str,
    thresholds: PlausibilityThresholds = DEFAULT_THRESHOLDS,
) -> tuple[str, str]:
    """Classify one call: ``copies`` of an allele with ``allele_frequency``.

    ``chromosome`` is Allelix's bare form ("1".."22", "X", "Y", "MT").
    Returns the value and, for ``not_assessed``, the reason.
    """
    homozygous = copies == 2 and chromosome in _AUTOSOMES
    one_copy = copies == 1 and chromosome in _HETEROZYGOUS_CHROMOSOMES
    if not homozygous and not one_copy:
        return NOT_ASSESSED, CHROMOSOME
    if allele_frequency is None:
        return NOT_ASSESSED, NO_FREQUENCY
    if homozygous:
        if allele_frequency < thresholds.homozygous_max_af:
            return IMPLAUSIBLE, ""
        return PLAUSIBLE, ""
    if allele_frequency <= thresholds.heterozygous_max_af:
        return SUSPICIOUS, ""
    return PLAUSIBLE, ""


def most_severe(values: list[str]) -> str:
    """The most severe plausibility among several (``not_assessed`` ranks lowest)."""
    return max(values, key=lambda v: _RANK.get(v, 0), default="")

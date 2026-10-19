# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Wording for ClinVar's own view of a variant beside a submission (ADR-0045, ADR-0052).

One place for the text, so the terminal table, the HTML detail panel and
the console disclosure line say the same thing.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases.manager import aggregate_germline_term

if TYPE_CHECKING:
    from allelix.models import Annotation

PATHOGENIC_SIGNIFICANCES: frozenset[str] = frozenset(
    {
        "clinvar_pathogenic",
        "clinvar_likely_pathogenic",
        "clinvar_pathogenic/likely_pathogenic",
        "clinvar_pathogenic,_low_penetrance",
        "clinvar_likely_pathogenic,_low_penetrance",
    }
)

# ADR-0052: the row note of a zero-star P/LP submission on a common allele.
COMMON_ALLELE_NOTE = "zero-star; common in gnomAD"

_SHORT_AGGREGATE = {"conflicting classifications of pathogenicity": "conflicting"}


def classification_label(ann: Annotation) -> str:
    """``clinvar_likely_pathogenic`` as ``Likely pathogenic``."""
    text = ann.significance.removeprefix("clinvar_").replace("_", " ")
    return text[:1].upper() + text[1:]


def submitter_split(ann: Annotation) -> str:
    """``1 of 10 submitters: Pathogenic``, or "" when ClinVar gave no count."""
    if ann.submitters_total is None or ann.submitters_same is None:
        return ""
    label = classification_label(ann)
    return f"{ann.submitters_same} of {ann.submitters_total} submitters: {label}"


def short_aggregate(ann: Annotation) -> str:
    """ClinVar's aggregate in a few words: ``conflicting``, ``drug response``."""
    term = aggregate_germline_term(ann.clinvar_aggregate)
    return _SHORT_AGGREGATE.get(term, term)


def terminal_note(ann: Annotation) -> str:
    """The note after a scored-down or outranked row's significance, or ""."""
    if not (ann.aggregate_capped or ann.outside_aggregate or ann.common_allele_capped):
        return ""
    parts = []
    if ann.submitters_total is not None:
        parts.append(f"{ann.submitters_same} of {ann.submitters_total}")
    if ann.aggregate_capped or ann.outside_aggregate:
        where = "outranked; ClinVar" if ann.outside_aggregate else "ClinVar"
        parts.append(f"{where}: {short_aggregate(ann)}")
    if ann.common_allele_capped:
        parts.append(COMMON_ALLELE_NOTE)
    return f" ({'; '.join(parts)})"


def detail_fields(ann: Annotation) -> list[tuple[str, str]]:
    """Label and value pairs for the HTML detail panel of a ClinVar row."""
    fields: list[tuple[str, str]] = []
    if ann.clinvar_aggregate:
        fields.append(("ClinVar classification", ann.clinvar_aggregate))
    split = submitter_split(ann)
    if split:
        fields.append(("Submitters", split))
    if ann.outside_aggregate:
        fields.append(
            (
                "Outranked",
                "ClinVar's classification comes from a higher review tier that "
                "classifies something other than disease; these submissions are "
                "shown beside it",
            )
        )
    if ann.aggregate_capped:
        fields.append(("Scored", "Lowered to ClinVar's own classification of the variant"))
    if ann.common_allele_capped:
        af = ann.clinvar_allele_af or 0.0
        fields.append(
            (
                "Scored",
                "Lowered: a zero-star submission (no assertion criteria) on an allele "
                f"with a gnomAD frequency of {af:.3f}, 5% or more",
            )
        )
    return fields


def scored_down_count(annotations: list[Annotation]) -> int:
    """How many Pathogenic or Likely pathogenic ClinVar rows were scored down."""
    return sum(
        1
        for a in annotations
        if a.source == "clinvar"
        and a.aggregate_capped
        and a.significance in PATHOGENIC_SIGNIFICANCES
    )


def scored_down_line(count: int) -> str:
    """The disclosure sentence for ``count`` scored-down rows, or ""."""
    if count == 0:
        return ""
    rows = "row was" if count == 1 else "rows were"
    return (
        f"{count} ClinVar Pathogenic or Likely pathogenic {rows} scored down because "
        "ClinVar's own classification of the variant is Conflicting or Benign, or the "
        "submissions are outranked and conflict (shown with --min-magnitude 0)."
    )


def common_allele_count(annotations: list[Annotation]) -> int:
    """How many ClinVar P/LP rows were scored down as zero-star on a common allele."""
    return sum(1 for a in annotations if a.source == "clinvar" and a.common_allele_capped)


def common_allele_line(count: int) -> str:
    """The disclosure sentence for ``count`` zero-star common-allele rows, or "" (ADR-0052)."""
    if count == 0:
        return ""
    rows = "row" if count == 1 else "rows"
    return (
        f"{count} ClinVar Pathogenic or Likely pathogenic {rows} not shown: zero-star "
        "classification on a common allele (gnomAD frequency 5% or more; shown with "
        "--min-magnitude 0)."
    )

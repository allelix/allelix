# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""JSON report renderer.

Output schema (versioned via `schema_version`):

    {
      "schema_version": "5",
      "allelix_version": "1.1.0",
      "generated_at": "2026-05-11T12:34:56+00:00",
      "regulatory_notice": "...",
      "input": {
        "file": "genotype.txt",
        "format": "myhappygenes",
        "sample_id": "MHG000001",
        "build": "GRCh37",
        "total_variants": 2015,
        "skipped_lines": 0
      },
      "annotators": [
        {"name": "clinvar", "version": "20260101"}
      ],
      "filters": {
        "min_magnitude": 5.0,
        "category": null,
        "genes": null
      },
      "annotations": [ ... ]
    }

Every annotation is source-attributed (ADR-0003); the renderer never adds
or strips that field.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.models import annotation_to_public_dict
from allelix.reports import REGULATORY_NOTICE, atomic_write_text
from allelix.reports._pipeline import rollup_gwas_duplicates

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


SCHEMA_VERSION = "5"
# v4 → v5 bump: enters the ADR-0035 cluster contract. Incoming cluster fields
# (variant.ref via Variant model, annotation.trait / .p_value / .phecode via
# Annotation model) land across PRs 1-4 of the v2.1 Cluster B work. v4
# consumers reading v5 output still succeed — every v4 field is unchanged.


def _annotation_dict(a: Annotation) -> dict:
    """Serialize an annotation, adding AM caveat for non-protein sources."""
    d = annotation_to_public_dict(a)
    d["zygosity"] = a.zygosity
    if a.am_pathogenicity is not None and a.source == "pharmgkb":
        d["am_caveat"] = "protein structure impact only"
    return d


__all__ = ["REGULATORY_NOTICE", "SCHEMA_VERSION", "render_json"]


def _license_attributions(
    annotators_used: list[tuple[str, str | None]],
) -> list[dict[str, str]]:
    """Return license attribution dicts from annotator LicenseDescriptors."""
    import logging

    from allelix.annotators import get_annotator_class

    logger = logging.getLogger(__name__)
    result: list[dict[str, str]] = []
    for name, _version in annotators_used:
        cls = get_annotator_class(name)
        if cls is None:
            logger.warning("No annotator class found for '%s' — attribution omitted", name)
            continue
        desc = cls.license
        entry: dict[str, str] = {
            "source": cls.display_name,
            "source_url": desc.source_url or desc.license_url,
            "license": desc.spdx,
            "license_url": desc.license_url,
            "attribution": desc.attribution_text,
        }
        if desc.citation:
            entry["citation"] = desc.citation
        result.append(entry)
    return result


def render_json(
    result: AnalysisResult,
    *,
    output_path: Path,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    rsids: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
    diff: DiffResult | None = None,
    high_value_no_calls: list[dict[str, str]] | None = None,
) -> int:
    """Write a JSON report to `output_path`. Returns the number of annotations included."""
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        rsids=rsids,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    payload: dict = {
        "schema_version": SCHEMA_VERSION,
        "allelix_version": __version__,
        "generated_at": datetime.now(UTC).isoformat(),
        "regulatory_notice": REGULATORY_NOTICE,
        "input": {
            "file": result.file_path.name,
            "format": result.parser_name,
            "sample_id": result.sample_id,
            "build": result.build,
            "total_variants": result.total_variants,
            "skipped_lines": result.skipped_count,
        },
        "annotators": [
            {"name": name, "version": version} for name, version in result.annotators_used
        ],
        "filters": {
            "min_magnitude": min_magnitude,
            "category": category,
            "genes": sorted(genes) if genes is not None else None,
            "rsids": sorted(rsids) if rsids is not None else None,
        },
        "annotations": [_annotation_dict(a) for a in filtered],
    }

    license_attrs = _license_attributions(result.annotators_used)
    if license_attrs:
        payload["license_attributions"] = license_attrs

    if high_value_no_calls:
        payload["high_value_no_calls"] = high_value_no_calls

    if diff is not None:
        from allelix.reports.diff import diff_annotation_to_dict, summarize_diff

        payload["diff"] = {
            "previous_report": diff.previous_generated_at,
            "summary": summarize_diff(diff),
            "new": [annotation_to_public_dict(a) for a in diff.new],
            "changed": [diff_annotation_to_dict(c) for c in diff.changed],
            "removed": diff.removed,
        }

    atomic_write_text(output_path, json.dumps(payload, indent=2, sort_keys=False) + "\n")
    return len(filtered)

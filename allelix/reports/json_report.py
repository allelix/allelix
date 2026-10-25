# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""JSON report renderer.

Output schema (versioned via `schema_version`):

    {
      "schema_version": "10",
      "allelix_version": "1.1.0",
      "generated_at": "2026-05-11T12:34:56+00:00",
      "regulatory_notice": "...",
      "input": {
        "file": "genotype.txt",
        "format": "myhappygenes",
        "method": null,
        "sample_id": "MHG000001",
        "build": "GRCh37",
        "build_source": "seeds",
        "build_evidence": {"file_label": null, "seeds_matched": 10, ...},
        "total_variants": 2015,
        "skipped_lines": 0,
        "gvcf": false,
        "vcf_quality": null,
        "array_quality": {"di_calls": 0, "unreadable_genotypes": 0,
                          "not_genotyped": 0, "wrong_column_rows": 0},
        "probe_classes": {"i_probes_at_plp_loci": 0, "di_calls_at_plp_rsids": 0}
      },
      "annotators": [
        {"name": "clinvar", "version": "20260101"}
      ],
      "sources": [
        {"name": "clinvar", "display_name": "ClinVar", "state": "used",
         "reason": null, "version": "20260101"},
        {"name": "cadd", "display_name": "CADD", "state": "skipped",
         "reason": "disabled", "version": null}
      ],
      "report_warnings": [],
      "filters": {
        "min_magnitude": 5.0,
        "category": null,
        "genes": null
      },
      "clinvar": {"pathogenic_rows_scored_down": 0, "unrecognised_rows_hidden": 0,
                  "common_allele_rows_scored_down": 0},
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
from allelix.utils.allele import STRAND_ABSTAIN_REASONS

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


SCHEMA_VERSION = "10"
# v9 -> v10 (2.5.1): ``input.vcf_quality`` adds the VCF genotype counts and the
# reference-source/ploidy conflict counts; each panel variant adds nullable
# ``reason`` for an explicit unusable VCF record (ADR-0055). Arrays emit null;
# annotation shape is unchanged.
# v8 -> v9 bump (one for every 2.5 report field, ADR-0045 first): every
# annotation carries ``clinvar_aggregate`` (ClinVar's own classification of
# the variant), ``submitters_same`` / ``submitters_total`` (the submitter
# split), ``aggregate_capped`` and ``outside_aggregate`` (empty / null /
# false on non-ClinVar rows). The top-level ``clinvar`` block counts the
# Pathogenic / Likely pathogenic rows scored down and the rows hidden for a
# classification Allelix does not score. The top-level ``sources`` list
# gives every registered source's part in the run (``used``,
# ``ready_not_used``, or ``skipped`` with a reason); ``annotators`` still
# lists only the sources that ran. ``report_warnings`` is a list of
# ``{code, message, count}`` (``clinvar_skipped`` when ClinVar did not run).
# Per annotation (ADR-0050): ``zygosity`` is a code (``homozygous``,
# ``heterozygous``, ``hemizygous``, ``half_call``, ``no_call``; it was the
# display string), ``plausibility`` is ``not_assessed`` where it was empty,
# and ``plausibility_reason`` says why. ``report_warnings`` may also carry
# ``multiple_homozygous_plp``, and ``forced_format`` when ``--format`` named a
# reader that does not recognise the file.
# ``panel_coverage`` (ADR-0051) is ``{requested, counts, variants}``: ``counts``
# has ``annotated``, ``typed_no_annotation``, ``no_call`` and
# ``not_on_chip`` (summing to ``requested``), and each ``variants`` entry
# is ``{rsid, gene, state, genotype, annotation_count, reason}``. It replaces
# ``found`` / ``missing`` / ``no_findings``.
# ``input.vcf_quality.clinvar_ref_mismatches`` counts one-base records at a
# ClinVar position whose REF is none of ClinVar's one-base REFs there.
# ``reference_carried_calls`` counts hom-ref calls at a site where the
# build's reference carries a ClinVar variant allele, assessed as homozygous
# for it; ``reference_carried_not_in_file`` the sites of that kind a plain VCF has
# no record at.
# ADR-0051: ``input.array_quality`` counts what the parse of an array file
# saw (``di_calls``, ``unreadable_genotypes``, ``not_genotyped``,
# ``wrong_column_rows``), and ``input.probe_classes`` the array calls
# ClinVar never matches at an aggregate-P/LP site
# (``i_probes_at_plp_loci``, ``di_calls_at_plp_rsids``). Both are null for
# VCF input; ``probe_classes`` is also null when ClinVar did not run.
# GWAS annotations carry ``effect_direction`` (``increases``,
# ``decreases``, ``unknown``; empty on other sources), and the top-level
# ``strand_abstentions`` gives, per strand-checked source that ran
# (``gwas``, ``snpedia``), the counts not assessed by reason
# (``not_in_identity_table``, ``palindromic``, ``no_fit``): GWAS
# association rows, SNPedia sites.
# v7 -> v8 bump (ADR-0033): ``input.build_source`` says where the build came
# from (ADR-0044 vocabulary: override, contig_lengths, seeds, assembly_header,
# file_label, default) and ``input.build_evidence`` records the evidence
# beside it (the file's own label, the contig-length verdict, the known-SNP
# counts). Both are null when the result carries no build diagnostics.
# v6 -> v7 bump (ADR-0033): every annotation carries ``plausibility``, the
# zygosity x gnomAD frequency quality-control flag (ADR-0041); the
# ``input`` block carries ``gvcf`` and ``vcf_quality`` (null for non-VCF
# input); ``high_value_no_calls`` entries carry ``reason``; ``input.method``
# names a declared genotyping method (null when none). All land in the
# same unreleased version, so one bump covers them.
# v4 → v5 bump: enters the ADR-0035 cluster contract. The cluster fields
# (variant.ref via Variant model, annotation.trait / .p_value / .phecode via
# Annotation model) all arrived under v5 in the v2.1 series. v4
# consumers reading v5 output still succeed — every v4 field is unchanged.


def _vcf_quality_dict(result: AnalysisResult) -> dict | None:
    """What the VCF quality gate did; None for non-VCF input."""
    stats = result.vcf_quality
    if stats is None:
        return None
    return {
        "filtered_rows": stats.filtered_rows,
        "filtered_by_label": dict(sorted(stats.filtered_by_label.items())),
        "kept_filtered": stats.kept_filtered,
        "low_quality_calls": stats.low_quality_calls,
        "min_gq": stats.min_gq,
        "spanning_deletion_calls": stats.spanning_deletion_calls,
        "spanning_deletion_sites": list(stats.spanning_deletion_sites),
        "symbolic_rows": stats.symbolic_rows,
        "half_calls": stats.half_calls,
        "malformed_rows": stats.malformed_rows,
        "no_genotype_rows": stats.no_genotype_rows,
        "unusable_genotype_calls": stats.unusable_genotype_calls,
        "symbolic_allele_calls": stats.symbolic_allele_calls,
        "trimmed_matches": stats.trimmed_matches,
        "clinvar_ref_mismatches": stats.clinvar_ref_mismatches,
        "reference_source_conflicts": stats.reference_source_conflicts,
        "reference_ploidy_conflicts": stats.reference_ploidy_conflicts,
        "reference_carried_calls": stats.reference_carried_calls,
        "reference_carried_not_in_file": stats.reference_carried_not_in_file,
        "genotype_records": stats.genotype_records,
        "hom_ref_calls": stats.hom_ref_calls,
        "non_reference_calls": stats.non_reference_calls,
    }


def _array_quality_dict(result: AnalysisResult) -> dict | None:
    """What the parse of an array file saw (ADR-0051); None for VCF input."""
    stats = result.array_quality
    if stats is None:
        return None
    return {
        "di_calls": stats.di_calls,
        "unreadable_genotypes": stats.unreadable_genotypes,
        "not_genotyped": stats.not_genotyped,
        "wrong_column_rows": stats.wrong_column_rows,
    }


def _annotation_dict(a: Annotation) -> dict:
    """Serialize an annotation, adding AM caveat for non-protein sources."""
    d = annotation_to_public_dict(a)
    d["zygosity"] = a.zygosity_code
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
    diag = result.build_diagnostics
    payload: dict = {
        "schema_version": SCHEMA_VERSION,
        "allelix_version": __version__,
        "generated_at": datetime.now(UTC).isoformat(),
        "regulatory_notice": REGULATORY_NOTICE,
        "input": {
            "file": result.file_path.name,
            "format": result.parser_name,
            "method": result.input_method or None,
            "sample_id": result.sample_id,
            "build": result.build,
            "build_source": diag.source if diag is not None else None,
            "build_evidence": diag.evidence() if diag is not None else None,
            "total_variants": result.total_variants,
            "skipped_lines": result.skipped_count,
            "gvcf": result.vcf_quality is not None and result.vcf_quality.is_gvcf,
            "vcf_quality": _vcf_quality_dict(result),
            "array_quality": _array_quality_dict(result),
            "probe_classes": dict(result.probe_classes)
            if result.probe_classes is not None
            else None,
        },
        "annotators": [
            {"name": name, "version": version} for name, version in result.annotators_used
        ],
        "sources": [record.to_dict() for record in result.sources],
        "report_warnings": [warning.to_dict() for warning in result.report_warnings],
        "filters": {
            "min_magnitude": min_magnitude,
            "category": category,
            "genes": sorted(genes) if genes is not None else None,
            "rsids": sorted(rsids) if rsids is not None else None,
        },
        "clinvar": {
            "pathogenic_rows_scored_down": result.clinvar_scored_down,
            "unrecognised_rows_hidden": result.clinvar_unrecognised,
            "common_allele_rows_scored_down": result.clinvar_common_allele,
        },
        "strand_abstentions": {
            source: {reason: counts.get(reason, 0) for reason in STRAND_ABSTAIN_REASONS}
            for source, counts in sorted(result.strand_abstentions.items())
        },
        "annotations": [_annotation_dict(a) for a in filtered],
    }

    license_attrs = _license_attributions(result.annotators_used)
    if license_attrs:
        payload["license_attributions"] = license_attrs

    if high_value_no_calls:
        payload["high_value_no_calls"] = high_value_no_calls

    # When --filter-file supplied a panel, every requested rsID is in one
    # of four states (ADR-0051), decided from the shown (post-filter)
    # annotations; ``variants`` is sorted by state, then rsID.
    panel_coverage = result.panel_coverage(filtered)
    if panel_coverage is not None:
        payload["panel_coverage"] = panel_coverage

    if diff is not None:
        from allelix.reports.diff import diff_annotation_to_dict, summarize_diff

        payload["diff"] = {
            "previous_report": diff.previous_generated_at,
            "summary": summarize_diff(diff),
            "new": [annotation_to_public_dict(a) for a in diff.new],
            "changed": [diff_annotation_to_dict(c) for c in diff.changed],
            "removed": diff.removed,
        }

    atomic_write_text(
        output_path,
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        input_path=result.file_path,
    )
    return len(filtered)

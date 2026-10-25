# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Annotated VCF output writer (ADR-0036).

Two-pass architecture: the pipeline runs pass 1 (annotation) via the
existing analyze path and produces an annotation dict keyed by
``(chrom, pos, ref, alt)``. This writer performs pass 2: re-reads the
raw input file line by line, injects provenance / attribution /
``##INFO`` declarations into the header, and stamps INFO fields onto
each data line using the per-allele ``Number=A`` semantics required by
the VCF spec for multi-allelic sites.

The writer does not modify the parser and does not touch
``AnalysisResult``. Everything the writer needs is passed in
explicitly: source path, ``annotators_used`` (for INFO / provenance /
attribution declarations), the annotations dict, and an optional
recovered-rsID map. ``run_analysis`` in ``_pipeline.py`` produces
those inputs when called with ``collect_vcf_write_inputs=True``.

Chromosome-normalization boundary: the pipeline keys
``annotations_by_key`` on ``Variant.chromosome`` which the parser has
already run through ``normalize_chromosome`` (``chr22`` → ``22``,
``chrM`` → ``MT``, standard-name uppercase). The writer must apply
the same normalization when hashing the raw input line back into the
dict; the output CHROM column preserves the file's original convention
unchanged. Real GRCh38 VCFs from GATK / DeepVariant / DRAGEN routinely
use ``chr``-prefixed contigs, so any drift on this boundary silently
zeroes annotation stamping on the primary user's file.

The ``##ALLELIX_Attribution=<>`` block mirrors ``license_attributions``
from the HTML and JSON renderers — three output formats, three matching
lists. Rule-of-Three-plus-one: extraction into a shared helper waits
for a fourth call site.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.annotators import get_annotator_class
from allelix.models import report_rank_key
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.vcf import _is_symbolic, _open_vcf, reference_row_key

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.parsers.vcf import Locus, RowKey

__all__ = ["VCF_LICENSE", "VCF_SCHEMA_VERSION", "VcfWriteInputs", "render_vcf"]


VCF_SCHEMA_VERSION = "0.5.0"
"""Semantic version of the annotated-VCF output contract.

Experimental — INFO field names and structure may change before v3.0.
Bump on any incompatible change to a declared INFO field, provenance
line, or attribution shape. ADR after real user feedback.

0.2.0: adds ``ALLELIX_FILTER`` (rows kept with ``--keep-filtered``), and
caller-rejected rows are no longer stamped by default.
0.3.0: ``##ALLELIX_Build=<ID=...,Source=...>`` records where the build came
from (ADR-0044); it was ``##ALLELIX_Build=<build>``.
0.4.0: ``##ALLELIX_SkippedSources=`` lists each registered source that did
not run, with its reason (``clinvar:outdated,cadd:disabled``), or ``none``,
so an annotated VCF never reads as complete when a source was skipped.
0.5.0: ``ALLELIX_CLINVAR`` also appears on an explicit ``0/0`` or haploid
``0`` record, where its value is a finding for the allele the reference
itself carries (ADR-0059). Re-annotation replaces earlier Allelix INFO.
"""

VCF_LICENSE = "AGPL-3.0-or-later"

_VariantKey = tuple[str, int, str, str]


@dataclass(frozen=True)
class VcfWriteInputs:
    """Everything the writer needs from the pipeline for pass 2.

    Populated by ``run_analysis`` when ``collect_vcf_write_inputs`` is
    True and attached to ``AnalysisResult.vcf_write_inputs``. Kept as a
    dedicated dataclass rather than two loose fields on ``AnalysisResult``
    so the writer's inputs stay together as one intentional surface.

    ``annotations_by_key`` values are references to the ``Annotation``
    objects on ``AnalysisResult.annotations`` — enrichment mutations
    (gnomAD allele frequency, AlphaMissense score, CADD phred) that
    happen after the streaming phase are visible through these refs
    without any per-annotation copy.
    """

    annotations_by_key: dict[_VariantKey, list[Annotation]]
    recovered_rsids: dict[_VariantKey, str]
    # Explicit reference records the analysis used, with the positions each
    # covers (``VcfQualityStats.reference_rows``).
    reference_rows: dict[RowKey, set[Locus]] = field(default_factory=dict)


_logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _InfoSpec:
    """Metadata for one Allelix-emitted INFO field.

    ``Number=A`` semantics apply to every field here — one value per
    ALT allele, comma-joined in ALT-column order.
    """

    field: str
    vcf_type: str
    description: str


# One entry per Allelix-emitted INFO field. Keys are internal source
# tags, not ``Annotation.source`` values directly.
#
# **Primary sources** (``_PRIMARY_SOURCES``: clinvar, gwas, pharmgkb,
# snpedia) map 1:1 to ``Annotation.source`` — one Annotation object per
# primary hit, formatted via ``_format_primary``.
#
# **Enrichment sources** (``_ENRICHMENT_SOURCES``: gnomad, alphamissense,
# cadd) never create Annotation objects; the pipeline's enrichment phase
# mutates FIELDS on existing primary Annotations (``allele_frequency``,
# ``am_pathogenicity``, ``am_class``, ``cadd_phred``). Their INFO fields
# are extracted via ``_extract_enrichment`` from any annotation at the
# variant key that carries a non-None enrichment value.
#
# Insertion order is the header emission order — alphabetical for
# stable diffs across runs.
_INFO_SPECS: dict[str, _InfoSpec] = {
    "alphamissense": _InfoSpec(
        field="ALLELIX_ALPHAMISSENSE",
        vcf_type="String",
        description="AlphaMissense class | pathogenicity score",
    ),
    "cadd": _InfoSpec(
        field="ALLELIX_CADD",
        vcf_type="Float",
        description="CADD PHRED score",
    ),
    "clinvar": _InfoSpec(
        field="ALLELIX_CLINVAR",
        vcf_type="String",
        description=(
            "ClinVar significance | rsID | gene. On a homozygous-reference or "
            "haploid reference record, the finding is for the allele the "
            "reference itself carries."
        ),
    ),
    "gnomad": _InfoSpec(
        field="ALLELIX_GNOMAD_AF",
        vcf_type="Float",
        description="gnomAD allele frequency",
    ),
    "gwas": _InfoSpec(
        field="ALLELIX_GWAS",
        vcf_type="String",
        description="GWAS trait | p_value | study",
    ),
    "pharmgkb": _InfoSpec(
        field="ALLELIX_CLINPGX",
        vcf_type="String",
        description="ClinPGx gene | significance",
    ),
    "snpedia": _InfoSpec(
        field="ALLELIX_SNPEDIA",
        vcf_type="String",
        description="SNPedia rsID | significance",
    ),
}

_PRIMARY_SOURCES: frozenset[str] = frozenset({"clinvar", "gwas", "pharmgkb", "snpedia"})
_ENRICHMENT_SOURCES: frozenset[str] = frozenset({"alphamissense", "cadd", "gnomad"})

_ORIGINAL_ID_FIELD = "ALLELIX_ORIGINAL_ID"
# the input row's FILTER, repeated on a stamped row that the caller
# rejected and the user kept with --keep-filtered.
_FILTER_FIELD = "ALLELIX_FILTER"
_UNFILTERED = frozenset({"PASS", ".", ""})
_CONTIG_DECLARATION_PREFIX = "##contig=<ID="

# Chars reserved by the VCF INFO grammar: ``;`` separates fields, ``,``
# separates per-allele values, ``=`` separates key from value, and
# whitespace is either a column boundary (tab) or forbidden inside a
# value (space). Newlines and carriage returns are record-structural —
# an embedded ``\n`` in an INFO value would split the data line into
# two, corrupting the file. Substituting underscore is the safe minimum
# for a v0.1 writer; full percent-encoding per VCF 4.3 §1.6 (including
# the escape char ``%`` itself) is a follow-up when we see a real user
# need for round-tripping arbitrary text.
_INFO_SUBSTITUTIONS = str.maketrans(
    {
        ";": "_",
        ",": "_",
        "=": "_",
        " ": "_",
        "\t": "_",
        "\n": "_",
        "\r": "_",
    }
)


def _sanitize_info_value(value: str) -> str:
    """Replace VCF-reserved chars in an INFO value with underscore."""
    return value.translate(_INFO_SUBSTITUTIONS)


def _format_number(value: float) -> str:
    """Format a float for VCF using ``%g``.

    Trims trailing zeros while preserving precision, matching how
    bcftools renders AF/CADD.
    """
    return f"{value:g}"


def _format_primary(ann: Annotation) -> tuple[str, str] | None:
    """Format a primary-source Annotation as (source_key, INFO_value).

    Returns None when ``ann.source`` is not a primary source known to
    this writer (enrichment sources go through ``_extract_enrichment``
    against the whole per-key annotation list, not source-dispatched).
    """
    source = ann.source.lower()
    if source == "clinvar":
        return (
            "clinvar",
            _sanitize_info_value(f"{ann.significance}|{ann.rsid or '.'}|{ann.gene or '.'}"),
        )
    if source == "gwas":
        pval = "." if ann.p_value is None else _format_number(ann.p_value)
        return (
            "gwas",
            _sanitize_info_value(f"{ann.trait or '.'}|{pval}|{ann.attribution}"),
        )
    if source == "pharmgkb":
        return (
            "pharmgkb",
            _sanitize_info_value(f"{ann.gene or '.'}|{ann.significance}"),
        )
    if source == "snpedia":
        return (
            "snpedia",
            _sanitize_info_value(f"{ann.rsid}|{ann.significance}"),
        )
    return None


def _extract_enrichment(anns: Iterable[Annotation]) -> dict[str, str]:
    """Extract enrichment INFO values from any annotation at one variant.

    Enrichment (gnomAD AF, AlphaMissense, CADD) is a property of the
    variant, not the primary source — the pipeline mutates
    ``allele_frequency`` / ``am_pathogenicity`` / ``am_class`` /
    ``cadd_phred`` on the existing primary Annotation objects. All
    annotations at the same ``(chrom, pos, ref, alt)`` share the same
    enrichment values, so this takes the first non-None per field.

    Returns ``{source_key: formatted_value}`` for populated fields;
    omitted keys mean the enrichment source had nothing to say at this
    variant.
    """
    result: dict[str, str] = {}
    af = next((a.allele_frequency for a in anns if a.allele_frequency is not None), None)
    if af is not None:
        result["gnomad"] = _format_number(af)
    am_score = next((a.am_pathogenicity for a in anns if a.am_pathogenicity is not None), None)
    am_cls = next((a.am_class for a in anns if a.am_class), None)
    if am_score is not None or am_cls:
        score_str = "." if am_score is None else _format_number(am_score)
        result["alphamissense"] = _sanitize_info_value(f"{am_cls or '.'}|{score_str}")
    cadd = next((a.cadd_phred for a in anns if a.cadd_phred is not None), None)
    if cadd is not None:
        result["cadd"] = _format_number(cadd)
    return result


def _build_info_field(per_source: dict[str, list[str]]) -> str:
    """Serialize per-allele INFO values into a single INFO-field string.

    ``per_source`` maps ``source_name → [value_per_ALT]``. Sources with
    no non-empty ALT slot are dropped; empty slots in a retained source
    render as ``.`` so ALT-column alignment holds. Emission order
    follows ``_INFO_SPECS`` insertion (alphabetical) for deterministic
    diffs.
    """
    fragments: list[str] = []
    for source in _INFO_SPECS:
        if source not in per_source:
            continue
        values = per_source[source]
        if not any(values):
            continue
        spec = _INFO_SPECS[source]
        joined = ",".join(v if v else "." for v in values)
        fragments.append(f"{spec.field}={joined}")
    return ";".join(fragments)


def _license_attributions(
    annotators_used: Iterable[tuple[str, str | None]],
) -> list[dict[str, str]]:
    """Attribution rows for licensed sources.

    VCF-shaped mirror of the ``_license_attributions`` helpers in
    ``html.py`` and ``json_report.py``. Three output formats, three
    matching lists.
    """
    rows: list[dict[str, str]] = []
    for name, _version in annotators_used:
        cls = get_annotator_class(name)
        if cls is None:
            _logger.warning("No annotator class found for '%s' — attribution omitted", name)
            continue
        desc = cls.license
        rows.append(
            {
                "name": cls.display_name,
                "spdx": desc.spdx,
                "url": desc.source_url or desc.license_url,
            }
        )
    return rows


def _emit_provenance_lines(
    *,
    build: str,
    build_source: str | None,
    annotators_used: list[tuple[str, str | None]],
    run_date: datetime,
    skipped_sources: list[tuple[str, str]] | None = None,
) -> list[str]:
    """Emit the provenance header block.

    Version, schema, license, run date, build, per-database version
    pins, and per-source attribution — everything a downstream tool
    needs to reproduce the annotation run from the file alone.
    """
    lines = [
        f"##ALLELIX_Version={__version__}",
        f"##ALLELIX_VCF_SchemaVersion={VCF_SCHEMA_VERSION}",
        f"##ALLELIX_License={VCF_LICENSE}",
        f"##ALLELIX_RunDate={run_date.isoformat()}",
        f"##ALLELIX_Build=<ID={build},Source={build_source or 'unknown'}>",
        "##ALLELIX_SkippedSources="
        + (",".join(f"{name}:{reason}" for name, reason in skipped_sources or []) or "none"),
    ]
    for name, version in annotators_used:
        version_str = version if version else "unknown"
        lines.append(f"##ALLELIX_Database=<Name={name},Version={version_str}>")
    for row in _license_attributions(annotators_used):
        lines.append(
            f"##ALLELIX_Attribution=<Name={row['name']},License={row['spdx']},URL={row['url']}>"
        )
    return lines


def _emit_info_declarations(
    sources_used: set[str], *, include_original_id: bool, keep_filtered: bool = False
) -> list[str]:
    """Emit ``##INFO=<>`` declarations for the sources used in this run.

    Sources not touched are omitted so the header only declares fields
    that actually appear in the file. ``ALLELIX_ORIGINAL_ID`` is
    declared once when the caller passed a non-empty ``recovered_rsids``
    map. ``ALLELIX_FILTER`` is declared when filtered rows were kept.
    """
    lines: list[str] = []
    for source, spec in _INFO_SPECS.items():
        if source in sources_used:
            lines.append(
                f"##INFO=<ID={spec.field},Number=A,Type={spec.vcf_type},"
                f'Description="{spec.description}">'
            )
    if include_original_id:
        lines.append(
            f"##INFO=<ID={_ORIGINAL_ID_FIELD},Number=1,Type=String,"
            'Description="Original VCF ID column value before Allelix rsID recovery">'
        )
    if keep_filtered:
        lines.append(
            f"##INFO=<ID={_FILTER_FIELD},Number=1,Type=String,"
            'Description="FILTER of a row the variant caller rejected, annotated by '
            'Allelix only because --keep-filtered was set">'
        )
    return lines


def _declared_contig_id(raw_line: str) -> str | None:
    """Return the raw contig ID from a ``##contig`` declaration."""
    if not raw_line.startswith(_CONTIG_DECLARATION_PREFIX):
        return None
    remainder = raw_line[len(_CONTIG_DECLARATION_PREFIX) :]
    contig_id = remainder.split(",", 1)[0].split(">", 1)[0]
    return contig_id or None


def _missing_contig_ids(input_path: Path) -> list[str]:
    """Return emitted contig IDs absent from the input header.

    IDs are compared and returned exactly as written in the VCF. The
    parser-normalized chromosome name is appropriate for annotation lookup,
    but a synthesized declaration must match the raw CHROM column consumed
    by downstream VCF tools. Missing IDs retain first-data-line order for
    deterministic output.
    """
    declared: set[str] = set()
    emitted: list[str] = []
    seen_emitted: set[str] = set()
    body_started = False

    with _open_vcf(input_path) as fin:
        for raw_line in fin:
            if not body_started:
                declared_id = _declared_contig_id(raw_line)
                if declared_id is not None:
                    declared.add(declared_id)
                if raw_line.startswith("#CHROM"):
                    body_started = True
                continue
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            contig_id = raw_line.split("\t", 1)[0]
            if contig_id and contig_id not in seen_emitted:
                emitted.append(contig_id)
                seen_emitted.add(contig_id)

    return [contig_id for contig_id in emitted if contig_id not in declared]


def _sources_used(
    annotations_by_key: Mapping[_VariantKey, list[Annotation]],
) -> set[str]:
    """Collect the internal source keys that actually appear in the run.

    Primary sources (``clinvar``, ``gwas``, ``pharmgkb``, ``snpedia``)
    are added when at least one Annotation object's ``source`` matches.
    Enrichment sources (``gnomad``, ``alphamissense``, ``cadd``) are
    added when at least one Annotation carries the corresponding
    populated field (``allele_frequency``, ``am_pathogenicity`` /
    ``am_class``, ``cadd_phred``) — the pipeline mutates these fields
    on primary annotations rather than emitting new source-tagged
    Annotation objects, so source-string dispatch alone would miss
    them entirely. Unknown primary sources are silently skipped so a
    future annotator can land without breaking the writer.
    """
    used: set[str] = set()
    for anns in annotations_by_key.values():
        for a in anns:
            src = a.source.lower()
            if src in _PRIMARY_SOURCES:
                used.add(src)
            if a.allele_frequency is not None:
                used.add("gnomad")
            if a.am_pathogenicity is not None or a.am_class:
                used.add("alphamissense")
            if a.cadd_phred is not None:
                used.add("cadd")
    return used


def render_vcf(
    *,
    input_path: Path,
    output_path: Path,
    build: str,
    annotators_used: list[tuple[str, str | None]],
    annotations_by_key: Mapping[_VariantKey, list[Annotation]],
    recovered_rsids: Mapping[_VariantKey, str] | None = None,
    run_date: datetime | None = None,
    keep_filtered: bool = False,
    build_source: str | None = None,
    skipped_sources: list[tuple[str, str]] | None = None,
    reference_rows: Mapping[RowKey, set[Locus]] | None = None,
) -> int:
    """Write an annotated VCF from ``input_path`` to ``output_path``.

    Pass 2 of the two-pass architecture. ``annotations_by_key`` and
    ``recovered_rsids`` are produced by ``run_analysis`` when called
    with ``collect_vcf_write_inputs=True``; the writer treats them as
    opaque lookups keyed by ``(chrom, pos, ref, alt)`` where ``chrom``
    is the parser-normalized bare form (``22``, ``X``, ``MT``).

    Multi-allelic sites use ``Number=A`` semantics — each Allelix
    INFO field holds one comma-joined value per ALT in column order,
    empty slots render as ``.``. Symbolic alleles (``<NON_REF>`` etc.)
    skip annotation lookup and preserve their column position.

    Header behaviour: every input ``##`` and blank line passes through
    verbatim, except the ``##INFO=<ID=ALLELIX_...`` declarations and
    ``##ALLELIX_...`` provenance lines an earlier run wrote; those are
    replaced by the current run's. ``ALLELIX_*`` INFO tokens on data rows
    are removed the same way before the current result is stamped.
    If the input omitted declarations for contigs present in its
    own data rows, minimal ``##contig=<ID=...>`` declarations are synthesized
    using the raw CHROM spelling. Those lines, injected ``##INFO``
    declarations, and the provenance / attribution block are inserted
    immediately before ``#CHROM``. Sample columns and FORMAT fields pass
    through unchanged.

    A row the caller rejected (FILTER not ``PASS`` / ``.``) is
    never freshly stamped, unless ``keep_filtered`` is set;
    then it is stamped like any row and also carries ``ALLELIX_FILTER``.
    Old ``ALLELIX_*`` tokens are removed from rejected rows and reference
    blocks too, but those rows get no new findings. An explicit reference
    call (``0/0`` or haploid ``0``) gets the pipeline's reference-carried
    ClinVar finding in every non-symbolic ALT slot. Only the records in
    ``reference_rows`` (the ones the analysis read that call from) get it.

    Returns the number of data lines stamped with at least one Allelix
    INFO field (either an annotation source or ``ALLELIX_ORIGINAL_ID``).
    """
    from allelix.reports import staged_output, validate_output_path

    validate_output_path(input_path, output_path, legacy_temp=False)
    if run_date is None:
        run_date = datetime.now(UTC)
    if recovered_rsids is None:
        recovered_rsids = {}

    sources_used = _sources_used(annotations_by_key)
    include_original_id = bool(recovered_rsids)
    missing_contig_ids = _missing_contig_ids(input_path)

    stamped = 0
    header_written = False
    reference_rows = reference_rows or {}

    with staged_output(output_path) as staged:
        with (
            _open_vcf(input_path) as fin,
            staged.open_text(encoding="utf-8", newline="\n") as fout,
        ):
            for raw_line in fin:
                if not header_written and raw_line.startswith("##"):
                    # Drop what an earlier Allelix run wrote; it is re-emitted fresh.
                    if not raw_line.startswith(("##INFO=<ID=ALLELIX_", "##ALLELIX_")):
                        fout.write(raw_line)
                    continue
                if not header_written and raw_line.startswith("#CHROM"):
                    for contig_id in missing_contig_ids:
                        fout.write(f"##contig=<ID={contig_id}>\n")
                    for h in _emit_info_declarations(
                        sources_used,
                        include_original_id=include_original_id,
                        keep_filtered=keep_filtered,
                    ):
                        fout.write(h + "\n")
                    for h in _emit_provenance_lines(
                        build=build,
                        build_source=build_source,
                        annotators_used=annotators_used,
                        run_date=run_date,
                        skipped_sources=skipped_sources,
                    ):
                        fout.write(h + "\n")
                    fout.write(raw_line)
                    header_written = True
                    continue
                if not raw_line.strip():
                    fout.write(raw_line)
                    continue
                if not header_written:
                    # Data line before #CHROM — malformed input, but the parser
                    # tolerates it and so do we. Pass through untouched.
                    fout.write(raw_line)
                    continue

                new_line, row_stamped = _stamp_data_line(
                    raw_line,
                    sources_used=sources_used,
                    annotations_by_key=annotations_by_key,
                    recovered_rsids=recovered_rsids,
                    keep_filtered=keep_filtered,
                    reference_rows=reference_rows,
                )
                fout.write(new_line)
                if row_stamped:
                    stamped += 1

        staged.publish()

    return stamped


def _stamp_data_line(
    raw_line: str,
    *,
    sources_used: set[str],
    annotations_by_key: Mapping[_VariantKey, list[Annotation]],
    recovered_rsids: Mapping[_VariantKey, str],
    keep_filtered: bool = False,
    reference_rows: Mapping[RowKey, set[Locus]] | None = None,
) -> tuple[str, bool]:
    """Rewrite one VCF data line with Allelix INFO fields.

    Returns ``(new_line_with_trailing_newline, was_stamped)``. Malformed
    rows (< 8 columns or non-integer POS) pass through unchanged and
    are reported as unstamped rather than raising — the writer never
    fails a whole file on one bad row.
    """
    line = raw_line.rstrip("\n").rstrip("\r")
    cols = line.split("\t")
    if len(cols) < 8:
        return raw_line, False
    chrom, pos_s, id_col, ref, alt_col, qual, filt, info = cols[:8]
    tail = cols[8:]
    try:
        pos = int(pos_s)
    except ValueError:
        return raw_line, False
    # Drop earlier Allelix tokens; other tokens keep their exact text.
    info = (
        ";".join(
            token for token in info.split(";") if not token.split("=", 1)[0].startswith("ALLELIX_")
        )
        or "."
    )
    cols[7] = info
    cleaned_line = "\t".join(cols) + "\n"
    caller_rejected = filt.strip() not in _UNFILTERED
    if caller_rejected and not keep_filtered:
        # Never put Allelix INFO on a row the caller rejected: other tools
        # would read it as an Allelix-endorsed call.
        return cleaned_line, False

    alts = [] if alt_col in ("", ".") else alt_col.split(",")
    per_source: dict[str, list[str]] = {src: ["" for _ in alts] for src in sources_used}
    row_stamped = False
    original_id: str | None = None
    recovered: str | None = None

    # The pipeline builds ``annotations_by_key`` from ``Variant.chromosome``
    # which the parser has already run through ``normalize_chromosome``
    # (``chr22`` → ``22``, ``chrM`` → ``MT``, standard-name uppercase).
    # The raw input line still carries the file's original convention
    # (real-world GRCh38 VCFs from GATK / DeepVariant / DRAGEN routinely
    # use ``chr``-prefixed contigs), so key lookups must use the
    # normalized form or every annotation on a ``chr``-prefixed file
    # silently misses. Only the lookup key is normalized — the output
    # column preserves the user's convention unchanged.
    key_chrom = normalize_chromosome(chrom)
    # A reference-carried ClinVar finding goes on the records it was read
    # from, at the base of REF that covers its position.
    reference_value = ""
    candidates = []
    for locus in (reference_rows or {}).get(reference_row_key(cols), ()):
        base = ref[locus[1] - pos].upper()
        for annotation in annotations_by_key.get((*locus, base, base), []):
            if (
                annotation.source == "clinvar"
                and annotation.reference_carries_alt
                and annotation.locus == locus
                and annotation.locus_ref == base
            ):
                candidates.append(annotation)
    for annotation in sorted(candidates, key=report_rank_key):
        primary = _format_primary(annotation)
        if primary is not None:
            reference_value = primary[1]
            break

    for i, alt in enumerate(alts):
        if _is_symbolic(alt):
            continue
        # The REF exception is ClinVar-only. REF enrichment values must
        # never be presented as properties of one of the listed ALTs.
        if reference_value:
            per_source["clinvar"][i] = reference_value
            row_stamped = True
        key: _VariantKey = (key_chrom, pos, ref.upper(), alt.upper())
        anns = annotations_by_key.get(key, [])
        # Primary sources: source-dispatched, one Annotation → one INFO
        # value slot. Where a source has several rows at the variant (two
        # ClinVar conditions, say), the slot takes the row the report
        # lists first, so the VCF never shows a different classification
        # from the report.
        for a in sorted(anns, key=report_rank_key):
            primary = _format_primary(a)
            if primary is None:
                continue
            src, value = primary
            if src not in per_source or per_source[src][i]:
                continue
            per_source[src][i] = value
            row_stamped = True
        # Enrichment sources: single per-variant lookup, not per-Annotation.
        # gnomAD AF / AlphaMissense / CADD are mutated onto primary
        # Annotation objects by the pipeline's enrichment phase, so they
        # aren't source-dispatched — extract once from any annotation at
        # this key that carries the field populated.
        if anns:
            enrichment = _extract_enrichment(anns)
            for src, value in enrichment.items():
                if src in per_source:
                    per_source[src][i] = value
                    row_stamped = True
        # First recovered rsID at this row wins. In practice the pipeline
        # emits at most one recovered rsID per position, but the loop is
        # defensive against future callers that key per-ALT.
        if recovered is None and key in recovered_rsids:
            recovered = recovered_rsids[key]
            original_id = id_col

    allelix_info = _build_info_field(per_source)
    if recovered is not None and original_id is not None:
        orig_frag = f"{_ORIGINAL_ID_FIELD}={_sanitize_info_value(original_id)}"
        allelix_info = f"{allelix_info};{orig_frag}" if allelix_info else orig_frag
    if allelix_info and caller_rejected:
        allelix_info += f";{_FILTER_FIELD}={_sanitize_info_value(filt.strip())}"

    if allelix_info:
        info = allelix_info if info in ("", ".") else f"{info};{allelix_info}"

    new_id = recovered if recovered is not None else id_col
    new_cols = [chrom, pos_s, new_id, ref, alt_col, qual, filt, info, *tail]
    return "\t".join(new_cols) + "\n", row_stamped or recovered is not None

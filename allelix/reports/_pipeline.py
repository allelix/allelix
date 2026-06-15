# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared analysis pipeline used by `analyze`, `methylation`, and `pharmacogenomics`.

The CLI builds an `AnalysisResult` once and hands it to a renderer
(terminal, JSON, HTML). Renderers never query the database or re-iterate
the parser — they receive a fully-populated value object.

ADR-0021: this pipeline owns build detection. Parsers report the
header-claimed build; the pipeline replaces each variant's `build`
with the build detected from position data (or the user's `--build`
override) before annotators see the variant.
"""

from __future__ import annotations

import contextlib
import logging
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    KNOWN_SNP_POSITIONS,
    detect_build,
    normalize_build_label,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.models import Annotation, Variant
    from allelix.parsers.base import GenotypeParser


# How many input variants to buffer while waiting for detection to
# converge. Detection completes once every entry in KNOWN_SNP_POSITIONS
# has been seen; typical files cover the table within the first ~5000
# probes. Cap so a file with no known SNPs doesn't buffer the whole
# input.
_DETECTION_BUFFER_LIMIT = 100_000

# How many variants to accumulate before flushing to each annotator's
# batch_annotate. Each annotator further chunks internally for SQL
# parameter limits, so this is a Python-side batch boundary, not an
# SQL limit. 5000 keeps memory bounded at ~1.5 MB per pipeline buffer
# while letting the annotators' IN clauses run in batches of several
# hundred. Tuneable; the pipeline is correct at any size including 1.
_BATCH_SIZE = 5000


@dataclass
class BuildDiagnostics:
    """What the pipeline learned about the file's genome build.

    `header_build` is the build claimed by the file header (normalized
    to GRCh37/GRCh38 via `normalize_build_label`; may be None if the
    header doesn't say or uses an unrecognized label).

    `detected_build` is what position data says (None if no known SNPs
    appeared in the input).

    `effective_build` is what was actually used for annotation — either
    a CLI `--build` override, the detected build, or a fallback. Always
    set when the pipeline ran on any data.

    `mismatch` is True when header_build and detected_build disagree
    AND no override was supplied. The CLI surfaces this as a warning.
    """

    header_build: str | None
    detected_build: str | None
    effective_build: str
    override: bool
    matched_count: int
    inspected_count: int

    @property
    def mismatch(self) -> bool:
        return (
            not self.override
            and self.header_build is not None
            and self.detected_build is not None
            and self.header_build != self.detected_build
        )


@dataclass
class AnalysisResult:
    """Everything a renderer needs to produce a report."""

    file_path: Path
    parser_name: str
    parser_display_name: str
    sample_id: str
    build: str
    total_variants: int
    skipped_count: int
    annotators_used: list[tuple[str, str | None]]
    annotations: list[Annotation] = field(default_factory=list)
    build_diagnostics: BuildDiagnostics | None = None
    hv_variants: list[Variant] = field(default_factory=list)

    def filter(
        self,
        *,
        min_magnitude: float = 0.0,
        category: str | None = None,
        genes: Iterable[str] | None = None,
        rsids: Iterable[str] | None = None,
        source_min_magnitudes: dict[str, float] | None = None,
    ) -> list[Annotation]:
        """Apply the standard filters and return a sorted list of annotations.

        Filters are independent and combine with AND. Sort is by magnitude
        descending, then rsid ascending (stable, deterministic).

        `source_min_magnitudes` overrides the floor for specific sources
        (e.g. ``{"gwas": 9.0, "snpedia": 2.0}``). When a source has an
        entry, that value IS the floor for that source — it can raise OR
        lower the global ``min_magnitude``. Sources without an entry use
        the global floor.

        `genes` and `rsids` combine with OR: when either is provided, an
        annotation passes if it matches the gene set OR the rsid set.
        Empty collections (vs None) mean "match nothing" — an empty
        filter file produces an empty report.
        """
        gene_set = {g.upper() for g in genes} if genes is not None else None
        rsid_set = {r.lower() for r in rsids} if rsids is not None else None
        out: list[Annotation] = []
        for a in self.annotations:
            if (
                source_min_magnitudes
                and a.source in source_min_magnitudes
                and not a.is_must_include
            ):
                floor = source_min_magnitudes[a.source]
            else:
                floor = min_magnitude
            if a.magnitude < floor:
                continue
            if category is not None and a.category != category:
                continue
            if gene_set is not None or rsid_set is not None:
                gene_match = gene_set is not None and (a.gene or "").upper() in gene_set
                rsid_match = rsid_set is not None and a.rsid.lower() in rsid_set
                if not gene_match and not rsid_match:
                    continue
            out.append(a)
        out.sort(key=lambda a: (-a.magnitude, a.rsid))
        return out


def _gwas_base_trait(description: str) -> str | None:
    """Extract trait text from a GWAS description, stripping MTAG suffix and PheCode label."""
    marker = "GWAS Catalog: "
    if marker not in description:
        return None
    s = description.split(marker, 1)[1]
    s = s.split(" (p=", 1)[0]
    if s.endswith(" (MTAG)"):
        s = s[: -len(" (MTAG)")]
    s = s.split(" (PheCode ", 1)[0]
    return s.strip().lower()


def _gwas_phecode_parent(description: str) -> str | None:
    """Extract PheCode parent (numeric prefix before the dot), or None."""
    idx = description.find("(PheCode ")
    if idx == -1:
        return None
    rest = description[idx + len("(PheCode ") :]
    end = rest.find(")")
    if end == -1:
        return None
    code = rest[:end].strip()
    parent = code.split(".", 1)[0]
    return parent if parent.isdigit() else None


def _gwas_p_value(description: str) -> float:
    """Extract p-value from a GWAS description. Returns inf if unparseable."""
    idx = description.find("(p=")
    if idx == -1:
        return float("inf")
    rest = description[idx + len("(p=") :]
    end = rest.find(",")
    if end == -1:
        end = rest.find(")")
    if end == -1:
        return float("inf")
    try:
        return float(rest[:end].strip())
    except ValueError:
        return float("inf")


def rollup_gwas_duplicates(annotations: list[Annotation]) -> list[Annotation]:
    """Collapse GWAS MTAG twins and PheCode parent/child hierarchies.

    Operates on the filtered annotation list (the output of
    AnalysisResult.filter). Non-GWAS rows pass through untouched.
    Must-include rows are never dropped.

    See ADR-0024 'MTAG and PheCode rollup' for rules.
    """
    survivors: list[Annotation] = []
    gwas_rows: list[Annotation] = []
    for a in annotations:
        (gwas_rows if a.source == "gwas" else survivors).append(a)

    if not gwas_rows:
        return annotations

    plain_keys = {
        (a.rsid, _gwas_base_trait(a.description))
        for a in gwas_rows
        if "(MTAG)" not in a.description
    }
    after_mtag = [
        a
        for a in gwas_rows
        if a.is_must_include
        or "(MTAG)" not in a.description
        or (a.rsid, _gwas_base_trait(a.description)) not in plain_keys
    ]

    by_parent: dict[tuple[str, str], list[Annotation]] = {}
    no_phecode: list[Annotation] = []
    for a in after_mtag:
        parent = _gwas_phecode_parent(a.description)
        if parent is None or a.is_must_include:
            no_phecode.append(a)
        else:
            by_parent.setdefault((a.rsid, parent), []).append(a)
    for group in by_parent.values():
        winner = min(group, key=lambda x: _gwas_p_value(x.description))
        no_phecode.append(winner)

    survivors.extend(no_phecode)
    survivors.sort(key=lambda a: (-a.magnitude, a.rsid))
    return survivors


def _lookup_user_allele(
    user_alt: str,
    coords: list[tuple[str, int, str, str]],
    scores: dict[tuple[str, int, str, str], float],
    resolve_strand: Callable[[str, str, str], str | None],
) -> float | None:
    """Find the CADD score for a specific user allele at a multi-allelic site.

    Prefers a direct allele match over a complement (minus-strand) match
    to avoid false positives where the complement of the user's allele
    coincidentally equals a different alt at the same position.
    """
    for chrom, pos, ref, alt in coords:
        if user_alt == alt:
            return scores.get((chrom, pos, ref, alt))
    for chrom, pos, ref, alt in coords:
        resolved = resolve_strand(user_alt, ref, alt)
        if resolved is not None and resolved == alt:
            return scores.get((chrom, pos, ref, alt))
    return None


def _enrich_cadd(
    annotations: list[Annotation],
    gnomad: GnomadAnnotator,
    cadd: CaddAnnotator,
) -> None:
    """Stamp annotations with CADD PHRED scores via coordinate resolution.

    Resolves rsIDs to genomic coordinates through gnomAD, normalizes
    alleles to reference-forward orientation, and looks up CADD scores.
    """
    from allelix.utils.allele import resolve_strand

    rsids = {a.rsid for a in annotations}
    coord_map = gnomad.bulk_resolve_coordinates(rsids)
    if not coord_map:
        return

    cadd_keys: set[tuple[str, int, str, str]] = set()
    for coords in coord_map.values():
        for chrom, pos, ref, alt in coords:
            cadd_keys.add((chrom, pos, ref, alt))
    scores = cadd.bulk_lookup(cadd_keys)
    if not scores:
        return

    for a in annotations:
        coords = coord_map.get(a.rsid)
        if not coords:
            continue
        if a.alt:
            score = _lookup_user_allele(a.alt, coords, scores, resolve_strand)
            a.cadd_phred = score
        else:
            best: float | None = None
            for chrom, pos, ref, alt in coords:
                score = scores.get((chrom, pos, ref, alt))
                if score is not None and (best is None or score > best):
                    best = score
            a.cadd_phred = best


def run_analysis(
    file_path: Path,
    parser: GenotypeParser,
    annotators: list[Annotator],
    skipped_count_provider: Callable[[], int] = lambda: 0,
    *,
    build_override: str | None = None,
    gnomad: GnomadAnnotator | None = None,
    alphamissense: AlphaMissenseAnnotator | None = None,
    cadd: CaddAnnotator | None = None,
    high_value_rsids: set[str] | None = None,
) -> AnalysisResult:
    """Stream the file once; batch-annotate; return a fully-populated result.

    Two phases over a single pass through the parser:

      Phase 1 — Build detection prelude. The first variants are
      buffered in `_BuildDetectionState` until the position-based
      detector converges (ADR-0021) or the detection buffer cap is
      hit. When the cap is hit without convergence, the prelude
      drains its buffer and the pipeline locks in the best-effort
      effective build.

      Phase 2 — Batched annotation. Each variant emerging from the
      build-detection phase (whether from the prelude drain or steady
      state) is appended to a batch buffer of `_BATCH_SIZE` variants.
      On every fill, the buffer is handed to each ready rsID-based
      annotator's `batch_annotate()` and cleared. EOF flushes the
      final partial batch.

    `high_value_rsids`, when supplied, collects variants whose rsID
    matches into `AnalysisResult.hv_variants` during the same pass —
    avoiding a second iteration of the parser for high-value no-call
    detection.

    Enrichment (gnomAD / AlphaMissense / CADD) runs after the streaming
    phase ends; those annotators already do bulk lookups against the
    collected annotation set.

    `build_override` short-circuits build detection: every variant
    gets that build and the position-data detector is skipped.

    Annotators are entered into a `contextlib.ExitStack` so their
    resources (e.g., SQLite connections) are deterministically closed.
    """
    metadata = parser.get_metadata(file_path)
    header_build = normalize_build_label(metadata.get("build"))

    annotations: list[Annotation] = []
    hv_variants: list[Variant] = []
    hv_set: set[str] = high_value_rsids or set()
    total = 0
    diag = _BuildDetectionState(override=build_override, header_build=header_build)
    # Coords for rsIDs the pipeline resolved on the fly (real-world VCFs from
    # variant callers emit ID=. — see GH #8). Lets the enrichment phase fall
    # back to position-keyed gnomAD / AlphaMissense lookups for resolved
    # rsIDs that don't appear in those caches' rsid index.
    resolved_coords: dict[str, tuple[str, int, str, str]] = {}

    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        clinvar_resolver = next((a for a in bound if a.name == "clinvar"), None)

        batch_buf: list[Variant] = []

        def _flush() -> None:
            if not batch_buf:
                return
            if clinvar_resolver is not None:
                # Resolve any variant whose ID column doesn't look like an
                # rsID: empty (GATK / DeepVariant emit "."), positional
                # synthetic IDs (1000 Genomes uses "22:10519265:CA:C"),
                # COSMIC IDs etc. None of those work as keys into the
                # rsID-indexed annotators. GH #8.
                rsidless = [v for v in batch_buf if not v.rsid.startswith("rs")]
                if rsidless:
                    resolution = clinvar_resolver.bulk_resolve_rsids(rsidless)
                    for (chrom, pos, ref, alt), rsid in resolution.items():
                        resolved_coords[rsid] = (chrom, pos, ref, alt)
            # High-value rsID match runs AFTER resolution so a variant
            # entering as ID=. that resolves to a high-value rsID is still
            # caught. GH #11.
            if hv_set:
                hv_variants.extend(v for v in batch_buf if v.rsid in hv_set)
            for annotator in bound:
                annotations.extend(annotator.batch_annotate(batch_buf))
            batch_buf.clear()

        def _accept(v: Variant) -> None:
            batch_buf.append(v)
            if len(batch_buf) >= _BATCH_SIZE:
                _flush()

        for variant in parser.parse(file_path):
            total += 1
            ready, batch = diag.feed(variant)
            if not ready:
                continue
            for v in batch:
                _accept(v)
        # End of stream: drain the build-detection buffer with the best
        # effective build we can resolve (detected → header → default),
        # then flush the final annotation batch.
        for v in diag.flush():
            _accept(v)
        _flush()

    if gnomad is not None and gnomad.is_ready():
        exact_keys = {(a.rsid, a.alt) for a in annotations if a.alt}
        max_rsids = {a.rsid for a in annotations if not a.alt}
        exact_freq = gnomad.bulk_lookup_by_alt(exact_keys)
        max_freq = gnomad.bulk_lookup(max_rsids)
        for a in annotations:
            if a.alt:
                a.allele_frequency = exact_freq.get((a.rsid, a.alt))
            else:
                a.allele_frequency = max_freq.get(a.rsid)
        # Position fallback for rsIDs resolved on the fly whose gnomAD rsid
        # index entry is missing/sparse. GH #8.
        if resolved_coords:
            missing_keys = {
                resolved_coords[a.rsid]
                for a in annotations
                if a.allele_frequency is None and a.rsid in resolved_coords
            }
            if missing_keys:
                pos_freq = gnomad.bulk_lookup_by_position(missing_keys)
                for a in annotations:
                    if a.allele_frequency is None and a.rsid in resolved_coords:
                        a.allele_frequency = pos_freq.get(resolved_coords[a.rsid])

    if alphamissense is not None and alphamissense.is_ready():
        exact_keys = {(a.rsid, a.alt) for a in annotations if a.alt}
        max_rsids = {a.rsid for a in annotations if not a.alt}
        exact_am = alphamissense.bulk_lookup_by_alt(exact_keys)
        max_am = alphamissense.bulk_lookup(max_rsids)
        for a in annotations:
            hit = exact_am.get((a.rsid, a.alt)) if a.alt else max_am.get(a.rsid)
            if hit is not None:
                a.am_pathogenicity, a.am_class = hit
        # Position fallback — see gnomAD block above. GH #8.
        if resolved_coords:
            missing_keys = {
                resolved_coords[a.rsid]
                for a in annotations
                if a.am_pathogenicity is None and a.rsid in resolved_coords
            }
            if missing_keys:
                pos_am = alphamissense.bulk_lookup_by_position(missing_keys)
                for a in annotations:
                    if a.am_pathogenicity is None and a.rsid in resolved_coords:
                        hit = pos_am.get(resolved_coords[a.rsid])
                        if hit is not None:
                            a.am_pathogenicity, a.am_class = hit

    if cadd is not None and cadd.is_ready() and gnomad is not None and gnomad.is_ready():
        if getattr(cadd, "_full_mode", False) and diag.effective_build != BUILD_GRCH38:
            logging.getLogger(__name__).warning(
                "CADD full mode requires GRCh38 coordinates; "
                "detected %s — skipping CADD enrichment",
                diag.effective_build,
            )
        else:
            _enrich_cadd(annotations, gnomad, cadd)

    annotators_used = [(a.name, a.version()) for a in annotators]
    if gnomad is not None and gnomad.is_ready():
        annotators_used.append((gnomad.name, gnomad.version()))
    if alphamissense is not None and alphamissense.is_ready():
        annotators_used.append((alphamissense.name, alphamissense.version()))
    if cadd is not None and cadd.is_ready():
        annotators_used.append((cadd.name, cadd.version()))

    return AnalysisResult(
        file_path=file_path,
        parser_name=parser.name,
        parser_display_name=parser.display_name,
        sample_id=metadata["sample_id"],
        build=diag.effective_build,
        total_variants=total,
        skipped_count=skipped_count_provider(),
        annotators_used=annotators_used,
        annotations=annotations,
        build_diagnostics=diag.diagnostics(),
        hv_variants=hv_variants,
    )


class _BuildDetectionState:
    """Buffer-and-flush state machine for build detection during streaming.

    `feed(variant)` returns (ready, batch). When `ready` is False, the
    variant has been buffered and the caller should keep streaming.
    When True, `batch` contains one or more variants with their build
    field set to the effective build, ready to be annotated.

    `flush()` is called at end of stream to drain anything still
    buffered (which only happens when detection never converged).
    """

    def __init__(self, *, override: str | None, header_build: str | None) -> None:
        self.header_build = header_build
        self.override = override
        # Effective build: starts as override (if given), else None until detection runs.
        self.effective: str | None = override
        self.detected: str | None = None
        self.matched_count = 0
        self.inspected_count = 0
        self._buffer: list[Variant] = []

    @property
    def effective_build(self) -> str:
        """Best-effort effective build at flush time."""
        return self.effective or self.header_build or BUILD_GRCH37

    def feed(self, variant: Variant) -> tuple[bool, list[Variant]]:
        if self.effective is not None:
            return True, [replace(variant, build=self.effective)]
        # Buffering until detection converges or we hit the cap.
        self._buffer.append(variant)
        if variant.rsid in KNOWN_SNP_POSITIONS:
            result = detect_build(self._buffer)
            if result.is_confident:
                self.detected = result.build
                self.matched_count = result.matched
                self.inspected_count = result.inspected
                self.effective = result.build
                batch = [replace(v, build=result.build) for v in self._buffer]
                self._buffer.clear()
                return True, batch
        if len(self._buffer) >= _DETECTION_BUFFER_LIMIT:
            # Buffer full before detection converged. Run partial detection
            # so the GRCh36 safety guard can fire (same logic as flush()).
            result = detect_build(self._buffer)
            if result.build is not None:
                self.detected = result.build
            self.matched_count = result.matched
            self.inspected_count = result.inspected
            if result.build == BUILD_GRCH36:
                self.effective = BUILD_GRCH36
            else:
                self.effective = self.header_build or BUILD_GRCH37
            batch = [replace(v, build=self.effective) for v in self._buffer]
            self._buffer.clear()
            return True, batch
        return False, []

    def flush(self) -> list[Variant]:
        if not self._buffer:
            return []
        # Detection never converged. Re-run on the full buffer to capture
        # partial counts even if not confident.
        result = detect_build(self._buffer)
        if result.is_confident:
            self.detected = result.build
            self.effective = result.build
        else:
            if result.build is not None:
                self.detected = result.build
            # GRCh36 must fail safe: there is no GRCh36 ClinVar cache,
            # so falling back to GRCh37 would silently query wrong
            # coordinates and bypass the GRCh36 safety guard.
            if result.build == BUILD_GRCH36:
                self.effective = BUILD_GRCH36
            else:
                self.effective = self.header_build or BUILD_GRCH37
        self.matched_count = result.matched
        self.inspected_count = result.inspected
        out = [replace(v, build=self.effective) for v in self._buffer]
        self._buffer.clear()
        return out

    def diagnostics(self) -> BuildDiagnostics:
        return BuildDiagnostics(
            header_build=self.header_build,
            detected_build=self.detected,
            effective_build=self.effective_build,
            override=self.override is not None,
            matched_count=self.matched_count,
            inspected_count=self.inspected_count,
        )


__all__ = [
    "BUILD_GRCH37",
    "BUILD_GRCH38",
    "AnalysisResult",
    "BuildDiagnostics",
    "rollup_gwas_duplicates",
    "run_analysis",
]

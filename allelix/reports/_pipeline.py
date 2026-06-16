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
from typing import TYPE_CHECKING, cast

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
    from allelix.annotators.clinvar import ClinVarAnnotator
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

    `chr_prefix_inferred` (GH #38): True when the effective build was
    picked using the ``chr``-prefixed contig heuristic (GRCh38
    convention). False whenever rsID detection or an explicit header
    build chose the answer, or when no chr-prefix signal was seen.
    Lets the CLI surface "inferred from chr-prefix" instead of the
    blind-default warning text.
    """

    header_build: str | None
    detected_build: str | None
    effective_build: str
    override: bool
    matched_count: int
    inspected_count: int
    chr_prefix_inferred: bool = False

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


# ADR-0035 PR 3 closed the GH #24 structural-half: GWAS rollup now reads
# trait / phecode / p_value from structured Annotation fields instead of
# regex-parsing the rendered ``description`` prose. The previous
# ``_gwas_base_trait`` / ``_gwas_phecode_parent`` / ``_gwas_p_value`` helpers
# and the ``_PHECODE_DELIM`` shared-delimiter constant (the v2.0.1 suppress-
# half) are gone — rollup keys directly off ``a.trait``, ``a.phecode``,
# ``a.p_value``. MTAG suffix still lives in ``description`` for display and
# is detected with ``"(MTAG)" in a.description``; promoting that to a
# structured flag is outside the ADR-0035 cluster manifest.


def _gwas_phecode_parent(phecode: str) -> str | None:
    """Return the integer parent of a structured PheCode (e.g., ``"411.4"`` → ``"411"``).

    Returns None when no PheCode is set or the prefix isn't an integer.
    """
    if not phecode:
        return None
    parent = phecode.split(".", 1)[0]
    return parent if parent.isdigit() else None


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

    plain_keys = {(a.rsid, a.trait.lower()) for a in gwas_rows if "(MTAG)" not in a.description}
    after_mtag = [
        a
        for a in gwas_rows
        if a.is_must_include
        or "(MTAG)" not in a.description
        or (a.rsid, a.trait.lower()) not in plain_keys
    ]

    by_parent: dict[tuple[str, str], list[Annotation]] = {}
    no_phecode: list[Annotation] = []
    for a in after_mtag:
        parent = _gwas_phecode_parent(a.phecode)
        if parent is None or a.is_must_include:
            no_phecode.append(a)
        else:
            by_parent.setdefault((a.rsid, parent), []).append(a)
    for group in by_parent.values():
        winner = min(group, key=lambda x: x.p_value if x.p_value is not None else float("inf"))
        no_phecode.append(winner)

    survivors.extend(no_phecode)
    survivors.sort(key=lambda a: (-a.magnitude, a.rsid))
    return survivors


def _lookup_user_allele(
    user_alt: str,
    coords: list[tuple[str, int, str, str]],
    scores: dict[tuple[str, int, str, str], float],
) -> float | None:
    """Find the CADD score for a specific user allele at a multi-allelic site.

    GH #18: only direct allele matches are accepted. The previous
    minus-strand fallback (``resolve_strand`` → complement match) could
    coincidentally hit a different alt at multi-allelic positions and
    stamp a wrong-allele CADD score. Until strand handling is properly
    plumbed (ADR-0010), enrichment is skipped rather than risked.
    """
    for chrom, pos, ref, alt in coords:
        if user_alt == alt:
            return scores.get((chrom, pos, ref, alt))
    return None


def _enrich_cadd(
    annotations: list[Annotation],
    gnomad: GnomadAnnotator,
    cadd: CaddAnnotator,
    resolved_coords: dict[str, tuple[str, int, str, str]] | None = None,
) -> None:
    """Stamp annotations with CADD PHRED scores via coordinate resolution.

    Resolves rsIDs to genomic coordinates through gnomAD, normalizes
    alleles to reference-forward orientation, and looks up CADD scores.

    GH #23 suppress-half: annotations without an explicit ``alt`` (raw
    GWAS rows) previously took a ``MAX(phred)`` fallback across every
    alt at the position — at multi-allelic sites that stamped the
    highest-CADD alt's score next to the annotation as if it described
    the user's variant. The MAX fallback is removed.

    ``resolved_coords`` (when supplied) carries the user's specific alt
    for rsIDs that the pipeline resolved on the fly via ClinVar's
    ``bulk_resolve_rsids`` (GH #8). That path is allele-specific —
    safe to fire for alt-less Annotation rows — and gives back the
    legitimate CADD enrichment that the MAX-fallback removal would
    otherwise lose for these rows. Symmetric with the gnomAD / AM
    position-fallback blocks in ``run_analysis``.
    """
    rsids = {a.rsid for a in annotations}
    coord_map = gnomad.bulk_resolve_coordinates(rsids)
    if not coord_map and not resolved_coords:
        return

    cadd_keys: set[tuple[str, int, str, str]] = set()
    for coords in coord_map.values():
        for chrom, pos, ref, alt in coords:
            cadd_keys.add((chrom, pos, ref, alt))
    if resolved_coords:
        cadd_keys.update(resolved_coords.values())
    scores = cadd.bulk_lookup(cadd_keys)
    if not scores:
        return

    for a in annotations:
        if a.alt:
            ann_coords = coord_map.get(a.rsid)
            if ann_coords:
                a.cadd_phred = _lookup_user_allele(a.alt, ann_coords, scores)
        elif resolved_coords and a.rsid in resolved_coords:
            # Position fallback: rsID was resolved on the fly via
            # ClinVar; the resolved tuple carries the user's actual
            # alt, so this is an allele-specific lookup, not a MAX.
            a.cadd_phred = scores.get(resolved_coords[a.rsid])


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
    # GH #38: chr-prefix on contigs is the strongest remaining heuristic
    # for the increasingly common case of rsID-less VCFs from modern
    # callers (DeepVariant / DRAGEN / GATK HC) that also lack
    # ``##contig assembly=`` tags. GRCh38 conventionally uses
    # ``chr1, chrX, chrM``; GRCh37 uses bare ``1, X, MT``. Only VCF
    # parsers populate this signal — consumer arrays always use bare
    # names regardless of build.
    chr_prefix_observed = bool(metadata.get("chr_prefix_observed", False))

    annotations: list[Annotation] = []
    hv_variants: list[Variant] = []
    hv_set: set[str] = high_value_rsids or set()
    total = 0
    diag = _BuildDetectionState(
        override=build_override,
        header_build=header_build,
        chr_prefix_observed=chr_prefix_observed,
    )
    # Coords for rsIDs the pipeline resolved on the fly (real-world VCFs from
    # variant callers emit ID=. — see GH #8). Lets the enrichment phase fall
    # back to position-keyed gnomAD / AlphaMissense lookups for resolved
    # rsIDs that don't appear in those caches' rsid index.
    resolved_coords: dict[str, tuple[str, int, str, str]] = {}

    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        # GH #36: the optional enrichment annotators were previously
        # constructed by callers and passed in by keyword without
        # context-management; their SQLite connections leaked at GC
        # time. Wire them into the same stack so cleanup is
        # deterministic alongside the primary annotators.
        for enrich in (gnomad, alphamissense, cadd):
            if enrich is not None:
                stack.enter_context(enrich)
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
                    resolution = cast("ClinVarAnnotator", clinvar_resolver).bulk_resolve_rsids(
                        rsidless
                    )
                    for (chrom, pos, ref, alt), rsid in resolution.items():
                        resolved_coords[rsid] = (chrom, pos, ref, alt)
            # ADR-0035 PR 4: populate Variant.ref for array data via gnomAD's
            # rsid → forward-REF map. VCF inputs already carry ref (PR 1).
            # Annotators downstream (ClinVar / PharmGKB strand-aware carrier
            # match; SNPedia / ClinPGx per-row alt threading from PR 2)
            # consume this. Variants without rsIDs (post-resolution) or
            # whose rsID isn't in gnomAD keep ref=None and degrade
            # gracefully — direct-match-only with no strand-flip.
            if gnomad is not None and gnomad.is_ready():
                need_ref = {v.rsid for v in batch_buf if v.ref is None and v.rsid.startswith("rs")}
                if need_ref:
                    coord_map = gnomad.bulk_resolve_coordinates(need_ref)
                    for v in batch_buf:
                        if v.ref is None and v.rsid in coord_map:
                            entries = coord_map[v.rsid]
                            # All multi-allelic rows at a position share the
                            # same REF; first entry's REF is canonical.
                            v.ref = entries[0][2] if entries else None
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

    # GH #23: enrichment is allele-specific. Annotations that carry an
    # explicit ``alt`` (ClinVar, ClinPGx, SNPedia, also GWAS rows whose
    # rsID resolved via ClinVar position lookup → recorded in
    # ``resolved_coords``) get an exact ``(rsid, alt)`` lookup. The old
    # code ran a ``MAX() GROUP BY rsid`` fallback for the alt-less case
    # (all original GWAS rows), which at multi-allelic sites stamps the
    # highest-frequency / highest-pathogenicity / highest-CADD alt's
    # value next to the user's annotation as if it described them. That
    # is the same wrong-allele hazard #18 fixed in the strand path.
    # Symmetric fix: skip enrichment rather than show a wrong-allele
    # number. The full fix — carrying the user's alt onto every
    # Annotation so GWAS rows can take the exact-alt path too — is
    # architectural and tracked for v2.1 (Variant.ref / per-annotation
    # allele tracking).
    if gnomad is not None and gnomad.is_ready():
        exact_keys = {(a.rsid, a.alt) for a in annotations if a.alt}
        exact_freq = gnomad.bulk_lookup_by_alt(exact_keys)
        for a in annotations:
            if a.alt:
                a.allele_frequency = exact_freq.get((a.rsid, a.alt))
        # Position fallback for rsIDs resolved on the fly via ClinVar
        # bulk_resolve_rsids (GH #8). The resolved-coords map already
        # carries the user's specific alt — this is an allele-specific
        # lookup, NOT a MAX, so it's safe to fire for alt-less rows
        # whose rsID was resolved this way.
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
        exact_am = alphamissense.bulk_lookup_by_alt(exact_keys)
        for a in annotations:
            if a.alt:
                hit = exact_am.get((a.rsid, a.alt))
                if hit is not None:
                    a.am_pathogenicity, a.am_class = hit
        # Position fallback — see gnomAD block above. Safe (allele-specific).
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
            _enrich_cadd(annotations, gnomad, cadd, resolved_coords=resolved_coords or None)

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

    def __init__(
        self,
        *,
        override: str | None,
        header_build: str | None,
        chr_prefix_observed: bool = False,
    ) -> None:
        self.header_build = header_build
        self.override = override
        # GH #38: ``chr``-prefixed contig names indicate GRCh38 in modern
        # variant callers. Tertiary signal — falls in priority after
        # override > rsID detection > header_build, ahead of the bare
        # GRCh37 fallback.
        self.chr_prefix_observed = chr_prefix_observed
        # Effective build: starts as override (if given), else None until detection runs.
        self.effective: str | None = override
        self.detected: str | None = None
        self.matched_count = 0
        self.inspected_count = 0
        self._buffer: list[Variant] = []

    @property
    def effective_build(self) -> str:
        """Best-effort effective build at flush time.

        Priority order:
        1. ``override`` (--build flag) — already applied to ``self.effective`` at init.
        2. Position-based ``detected`` (set by ``feed()`` / ``flush()``).
        3. ``header_build`` (``##contig assembly=...`` tag normalized).
        4. ``chr_prefix_observed`` → GRCh38 (GH #38).
        5. ``BUILD_GRCH37`` fallback.
        """
        if self.effective:
            return self.effective
        if self.header_build:
            return self.header_build
        if self.chr_prefix_observed:
            return BUILD_GRCH38
        return BUILD_GRCH37

    def feed(self, variant: Variant) -> tuple[bool, list[Variant]]:
        if self.effective is not None:
            return True, [replace(variant, build=self.effective)]
        # Buffering until detection converges or we hit the cap.
        self._buffer.append(variant)
        if variant.rsid in KNOWN_SNP_POSITIONS:
            result = detect_build(self._buffer)
            if result.is_confident:
                assert result.build is not None  # is_confident ⇒ build set
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
                # Fallback priority matches `effective_build` property:
                # header_build > chr_prefix_observed (GRCh38) > GRCh37.
                self.effective = self.header_build or (
                    BUILD_GRCH38 if self.chr_prefix_observed else BUILD_GRCH37
                )
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
                # Fallback priority matches `effective_build` property:
                # header_build > chr_prefix_observed (GRCh38) > GRCh37.
                self.effective = self.header_build or (
                    BUILD_GRCH38 if self.chr_prefix_observed else BUILD_GRCH37
                )
        self.matched_count = result.matched
        self.inspected_count = result.inspected
        assert self.effective is not None  # all branches above set it
        out = [replace(v, build=self.effective) for v in self._buffer]
        self._buffer.clear()
        return out

    def diagnostics(self) -> BuildDiagnostics:
        # GH #38: chr_prefix_inferred is True only when the
        # chr-prefix signal is what actually picked the effective
        # build — i.e., no override, no rsID detection, no header
        # build, and the chr-prefix signal flipped the fallback from
        # GRCh37 to GRCh38. Matches the priority order in the
        # ``effective_build`` property.
        chr_prefix_inferred = (
            self.override is None
            and self.detected is None
            and self.header_build is None
            and self.chr_prefix_observed
            and self.effective_build == BUILD_GRCH38
        )
        return BuildDiagnostics(
            header_build=self.header_build,
            detected_build=self.detected,
            effective_build=self.effective_build,
            override=self.override is not None,
            matched_count=self.matched_count,
            inspected_count=self.inspected_count,
            chr_prefix_inferred=chr_prefix_inferred,
        )


__all__ = [
    "BUILD_GRCH37",
    "BUILD_GRCH38",
    "AnalysisResult",
    "BuildDiagnostics",
    "rollup_gwas_duplicates",
    "run_analysis",
]

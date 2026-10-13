# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared analysis pipeline used by `analyze`, `methylation`, and `pharmacogenomics`.

The CLI builds an `AnalysisResult` once and hands it to a renderer
(terminal, JSON, HTML). Renderers never query the database or re-iterate
the parser — they receive a fully-populated value object.

ADR-0044: the build is decided once, before any annotation, from the
evidence the file carries (``decide_file_build``). Parsers report the
build the file states, if any; the pipeline sets each variant's `build`
to the decided build before annotators see it, or stops.
"""

from __future__ import annotations

import contextlib
import logging
import threading
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, TypedDict, cast

from allelix.annotators._resolution import resolve_batch
from allelix.models import DEFAULT_BUILD
from allelix.utils.build_detect import (
    BUILD_GRCH37,
    BUILD_GRCH38,
    KIND_ARRAY,
    KIND_VCF,
    NO_SEEDS,
    SOURCE_CONTIG_LENGTHS,
    SOURCE_OVERRIDE,
    BuildDecision,
    BuildDecisionError,
    BuildDetectionResult,
    check_late_seeds,
    contig_length_verdict,
    decide_build,
    detect_build,
    known_snp_seed,
    normalize_build_label,
)
from allelix.utils.plausibility import (
    DEFAULT_THRESHOLDS,
    PlausibilityThresholds,
    classify,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.models import Annotation, Variant
    from allelix.parsers.base import GenotypeMetadata, GenotypeParser
    from allelix.parsers.vcf import Locus, LocusState, VcfParser, VcfQualityStats
    from allelix.reports.vcf import VcfWriteInputs


# How many variants to accumulate before flushing to each annotator's
# batch_annotate. Each annotator further chunks internally for SQL
# parameter limits, so this is a Python-side batch boundary, not an
# SQL limit. 5000 keeps memory bounded at ~1.5 MB per pipeline buffer
# while letting the annotators' IN clauses run in batches of several
# hundred. Tuneable; the pipeline is correct at any size including 1.
_BATCH_SIZE = 5000


@dataclass
class BuildDiagnostics:
    """How the run's genome build was decided (ADR-0044).

    `effective_build` is the build used for annotation. `source` is one of
    ``allelix.utils.build_detect.BUILD_SOURCES``. `file_label` is the build
    the file states (None when it states none). `detected_build`,
    `matched_count`, `inspected_count` and `position_confident` describe the
    known-SNP evidence (from the whole-file pre-scan, or from the annotation
    pass when contig lengths decided). `lengths_status` / `lengths_detail`
    describe a VCF's contig lengths (None for other inputs).
    `label_overridden` is True when confident position evidence replaced a
    label the file stated; `tentative_grch36` marks the ADR-0025 fail-safe;
    `mt_skipped` counts MT records left out because chrM has the hg19
    length.
    """

    effective_build: str
    source: str
    file_label: str | None
    detected_build: str | None
    matched_count: int
    inspected_count: int
    position_confident: bool
    lengths_status: str | None = None
    lengths_detail: str | None = None
    label_overridden: bool = False
    tentative_grch36: bool = False
    mt_skipped: int = 0
    tied_builds: tuple[str, ...] = ()
    tied_votes: int = 0

    @property
    def inconclusive(self) -> bool:
        """True when the known-SNP positions tied and the build is one of the tied builds."""
        return bool(self.tied_builds) and not self.override

    @property
    def override(self) -> bool:
        """True when ``--build`` forced the build."""
        return self.source == SOURCE_OVERRIDE

    @property
    def mismatch(self) -> bool:
        """True when position evidence replaced the build the file states."""
        return self.label_overridden

    @classmethod
    def from_decision(
        cls,
        decision: BuildDecision,
        seeds: BuildDetectionResult | None = None,
        mt_skipped: int = 0,
    ) -> BuildDiagnostics:
        """Diagnostics for ``decision``; ``seeds`` replaces its seed evidence when given."""
        evidence = decision.seeds if seeds is None else seeds
        lengths = decision.lengths
        return cls(
            effective_build=decision.build,
            source=decision.source,
            file_label=decision.file_label,
            detected_build=evidence.build,
            matched_count=evidence.matched,
            inspected_count=evidence.inspected,
            position_confident=evidence.is_confident,
            lengths_status=None if lengths is None else lengths.status,
            lengths_detail=None if lengths is None else lengths.detail,
            label_overridden=decision.label_overridden,
            tentative_grch36=decision.tentative_grch36,
            mt_skipped=mt_skipped,
            tied_builds=evidence.tied,
            tied_votes=evidence.tied_votes,
        )

    def evidence(self) -> dict[str, object]:
        """JSON-ready evidence record (``input.build_evidence``)."""
        return {
            "file_label": self.file_label,
            "contig_lengths": self.lengths_status,
            "contig_lengths_detail": self.lengths_detail,
            "seeds_build": self.detected_build,
            "seeds_matched": self.matched_count,
            "seeds_inspected": self.inspected_count,
            "seeds_confident": self.position_confident,
            "label_overridden": self.label_overridden,
            "mt_records_skipped": self.mt_skipped,
            "seeds_tied": list(self.tied_builds),
        }


class PanelCoverage(TypedDict):
    """Shape returned by AnalysisResult.panel_coverage().

    Same shape lands in the JSON report under ``panel_coverage`` and
    drives the HTML Panel Coverage section.
    """

    requested: int
    found: int
    missing: list[str]
    no_findings: list[str]


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
    # when --filter-file supplies a rsid panel, the pipeline
    # records which panel rsIDs the input file actually carries so
    # reporters can surface coverage (state 2: genotyped but no
    # annotations; state 3: not on the user's chip). Both None when
    # no panel was supplied — that path is the default analyze run
    # with no behavior change.
    panel_rsids: frozenset[str] | None = None
    genotyped_panel_rsids: frozenset[str] | None = None
    panel_annotated_rsids: frozenset[str] | None = None
    # per-rsID genotype string (e.g. "C/T", "A/-", "-/-")
    # for each panel rsID present in the input file. Populated only
    # when --filter-file supplied a panel. None on default analyze
    # runs (zero behavior change on the default path). The HTML
    # Panel coverage section renders this in a GT column; JSON shape
    # stays v6 untouched — display-only addition. The genotype-aware
    # MECE redesign that lands this in JSON + terminal is planned for v2.3.
    panel_genotypes: dict[str, str] | None = None
    # variants reaching the annotator phase with ref=None.
    # Used by the CLI to surface the strand-aware-matching-inactive
    # info line once per run for array inputs.
    no_ref_variant_count: int = 0
    # variants reaching the annotator phase whose rsID isn't
    # an "rs..." string (variant-caller VCFs emit ID=. — see the rsID-less resolution step).
    # Used by the CLI to surface the GRCh38 rsID-less undercount
    # warning until the dbSNP resolution lands.
    rsidless_variant_count: int = 0
    # pipeline-produced inputs for the annotated VCF writer.
    # Populated only when ``run_analysis`` is called with
    # ``collect_vcf_write_inputs=True`` — the ``--vcf-out`` code path.
    # None on every other analyze run (zero memory / CPU overhead on the
    # default path). Consumed by ``allelix.reports.vcf.render_vcf``.
    vcf_write_inputs: VcfWriteInputs | None = None
    # what the VCF quality gate did (caller-filtered rows skipped or
    # kept, low-GQ calls turned into no-calls). None for non-VCF input.
    vcf_quality: VcfQualityStats | None = None
    # How the input's genotypes were produced, when the file declares it
    # (GenotypeMetadata ``method``, e.g. "low-pass-wgs-imputed" for
    # MyHeritage low-pass WGS). Empty when not declared.
    input_method: str = ""
    # per high-value rsID, what a VCF says at its position on the
    # effective build: (called | reference | no_call, reason). A gVCF
    # position no record covers is (no_call, no_coverage). Absent rsIDs
    # are unknown (plain VCF, or no placement). None for non-VCF input.
    high_value_vcf_states: dict[str, LocusState] | None = None

    def panel_coverage(
        self, filtered_annotations: list[Annotation] | None = None
    ) -> PanelCoverage | None:
        """Return panel-coverage info when a rsid panel was supplied.

        States,

        - **state 1** (genotyped + has annotations): rsid is in
          ``genotyped_panel_rsids`` AND ``panel_annotated_rsids``;
          rendered normally in the main annotation table.
        - **state 2** (genotyped + no annotations): rsid is in
          ``genotyped_panel_rsids`` but NOT in
          ``panel_annotated_rsids``; rendered as a "no findings" entry
          in panel coverage so the user knows it was on their chip but
          carries no risk allele.
        - **state 3** (not genotyped): rsid is in ``panel_rsids`` but
          NOT in ``genotyped_panel_rsids``; rendered as "not on chip".

        Returns None when no panel was supplied (no coverage section
        should appear). Returns the dict shape rendered into the JSON
        report and the HTML Panel Coverage section.

        A no-call (``-/-``) panel variant that IS in the input file
        counts as ``found`` and lands in ``no_findings`` (state 2),
        not ``missing`` (state 3). This matches the issue's literal
        "not in input file" definition for state 3 — "found" reports
        chip presence, not call success. High-value no-calls are
        separately surfaced by the ``hv_warnings`` path in the CLI;
        readers should not interpret ``found`` as "called."

        (v2.2.1 patch): when ``filtered_annotations`` is
        supplied — the post-magnitude-filter list each renderer
        already produces for the main annotation table — the
        "annotated" set is derived from THAT list, not from
        ``panel_annotated_rsids`` (which is populated from the
        unfiltered annotation set in ``run_analysis``). Without
        this, panel rsIDs whose only annotations score below the
        analyze display floor land in ``found`` but render nowhere
        — they're counted as "found" then filtered out of the
        annotation table, then NOT folded into ``no_findings``
        (because the pipeline thinks they had annotations). Pure
        limbo, a flat lie in ``found``. Passing the filtered list
        in makes "annotated" equal "actually rendered" and the
        below-floor variants drop into ``no_findings`` honestly.

        The redesign (genotype-aware MECE 4-state with `no_call`
        distinct from `typed_no_annotation`, schema 6→7) lands in
        v2.3.
        """
        if self.panel_rsids is None:
            return None
        genotyped = self.genotyped_panel_rsids or frozenset()
        panel_rsids = self.panel_rsids
        if filtered_annotations is not None:
            # patch: "annotated" is the set of panel rsIDs that
            # actually got rendered — derived from the post-filter
            # annotation list the renderer passes in. Sub-floor-only
            # panel rsids drop out of `annotated`, so they fall into
            # `no_findings` instead of limbo. The `found` field still
            # means "in input file" (genotyped count) — same v6
            # contract — only `no_findings` widens to be honest.
            annotated: frozenset[str] = frozenset(
                a.rsid for a in filtered_annotations if a.rsid in panel_rsids
            )
        else:
            # Back-compat path for callers that don't thread the
            # filtered list. `annotated` stays the raw pre-filter set,
            # reproducing the v6.0 behavior (sub-floor panel rsids
            # land in `found` and vanish from every rendered surface —
            # the limbo bug the fix retires when the filtered
            # list IS supplied).
            annotated = self.panel_annotated_rsids or frozenset()
        missing = sorted(panel_rsids - genotyped)
        # Below-floor panel rsids land here on the filtered path: in
        # genotyped but not in the post-filter annotated set. They had
        # annotations in the cache but none survived the display floor.
        no_findings = sorted(genotyped - annotated)
        return {
            "requested": len(panel_rsids),
            "found": len(genotyped),
            "missing": missing,
            "no_findings": no_findings,
        }

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


# ADR-0035: GWAS rollup reads trait / phecode / p_value from structured
# Annotation fields instead of regex-parsing the rendered ``description``
# prose. The prose helpers ``_gwas_base_trait`` and ``_gwas_p_value`` and the
# ``_PHECODE_DELIM`` shared-delimiter constant (the v2.0.1 suppress-half) are
# gone; rollup keys directly off ``a.trait``, ``a.phecode`` and ``a.p_value``,
# and ``_gwas_phecode_parent`` below takes the structured ``a.phecode``. The
# MTAG suffix still lives in ``description`` for display and is detected with
# ``"(MTAG)" in a.description``; promoting that to a structured flag is outside
# the ADR-0035 cluster manifest.


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


_IDENTITY_BUILDS: frozenset[str] = frozenset({BUILD_GRCH37, BUILD_GRCH38})

EnrichmentKey = tuple[str, int, str, str]


def watch_high_value(
    parser: GenotypeParser,
    hv_set: set[str],
    identity: RsidIdentityAnnotator | None,
) -> dict[str, dict[str, Locus]] | None:
    """Ask a VCF parser to record what the file says at high-value sites.

    Positions come from the identity table on both builds (the effective
    build is not known until streaming ends); rsIDs are watched directly
    as well. Returns the placements, or None when nothing is watched
    (non-VCF input or no high-value list).
    """
    from allelix.parsers.vcf import VcfParser

    if not hv_set or not isinstance(parser, VcfParser):
        return None
    placements: dict[str, dict[str, Locus]] = {}
    if identity is not None and identity.is_ready():
        placements = identity.loci_of_rsids(hv_set)
    parser.watch({locus for p in placements.values() for locus in p.values()}, hv_set)
    return placements


def high_value_states(
    stats: VcfQualityStats | None,
    placements: dict[str, dict[str, Locus]] | None,
    hv_set: set[str],
    build: str,
) -> dict[str, LocusState] | None:
    """Each high-value rsID's state on the effective build, from the watch results.

    A record carrying the rsID and a record covering its position both
    count; the stronger state wins. On a gVCF, a placed rsID that no
    record covers is a no-call (``no_coverage``). On a plain VCF it is
    unknown and left out.
    """
    from allelix.parsers.vcf import NO_CALL, REASON_NO_COVERAGE, strongest_state

    if stats is None or placements is None:
        return None
    states: dict[str, LocusState] = {}
    for rsid in hv_set:
        locus = placements.get(rsid, {}).get(build)
        state = stats.rsid_states.get(rsid)
        if locus is not None:
            state = strongest_state(state, stats.locus_states.get(locus))
        if state is None and locus is not None and stats.is_gvcf:
            state = (NO_CALL, REASON_NO_COVERAGE)
        if state is not None:
            states[rsid] = state
    return states


def _enrichment_allele(a: Annotation) -> str | None:
    """The non-reference allele an annotation's enrichment describes, or None.

    The annotation's ``alt``; when it has none (a GWAS row with no reported
    risk allele), the user's single carried non-reference allele. None when
    the locus or its REF is unknown, when there is no such allele or more
    than one, or when the allele is the reference, unless the annotation
    says the build's reference carries the variant allele
    (``reference_carries_alt``, F5 c.1601G>A on GRCh37): that allele is
    the variant, and its GRCh38 key is the real variant there.
    """
    if a.locus is None or a.locus_ref is None:
        return None
    allele = a.alt
    if not allele:
        carried = {x for x in a.carried_alleles if x != a.locus_ref}
        if len(carried) != 1:
            return None
        allele = carried.pop()
    if allele == a.locus_ref and not a.reference_carries_alt:
        return None
    return allele


def _flag_plausibility(
    keyed: list[tuple[Annotation, EnrichmentKey]],
    thresholds: PlausibilityThresholds,
) -> None:
    """Stamp each enriched annotation's plausibility (ADR-0041).

    Runs after gnomAD has stamped ``allele_frequency``. The copy count is
    of the same allele the frequency describes, counted in the user's
    alleles in the annotation's orientation.
    """
    for a, _key in keyed:
        allele = _enrichment_allele(a)
        if allele is None or a.locus is None:
            continue
        a.plausibility = classify(
            a.allele_frequency, a.carried_alleles.count(allele), a.locus[0], thresholds
        )


def _enrichment_keys(
    annotations: list[Annotation],
    build: str,
    identity: RsidIdentityAnnotator | None,
) -> list[tuple[Annotation, EnrichmentKey]]:
    """Pair each annotation with the GRCh38 key its enrichment is looked up by.

    The input-build key is ``(locus chrom, locus pos, locus_ref, allele)``.
    The allele is the annotation's ``alt``; an annotation without one (a
    GWAS row with no reported risk allele) uses the user's single carried
    non-reference allele, as the rsID-less position fallback always did,
    and gets no key when there is none or more than one. An allele equal
    to the reference has no variant to look up, except where the build's
    reference carries the variant allele (``reference_carries_alt``). GRCh38 keys are used as
    they are; GRCh37 keys go through the identity table; any other build
    gets no keys.
    """
    keyed: list[tuple[Annotation, EnrichmentKey]] = []
    for a in annotations:
        allele = _enrichment_allele(a)
        if allele is None or a.locus is None or a.locus_ref is None:
            continue
        keyed.append((a, (a.locus[0], a.locus[1], a.locus_ref, allele)))
    if build == BUILD_GRCH38:
        return keyed
    if build == BUILD_GRCH37 and identity is not None and identity.is_ready():
        translated = identity.translate_to_grch38({key for _a, key in keyed})
        return [(a, translated[key]) for a, key in keyed if key in translated]
    return []


def _enrich(
    keyed: list[tuple[Annotation, EnrichmentKey]],
    gnomad: GnomadAnnotator | None,
    alphamissense: AlphaMissenseAnnotator | None,
    cadd: CaddAnnotator | None,
) -> list[Annotator]:
    """Stamp gnomAD AF, AlphaMissense and CADD by GRCh38 key.

    Returns the enrichment sources that were queried, for provenance.
    """
    used: list[Annotator] = []
    if not keyed:
        return used
    keys = {key for _a, key in keyed}
    if gnomad is not None and gnomad.is_ready():
        used.append(gnomad)
        freqs = gnomad.bulk_lookup_by_position(keys)
        for a, key in keyed:
            a.allele_frequency = freqs.get(key)
    if alphamissense is not None and alphamissense.is_ready():
        used.append(alphamissense)
        scores = alphamissense.bulk_lookup_by_position(keys)
        for a, key in keyed:
            hit = scores.get(key)
            if hit is not None:
                a.am_pathogenicity, a.am_class = hit
    if cadd is not None and cadd.is_ready():
        used.append(cadd)
        phreds = cadd.bulk_lookup(keys)
        for a, key in keyed:
            a.cadd_phred = phreds.get(key)
    return used


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
    identity: RsidIdentityAnnotator | None = None,
    plausibility: PlausibilityThresholds = DEFAULT_THRESHOLDS,
    high_value_rsids: set[str] | None = None,
    panel_rsids: frozenset[str] | None = None,
    collect_vcf_write_inputs: bool = False,
) -> AnalysisResult:
    """Stream the file once; batch-annotate; return a fully-populated result.

    Two phases over a single pass through the parser:

      Phase 1 — Build decision (ADR-0044), before any annotation:
      ``decide_file_build`` reads the evidence the file carries (contig
      lengths, or known-SNP positions from a whole-file pre-scan) and
      either decides the build or raises ``BuildDecisionError``.

      Phase 2 — Batched annotation. Each variant, with its build set to
      the decided build, is appended to a batch buffer of `_BATCH_SIZE`
      variants.
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

    `plausibility` holds the allele-frequency cut-offs for the
    zygosity x gnomAD plausibility flag (ADR-0041), applied after gnomAD
    enrichment.

    `build_override` short-circuits the decision: every variant gets
    that build. When contig lengths decided, the known-SNP positions met
    during the pass are checked afterwards; a contradiction raises
    ``BuildDecisionError`` before anything is returned or written.

    Annotators are entered into a `contextlib.ExitStack` so their
    resources (e.g., SQLite connections) are deterministically closed.
    """
    metadata = parser.get_metadata(file_path)
    decision = decide_file_build(parser, file_path, metadata, override=build_override)

    annotations: list[Annotation] = []
    hv_variants: list[Variant] = []
    hv_set: set[str] = high_value_rsids or set()
    # rsIDs the input file actually carries that are in the
    # supplied panel. Only collected when ``panel_rsids`` is non-empty
    # (so a normal analyze run has zero overhead). Filtered at insert
    # so a 600k-variant input contributes ~panel-size memory, not
    # 600k rsids.
    genotyped_panel: set[str] = set()
    # per-rsID genotype for panel members. First-occurrence
    # wins (a 23andMe-style chip emits one row per rsID anyway; multi-
    # row input formats keep the earliest reading). HTML-only surface
    # in v2.2.1; JSON + terminal addition planned for v2.3.
    panel_genotypes_map: dict[str, str] = {}
    # Per-run diagnostic totals for ref-None and
    # rsidless variants reaching the annotator phase. Used by the CLI
    # to emit the "strand-aware inactive" info line and the
    # "GRCh38 rsID-less undercount" warning.
    no_ref_total = 0
    rsidless_total = 0
    total = 0
    effective_build = decision.build
    # Known-SNP records met during the pass; only kept when contig lengths
    # decided, so check_late_seeds can compare them with the lengths.
    late_seeds: list[Variant] = []
    track_seeds = decision.source == SOURCE_CONTIG_LENGTHS
    mt_skipped = 0
    # Coords for rsIDs the pipeline resolved on the fly (real-world VCFs from
    # variant callers emit ID=.). Lets the enrichment phase fall
    # back to position-keyed gnomAD / AlphaMissense lookups for resolved
    # rsIDs that don't appear in those caches' rsid index.
    resolved_coords: dict[str, tuple[str, int, str, str]] = {}
    # annotation dict keyed by (chrom, pos, ref, alt) for the
    # annotated-VCF writer. Populated per batch inside ``_flush`` when
    # ``collect_vcf_write_inputs`` is True. Values are references to the
    # Annotation objects in ``annotations`` — enrichment mutations
    # (gnomAD AF, AlphaMissense, CADD) that happen after streaming are
    # visible through these refs without a per-annotation copy.
    vcf_annotations_by_key: dict[tuple[str, int, str, str], list[Annotation]] = {}

    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        # the optional enrichment annotators were previously
        # constructed by callers and passed in by keyword without
        # context-management; their SQLite connections leaked at GC
        # time. Wire them into the same stack so cleanup is
        # deterministic alongside the primary annotators.
        for enrich in (gnomad, alphamissense, cadd, identity):
            if enrich is not None:
                stack.enter_context(enrich)
        clinvar_resolver = next((a for a in bound if a.name == "clinvar"), None)

        batch_buf: list[Variant] = []

        def _flush() -> None:
            if not batch_buf:
                return
            # REF population and rsID-less resolution are one shared
            # step (ClinVar, identity table, gnomAD), the same for every
            # command that needs rsIDs. See allelix.annotators._resolution.
            resolve_batch(
                batch_buf,
                clinvar=cast("ClinVarAnnotator | None", clinvar_resolver),
                gnomad=gnomad,
                identity=identity,
                resolved_coords=resolved_coords,
            )
            # High-value rsID match runs AFTER resolution so a variant
            # entering as ID=. that resolves to a high-value rsID is still
            # caught.
            if hv_set:
                hv_variants.extend(v for v in batch_buf if v.rsid in hv_set)
            if panel_rsids:
                for v in batch_buf:
                    if v.rsid in panel_rsids:
                        genotyped_panel.add(v.rsid)
                        # First-occurrence wins; renders as "A/G",
                        # "C/C", "A/-", "-/-" honestly — uses the
                        # same Variant.genotype property the rest of
                        # the codebase uses.
                        panel_genotypes_map.setdefault(v.rsid, v.genotype)
            # Per-batch diagnostic counters. Tallied AFTER
            # the ref-population and resolver passes so the counts
            # reflect the final state annotators see — not the raw input.
            nonlocal no_ref_total, rsidless_total
            no_ref_total += sum(1 for v in batch_buf if v.ref is None)
            rsidless_total += sum(1 for v in batch_buf if not v.rsid.startswith("rs"))
            batch_new: list[Annotation] = []
            for annotator in bound:
                batch_new.extend(annotator.batch_annotate(batch_buf))
            annotations.extend(batch_new)
            # collect writer inputs for this batch. Only touches
            # variants whose ``ref`` is populated (VCF inputs; array
            # data with ref=None can't produce a well-formed VCF key).
            # Skips annotations with empty ``alt`` — GWAS rows keyed by
            # rsID alone can't be VCF-stamped without allele context.
            if collect_vcf_write_inputs:
                rsid_to_coords: dict[str, tuple[str, int, str]] = {
                    v.rsid: (v.chromosome, v.position, v.ref)
                    for v in batch_buf
                    if v.ref is not None
                }
                for a in batch_new:
                    if not a.alt:
                        continue
                    coords = rsid_to_coords.get(a.rsid)
                    if coords is None:
                        continue
                    chrom, pos, ref = coords
                    vcf_annotations_by_key.setdefault((chrom, pos, ref, a.alt), []).append(a)
            batch_buf.clear()

        def _accept(v: Variant) -> None:
            batch_buf.append(v)
            if len(batch_buf) >= _BATCH_SIZE:
                _flush()

        # watch high-value positions so gVCF reference blocks and
        # no-call records, which carry no rsIDs, can be classified.
        hv_placements = watch_high_value(parser, hv_set, identity)
        try:
            for variant in parser.parse(file_path):
                total += 1
                if decision.skip_mt and variant.chromosome == "MT":
                    mt_skipped += 1
                    continue
                if track_seeds and known_snp_seed(variant) is not None:
                    late_seeds.append(variant)
                _accept(replace(variant, build=effective_build))
        finally:
            if hv_placements is not None:
                cast("VcfParser", parser).watch(set(), set())
        _flush()
        late = detect_build(late_seeds) if track_seeds else None
        if late is not None:
            check_late_seeds(decision, late, file_path.name)

    # ADR-0040: one enrichment step for every build. Each annotation gets a
    # GRCh38 ``(chrom, pos, ref, alt)`` key from its own locus (GRCh38
    # input) or through the identity table (GRCh37 input; palindromic and
    # reference-swapped alleles abstain), and gnomAD / AlphaMissense / CADD
    # are all queried by that key. No rsID-keyed or MAX lookups:
    # an annotation without an allele-specific key gets no value.
    enrichment_keys = _enrichment_keys(annotations, effective_build, identity)
    enriched_by = _enrich(enrichment_keys, gnomad, alphamissense, cadd)
    # ADR-0041: zygosity x gnomAD frequency, only where gnomAD was queried.
    if gnomad in enriched_by:
        _flag_plausibility(enrichment_keys, plausibility)

    annotators_used = [(a.name, a.version()) for a in annotators]
    if identity is not None and identity.is_ready() and effective_build in _IDENTITY_BUILDS:
        annotators_used.append((identity.name, identity.version()))
    annotators_used.extend((e.name, e.version()) for e in enriched_by)

    # panel coverage derives genotyped rsIDs above (state 2+3
    # boundary) and annotated rsIDs from the post-enrichment annotation
    # set (state 1 boundary). Empty frozensets vs None disambiguate
    # "no panel supplied" from "panel supplied, zero hits".
    if panel_rsids:
        result_panel = frozenset(panel_rsids)
        result_genotyped = frozenset(genotyped_panel)
        result_annotated = frozenset(a.rsid for a in annotations if a.rsid in panel_rsids)
        result_genotypes: dict[str, str] | None = dict(panel_genotypes_map)
    else:
        result_panel = None
        result_genotyped = None
        result_annotated = None
        result_genotypes = None

    vcf_write_inputs: VcfWriteInputs | None = None
    if collect_vcf_write_inputs:
        # Import here to avoid a top-level cycle with allelix.reports.vcf
        # (which reaches into the annotator registry for attribution).
        from allelix.reports.vcf import VcfWriteInputs as _VcfWriteInputs

        vcf_recovered_rsids = {
            (chrom, pos, ref, alt): rsid
            for rsid, (chrom, pos, ref, alt) in resolved_coords.items()
        }
        vcf_write_inputs = _VcfWriteInputs(
            annotations_by_key=vcf_annotations_by_key,
            recovered_rsids=vcf_recovered_rsids,
        )

    return AnalysisResult(
        file_path=file_path,
        parser_name=parser.name,
        parser_display_name=parser.display_name,
        sample_id=metadata["sample_id"],
        input_method=metadata.get("method", ""),
        build=effective_build,
        total_variants=total,
        skipped_count=skipped_count_provider(),
        annotators_used=annotators_used,
        annotations=annotations,
        build_diagnostics=BuildDiagnostics.from_decision(decision, late, mt_skipped),
        hv_variants=hv_variants,
        panel_rsids=result_panel,
        genotyped_panel_rsids=result_genotyped,
        panel_annotated_rsids=result_annotated,
        panel_genotypes=result_genotypes,
        no_ref_variant_count=no_ref_total,
        rsidless_variant_count=rsidless_total,
        vcf_write_inputs=vcf_write_inputs,
        vcf_quality=getattr(parser, "quality_stats", None),
        high_value_vcf_states=high_value_states(
            getattr(parser, "quality_stats", None), hv_placements, hv_set, effective_build
        ),
    )


# The pre-scan parses the file once more; its malformed-line warnings would
# repeat the ones the real pass reports. A thread-local flag (not a global
# logging switch) keeps another thread's parser warnings intact.
_PRESCAN = threading.local()


class _PrescanFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return not getattr(_PRESCAN, "active", False)


_PRESCAN_FILTER = _PrescanFilter()


@contextlib.contextmanager
def _quiet_prescan() -> Iterator[None]:
    """Drop ``allelix.parsers`` log records raised by this thread inside the block."""
    for name in list(logging.Logger.manager.loggerDict):
        if name == "allelix.parsers" or name.startswith("allelix.parsers."):
            parser_logger = logging.getLogger(name)
            if _PRESCAN_FILTER not in parser_logger.filters:
                parser_logger.addFilter(_PRESCAN_FILTER)
    _PRESCAN.active = True
    try:
        yield
    finally:
        _PRESCAN.active = False


class UnusableInputError(ValueError):
    """The file cannot be used whatever its build; raised before the build decision."""


def structural_only_message(stats: VcfQualityStats) -> str:
    """The error for a VCF of structural-variant / copy-number records only."""
    kinds = ", ".join(f"<{kind}>" for kind, _n in stats.symbolic_kinds.most_common(4))
    return (
        f"This VCF contains only structural-variant / copy-number records ({kinds}; "
        f"{stats.symbolic_rows:,} record(s)). Allelix annotates SNVs and indels; "
        "for Sequencing.com, upload the snp-indel.genome.vcf.gz file."
    )


def reject_unusable_input(parser: GenotypeParser, file_path: Path, records: int) -> None:
    """Raise ``UnusableInputError`` when a read that yielded ``records`` found nothing usable.

    That is a VCF of structural-variant / copy-number records only, or a
    file with no genotype records. A VCF whose records were all dropped by
    the quality gate is not empty and passes. Used before the build
    decision (pre-scan) and after a read the build did not need.
    """
    if records:
        return
    stats = getattr(parser, "quality_stats", None)
    if stats is not None and stats.symbolic_rows:
        raise UnusableInputError(structural_only_message(stats))
    if stats is not None and (stats.filtered_rows or stats.low_quality_calls):
        return  # records exist but the quality gate dropped them; not an empty file
    raise UnusableInputError(f"{file_path.name} contains no genotype records.")


def decide_file_build(
    parser: GenotypeParser,
    file_path: Path,
    metadata: GenotypeMetadata,
    *,
    override: str | None,
    variants: list[Variant] | None = None,
) -> BuildDecision:
    """Decide ``file_path``'s build (ADR-0044); the one entry point for every command.

    Contig lengths decide a VCF without reading its records. Otherwise the
    known-SNP positions of the whole file are collected first (a pre-scan
    through the parser, or ``variants`` when the caller already holds
    them), so evidence far into a chromosome-sorted file counts. Raises
    ``BuildDecisionError`` when the run must stop.

    "Wrong file" errors come first (ADR-0044): when the scan finds no
    genotype record at all, ``UnusableInputError`` says so (a VCF of
    structural-variant / copy-number records only, or a file with no data
    rows) instead of asking for a build the file could not use anyway.
    When the header alone stops the decision (contig lengths that are
    mixed, unsupported or against ``assembly=``), the records are read up
    to the first usable one before the stop is raised.
    """
    kind = KIND_VCF if parser.name == "vcf" else KIND_ARRAY
    lengths = (
        contig_length_verdict(dict(metadata.get("contig_lengths") or {}))
        if kind == KIND_VCF
        else None
    )
    file_label = normalize_build_label(metadata.get("build"))
    checked = False  # whether the records were already checked for a wrong file
    if override or (lengths is not None and lengths.status != "none"):
        seeds = NO_SEEDS
    elif variants is not None:
        reject_unusable_input(parser, file_path, len(variants))
        checked = True
        seeds = detect_build(variants)
    else:
        count = [0]

        def _counted() -> Iterator[Variant]:
            for variant in parser.parse(file_path):
                count[0] += 1
                yield variant

        with _quiet_prescan():
            seeds = detect_build(_counted())
        reject_unusable_input(parser, file_path, count[0])
        checked = True
    try:
        return decide_build(
            file_name=file_path.name,
            kind=kind,
            override=override,
            file_label=file_label,
            default=DEFAULT_BUILD if kind == KIND_ARRAY else None,
            lengths=lengths,
            seeds=seeds,
        )
    except BuildDecisionError:
        if not checked:
            _reject_without_a_record(parser, file_path, variants)
        raise


def _reject_without_a_record(
    parser: GenotypeParser, file_path: Path, variants: list[Variant] | None
) -> None:
    """Before a stop on header evidence: a wrong file gets its own error.

    Reads up to the first usable genotype record (normally the first
    record); only a file without one is read to the end.
    """
    if variants is not None:
        reject_unusable_input(parser, file_path, len(variants))
        return
    stream = parser.parse(file_path)
    try:
        with _quiet_prescan():
            first = next(stream, None)
    finally:
        close = getattr(stream, "close", None)
        if close is not None:
            close()
    reject_unusable_input(parser, file_path, 0 if first is None else 1)


@dataclass
class StreamBuild:
    """What ``with_effective_build`` learnt while streaming, for the build notices.

    ``late`` is the known-SNP evidence met during the stream when contig
    lengths decided (None otherwise); ``mt_skipped`` counts dropped MT
    records. ``diagnostics()`` gives the same record ``analyze`` reports.
    """

    decision: BuildDecision
    late: BuildDetectionResult | None = None
    mt_skipped: int = 0

    def diagnostics(self) -> BuildDiagnostics:
        return BuildDiagnostics.from_decision(self.decision, self.late, self.mt_skipped)


def with_effective_build(
    variants: Iterable[Variant],
    *,
    decision: BuildDecision,
    file_name: str,
    record: StreamBuild | None = None,
) -> Iterator[Variant]:
    """Yield ``variants`` with ``build`` set to ``decision.build`` (ADR-0044).

    For command paths that need the build before rsID-less resolution
    (``stats``, ``extract``, ``export plink``). Order is kept. MT records
    are dropped when chrM has the hg19 length. When contig lengths decided,
    the known-SNP records met on the way are checked after the last
    variant; a contradiction raises ``BuildDecisionError``, so callers must
    not finalise any output before the iterator is exhausted. ``record``, when
    given, receives the late evidence (including a tie that let the run
    proceed) and the MT count, so the caller can show the build notices.
    """
    track = decision.source == SOURCE_CONTIG_LENGTHS
    seen: list[Variant] = []
    for variant in variants:
        if decision.skip_mt and variant.chromosome == "MT":
            if record is not None:
                record.mt_skipped += 1
            continue
        if track and known_snp_seed(variant) is not None:
            seen.append(variant)
        yield replace(variant, build=decision.build)
    if track:
        late = detect_build(seen)
        check_late_seeds(decision, late, file_name)
        if record is not None:
            record.late = late


__all__ = [
    "BUILD_GRCH37",
    "BUILD_GRCH38",
    "AnalysisResult",
    "BuildDiagnostics",
    "decide_file_build",
    "rollup_gwas_duplicates",
    "run_analysis",
    "with_effective_build",
]

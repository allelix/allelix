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
import re
import threading
from dataclasses import dataclass, field, replace
from itertools import chain
from typing import TYPE_CHECKING, TypedDict, cast

from allelix.annotators._resolution import resolve_batch
from allelix.annotators.base import variant_locus
from allelix.databases.manager import is_pathogenic_aggregate
from allelix.models import (
    DEFAULT_BUILD,
    DI_ALLELES,
    HALF_CALL_CODE,
    HOMOZYGOUS,
    Variant,
    half_call_genotype,
    report_rank_key,
)
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.base import ArrayInputStats
from allelix.reports._sources import ReportWarning
from allelix.reports.association_text import STRAND_CHECKED_SOURCES
from allelix.reports.clinvar_text import (
    PATHOGENIC_SIGNIFICANCES,
    common_allele_count,
    scored_down_count,
)
from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    KIND_ARRAY,
    KIND_VCF,
    NO_SEEDS,
    RCRS_LENGTH,
    SOURCE_OVERRIDE,
    BuildDecision,
    BuildDecisionError,
    BuildDetectionResult,
    check_late_seeds,
    contig_length_verdict,
    decide_build,
    detect_build,
    known_snp_seed,
    needs_late_seeds,
    normalize_build_label,
)
from allelix.utils.plausibility import (
    DEFAULT_THRESHOLDS,
    GNOMAD_NOT_RUN,
    MULTI_HOMOZYGOUS_PLP_THRESHOLD,
    NO_ALLELE,
    NOT_ASSESSED,
    PlausibilityThresholds,
    classify,
)
from allelix.utils.plausibility import HALF_CALL as REASON_HALF_CALL

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.models import Annotation
    from allelix.parsers.base import GenotypeMetadata, GenotypeParser
    from allelix.parsers.vcf import Locus, LocusState, VcfParser, VcfQualityStats
    from allelix.reports._sources import SourceRecord
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
    `mt_skipped` counts MT records left out because chrM is declared with
    the hg19 length, or (`mt_ucsc_name`) because on GRCh37 a record names
    the mitochondrion chrM or M and that name is not declared with the rCRS
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
    mt_ucsc_name: bool = False

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
        ucsc_mito_dropped: int = 0,
    ) -> BuildDiagnostics:
        """Diagnostics for ``decision``; ``seeds`` replaces its seed evidence when given.

        ``ucsc_mito_dropped`` counts the MT records the VCF parser left out by
        their UCSC name (``chrM``, ``M``); they join ``mt_skipped``.
        """
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
            mt_skipped=mt_skipped + ucsc_mito_dropped,
            tied_builds=evidence.tied,
            tied_votes=evidence.tied_votes,
            mt_ucsc_name=ucsc_mito_dropped > 0,
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


# Panel coverage states (ADR-0051), in report order. Each requested rsID is in
# exactly one.
PANEL_ANNOTATED = "annotated"
PANEL_TYPED_NO_ANNOTATION = "typed_no_annotation"
PANEL_NO_CALL = "no_call"
PANEL_NOT_ON_CHIP = "not_on_chip"
PANEL_STATES = (PANEL_ANNOTATED, PANEL_TYPED_NO_ANNOTATION, PANEL_NO_CALL, PANEL_NOT_ON_CHIP)
_PANEL_STATE_LABELS = {
    PANEL_ANNOTATED: "annotated",
    PANEL_TYPED_NO_ANNOTATION: "typed (no annotation)",
    PANEL_NO_CALL: "no call",
    PANEL_NOT_ON_CHIP: "not on chip",
}


def is_vcf_parser(parser_name: str) -> bool:
    """True when the input was read as a VCF or gVCF (it has no chip)."""
    return parser_name.lower() in {"vcf", "gvcf", "vcf parser"}


def panel_state_label(state: str, *, vcf: bool) -> str:
    """A panel state for display; ``not_on_chip`` is "not in file" for VCF input."""
    if vcf and state == PANEL_NOT_ON_CHIP:
        return "not in file"
    return _PANEL_STATE_LABELS[state]


def panel_reason_label(reason: str | None) -> str:
    """Plain-language explanation of an unusable explicit VCF panel record."""
    return {
        "no_call": "Caller no-call",
        "filtered": "Caller FILTER rejected",
        "low_gq_call": "Low genotype quality or zero depth",
        "low_gq_block": "Low genotype quality or zero depth",
        "spanning_deletion": "Spanning deletion not assessed",
        "symbolic_allele": "Symbolic allele not assessed",
        "excluded_mt": "MT excluded by existing build eligibility",
        "excluded_ucsc_mito": "UCSC mitochondrial reference excluded",
        "unusable_record": "No usable sample genotype",
        "reference_block_or_symbolic": "Reference block or symbolic record not assessed",
    }.get(reason or "", reason or "")


def gene_tokens(value: str, *, source: str = "") -> frozenset[str]:
    """Exact symbols in list fields; spaced GWAS flanking/interaction separators only."""
    separator = r"[;,]|\s+-\s+|\s+[xX]\s+" if source == "gwas" else r"[;,]"
    return frozenset(
        token.strip().upper() for token in re.split(separator, value) if token.strip()
    )


class PanelVariant(TypedDict):
    """One requested panel rsID and its state."""

    rsid: str
    gene: str | None
    state: str
    genotype: str | None
    annotation_count: int
    reason: str | None


class PanelCoverage(TypedDict):
    """Shape returned by AnalysisResult.panel_coverage().

    Same shape lands in the JSON report under ``panel_coverage`` and
    drives the HTML and terminal panel coverage. ``counts`` has one entry
    per state in ``PANEL_STATES``, and they sum to ``requested``.
    """

    requested: int
    counts: dict[str, int]
    variants: list[PanelVariant]


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
    # Per-rsID genotype as the reports display it, collected only for a
    # supplied panel. VCF observations also include explicit reference
    # calls omitted from annotation. None on the default analyze path.
    panel_genotypes: dict[str, str] | None = None
    # Panel rsIDs whose call in the input is a no-call (either allele
    # missing), or an explicit VCF record rejected by quality/eligibility.
    # None when no panel was supplied.
    no_call_panel_rsids: frozenset[str] | None = None
    # Bounded, explicit VCF rejection/exclusion reasons; arrays keep their existing calls.
    panel_reasons: dict[str, str] | None = None
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
    # ADR-0051: what the parse of an array file saw, and the array calls
    # ClinVar can never match at an aggregate-P/LP site
    # (``i_probes_at_plp_loci``, ``di_calls_at_plp_rsids``). None for VCF
    # input; ``probe_classes`` is also None when ClinVar did not run.
    array_quality: ArrayInputStats | None = None
    probe_classes: dict[str, int] | None = None
    # How the input's genotypes were produced, when the file declares it
    # (GenotypeMetadata ``method``, e.g. "low-pass-wgs-imputed" for
    # MyHeritage low-pass WGS). Empty when not declared.
    input_method: str = ""
    # per high-value rsID, what a VCF says at its position on the
    # effective build: (called | reference | no_call, reason). A gVCF
    # position no record covers is (no_call, no_coverage). Absent rsIDs
    # are unknown (plain VCF, or no placement). None for non-VCF input.
    high_value_vcf_states: dict[str, LocusState] | None = None
    # ADR-0045 disclosure: ClinVar Pathogenic / Likely pathogenic rows
    # scored down to ClinVar's own view of the variant, and ClinVar
    # rows not shown because their classification has no score.
    clinvar_scored_down: int = 0
    clinvar_unrecognised: int = 0
    # ADR-0052 disclosure: zero-star ClinVar P/LP rows on a common allele,
    # scored below the default floor.
    clinvar_common_allele: int = 0
    # Per strand-checked source that ran (GWAS, SNPedia), what was
    # not assessed because its strand was not settled, by reason. GWAS
    # counts association rows, SNPedia counts sites.
    strand_abstentions: dict[str, dict[str, int]] = field(default_factory=dict)
    # Every registered source with its part in the run (used,
    # ready_not_used, or skipped with a reason), and report-level
    # warnings. Set by the CLI after the run; the renderers only read
    # them.
    sources: list[SourceRecord] = field(default_factory=list)
    report_warnings: list[ReportWarning] = field(default_factory=list)

    def panel_coverage(
        self, filtered_annotations: list[Annotation] | None = None
    ) -> PanelCoverage | None:
        """Every requested panel rsID in one of four states (ADR-0051).

        - ``annotated``: at least one shown annotation. It wins over a
          no-call, so a half call with a finding is annotated.
        - ``typed_no_annotation``: in the input and called, nothing shown.
        - ``no_call``: present but unusable (missing, rejected or excluded).
        - ``not_on_chip``: not in the input ("not in file" for VCF input,
          where an absent rsID may be reference at that position).

        "Shown" is ``filtered_annotations`` when the renderer passes its
        post-filter list, so an rsID whose annotations are all below the
        display floor is ``typed_no_annotation``, not annotated. Without
        it, ``panel_annotated_rsids`` (every annotation of the run) does.

        Returns None when no panel was supplied.
        """
        if self.panel_rsids is None:
            return None
        panel_rsids = self.panel_rsids
        genotyped = self.genotyped_panel_rsids or frozenset()
        no_calls = self.no_call_panel_rsids or frozenset()
        genotypes = self.panel_genotypes or {}
        if filtered_annotations is not None:
            shown = [a for a in filtered_annotations if a.rsid in panel_rsids]
            annotated = frozenset(a.rsid for a in shown)
        else:
            annotated = self.panel_annotated_rsids or frozenset()
            shown = [a for a in self.annotations if a.rsid in annotated]
        by_rsid: dict[str, list[Annotation]] = {}
        for a in shown:
            by_rsid.setdefault(a.rsid, []).append(a)
        variants: list[PanelVariant] = []
        for rsid in panel_rsids:
            hits = by_rsid.get(rsid, [])
            if rsid in annotated:
                state = PANEL_ANNOTATED
            elif rsid not in genotyped:
                state = PANEL_NOT_ON_CHIP
            elif rsid in no_calls:
                state = PANEL_NO_CALL
            else:
                state = PANEL_TYPED_NO_ANNOTATION
            variants.append(
                {
                    "rsid": rsid,
                    "gene": next((a.gene for a in hits if a.gene), None),
                    "state": state,
                    "genotype": genotypes.get(rsid) if rsid in genotyped else None,
                    "annotation_count": len(hits),
                    "reason": (self.panel_reasons or {}).get(rsid),
                }
            )
        variants.sort(key=lambda v: (PANEL_STATES.index(v["state"]), v["rsid"]))
        counts = dict.fromkeys(PANEL_STATES, 0)
        for v in variants:
            counts[v["state"]] += 1
        return {"requested": len(panel_rsids), "counts": counts, "variants": variants}

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

        `source_min_magnitudes` sets the floor for specific sources
        (e.g. ``{"gwas": 9.0, "snpedia": 2.0}``). GWAS rows must meet the
        larger of it and ``min_magnitude`` (ADR-0024); must-include GWAS
        rsIDs are exempt from the GWAS floor and meet ``min_magnitude``.
        For other sources the entry IS the floor: it can raise OR lower
        the global ``min_magnitude``. Sources without an entry use the
        global floor.

        `genes` and `rsids` combine with OR: when either is provided, an
        annotation passes if it matches the gene set OR the rsid set.
        Empty collections (vs None) mean "match nothing" — an empty
        filter file produces an empty report.
        """
        # Materialize once so generator queries work across all source fields.
        gene_queries = tuple(genes) if genes is not None else None
        gene_sets: dict[str, set[str]] = {}
        rsid_set = {r.lower() for r in rsids} if rsids is not None else None
        out: list[Annotation] = []
        for a in self.annotations:
            floor = min_magnitude
            if source_min_magnitudes and a.source in source_min_magnitudes:
                source_floor = source_min_magnitudes[a.source]
                if a.source != "gwas":
                    floor = source_floor
                elif not a.is_must_include:
                    floor = max(min_magnitude, source_floor)
            if a.magnitude < floor:
                continue
            if category is not None and a.category != category:
                continue
            if gene_queries is not None or rsid_set is not None:
                gene_match = False
                if gene_queries is not None:
                    if a.source not in gene_sets:
                        gene_sets[a.source] = {
                            token
                            for query in gene_queries
                            for token in gene_tokens(query, source=a.source)
                        }
                    gene_match = bool(
                        gene_sets[a.source] & gene_tokens(a.gene or "", source=a.source)
                    )
                rsid_match = rsid_set is not None and a.rsid.lower() in rsid_set
                if not gene_match and not rsid_match:
                    continue
            out.append(a)
        out.sort(key=report_rank_key)
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
    survivors.sort(key=report_rank_key)
    return survivors


_IDENTITY_BUILDS: frozenset[str] = frozenset({BUILD_GRCH37, BUILD_GRCH38})

EnrichmentKey = tuple[str, int, str, str]


def watch_high_value(
    parser: GenotypeParser,
    hv_set: set[str],
    identity: RsidIdentityAnnotator | None,
    extra_loci: set[Locus] | None = None,
    *,
    rows: bool = False,
) -> dict[str, dict[str, Locus]] | None:
    """Ask a VCF parser to record what the file says at high-value sites.

    Positions come from the identity table on both builds (the effective
    build is not known until streaming ends); rsIDs are watched directly
    as well, and so are ``extra_loci`` (ClinVar's reference-carried
    sites). ``rows`` keeps the reference records for the annotated VCF
    (``VcfParser.watch``). Returns the high-value placements, or None when there is no
    high-value list or the input is not a VCF.
    """
    from allelix.parsers.vcf import VcfParser

    if not isinstance(parser, VcfParser) or not (hv_set or extra_loci):
        return None
    placements: dict[str, dict[str, Locus]] = {}
    if hv_set and identity is not None and identity.is_ready():
        placements = identity.loci_of_rsids(hv_set)
    loci = {locus for p in placements.values() for locus in p.values()} | (extra_loci or set())
    parser.watch(loci, hv_set, rows=rows)
    return placements if hv_set else None


def reference_carried_calls(
    stats: VcfQualityStats | None,
    sites: dict[Locus, tuple[str, str, bool]],
    build: str,
    identity: RsidIdentityAnnotator | None = None,
) -> list[Variant]:
    """Hom-ref calls at ClinVar's reference-carried sites, as that homozygous genotype.

    At a site where ``build``'s reference carries a ClinVar row's allele (T
    at F5 rs6025 on GRCh37; the normal G at CASP12 rs497116), a ``0/0`` call
    or a reference block that passes the GQ gate is homozygous for it, which
    no record states. Each becomes a REF/REF call for the annotators. Sets
    the two counts on ``stats``: those calls, and on a plain VCF the
    variant sites (not normal-allele ones) no record covers: a homozygous
    carrier cannot be told from no coverage there.
    """
    from allelix.parsers.vcf import REFERENCE

    if stats is None or not sites:
        return []
    calls: list[Variant] = []
    unseen = 0
    identity_rows = (
        identity.alleles_at_loci(build, set(sites))
        if identity is not None and identity.is_ready()
        else {}
    )
    for (chrom, pos), (allele, rsid, normal_allele) in sorted(sites.items()):
        state = stats.locus_states.get((chrom, pos))
        if state is not None and state[0] == REFERENCE:
            matching = []
            matching_count = 0
            known_refs = {
                ref.upper()
                for ref, _alt, _placed_rsid in identity_rows.get((chrom, pos), [])
                if len(ref) == 1 and ref.upper() in "ACGT"
            }
            conflict = len(known_refs | {allele}) > 1
            for source, count in stats.reference_calls.get((chrom, pos), {}).items():
                if conflict:
                    stats.reference_source_conflicts += count
                elif source.ref is not None and source.ref != allele:
                    stats.clinvar_ref_mismatches += count
                else:
                    matching.append(source)
                    matching_count += count
            if len({source.haploid for source in matching}) > 1:
                stats.reference_ploidy_conflicts += matching_count
                continue
            if matching:
                calls.append(
                    Variant(
                        rsid,
                        chrom,
                        pos,
                        allele,
                        allele,
                        build=build,
                        ref=allele,
                        ref_from_input=any(source.ref is not None for source in matching),
                        haploid=matching[0].haploid,
                    )
                )
        elif state is None and not stats.is_gvcf and not normal_allele:
            unseen += 1
    stats.reference_carried_calls = len(calls)
    stats.reference_carried_not_in_file = unseen
    return calls


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


_AUTOSOMES = frozenset(str(n) for n in range(1, 23))


def _mark_haploid(annotations: list[Annotation], variants: list[Variant]) -> None:
    """Mark the annotations of calls the file declares haploid (ADR-0050).

    Matching used the called allele twice; the report shows it once, as a
    hemizygous call. An annotation finds its call by the input record it
    names (``record_key``: a VCF allele matched after trimming, whose locus
    is the trimmed position, ADR-0047), else by locus, else by rsID when it
    has no locus.
    """
    haploid = [v for v in variants if v.haploid]
    if not haploid:
        return
    loci = {variant_locus(v) for v in haploid}
    records = {(v.chromosome, v.position, v.ref) for v in haploid if v.ref is not None}
    rsids = {v.rsid for v in haploid if v.rsid}
    for a in annotations:
        if a.record_key is not None:
            matched = a.record_key[:3] in records
        elif a.locus is not None:
            matched = a.locus in loci
        else:
            matched = a.rsid in rsids
        if matched:
            a.haploid = True
            gm = a.genotype_match
            half = len(gm) // 2
            first, second = gm.split("/", 1) if "/" in gm else (gm[:half], gm[half:])
            if first == second:
                a.genotype_match = first


def _flag_plausibility(
    annotations: list[Annotation],
    keyed: list[tuple[Annotation, EnrichmentKey]],
    thresholds: PlausibilityThresholds,
    *,
    gnomad_ran: bool,
) -> None:
    """Stamp every annotation's plausibility and, when not assessed, why (ADR-0041, ADR-0050).

    Runs after gnomAD has stamped ``allele_frequency``. The copy count is
    of the same allele the frequency describes, counted in the user's
    alleles in the annotation's orientation; a call the file declares
    haploid is one copy.
    """
    keyed_ids = {id(a) for a, _key in keyed}
    for a in annotations:
        allele = _enrichment_allele(a)
        if not gnomad_ran:
            a.plausibility, a.plausibility_reason = NOT_ASSESSED, GNOMAD_NOT_RUN
        elif a.zygosity_code == HALF_CALL_CODE:
            # a half call's copy count is unknown (ADR-0047).
            a.plausibility, a.plausibility_reason = NOT_ASSESSED, REASON_HALF_CALL
        elif allele is None or a.locus is None or id(a) not in keyed_ids:
            a.plausibility, a.plausibility_reason = NOT_ASSESSED, NO_ALLELE
        else:
            copies = 1 if a.haploid else a.carried_alleles.count(allele)
            a.plausibility, a.plausibility_reason = classify(
                a.allele_frequency, copies, a.locus[0], thresholds
            )


# ADR-0052: the per-submission review statuses ClinVar shows as zero stars.
ZERO_STAR_REVIEW_STATUSES: frozenset[str] = frozenset(
    {"no assertion criteria provided", "no classification provided"}
)
# ACMG/AMP BA1: a frequency of 5% or more is stand-alone benign evidence.
COMMON_ALLELE_AF = 0.05
# Below the default display floor of 5.0, with Conflicting's cap (ADR-0045).
COMMON_ALLELE_MAGNITUDE = 4.0


def _is_common_allele_candidate(a: Annotation) -> bool:
    """A zero-star ClinVar P/LP row above the cap, named by an rsID (ADR-0052)."""
    return (
        a.source == "clinvar"
        and a.significance in PATHOGENIC_SIGNIFICANCES
        and a.review_status.strip().lower() in ZERO_STAR_REVIEW_STATUSES
        and a.magnitude > COMMON_ALLELE_MAGNITUDE
        and a.rsid.startswith("rs")
        and a.locus is not None
        and a.locus_ref is not None
        and bool(a.alt)
        and a.alt != a.locus_ref
    )


def _clinvar_allele_keys(
    candidates: list[Annotation],
    build: str,
    clinvar: ClinVarAnnotator,
    identity: RsidIdentityAnnotator,
) -> dict[int, EnrichmentKey]:
    """The GRCh38 gnomAD key of each candidate's ClinVar allele, by ``id()`` (ADR-0052).

    The key is ClinVar's own forward-strand REF and ALT at the rsID's
    GRCh38 locus in the identity table, whatever the user's call. The
    palindromic flag does not apply: it guards the orientation of a call,
    and here no call is translated. A candidate gets no key on any doubt:
    the table places the rsID at more than one locus on either build; one
    ClinVar build has no readable cache, or lacks the same REF and ALT at
    that locus; the table has no row with that REF and ALT on both builds;
    or the row is not at the table's locus on the analysis build.
    """
    rsids = {a.rsid for a in candidates}
    records: dict[str, dict[str, set[tuple[str, int, str, str]]]] = {}
    for b in (BUILD_GRCH37, BUILD_GRCH38):
        found = clinvar.records_of_rsids(rsids, b)
        if found is None:
            return {}
        records[b] = found
    if build not in records:
        return {}
    placements = identity.placements(rsids)
    keys: dict[int, EnrichmentKey] = {}
    for a in candidates:
        rows = placements.get(a.rsid, [])
        loci = {
            BUILD_GRCH37: {(chrom, pos37) for chrom, pos37, *_rest in rows},
            BUILD_GRCH38: {(row[0], row[4]) for row in rows},
        }
        if len(loci[BUILD_GRCH37]) != 1 or len(loci[BUILD_GRCH38]) != 1:
            continue
        on_build = {b: next(iter(places)) for b, places in loci.items()}
        ref, alt = cast("str", a.locus_ref), a.alt
        if a.locus != on_build[build]:
            continue
        if not all((*on_build[b], ref, alt) in records[b].get(a.rsid, set()) for b in records):
            continue
        if (ref, alt, ref, alt) not in {
            (r37, a37, r38, a38) for _c, _p, r37, a37, _q, r38, a38 in rows
        }:
            continue
        chrom, pos38 = on_build[BUILD_GRCH38]
        keys[id(a)] = (chrom, pos38, ref, alt)
    return keys


def _score_down_common_zero_star(
    annotations: list[Annotation],
    build: str,
    clinvar: ClinVarAnnotator | None,
    identity: RsidIdentityAnnotator | None,
    gnomad: GnomadAnnotator | None,
) -> None:
    """Score a zero-star ClinVar P/LP row on a common allele below the default floor (ADR-0052).

    The frequency is gnomAD's for ClinVar's allele itself, looked up at the
    rsID's GRCh38 locus (``_clinvar_allele_keys``), not the frequency
    enrichment gave the user's call. Without ClinVar on both builds, the
    identity table or gnomAD, or on any doubt about the allele, a row is
    unchanged; so is a row whose submission has one star or more.
    """
    candidates = [a for a in annotations if _is_common_allele_candidate(a)]
    if not candidates or clinvar is None or identity is None or gnomad is None:
        return
    if not identity.is_ready() or not gnomad.is_ready():
        return
    keys = _clinvar_allele_keys(candidates, build, clinvar, identity)
    freqs = gnomad.bulk_lookup_by_position(set(keys.values())) if keys else {}
    for a in candidates:
        af = freqs.get(keys[id(a)]) if id(a) in keys else None
        if af is None:
            continue
        a.clinvar_allele_af = af
        if af >= COMMON_ALLELE_AF:
            a.magnitude = COMMON_ALLELE_MAGNITUDE
            a.common_allele_capped = True


def _multiple_homozygous_plp(
    annotations: list[Annotation], thresholds: PlausibilityThresholds
) -> ReportWarning | None:
    """The ADR-0050 warning: several rare homozygous ClinVar P/LP calls.

    Counts distinct autosomal variants, before any display filter, called
    homozygous for an allele whose ClinVar aggregate is Pathogenic or
    Likely pathogenic, and whose gnomAD frequency is below
    ``homozygous_max_af`` or missing.
    """
    variants = {
        a.locus
        for a in annotations
        if a.source == "clinvar"
        and a.locus is not None
        and a.locus[0] in _AUTOSOMES
        and a.zygosity_code == HOMOZYGOUS
        and is_pathogenic_aggregate(a.clinvar_aggregate)
        and (a.allele_frequency is None or a.allele_frequency < thresholds.homozygous_max_af)
    }
    if len(variants) < MULTI_HOMOZYGOUS_PLP_THRESHOLD:
        return None
    return ReportWarning(
        "multiple_homozygous_plp",
        f"{len(variants)} variants are called homozygous for an allele ClinVar classifies "
        "as Pathogenic or Likely pathogenic, and are rare in gnomAD or have no gnomAD "
        "frequency. Several such calls in one file usually point to a genotyping or file "
        "problem rather than to findings; confirm any of them with an independent test.",
        len(variants),
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


def _with_trimmed_matches(
    stats: VcfQualityStats | None, clinvar_counts: dict[str, int]
) -> VcfQualityStats | None:
    """The VCF stats with the ClinVar annotator's trimmed matches and REF mismatches."""
    if stats is not None:
        stats.trimmed_matches = clinvar_counts.get("trimmed_matches", 0)
        stats.clinvar_ref_mismatches += clinvar_counts.get("ref_mismatches", 0)
    return stats


def _probe_classes(
    clinvar: Annotator | None, clinvar_counts: dict[str, int], build: str
) -> dict[str, int] | None:
    """ADR-0051: the array calls ClinVar never matches at aggregate-P/LP sites."""
    if clinvar is None or not clinvar.serves_build(build):
        return None
    return {
        "i_probes_at_plp_loci": clinvar_counts.get("i_probes_at_plp_loci", 0),
        "di_calls_at_plp_rsids": clinvar_counts.get("di_calls_at_plp_rsids", 0),
    }


def _record_for(
    a: Annotation,
    by_locus: dict[tuple[str, int], list[Variant]],
    by_rsid: dict[str, Variant],
) -> Variant | None:
    """The input record an annotation belongs to, for ``--vcf-out``.

    Found by the annotation's locus: among the records there, those that
    carry the annotation's allele (all of them when none does, as a GWAS
    effect allele need not be carried), then the one with the annotation's
    rsID when several remain. An annotation without a locus falls back to
    its rsID.
    """
    if a.locus is None:
        return by_rsid.get(a.rsid)
    at_locus = by_locus.get(a.locus, [])
    hits = [v for v in at_locus if a.alt in (v.allele1, v.allele2)] or at_locus
    if len(hits) > 1:
        hits = [v for v in hits if v.rsid and v.rsid == a.rsid] or hits
    return hits[0] if hits else None


@contextlib.contextmanager
def _cleanup_on_exit() -> Iterator[contextlib.ExitStack]:
    """Attempt every registered cleanup without replacing a primary error."""
    stack = contextlib.ExitStack()
    try:
        yield stack
    except BaseException:
        # ExitStack still attempts every close. Cleanup must not replace the
        # computation failure (including an interrupt) already being raised.
        with contextlib.suppress(BaseException):
            stack.close()
        raise
    else:
        stack.close()


@contextlib.contextmanager
def _owned_annotators(owners: Iterable[Annotator | None]) -> Iterator[dict[int, Annotator]]:
    """Close each per-run owner once, preserving a primary computation error."""
    with _cleanup_on_exit() as stack:
        entered: dict[int, Annotator] = {}
        for owner in owners:
            if owner is not None and id(owner) not in entered:
                entered[id(owner)] = stack.enter_context(owner)
        yield entered


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
    resource_owners: Iterable[Annotator] = (),
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

    All active annotators, optional providers and ``resource_owners`` are
    owned from before metadata/build work through result construction.
    Each distinct object closes once per run, including on failure.
    Ownership does not activate an annotator or acquire its resources.
    Direct callers must supply callback-only providers in ``resource_owners``;
    arbitrary callback closures are not inspected. The CLI supplies its
    complete registry, including skipped providers used by other sources.
    """
    from allelix.parsers.vcf import VcfParser

    with _owned_annotators(
        chain(annotators, (gnomad, alphamissense, cadd, identity), resource_owners)
    ) as owners:
        bound = [owners[id(a)] for a in annotators]
        metadata = parser.get_metadata(file_path)
        decision = decide_file_build(parser, file_path, metadata, override=build_override)
        is_array = not isinstance(parser, VcfParser)
        di_calls = 0

        annotations: list[Annotation] = []
        hv_variants: list[Variant] = []
        hv_set: set[str] = high_value_rsids or set()
        # rsIDs the input file actually carries that are in the
        # supplied panel. Only collected when ``panel_rsids`` is non-empty
        # (so a normal analyze run has zero overhead). Filtered at insert
        # so a 600k-variant input contributes ~panel-size memory, not
        # 600k rsids.
        genotyped_panel: set[str] = set()
        # First-occurrence genotype for panel members. VCF explicit observations
        # replace the whole tuple below, including non-yielded reference/FILTER rows.
        panel_genotypes_map: dict[str, str] = {}
        no_call_panel: set[str] = set()
        # Per-run diagnostic totals for ref-None and
        # rsidless variants reaching the annotator phase. Used by the CLI
        # to emit the "strand-aware inactive" info line and the
        # "GRCh38 rsID-less undercount" warning.
        no_ref_total = 0
        rsidless_total = 0
        total = 0
        effective_build = decision.build
        # Known-SNP records met during the pass; only kept when contig lengths
        # decided (with or without --build), so check_late_seeds can compare
        # them with the decision.
        late_seeds: list[Variant] = []
        track_seeds = needs_late_seeds(decision)
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
                        # The first occurrence wins: "A/G", "C/C", "A/-"
                        # or "-/-" as the file wrote it, and a VCF half
                        # call as the annotations show it ("T/.").
                        if v.rsid not in panel_genotypes_map:
                            panel_genotypes_map[v.rsid] = half_call_genotype(v) or v.genotype
                            if v.is_no_call:
                                no_call_panel.add(v.rsid)
            # Per-batch diagnostic counters. Tallied AFTER
            # the ref-population and resolver passes so the counts
            # reflect the final state annotators see — not the raw input.
            nonlocal no_ref_total, rsidless_total
            no_ref_total += sum(1 for v in batch_buf if v.ref is None)
            rsidless_total += sum(1 for v in batch_buf if not v.rsid.startswith("rs"))
            batch_new: list[Annotation] = []
            for annotator in bound:
                batch_new.extend(annotator.batch_annotate(batch_buf))
            _mark_haploid(batch_new, batch_buf)
            annotations.extend(batch_new)
            # collect writer inputs for this batch. Only touches
            # variants whose ``ref`` is populated (VCF inputs; array
            # data with ref=None can't produce a well-formed VCF key).
            # Skips annotations with empty ``alt`` — GWAS rows keyed by
            # rsID alone can't be VCF-stamped without allele context.
            # An annotation finds its record by locus and carried ALT, so
            # records without an rsID (ADR-0046) do not collide.
            if collect_vcf_write_inputs:
                with_ref = [v for v in batch_buf if v.ref is not None]
                by_locus: dict[tuple[str, int], list[Variant]] = {}
                for v in with_ref:
                    by_locus.setdefault(variant_locus(v), []).append(v)
                by_rsid = {v.rsid: v for v in with_ref if v.rsid}
                for a in batch_new:
                    if not a.alt:
                        continue
                    if a.record_key is not None:
                        # ADR-0047: matched after trimming; stamp the
                        # allele as the record wrote it.
                        vcf_annotations_by_key.setdefault(a.record_key, []).append(a)
                        continue
                    record = _record_for(a, by_locus, by_rsid)
                    if record is None or record.ref is None:
                        continue
                    key = (record.chromosome, record.position, record.ref, a.alt)
                    vcf_annotations_by_key.setdefault(key, []).append(a)
            batch_buf.clear()

        def _accept(v: Variant) -> None:
            batch_buf.append(v)
            if len(batch_buf) >= _BATCH_SIZE:
                _flush()

        # watch high-value positions so gVCF reference blocks and
        # no-call records, which carry no rsIDs, can be classified.
        # A VCF says nothing at ClinVar's reference-carried sites for a
        # person homozygous for the variant allele the reference carries:
        # watch them too, and assess a hom-ref call there.
        rc_sites = (
            cast("ClinVarAnnotator", clinvar_resolver).reference_carried_sites(effective_build)
            if clinvar_resolver is not None and parser.name == "vcf"
            else {}
        )
        hv_placements = None
        with _cleanup_on_exit() as cleanup:
            if isinstance(parser, VcfParser):
                # Register cleanup BEFORE invoking fallible watch methods:
                # they may register state and then fail or be interrupted.
                if hv_set or rc_sites:
                    cleanup.callback(parser.watch, set(), set())
                cleanup.callback(parser.watch_panel, frozenset())
            hv_placements = watch_high_value(
                parser, hv_set, identity, set(rc_sites), rows=collect_vcf_write_inputs
            )
            if isinstance(parser, VcfParser):
                panel_loci: dict[Locus, set[str]] = {}
                if (
                    panel_rsids
                    and identity is not None
                    and identity.is_ready()
                    and effective_build in _IDENTITY_BUILDS
                ):
                    # The existing per-build API abstains when identity rows disagree
                    # on the locus; multiple alleles at one locus remain supported.
                    for rsid, (chrom, pos, _alleles) in identity.forward_sites(
                        set(panel_rsids), effective_build
                    ).items():
                        panel_loci.setdefault((chrom, pos), set()).add(rsid)
                parser.watch_panel(
                    panel_rsids or frozenset(),
                    {locus: frozenset(rsids) for locus, rsids in panel_loci.items()},
                    skip_mt=decision.skip_mt,
                )
            for variant in parser.parse(file_path):
                total += 1
                if is_array and not {variant.allele1, variant.allele2}.isdisjoint(DI_ALLELES):
                    di_calls += 1
                if decision.skip_mt and variant.chromosome == "MT":
                    mt_skipped += 1
                    continue
                if track_seeds and known_snp_seed(variant) is not None:
                    late_seeds.append(variant)
                _accept(replace(variant, build=effective_build))
        for call in reference_carried_calls(
            getattr(parser, "quality_stats", None), rc_sites, effective_build, identity
        ):
            if decision.skip_mt and call.chromosome == "MT":
                mt_skipped += 1
                cast("VcfParser", parser).quality_stats.reference_carried_calls -= 1
                continue
            _accept(call)
        _flush()
        if isinstance(parser, VcfParser) and panel_rsids:
            # Explicit observations include rows the annotation stream intentionally
            # omits. Replace the complete first-record tuple, never just its state.
            from allelix.parsers.vcf import NO_CALL

            for rsid, observation in parser.panel_observations.items():
                genotyped_panel.add(rsid)
                if observation.genotype is None:
                    panel_genotypes_map.pop(rsid, None)
                else:
                    panel_genotypes_map[rsid] = observation.genotype
                if observation.state == NO_CALL:
                    no_call_panel.add(rsid)
                else:
                    no_call_panel.discard(rsid)
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
        # ADR-0041: zygosity x gnomAD frequency where gnomAD was queried;
        # every other annotation says why it was not assessed (ADR-0050).
        _flag_plausibility(
            annotations, enrichment_keys, plausibility, gnomad_ran=gnomad in enriched_by
        )
        _score_down_common_zero_star(
            annotations,
            effective_build,
            cast("ClinVarAnnotator | None", clinvar_resolver),
            identity,
            gnomad,
        )
        warnings = [w for w in [_multiple_homozygous_plp(annotations, plausibility)] if w]

        annotators_used = [
            (a.name, a.version()) for a in annotators if a.serves_build(effective_build)
        ]
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
            result_no_calls: frozenset[str] | None = frozenset(no_call_panel)
        else:
            result_panel = None
            result_genotyped = None
            result_annotated = None
            result_genotypes = None
            result_no_calls = None

        vcf_write_inputs: VcfWriteInputs | None = None
        if collect_vcf_write_inputs:
            # Import here to avoid a top-level cycle with allelix.reports.vcf
            # (which reaches into the annotator registry for attribution).
            from allelix.reports.vcf import VcfWriteInputs as _VcfWriteInputs

            vcf_recovered_rsids = {
                (chrom, pos, ref, alt): rsid
                for rsid, (chrom, pos, ref, alt) in resolved_coords.items()
            }
            vcf_stats = getattr(parser, "quality_stats", None)
            vcf_write_inputs = _VcfWriteInputs(
                annotations_by_key=vcf_annotations_by_key,
                recovered_rsids=vcf_recovered_rsids,
                reference_rows=vcf_stats.reference_rows if vcf_stats is not None else {},
            )

        clinvar_counts = clinvar_resolver.run_counts() if clinvar_resolver is not None else {}
        # GWAS association rows and SNPedia sites whose strand was not settled.
        strand_abstentions = {
            a.name: a.run_counts()
            for a in annotators
            if a.name in STRAND_CHECKED_SOURCES and a.serves_build(effective_build)
        }
        array_quality: ArrayInputStats | None = None
        if is_array:
            # The parser's counts from this, its last, parse (ADR-0051).
            array_quality = parser.input_stats() or ArrayInputStats()
            array_quality.di_calls = di_calls

        return AnalysisResult(
            report_warnings=warnings,
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
            build_diagnostics=BuildDiagnostics.from_decision(
                decision, late, mt_skipped, getattr(parser, "ucsc_mito_dropped", 0)
            ),
            hv_variants=hv_variants,
            panel_rsids=result_panel,
            genotyped_panel_rsids=result_genotyped,
            panel_annotated_rsids=result_annotated,
            panel_genotypes=result_genotypes,
            no_call_panel_rsids=result_no_calls,
            panel_reasons={
                rsid: observation.reason
                for rsid, observation in parser.panel_observations.items()
                if observation.reason is not None
            }
            if isinstance(parser, VcfParser) and panel_rsids
            else None,
            no_ref_variant_count=no_ref_total,
            rsidless_variant_count=rsidless_total,
            vcf_write_inputs=vcf_write_inputs,
            vcf_quality=_with_trimmed_matches(
                getattr(parser, "quality_stats", None), clinvar_counts
            ),
            array_quality=array_quality,
            probe_classes=_probe_classes(clinvar_resolver, clinvar_counts, effective_build)
            if is_array
            else None,
            high_value_vcf_states=high_value_states(
                getattr(parser, "quality_stats", None), hv_placements, hv_set, effective_build
            ),
            clinvar_scored_down=scored_down_count(annotations),
            clinvar_unrecognised=clinvar_counts.get("unrecognised_rows", 0),
            clinvar_common_allele=common_allele_count(annotations),
            strand_abstentions=strand_abstentions,
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
    if stats is not None and stats.hom_ref_calls:
        return
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
    if hasattr(parser, "drop_ucsc_mito"):
        parser.drop_ucsc_mito = False  # set again below, after any pre-scan
    lengths = (
        contig_length_verdict(dict(metadata.get("contig_lengths") or {}))
        if kind == KIND_VCF
        else None
    )
    file_label = normalize_build_label(metadata.get("build"))
    checked = False  # whether the records were already checked for a wrong file
    if override == BUILD_GRCH36 or (lengths is not None and lengths.status != "none"):
        # Contig lengths decide without a pre-scan; their seeds are checked
        # during the pass (needs_late_seeds). --build grch36 is never checked
        # (ADR-0025). Any other --build is checked against the seeds.
        seeds = NO_SEEDS
    elif variants is not None:
        if not override:
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
        if not override:
            # Under --build the pre-scan only checks the forced build; a file
            # without a usable record runs as it did before the check.
            reject_unusable_input(parser, file_path, count[0])
            checked = True
    try:
        decision = decide_build(
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
    if hasattr(parser, "drop_ucsc_mito") and hasattr(parser, "rcrs_mito_names"):
        # On GRCh37, a record named as UCSC hg19 names the mitochondrion
        # (chrM, M) is the 16,571 bp sequence unless that same name is
        # declared with the rCRS length: on the pass it is left out and
        # counted. MT and chrMT records are GRCh37's own and are read. Not
        # needed when a declared hg19 length leaves out every MT record.
        parser.drop_ucsc_mito = (
            kind == KIND_VCF and decision.build == BUILD_GRCH37 and not decision.skip_mt
        )
        parser.rcrs_mito_names = frozenset(
            name
            for name, length in (metadata.get("contig_lengths") or {}).items()
            if normalize_chromosome(name) == "MT" and length == RCRS_LENGTH
        )
    return decision


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
    records; ``parser``, when given, adds the MT records a VCF parser left
    out by their UCSC name. ``diagnostics()`` gives the same record
    ``analyze`` reports.
    """

    decision: BuildDecision
    late: BuildDetectionResult | None = None
    mt_skipped: int = 0
    parser: GenotypeParser | None = None

    @property
    def ucsc_mito_dropped(self) -> int:
        """MT records the VCF parser left out by their UCSC name."""
        return getattr(self.parser, "ucsc_mito_dropped", 0)

    def diagnostics(self) -> BuildDiagnostics:
        return BuildDiagnostics.from_decision(
            self.decision, self.late, self.mt_skipped, self.ucsc_mito_dropped
        )


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
    track = needs_late_seeds(decision)
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

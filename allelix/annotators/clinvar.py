# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar annotator. Source-attributed pathogenicity calls (ADR-0003).

ADR-0021: per-build SQLite caches. ClinVar publishes separate VCFs for
GRCh37 and GRCh38, and the strand orientation of REF/ALT can invert
between builds for the ~0.4% of the genome where the reference
assembly was rebuilt. Carrier-rule matches (ADR-0007) MUST be done
against the same build the user's data is on. The annotator holds one
SQLite cache per build (`clinvar.GRCh37.sqlite`, `clinvar.GRCh38.sqlite`)
and dispatches per-variant by `variant.build`.
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor, variant_locus
from allelix.databases import manager as _manager_module
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.cache_state import (
    READY,
    CacheState,
    CacheStatus,
    cache_state,
    connect_readonly,
)
from allelix.databases.manager import (
    CONTRIBUTES_TO_AGGREGATE,
    OUTSIDE_AGGREGATE_CONFLICTING,
    aggregate_germline_term,
    download,
    fetch_remote_text,
    get_database_info,
    is_disease_classification,
    is_pathogenic_aggregate,
    load_clinvar_tsv,
    stamp_existing_clinvar_cache,
    verify_file_hash,
)
from allelix.models import DI_ALLELES, NO_CALL_MARKER, Annotation, diploid_label
from allelix.utils.allele import (
    resolve_rsid_at_locus,
    trim_allele_pair,
)
from allelix.utils.build_detect import BUILD_GRCH37, BUILD_GRCH38

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Callable, Iterable, Iterator

    from allelix.models import Variant

logger = logging.getLogger(__name__)

CLINVAR_SUPPORTED_BUILDS: tuple[str, ...] = ("GRCh37", "GRCh38")

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999

_MD5_HEX_RE = re.compile(r"^[0-9a-fA-F]{32}$")
_CHECKSUM_MAX_BYTES = 4096


def _fetch_checksum(source_url: str) -> str | None:
    """Read a bounded MD5 sidecar associated with this source's basename."""
    body = fetch_remote_text(source_url + ".md5", max_bytes=_CHECKSUM_MAX_BYTES)
    if not body or len(body.encode("utf-8")) > _CHECKSUM_MAX_BYTES:
        return None
    expected_name = source_url.rsplit("/", 1)[-1]
    digests: set[str] = set()
    for line in body.splitlines():
        parts = line.split(None, 1)
        if not parts:
            continue
        if not _MD5_HEX_RE.fullmatch(parts[0]):
            return None
        if len(parts) == 2:
            filename = parts[1].strip().removeprefix("*")
            if filename.rsplit("/", 1)[-1] != expected_name:
                return None
        digests.add(parts[0].lower())
    return digests.pop() if len(digests) == 1 else None


def clinvar_db_filename(build: str) -> str:
    """Per-build cache filename. Two coexisting SQLite files per data_dir."""
    return f"clinvar.{build}.sqlite"


def clinvar_record_name(build: str) -> str:
    """`database_versions` row name for a given build."""
    return f"clinvar.{build}"


# Allelix-derived magnitude scoring from ClinVar's CLNSIG. See ADR-0008.
_CLNSIG_MAGNITUDE: dict[str, float] = {
    "pathogenic": 9.0,
    "pathogenic/likely_pathogenic": 8.5,
    "likely_pathogenic": 7.0,
    "drug_response": 6.5,
    "risk_factor": 6.0,
    "uncertain_significance": 4.0,
    "conflicting_interpretations_of_pathogenicity": 4.0,
    "conflicting_classifications_of_pathogenicity": 4.0,
    # ADR-0045: the per-submission terms ClinVar uses beside the five-tier
    # scale.
    "pathogenic,_low_penetrance": 7.5,
    "likely_pathogenic,_low_penetrance": 6.5,
    "established_risk_allele": 6.0,
    "confers_sensitivity": 6.0,
    "likely_risk_allele": 5.0,
    "affects": 5.0,
    "protective": 5.0,
    "vus-high": 4.5,
    "uncertain_risk_allele": 4.0,
    "vus-mid": 4.0,
    "vus-low": 3.5,
    "not_provided": 2.0,
    "no_classification_for_the_single_variant": 2.0,
    "likely_benign": 2.0,
    "benign/likely_benign": 1.5,
    "benign": 1.0,
}


_BENIGN_CLNSIGS = frozenset({"benign", "likely_benign", "benign/likely_benign"})

# ADR-0045: when ClinVar's own aggregate classification of the variant is
# one of these, a submission is scored no higher than the aggregate.
_AGGREGATE_CAP: dict[str, float] = {
    "conflicting classifications of pathogenicity": 4.0,
    "benign": 1.0,
    "likely benign": 2.0,
    "benign/likely benign": 1.5,
}
# ADR-0045: an outranked P/LP submission whose review tier holds other
# disease classifications is scored as Conflicting.
_OUTSIDE_CONFLICT_CAP = 4.0


def _scored_magnitude(own: float, clnsig: str, aggregate: str, outside_aggregate: int) -> float:
    """A submission's magnitude, capped by ClinVar's view of the variant (ADR-0045).

    The aggregate's germline term judges disease classifications only. A
    drug response, risk factor or similar submission under "Benign; risk
    factor" is a different kind of claim and keeps its own magnitude.
    """
    if not is_disease_classification(clnsig):
        return own
    cap = _AGGREGATE_CAP.get(aggregate_germline_term(aggregate))
    if outside_aggregate == OUTSIDE_AGGREGATE_CONFLICTING:
        cap = _OUTSIDE_CONFLICT_CAP if cap is None else min(cap, _OUTSIDE_CONFLICT_CAP)
    return own if cap is None else min(own, cap)


def _normalize_clnsig(value: str) -> str:
    return value.strip().lower().replace(" ", "_")


def _magnitude(clnsig: str) -> float | None:
    """The magnitude for a classification, or None when it has no score (ADR-0045)."""
    return _CLNSIG_MAGNITUDE.get(_normalize_clnsig(clnsig))


# ClinVar review-status ranking — lower index = stronger
# evidence per ClinVar's star-rating documentation. Used to pick the
# winning row when N SCVs with identical (significance, condition,
# gene) tuples collapse during per-SCV dedup; the surviving Annotation
# carries the highest-rated review_status across the group instead of
# silently dropping provenance.
#
# Real cache rows use the comma-and-underscore form
# (`criteria_provided,_multiple_submitters,_no_conflicts`). The rank
# helper normalizes via _normalize_clnsig so it's robust to spacing
# variation. Unknown values fall through to "weakest known."
_REVIEW_STATUS_RANK: tuple[str, ...] = (
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "criteria_provided,_conflicting_classifications",
    "criteria_provided,_conflicting_interpretations",
    "criteria_provided,_single_submitter",
    "no_assertion_criteria_provided",
    "no_assertion_provided",
    "no_classification_provided",
)


def _review_status_rank(status: str) -> int:
    """Return the rank index for a review_status string.

    Lower index = stronger evidence. Unknown / unrecognized values get
    the max rank (weakest), so a known SCV always beats an unknown one
    when collapsing.
    """
    if not status:
        return len(_REVIEW_STATUS_RANK)
    normalized = _normalize_clnsig(status)
    for i, candidate in enumerate(_REVIEW_STATUS_RANK):
        if candidate == normalized:
            return i
    return len(_REVIEW_STATUS_RANK)


def _collapse_identical_scvs(annotations: list[Annotation]) -> list[Annotation]:
    """Collapse identical per-SCV rows from the TSV loader.

    Two Annotations are identical for dedup purposes if they share
    (rsid, significance, condition, gene, alt) and the row's locus and REF.
    The surviving row carries:

      - the strongest ``review_status`` across the collapsed group
        (per ``_REVIEW_STATUS_RANK``)
      - the union of ``references`` across the collapsed group
        (preserving first-seen order — never silently drops provenance)

    Distinct (significance, condition) pairs on the same variant are
    preserved as separate rows — that's exactly what the per-SCV loader was built to
    surface, and collapsing them would be the Frankenstein-pair
    regression the per-SCV loader retired.

    ``alt`` is in the key for the
    wrong-allele safety invariant. Two SCVs at a multi-allelic site
    can share (significance, condition, gene) but differ in alt — e.g.
    rs1065852 (G→A drug_response/Codeine_response + G→C drug_response/
    Codeine_response would be a hypothetical). Collapsing across alts
    would (a) drop a distinct-allele annotation and (b) mis-attach the
    surviving row's ``alt`` field to the merged result, breaking the
    downstream exact-(rsid, alt) lookup that gnomAD / AlphaMissense /
    CADD use for allele-specific enrichment. The previous v2.2.0 code
    emitted one Annotation per cache row, so this case was always
    handled correctly; the current dedup must preserve that property.

    The row's locus and REF are in the key too (ADR-0047): a call written
    untrimmed can match two ClinVar variants with one ALT, such as
    ``TTCTC>TTC,T`` 2/1 matching ``TTC>T`` (trimmed) and ``TTCTC>T``.
    ClinVar can give both one rsID, classification and condition.

    Input order is preserved for the kept rows: each group is replaced
    by its winner at the position of the group's first occurrence.
    """
    from dataclasses import replace

    groups: dict[tuple, list[Annotation]] = {}
    for a in annotations:
        key = (a.rsid, a.significance, a.condition, a.gene, a.alt, a.locus, a.locus_ref)
        groups.setdefault(key, []).append(a)

    out: list[Annotation] = []
    for group in groups.values():
        if len(group) == 1:
            out.append(group[0])
            continue
        winner = min(group, key=lambda a: _review_status_rank(a.review_status))
        seen_refs: set[str] = set()
        merged_refs: list[str] = []
        for a in group:
            for ref in a.references:
                if ref not in seen_refs:
                    seen_refs.add(ref)
                    merged_refs.append(ref)
        out.append(replace(winner, references=merged_refs))
    return out


_READABLE = frozenset({CacheState.READY, CacheState.OUTDATED})


class ClinVarAnnotator(Annotator):
    """Annotates variants with ClinVar's clinical significance classifications.

    Per-build aware (ADR-0021). At `setup()` time, downloads each
    requested build's VCF (default: both). At `annotate()` time,
    dispatches to the cache matching `variant.build`. If the matching
    cache is missing, the variant is skipped and a warning logged
    (db update needed).
    """

    name: ClassVar[str] = "clinvar"
    display_name: ClassVar[str] = "ClinVar"
    attribution: ClassVar[str] = "ClinVar"
    requires_download: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="custom-clinvar",
        license_url="https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",
        attribution_text="ClinVar variant classifications from NCBI.",
        source_url="https://www.ncbi.nlm.nih.gov/clinvar/",
        commercial_ok=True,
    )

    def __init__(
        self,
        data_dir: Path,
        builds: tuple[str, ...] = CLINVAR_SUPPORTED_BUILDS,
        *,
        include_benign: bool = False,
        other_build_refs: Callable[[dict[str, str], str], dict[str, str]] | None = None,
    ) -> None:
        """Resolve per-build SQLite cache paths within `data_dir`.

        `builds` selects which builds this annotator instance manages.
        Default is both GRCh37 and GRCh38. Passing a single-element
        tuple (e.g. `("GRCh38",)`) restricts setup/refresh to that
        build — used by the CLI's `--build` flag.

        `include_benign` controls whether Benign/Likely_benign annotations
        are emitted. Default False suppresses them (ADR-0008 amendment).

        `other_build_refs` (``RsidIdentityAnnotator.normal_alleles_on_other_build``)
        gives the normal allele at a reference-carried rsID when the other
        build's ClinVar cache cannot (it is not installed, or not current).
        """
        super().__init__(data_dir)
        self._other_build_refs = other_build_refs
        self._builds = tuple(builds)
        self._include_benign = include_benign
        # ADR-0045: rows whose classification has no score here; not shown.
        self._unrecognised_rows = 0
        # VCF records at a ClinVar position whose REF is none of the rows'.
        self._ref_mismatches = 0
        # ADR-0047: VCF calls that matched only after trimming an allele.
        self._trimmed_matches = 0
        # ADR-0051: array calls never matched against an aggregate-P/LP
        # row: i-probes at its locus, and D/I calls at its rsID.
        self._i_probes_at_plp_loci = 0
        self._di_calls_at_plp_rsids = 0
        for build in self._builds:
            if build not in CLINVAR_SUPPORTED_BUILDS:
                raise ValueError(
                    f"Unsupported ClinVar build {build!r}; expected one of "
                    f"{CLINVAR_SUPPORTED_BUILDS}"
                )
        self._db_paths: dict[str, Path] = {
            build: data_dir / clinvar_db_filename(build) for build in self._builds
        }
        self._conns: dict[str, sqlite3.Connection] = {}
        # ADR-0023: per-build (rsid -> single-base REF) cache. ClinPGx
        # consults this as its primary non-finding filter. Built lazily
        # on first lookup per build.
        self._ref_lookups: dict[str, dict[str, str]] = {}
        # Per build, the rsIDs where ``reference_for`` gives the normal
        # allele from the other build (the build's reference carries a
        # ClinVar variant allele), with that allele.
        self._normal_alleles: dict[str, dict[str, str]] = {}
        # Per build, whether its cache is ready for this version: rsID
        # resolution outside annotation uses only a ready cache.
        self._ready_builds: dict[str, bool] = {}

    def serves_build(self, build: str) -> bool:
        """Only the builds this run manages (the ready caches, when any are)."""
        return build in self._db_paths

    def run_counts(self) -> dict[str, int]:
        """ADR-0045, ADR-0047 and ADR-0051 counts for the run so far."""
        return {
            "unrecognised_rows": self._unrecognised_rows,
            "trimmed_matches": self._trimmed_matches,
            "i_probes_at_plp_loci": self._i_probes_at_plp_loci,
            "di_calls_at_plp_rsids": self._di_calls_at_plp_rsids,
            "ref_mismatches": self._ref_mismatches,
        }

    def _connection(self, build: str) -> sqlite3.Connection | None:
        """Return a lazy read-only connection to the per-build cache, or None if missing.

        ``idx_clinvar_position`` comes from the loader schema (every cache
        since v2.0.0), so nothing is built here (ADR-0048).
        """
        if build not in self._db_paths:
            return None
        if build not in self._conns:
            db_path = self._db_paths[build]
            # An outdated cache is still read (ClinPGx's REF filter uses it,
            # as before); one that cannot be read is not opened.
            if self._build_status(build).state not in _READABLE:
                return None
            self._conns[build] = connect_readonly(db_path)
        return self._conns[build]

    def setup(self) -> None:
        """Download the per-SCV TSV pair once and ingest both builds from it.

        Since v2.2 the VCF-per-build path is replaced by the
        variant_summary + submission_summary TSV pair. variant_summary
        carries per-(VariationID, Assembly) rows so the same downloaded
        file populates both GRCh37 and GRCh38 caches with no second
        download. submission_summary is build-agnostic; one fetch
        applies across builds.

        Freshness signal is derived from submission_summary's md5
        endpoint — that file drives the meaningful change cadence
        (per-SCV submissions land daily) and variant_summary turns over
        weekly on a slower clock.

        Both compressed sources are checked before ingest. All selected
        builds meet fixed/relative row floors and have final metadata
        before publication. Replacements are individually atomic; a later
        replacement failure reports any build already updated.

        The freshness signal is NCBI's submission_summary checksum, which
        also verifies that download, so a refresh is never done without it.
        A failed fetch is tried once more (a transient failure may clear).
        """
        signal = self._fetch_remote_signal_for_tsv() or self._fetch_remote_signal_for_tsv()
        if signal is None:
            raise RuntimeError(
                "clinvar: ClinVar not refreshed: NCBI's checksum for submission_summary "
                "is unavailable, so the download cannot be verified. Existing caches "
                "kept; try again later."
            )

        vs_url = _manager_module.CLINVAR_VARIANT_SUMMARY_URL
        ss_url = _manager_module.CLINVAR_SUBMISSION_SUMMARY_URL
        previous_counts = {build: self._previous_record_count(build) for build in self._builds}
        variant_md5 = _fetch_checksum(vs_url)
        if variant_md5 is None:
            raise RuntimeError(
                "clinvar: cannot verify variant_summary checksum; refresh aborted, "
                "existing caches kept"
            )
        self.data_dir.mkdir(parents=True, exist_ok=True)
        # All owned downloads and candidates live on the cache filesystem.
        # Validation completes for every selected build before publication.
        with tempfile.TemporaryDirectory(prefix=".clinvar-", dir=self.data_dir) as temporary:
            staging = Path(temporary)
            vs_path = staging / "variant_summary.txt.gz"
            ss_path = staging / "submission_summary.txt.gz"
            download(vs_url, vs_path)
            download(ss_url, ss_path)
            verify_file_hash(vs_path, "md5", variant_md5)
            verify_file_hash(ss_path, "md5", signal.removeprefix("md5:"))
            candidates = {}
            for build in self._builds:
                candidate = staging / clinvar_db_filename(build)
                self._setup_one_from_tsv(
                    vs_path, ss_path, candidate, build, signal, vs_url, previous_counts[build]
                )
                candidates[build] = candidate
            self.close()
            published: list[str] = []
            for build, candidate in candidates.items():
                try:
                    os.replace(candidate, self._db_paths[build])
                except OSError as exc:
                    kept = (
                        "old cache kept"
                        if self._db_paths[build].exists()
                        else "no cache installed"
                    )
                    raise RuntimeError(
                        f"ClinVar {build}: publication failed; {kept}; "
                        f"already updated: {', '.join(published) or 'none'}. {exc}"
                    ) from exc
                published.append(build)

    def _previous_record_count(self, build: str) -> int | None:
        """Use existing metadata, refusing an unknown baseline without a full scan."""
        path = self._db_paths[build]
        if not path.exists():
            return None
        info = get_database_info(path, clinvar_record_name(build))
        if info is None:
            info = get_database_info(path, "clinvar")
        count = info["record_count"] if info is not None else None
        if type(count) is not int or count <= 0:
            raise RuntimeError(
                f"ClinVar {build}: existing cache record count is missing or invalid; "
                "cannot establish the relative floor, existing caches kept"
            )
        return count

    def _setup_one_from_tsv(
        self,
        vs_path: Path,
        ss_path: Path,
        candidate: Path,
        build: str,
        signal: str,
        source_url: str,
        previous_count: int | None,
    ) -> None:
        """Finish and validate a candidate against the downloaded source pair."""
        load_clinvar_tsv(
            vs_path,
            ss_path,
            candidate,
            build,
            source_url=source_url,
            remote_signal=signal,
            record_name=clinvar_record_name(build),
            min_rows=_manager_module.CLINVAR_MIN_ROWS[build],
            previous_count=previous_count,
        )

    def is_ready(self) -> bool:
        """True iff EVERY managed build has a populated cache with the current tag."""
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        """The first managed build's cache that is not ready, or ready (ADR-0048)."""
        for build in self._builds:
            status = self._build_status(build)
            if not status.ready:
                return status
        return READY

    def _build_status(self, build: str) -> CacheStatus:
        return cache_state(
            self._db_paths[build],
            clinvar_record_name(build),
            tag=f"iv:{CLINVAR_INTERPRETER_VERSION}",
            tables=("clinvar_variants",),
        )

    def prepare_cache(self) -> bool:
        """Self-heal a pre-mechanism cache whose stamp says the current version.

        Caches with the tag baked into ``remote_signal`` are stamped once
        via ``stamp_existing_clinvar_cache``; any other outdated cache is
        left for ``db update`` to rebuild.
        """
        changed = False
        for build in self._builds:
            if self._build_status(build).state is not CacheState.OUTDATED:
                continue
            changed = stamp_existing_clinvar_cache(self._db_paths[build]) or changed
        return changed

    def _cache_is_current(self, build: str) -> bool:
        """True when ``build``'s cache carries the current interpreter-version tag."""
        if build not in self._db_paths:
            return False
        info = get_database_info(self._db_paths[build], clinvar_record_name(build))
        return info is not None and info.get("local_version_tag") == (
            f"iv:{CLINVAR_INTERPRETER_VERSION}"
        )

    def version(self) -> str | None:
        """Composite version string across managed builds.

        Format: `"GRCh37:<v>; GRCh38:<v>"` when both present, or a
        single `<build>:<v>` when only one is managed. None if none.
        """
        parts: list[str] = []
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is not None:
                parts.append(f"{build}:{info['version']}")
        return "; ".join(parts) if parts else None

    def record_count(self) -> int | None:
        """Total record count across managed build caches, or None if none cached."""
        total = 0
        any_present = False
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is not None:
                any_present = True
                total += info["record_count"]
        return total if any_present else None

    def close(self) -> None:
        """Close all open per-build connections. Safe to call repeatedly."""
        for conn in self._conns.values():
            conn.close()
        self._conns.clear()
        self._ref_lookups.clear()
        self._normal_alleles.clear()
        self._ready_builds.clear()

    def reference_for(self, rsid: str, build: str) -> str | None:
        """Return ClinVar's single-base REF allele for `rsid` in `build`, or None.

        ADR-0023: ClinPGx's primary non-finding filter calls this. If the
        return value matches both of the user's alleles, the user is
        homozygous reference and the ClinPGx annotation is a non-finding.

        Lazily builds an in-memory `(rsid -> REF)` map per build on first
        call so subsequent lookups are O(1). Multi-base REFs (indels) are
        skipped — array-based parsers can't call indels, so a multi-base
        REF can't validly suppress a single-base genotype.

        Where this build's reference carries a ClinVar variant allele (a
        REF == ALT row, F5 c.1601G>A on GRCh37), a homozygous-reference
        call is the variant genotype. There the normal allele is returned
        instead: the REF of the same ClinVar variant on the other build,
        when that row is a single-base REF != ALT (C for F5), so the same
        call is a non-finding on both builds. Without such a row (only one
        build's cache installed, or REF == ALT on both builds) it returns
        None and nothing is suppressed.

        Returns None when ClinVar has no data for the rsid in this build
        (or has only indel REFs). Callers fall through to secondary tiers.
        """
        if build not in self._db_paths:
            return None
        if build not in self._ref_lookups:
            self._ref_lookups[build] = self._load_ref_lookup(build)
        return self._ref_lookups[build].get(rsid)

    def records_of_rsids(
        self, rsids: set[str], build: str
    ) -> dict[str, set[tuple[str, int, str, str]]] | None:
        """Each rsID's ``(chrom, pos, ref, alt)`` records in ``build``'s cache (ADR-0052).

        None when this run has no readable cache for ``build``; rsIDs the
        cache does not hold are absent.
        """
        conn = self._connection(build)
        if conn is None:
            return None
        wanted = sorted(rsids)
        out: dict[str, set[tuple[str, int, str, str]]] = {}
        for start in range(0, len(wanted), _BATCH_CHUNK):
            chunk = wanted[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            for rsid, chrom, pos, ref, alt in conn.execute(
                f"SELECT DISTINCT rsid, chromosome, position, ref, alt FROM clinvar_variants "
                f"WHERE rsid IN ({placeholders})",
                chunk,
            ):
                out.setdefault(rsid, set()).add((chrom, pos, ref, alt))
        return out

    def reference_at(self, rsid: str, build: str, locus: tuple[str, int]) -> str | None:
        """``reference_for``, only where ClinVar places ``rsid`` at ``locus`` on ``build``.

        The ADR-0023 hom-ref filter of ClinPGx and SNPedia uses this, so a
        call is suppressed only by ClinVar's REF at the call's own position.
        An rsID ClinVar places elsewhere, or a call with no position, gives
        None and nothing is suppressed.
        """
        ref = self.reference_for(rsid, build)
        if ref is None:
            return None
        # One indexed query per call: a map of every rsID's loci would
        # hold millions of entries for the few thousand calls made.
        conn = self._connection(build)
        if conn is None:
            return None
        chrom, pos = locus
        row = conn.execute(
            "SELECT 1 FROM clinvar_variants WHERE rsid = ? AND chromosome = ? "
            "AND position = ? AND length(ref) = 1 LIMIT 1",
            (rsid, chrom, pos),
        ).fetchone()
        return ref if row is not None else None

    def _load_ref_lookup(self, build: str) -> dict[str, str]:
        """Read the per-build cache once and build the `(rsid -> REF)` map."""
        conn = self._connection(build)
        if conn is None:
            return {}
        # Single-base REFs only: indel anchor-base encoding (REF=CTT, etc.)
        # can't suppress a single-base array readout. The per-build cache
        # may have BOTH SNV and indel rows for the same rsid; the WHERE
        # filters those out so we keep only the SNV REF.
        #
        # Iterate the cursor directly instead of fetchall()
        # so peak memory is one row + the dict, not the entire result set
        # twice.
        # Rows without an rsID (ADR-0046) have nothing to look up by.
        # An outdated cache has no ``normal_allele`` column (before iv:9).
        normal_col = "normal_allele" if self._cache_is_current(build) else "0"
        cursor = conn.execute(
            f"SELECT DISTINCT rsid, ref, alt, allele_id, {normal_col} "
            "FROM clinvar_variants WHERE length(ref) = 1 AND rsid <> ''"
        )
        out: dict[str, str] = {}
        # rsID -> VariationIDs whose row is REF == ALT: this build's
        # reference carries their variant allele. A normal-allele record
        # written REF == ALT names the reference's own (normal) allele and
        # gives no REF (below).
        reference_carried: dict[str, set[int]] = {}
        # rsID -> the variant allele this build's reference carries.
        carried_ref: dict[str, str] = {}
        for rsid, ref, alt, variation_id, normal_allele in cursor:
            if ref == alt and not normal_allele:
                reference_carried.setdefault(rsid, set()).add(variation_id)
                carried_ref.setdefault(rsid, ref)
            if ref == alt and normal_allele:
                # A normal-allele record written REF == ALT gives no REF: its
                # normal allele can be the classified one (CASP12), so it must
                # not make a homozygous-reference call a ClinPGx or SNPedia
                # non-finding. Caches before iv:9 did not hold these records.
                continue
            # If a rsid has multiple single-base REFs (shouldn't happen at
            # one position but defending against future data shapes), keep
            # the first.
            if rsid not in out:
                out[rsid] = ref
        # Only a current cache gives REF == ALT its one meaning; an older
        # one still holds normal-allele records written REF == ALT, and
        # keeps the locus REF (the 2.4.1 behaviour) until `db update`.
        normal: dict[str, str] = {}
        if reference_carried and self._cache_is_current(build):
            normal = self._normal_alleles_on_other_build(build, reference_carried, carried_ref)
            for rsid in reference_carried:
                if rsid in normal:
                    out[rsid] = normal[rsid]
                else:
                    out.pop(rsid, None)
        self._normal_alleles[build] = normal
        return out

    def reference_carried_sites(self, build: str) -> dict[tuple[str, int], tuple[str, str, bool]]:
        """``{(chrom, pos): (allele, rsid)}`` where ``build``'s reference carries a row's allele.

        The one-base REF == ALT rows of a current cache: a variant allele
        the reference carries (T at F5 rs6025 on GRCh37), or the normal
        allele of a normal-allele record (G at CASP12 rs497116 on GRCh37).
        A VCF says nothing at such a site for a person homozygous for that
        allele: no record, a ``0/0`` call or a reference block.
        The pipeline watches these positions and assesses a hom-ref call
        there as the homozygous genotype. The flag is True for a
        normal-allele record. Rows without an rsID, and positions with more
        than one such allele or rsID, are left out.
        """
        conn = self._connection(build)
        if conn is None or not self._cache_is_current(build):
            return {}
        found: dict[tuple[str, int], set[tuple[str, str, bool]]] = {}
        for chrom, pos, ref, rsid, normal_allele in conn.execute(
            "SELECT DISTINCT chromosome, position, ref, rsid, normal_allele "
            "FROM clinvar_variants WHERE ref = alt AND length(ref) = 1 AND rsid <> ''"
        ):
            found.setdefault((chrom, pos), set()).add((ref, rsid, bool(normal_allele)))
        return {locus: next(iter(v)) for locus, v in found.items() if len(v) == 1}

    def normal_allele_for(self, rsid: str, build: str) -> str | None:
        """The normal allele at `rsid` where `build`'s reference carries a variant, or None.

        Set only where ``reference_for`` returns the other build's REF of
        the same ClinVar variant (C at F5 rs6025 on GRCh37). ClinPGx and
        SNPedia derive their ``alt`` against it there; everywhere else
        they keep ``Variant.ref``, even where ClinVar's REF differs from it.
        """
        self.reference_for(rsid, build)
        return self._normal_alleles.get(build, {}).get(rsid)

    def _normal_alleles_on_other_build(
        self,
        build: str,
        reference_carried: dict[str, set[int]],
        carried_ref: dict[str, str],
    ) -> dict[str, str]:
        """The normal allele at each reference-carried rsID, from the other build.

        For each rsID, the single-base REF of the same VariationIDs on the
        other build where that row's REF differs from its ALT. An rsID with
        no such row, or more than one candidate REF, is left out.

        Without a current ClinVar cache for the other build (one build
        installed), the rsID identity table's REF on the other build is used
        instead (``other_build_refs``), when it is one base that differs
        from the allele this build's reference carries.
        """
        other = BUILD_GRCH38 if build == BUILD_GRCH37 else BUILD_GRCH37
        conn = self._connection(other)
        if conn is None or not self._cache_is_current(other):
            if self._other_build_refs is None:
                return {}
            return self._other_build_refs(carried_ref, build)
        rsids = sorted(reference_carried)
        candidates: dict[str, set[str]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            for rsid, ref, alt, variation_id in conn.execute(
                f"SELECT DISTINCT rsid, ref, alt, allele_id FROM clinvar_variants "
                f"WHERE rsid IN ({placeholders})",
                chunk,
            ):
                if variation_id in reference_carried[rsid] and len(ref) == 1 and ref != alt:
                    candidates.setdefault(rsid, set()).add(ref)
        return {rsid: refs.pop() for rsid, refs in candidates.items() if len(refs) == 1}

    def fetch_remote_signal(self) -> str | None:
        r"""Freshness signal derived from the per-SCV submission_summary md5.

        Before v2.2 this returned a composite VCF-md5 per build. Since
        v2.2, the cache is built from a build-agnostic
        TSV pair, and submission_summary's md5 endpoint is the change
        cadence we want to track (per-SCV submissions land daily; the
        slower-moving variant_summary only adds new VariationIDs).
        Returns the signal as ``"md5:<hex>"`` to keep ``setup()``'s
        ``signal.removeprefix("md5:")`` contract unchanged.
        """
        return self._fetch_remote_signal_for_tsv()

    @staticmethod
    def _fetch_remote_signal_for_tsv() -> str | None:
        digest = _fetch_checksum(_manager_module.CLINVAR_SUBMISSION_SUMMARY_URL)
        return f"md5:{digest}" if digest is not None else None

    def cached_remote_signal(self) -> str | None:
        """Return the signal stamped at the last successful setup, or None.

        None if any managed build's cache is absent / unstamped or if the
        builds' stamps diverge (an interrupted partial reingest).

        Since v2.2 both builds share the same source TSVs and therefore
        the same remote_signal. We still cross-check that EVERY managed
        build's cache carries the stamp; a missing build means the
        last-setup didn't complete, treat as "no signal" and let
        ``db update`` decide.
        """
        signals: set[str] = set()
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is None or info["remote_signal"] is None:
                return None
            sig = info["remote_signal"]
            if not sig:
                return None
            signals.add(sig)
        if not signals:
            return None
        # All builds should carry the same signal (single download, single
        # ingest run). If they diverge — e.g. someone reingested one build
        # against a newer TSV — treat as "no signal" so the next
        # `db update` re-syncs both.
        if len(signals) > 1:
            return None
        return signals.pop()

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return ClinVar annotations for the alleles the user carries.

        Rows are selected by locus, never by rsID alone:
        the cache is queried at ``(chromosome, position)`` and a row
        annotates only when its REF agrees with ``variant.ref`` (when
        known) and the user carries its ALT. See ``_annotations_for``.

        ADR-0007 carrier rule, ADR-0011 indel-anchor protection and
        ADR-0035 strand-aware matching apply per row. ADR-0021: dispatch
        by ``variant.build``; a missing cache yields no annotations (the
        analyze-time build warning already told the user).
        """
        if variant.is_no_call and not variant.is_half_call:
            return []
        conn = self._connection(variant.build)
        if conn is None:
            return []
        self._count_di_calls_at_plp(conn, [variant])
        rows = _candidate_rows(variant, self._rows_at_loci(conn, [variant]))
        known = self._loci_of_rsids(conn, [variant])
        return self._annotations_for(variant, rows, known.get(variant.rsid))

    def bulk_resolve_rsids(self, variants: list[Variant]) -> dict[tuple[str, int, str, str], str]:
        """Resolve rsIDs by locus for variants with no usable ID.

        Real-world VCFs from variant callers (GATK HaplotypeCaller,
        DeepVariant) emit ``ID=.``. This queries ClinVar at
        ``(chromosome, position)`` and stamps ``variant.rsid`` in place
        when ``resolve_rsid_at_locus`` finds exactly one rsID whose rows
        explain every allele the user carries (duplicate
        per-SCV rows and ``1/2`` calls carrying two alts both resolve).

        Groups by ``variant.build`` and dispatches per-build (ADR-0021).
        A build whose cache is missing or not ready for this version (an
        older cache lacks columns this version reads) is skipped, as
        annotation skips it; the cache is never changed here.

        The returned ``{(chrom, pos, ref, alt): rsid}`` mapping carries
        the resolved coordinates so the pipeline's enrichment phase can
        do position-keyed gnomAD / AlphaMissense lookups. Only variants
        carrying a single alt contribute an entry: a ``1/2`` call has no
        single allele to stand for the rsID.
        """
        by_build: dict[str, list[Variant]] = {}
        for v in variants:
            if (v.is_no_call and not v.is_half_call) or not v.chromosome or v.position <= 0:
                continue
            if v.ref is None:
                # A call with neither an rsID nor a reference allele (array
                # vendor probes such as 23andMe ``i5000123``) has nothing to
                # confirm that the ClinVar row at its position is the same
                # variant; stay silent rather than guess.
                continue
            by_build.setdefault(v.build, []).append(v)

        resolved: dict[tuple[str, int, str, str], str] = {}
        for build, build_variants in by_build.items():
            if build not in self._db_paths:
                continue
            if build not in self._ready_builds:
                self._ready_builds[build] = self._build_status(build).ready
            conn = self._connection(build) if self._ready_builds[build] else None
            if conn is None:
                continue
            rows_at = self._rows_at_loci(conn, build_variants)
            for v in build_variants:
                rows = rows_at.get(variant_locus(v))
                if not rows:
                    continue
                hit = resolve_rsid_at_locus(
                    v.ref,
                    *v.carried_pair,
                    [(row[_R_REF], row[_R_ALT], row[_R_RSID]) for row in rows],
                )
                if hit is None:
                    continue
                rsid, pair = hit
                # stash the pre-resolution pseudo-rsID for audit /
                # debug before stamping. Debug-only — see
                # Variant.original_rsid docstring; the serialization-guard
                # test asserts this never reaches public output.
                v.original_rsid = v.rsid
                v.rsid = rsid
                if pair is not None:
                    resolved[(v.chromosome, v.position, *pair)] = rsid
        return resolved

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Bulk-annotate via chunked locus queries per build.

        ClinVar holds one SQLite cache per build (ADR-0021). Variants are
        grouped by ``variant.build`` and each group's loci are fetched
        with chunked ``WHERE chromosome = ? AND position IN (...)``
        queries. Per-row filtering is shared with ``annotate`` through
        ``_annotations_for``, so the two paths cannot diverge.

        Variants whose build has no cache loaded are skipped silently —
        the build warning was already emitted at analyze time.
        """
        by_build: dict[str, list[Variant]] = {}
        for v in variants:
            if not v.is_no_call or v.is_half_call:
                by_build.setdefault(v.build, []).append(v)
        for build, build_variants in by_build.items():
            conn = self._connection(build)
            if conn is None:
                continue
            self._count_di_calls_at_plp(conn, build_variants)
            rows_at = self._rows_at_loci(conn, build_variants)
            known = self._loci_of_rsids(conn, build_variants)
            for variant in build_variants:
                rows = _candidate_rows(variant, rows_at)
                if rows:
                    yield from self._annotations_for(variant, rows, known.get(variant.rsid))

    def _rows_at_loci(
        self, conn: sqlite3.Connection, variants: Iterable[Variant]
    ) -> dict[tuple[str, int], list[tuple]]:
        """Fetch every cache row at the variants' loci, keyed by ``(chrom, pos)``.

        Rows come back in insertion order (``ORDER BY id``) so the per-SCV
        collapse and the resolver see a deterministic sequence regardless
        of SQLite's physical row order.
        """
        by_chrom: dict[str, set[int]] = {}
        for v in variants:
            if v.chromosome and v.position > 0:
                chrom, pos = variant_locus(v)
                by_chrom.setdefault(chrom, set()).add(pos)
                by_chrom[chrom].update(t[0] for t in _trimmed_alleles(v))
        out: dict[tuple[str, int], list[tuple]] = {}
        for chrom, positions in by_chrom.items():
            position_list = sorted(positions)
            for start in range(0, len(position_list), _BATCH_CHUNK):
                chunk = position_list[start : start + _BATCH_CHUNK]
                placeholders = ",".join("?" * len(chunk))
                cursor = conn.execute(
                    f"SELECT chromosome, position, ref, alt, clinical_significance, "
                    f"condition, gene, review_status, allele_id, rsid, aggregate, "
                    f"submitters_same, submitters_total, outside_aggregate "
                    f"FROM clinvar_variants "
                    f"WHERE chromosome = ? AND position IN ({placeholders}) ORDER BY id",
                    (chrom, *chunk),
                )
                for row in cursor:
                    out.setdefault((row[0], row[1]), []).append(row)
        return out

    def _count_di_calls_at_plp(self, conn: sqlite3.Connection, variants: list[Variant]) -> None:
        """ADR-0051: count D/I array calls whose rsID has an aggregate-P/LP row.

        ``D`` and ``I`` are an array's codes for a deletion and an insertion
        allele, so they never equal a ClinVar ALT and are never matched.
        """
        di_calls = [
            v
            for v in variants
            if v.rsid.startswith("rs") and not {v.allele1, v.allele2}.isdisjoint(DI_ALLELES)
        ]
        rsids = sorted({v.rsid for v in di_calls})
        plp: set[str] = set()
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            for rsid, aggregate in conn.execute(
                f"SELECT DISTINCT rsid, aggregate FROM clinvar_variants "
                f"WHERE rsid IN ({placeholders})",
                chunk,
            ):
                if is_pathogenic_aggregate(aggregate or ""):
                    plp.add(rsid)
        self._di_calls_at_plp_rsids += sum(1 for v in di_calls if v.rsid in plp)

    def _loci_of_rsids(
        self, conn: sqlite3.Connection, variants: Iterable[Variant]
    ) -> dict[str, set[tuple[str, int]]]:
        """Where ClinVar places each input rsID on this build, keyed by rsID."""
        rsids = sorted({v.rsid for v in variants if v.rsid.startswith("rs")})
        out: dict[str, set[tuple[str, int]]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            for rsid, chrom, pos in conn.execute(
                f"SELECT DISTINCT rsid, chromosome, position FROM clinvar_variants "
                f"WHERE rsid IN ({placeholders})",
                chunk,
            ):
                out.setdefault(rsid, set()).add((chrom, pos))
        return out

    def _annotations_for(
        self,
        variant: Variant,
        rows: list[tuple],
        rsid_loci: set[tuple[str, int]] | None = None,
    ) -> list[Annotation]:
        """Turn the cache rows at a variant's locus into Annotations.

        A row annotates only when all of these hold:

        - REF == ALT means this build's reference carries the row's allele:
          a variant allele (F5 c.1601G>A is T>T on GRCh37, whose reference
          base is the Leiden allele) or, for a normal-allele record (HGVS
          ``c.N=``, ``normal_allele``), the normal allele (CASP12 rs497116
          is G>G on GRCh37 and A>G on GRCh38; cache ``iv:9``). The loader
          drops REF == ALT rows that are not SNVs. Such a row annotates a
          user who carries its ALT, counting copies of it, like any other
          row; P/LP classifications of normal-allele records are dropped
          at load.
        - The row's REF equals ``variant.ref`` when the input supplies one.
          A VCF is forward strand by specification and an array call's REF
          comes from the identity table at the same forward locus, so a
          complemented REF is not the same variant: a one-base VCF REF that
          is none of the rows' one-base REFs is counted (``ref_mismatches``).
          Together with the locus this identifies the variant by
          ``(chrom, pos, ref, alt)`` rather than by rsID.
        - ADR-0011: ClinVar indels never match single-base array calls.
          A call whose REF came from the input (a VCF record) is exempt:
          REF equality already pins the variant, and a homozygous deletion
          to one base carries no multi-base allele (ADR-0047).
        - ADR-0035 / ADR-0007: the user carries the ALT.
        - ADR-0047: a VCF allele written untrimmed (``GTT>GT,G``,
          ``CA>TA``) also matches the row equal to its trimmed form, at
          the trimmed position. The annotation then carries the row's
          position and alleles; the genotype is shown as the file wrote it.
        - ADR-0047: a VCF half call (``./1``) is assessed on its called
          allele and shown as ``T/.``.
        - Benign rows are suppressed unless ``include_benign``.

        Two guards decide whether the locus may be matched at all:

        - A call with no reference allele is locus-matched only when
          ClinVar holds its rsID on this build (``rsid_loci``, every
          position ClinVar records for that rsID). That excludes array
          vendor probes without an rsID (23andMe ``i5000123``) and array
          calls whose rsID neither ClinVar nor the identity table holds (a
          withdrawn or merged ID): nothing confirms that the row at the
          position describes the same variant. The pipeline fills REF only
          from the identity rows of the call's own rsID, so a call with a
          REF is either one of those or a VCF record whose REF came from
          the file, and the REF check applies to it.
        - When the input carries an rsID that ClinVar places elsewhere
          and not at the input's locus, the rows at the locus belong to a
          different variant and are not reported under the input's rsID.

        - ADR-0046: a row without an rsID matches only a call whose REF
          came from the input file (``Variant.ref_from_input``, VCF).

        The Annotation's rsID is the matching row's when the input names
        another rsID ClinVar holds among these rows (another ALT of a
        multi-allelic record, an untrimmed allele's other form). Otherwise
        it is the input's when it has one (a merged or alias rsID is kept),
        then the row's, then the input's own ID (``""`` for a VCF ID of
        "."). Identical per-SCV rows collapse.
        """
        if variant.is_no_call and not variant.is_half_call:
            return []
        if variant.ref is None and not rsid_loci:
            if variant.rsid.startswith("i") and any(
                row[_R_POS] == variant_locus(variant)[1]
                and is_pathogenic_aggregate(row[_R_AGGREGATE] or "")
                for row in rows
            ):
                # ADR-0051: a vendor probe at a locus with an aggregate-P/LP row.
                self._i_probes_at_plp_loci += 1
            return []
        record_locus = variant_locus(variant)
        trimmed = _trimmed_alleles(variant)
        if rsid_loci and not rsid_loci & (
            {record_locus} | {(record_locus[0], t[0]) for t in trimmed}
        ):
            return []
        allele1, allele2 = variant.carried_pair
        user_is_multibase = len(allele1) > 1 or len(allele2) > 1
        # ADR-0023: report the user's actual diploid call consistently
        # across annotators, not the matched ALT base alone.
        user_diploid = diploid_label(variant)
        matched_trimmed = False
        ref_agrees = False
        # The rsIDs ClinVar holds among these rows. A record naming one of
        # them can match the rows of the others (another ALT of a
        # multi-allelic record, an untrimmed allele's other form): each
        # such row carries its own rsID.
        locus_rsids = {row[_R_RSID] for row in rows if row[_R_RSID]}
        annotations: list[Annotation] = []
        for row in rows:
            (
                _chrom,
                row_pos,
                ref,
                alt,
                clnsig,
                condition,
                gene,
                review_status,
                allele_id,
                row_rsid,
                aggregate,
                submitters_same,
                submitters_total,
                outside_aggregate,
            ) = row
            if not row_rsid and not variant.ref_from_input:
                # ADR-0046: a variant without an rsID is identified by its
                # locus and alleles alone, so only a REF read from the
                # input (VCF) can confirm it. An array call's REF comes
                # from the rsID identity table of its own rsID.
                continue
            written = trimmed.get((row_pos, ref, alt))
            if written is not None:
                ref_agrees = True
                # ADR-0047: the row is a carried allele's trimmed form.
                carried = tuple(
                    alt if a == written else ref if a == variant.ref else a
                    for a in (allele1, allele2)
                )
                locus = (record_locus[0], row_pos)
            else:
                if row_pos != record_locus[1]:
                    continue
                # A REF from a VCF is forward strand by specification, and
                # an array call's REF comes from the identity table at the
                # same forward locus as ClinVar's: only an equal REF is the
                # same variant. A complemented REF is a mismatch.
                if variant.ref is not None and variant.ref != ref:
                    continue
                ref_agrees = True
                clinvar_is_indel = len(ref) > 1 or len(alt) > 1
                if clinvar_is_indel and not user_is_multibase and not variant.ref_from_input:
                    continue
                if alt not in (allele1, allele2):
                    continue
                carried = (allele1, allele2)
                locus = record_locus
            sig_label = _normalize_clnsig(clnsig) if clnsig else "unknown"
            if not self._include_benign and sig_label in _BENIGN_CLNSIGS:
                continue
            own_magnitude = _CLNSIG_MAGNITUDE.get(sig_label)
            if own_magnitude is None:
                # ADR-0045: a classification Allelix has no score for is not
                # shown, and the run counts it.
                self._unrecognised_rows += 1
                continue
            magnitude = _scored_magnitude(
                own_magnitude, clnsig or "", aggregate or "", outside_aggregate
            )
            description = (
                f"ClinVar classifies this allele as "
                f"{clnsig.replace('_', ' ') if clnsig else 'unknown significance'}"
            )
            references = [f"clinvar:allele/{allele_id}"] if allele_id else []
            matched_trimmed = matched_trimmed or written is not None
            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=(
                        row_rsid
                        if row_rsid and variant.rsid in locus_rsids
                        else variant.rsid
                        if variant.rsid.startswith("rs")
                        else row_rsid or variant.rsid
                    ),
                    significance=f"clinvar_{sig_label}",
                    category="clinical",
                    magnitude=magnitude,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition="" if not condition or condition == "." else condition,
                    gene=gene or "",
                    review_status=review_status or "",
                    alt=alt,
                    locus_ref=ref,
                    carried_alleles=carried,
                    locus=locus,
                    record_key=None
                    if written is None or variant.ref is None
                    else (variant.chromosome, variant.position, variant.ref, written),
                    reference_carries_alt=ref == alt,
                    clinvar_aggregate=aggregate or "",
                    submitters_same=submitters_same,
                    submitters_total=submitters_total,
                    aggregate_capped=magnitude < own_magnitude,
                    outside_aggregate=outside_aggregate != CONTRIBUTES_TO_AGGREGATE,
                )
            )
        if matched_trimmed:
            self._trimmed_matches += 1
        if (
            variant.ref_from_input
            and not ref_agrees
            and len(variant.ref or "") == 1
            and any(row[_R_POS] == record_locus[1] and len(row[_R_REF]) == 1 for row in rows)
        ):
            # A one-base VCF REF that is none of ClinVar's one-base REFs at
            # the position: mis-coordinated input, the complement included.
            self._ref_mismatches += 1
        return _collapse_identical_scvs(annotations)


# Column positions in ``_rows_at_loci`` rows used outside the unpacking.
_R_POS = 1
_R_REF = 2
_R_ALT = 3
_R_RSID = 9
_R_AGGREGATE = 10


def _trimmed_alleles(variant: Variant) -> dict[tuple[int, str, str], str]:
    """A VCF call's carried alleles that its record wrote untrimmed (ADR-0047).

    Maps each one's trimmed ``(position, ref, alt)`` to the allele as
    written. Empty for a call whose REF did not come from the input:
    trimming against a REF from a reference table would be a guess.
    """
    if not variant.ref_from_input or variant.ref is None:
        return {}
    out: dict[tuple[int, str, str], str] = {}
    for allele in set(variant.carried_pair):
        if allele in (variant.ref, NO_CALL_MARKER):
            continue
        key = trim_allele_pair(variant.position, variant.ref, allele)
        if key != (variant.position, variant.ref, allele):
            out[key] = allele
    return out


def _candidate_rows(variant: Variant, rows_at: dict[tuple[str, int], list[tuple]]) -> list[tuple]:
    """The cache rows at a call's locus, plus those at its trimmed alleles' positions."""
    chrom, pos = variant_locus(variant)
    rows = list(rows_at.get((chrom, pos), []))
    for trimmed_pos in sorted({t[0] for t in _trimmed_alleles(variant)} - {pos}):
        rows.extend(rows_at.get((chrom, trimmed_pos), []))
    return rows

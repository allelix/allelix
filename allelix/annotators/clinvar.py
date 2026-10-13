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
import re
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor, variant_locus
from allelix.databases import manager as _manager_module
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import (
    download,
    fetch_remote_text,
    get_database_info,
    load_clinvar_tsv,
    stamp_existing_clinvar_cache,
    verify_file_hash,
)
from allelix.models import Annotation
from allelix.utils.allele import (
    carried_in_source_orientation,
    complement,
    resolve_rsid_at_locus,
    strand_aware_carrier_match,
)
from allelix.utils.build_detect import BUILD_GRCH37, BUILD_GRCH38

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

CLINVAR_SUPPORTED_BUILDS: tuple[str, ...] = ("GRCh37", "GRCh38")

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999

# a remote .md5 endpoint can return an HTML error page on a
# transient blip. The first whitespace-separated token of the body is
# what we treat as the hash, so without this gate `<!DOCTYPE` would be
# accepted as the "signal" and later passed to `verify_file_hash`, which
# would then delete the freshly downloaded VCF. MD5 is exactly 32 hex
# digits; reject anything else.
_MD5_HEX_RE = re.compile(r"^[0-9a-fA-F]{32}$")


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
    "not_provided": 2.0,
    "no_classification_for_the_single_variant": 2.0,
    "likely_benign": 2.0,
    "benign/likely_benign": 1.5,
    "benign": 1.0,
}


_BENIGN_CLNSIGS = frozenset({"benign", "likely_benign", "benign/likely_benign"})


def _normalize_clnsig(value: str) -> str:
    return value.strip().lower().replace(" ", "_")


def _magnitude(clnsig: str) -> float:
    return _CLNSIG_MAGNITUDE.get(_normalize_clnsig(clnsig), 5.0)


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
    (rsid, significance, condition, gene, alt). The surviving row carries:

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

    Input order is preserved for the kept rows: each group is replaced
    by its winner at the position of the group's first occurrence.
    """
    from dataclasses import replace

    groups: dict[tuple[str, str, str, str, str], list[Annotation]] = {}
    for a in annotations:
        key = (a.rsid, a.significance, a.condition, a.gene, a.alt)
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
    ) -> None:
        """Resolve per-build SQLite cache paths within `data_dir`.

        `builds` selects which builds this annotator instance manages.
        Default is both GRCh37 and GRCh38. Passing a single-element
        tuple (e.g. `("GRCh38",)`) restricts setup/refresh to that
        build — used by the CLI's `--build` flag.

        `include_benign` controls whether Benign/Likely_benign annotations
        are emitted. Default False suppresses them (ADR-0008 amendment).
        """
        super().__init__(data_dir)
        self._builds = tuple(builds)
        self._include_benign = include_benign
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

    def _connection(self, build: str) -> sqlite3.Connection | None:
        """Return a lazy connection to the per-build cache, or None if missing.

        First-touch ensures ``idx_clinvar_position`` exists. ``bulk_resolve_rsids``
        (introduced for rsID-less VCFs from variant callers) joins on
        ``(chromosome, position)``; without the index that would table-scan
        a ~5M-row cache per query. CREATE INDEX IF NOT EXISTS is cheap when
        the index already exists, and a one-time several-second cost when
        migrating a pre-existing cache.
        """
        if build not in self._db_paths:
            return None
        if build not in self._conns:
            db_path = self._db_paths[build]
            if not db_path.exists():
                return None
            conn = sqlite3.connect(db_path)
            try:
                conn.execute(
                    "CREATE INDEX IF NOT EXISTS idx_clinvar_position "
                    "ON clinvar_variants(chromosome, position)"
                )
            except sqlite3.OperationalError as exc:
                # Concurrent process holds a write lock. The next opener
                # retries the IF NOT EXISTS. Position queries still work
                # via table scan in the meantime — degraded perf, not
                # incorrectness.
                logger.debug("clinvar(%s) position-index migration deferred: %s", build, exc)
            self._conns[build] = conn
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
        """
        signal = self._fetch_remote_signal_for_tsv()
        if signal is None:
            msg = (
                "clinvar: cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)

        vs_url = _manager_module.CLINVAR_VARIANT_SUMMARY_URL
        ss_url = _manager_module.CLINVAR_SUBMISSION_SUMMARY_URL
        vs_path = self.data_dir / "variant_summary.txt.gz"
        ss_path = self.data_dir / "submission_summary.txt.gz"

        download(vs_url, vs_path)
        download(ss_url, ss_path)
        # We don't md5-verify variant_summary; signal is calibrated against
        # submission_summary (the more frequently-updated file) and that
        # endpoint's md5 is what we hashed below. variant_summary integrity
        # is asserted structurally by the loader (rows that fail to parse
        # are skipped, mismatched assemblies are dropped).
        try:
            verify_file_hash(ss_path, "md5", signal.removeprefix("md5:"))
            for build in self._builds:
                self._setup_one_from_tsv(vs_path, ss_path, build, signal, vs_url)
        finally:
            for tmp in (vs_path, ss_path):
                try:
                    tmp.unlink()
                except FileNotFoundError:
                    pass
                except OSError:
                    logger.warning("Could not remove staged TSV at %s", tmp)

    def _setup_one_from_tsv(
        self,
        vs_path: Path,
        ss_path: Path,
        build: str,
        signal: str,
        source_url: str,
    ) -> None:
        """Run the per-build cache build against an already-downloaded TSV pair."""
        load_clinvar_tsv(
            vs_path,
            ss_path,
            self._db_paths[build],
            build,
            source_url=source_url,
            remote_signal=signal,
        )
        # load_clinvar_tsv stamps the record name as "clinvar"; for
        # backward compatibility with the per-build record-name scheme
        # (clinvar.GRCh37 / clinvar.GRCh38), rename the row in place.
        # This keeps is_ready() / version() per-build dispatch working
        # against the existing database_versions table layout.
        import contextlib
        import sqlite3

        per_build_name = clinvar_record_name(build)
        if per_build_name != "clinvar":
            with contextlib.closing(sqlite3.connect(self._db_paths[build])) as conn:
                conn.execute(
                    "UPDATE database_versions SET name = ? WHERE name = 'clinvar'",
                    (per_build_name,),
                )
                conn.commit()

    def is_ready(self) -> bool:
        """True iff EVERY managed build has a populated, version-stamped cache.

        Checks ``local_version_tag`` for the current interpreter version.
        Pre-mechanism caches (tag missing or baked into ``remote_signal``)
        are self-healed once via ``stamp_existing_clinvar_cache``.
        """
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is None:
                return False
            tag = info.get("local_version_tag") or ""
            if tag == f"iv:{CLINVAR_INTERPRETER_VERSION}":
                continue
            if stamp_existing_clinvar_cache(self._db_paths[build]):
                continue
            return False
        return True

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
        cursor = conn.execute(
            "SELECT DISTINCT rsid, ref, alt, allele_id FROM clinvar_variants WHERE length(ref) = 1"
        )
        out: dict[str, str] = {}
        # rsID -> VariationIDs whose row is REF == ALT: this build's
        # reference carries their variant allele.
        reference_carried: dict[str, set[int]] = {}
        for rsid, ref, alt, variation_id in cursor:
            if ref == alt:
                reference_carried.setdefault(rsid, set()).add(variation_id)
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
            normal = self._normal_alleles_on_other_build(build, reference_carried)
            for rsid in reference_carried:
                if rsid in normal:
                    out[rsid] = normal[rsid]
                else:
                    out.pop(rsid, None)
        self._normal_alleles[build] = normal
        return out

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
        self, build: str, reference_carried: dict[str, set[int]]
    ) -> dict[str, str]:
        """The normal allele at each reference-carried rsID, from the other build.

        For each rsID, the single-base REF of the same VariationIDs on the
        other build where that row's REF differs from its ALT. An rsID with
        no such row, or more than one candidate REF, is left out.
        """
        other = BUILD_GRCH38 if build == BUILD_GRCH37 else BUILD_GRCH37
        conn = self._connection(other)
        if conn is None or not self._cache_is_current(other):
            return {}
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
        body = fetch_remote_text(_manager_module.CLINVAR_SUBMISSION_SUMMARY_URL + ".md5")
        if not body:
            return None
        first_token = body.strip().split(None, 1)[0] if body.strip() else ""
        if not _MD5_HEX_RE.fullmatch(first_token):
            # CDN error page, redirect interstitial, or empty body. Treat
            # as a transient signal failure rather than poisoning the
            # cache: callers handle `None` as "freshness unknown, skip"
            # in `db update`, and `setup()` raises rather than passing
            # garbage to `verify_file_hash` (which would delete the TSV).
            logger.warning(
                "clinvar: submission_summary.md5 endpoint returned a body whose "
                "first token is not a 32-char hex digest (got %r); treating as "
                "no signal",
                first_token[:32],
            )
            return None
        return f"md5:{first_token}"

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
        if variant.is_no_call:
            return []
        conn = self._connection(variant.build)
        if conn is None:
            return []
        rows = self._rows_at_loci(conn, [variant]).get(variant_locus(variant), [])
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

        The returned ``{(chrom, pos, ref, alt): rsid}`` mapping carries
        the resolved coordinates so the pipeline's enrichment phase can
        do position-keyed gnomAD / AlphaMissense lookups. Only variants
        carrying a single alt contribute an entry: a ``1/2`` call has no
        single allele to stand for the rsID.
        """
        by_build: dict[str, list[Variant]] = {}
        for v in variants:
            if v.is_no_call or not v.chromosome or v.position <= 0:
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
            conn = self._connection(build)
            if conn is None:
                continue
            rows_at = self._rows_at_loci(conn, build_variants)
            for v in build_variants:
                rows = rows_at.get(variant_locus(v))
                if not rows:
                    continue
                hit = resolve_rsid_at_locus(
                    v.ref,
                    v.allele1,
                    v.allele2,
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
            if not v.is_no_call:
                by_build.setdefault(v.build, []).append(v)
        for build, build_variants in by_build.items():
            conn = self._connection(build)
            if conn is None:
                continue
            rows_at = self._rows_at_loci(conn, build_variants)
            known = self._loci_of_rsids(conn, build_variants)
            for variant in build_variants:
                rows = rows_at.get(variant_locus(variant))
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
        out: dict[tuple[str, int], list[tuple]] = {}
        for chrom, positions in by_chrom.items():
            position_list = sorted(positions)
            for start in range(0, len(position_list), _BATCH_CHUNK):
                chunk = position_list[start : start + _BATCH_CHUNK]
                placeholders = ",".join("?" * len(chunk))
                cursor = conn.execute(
                    f"SELECT chromosome, position, ref, alt, clinical_significance, "
                    f"condition, gene, review_status, allele_id, rsid "
                    f"FROM clinvar_variants "
                    f"WHERE chromosome = ? AND position IN ({placeholders}) ORDER BY id",
                    (chrom, *chunk),
                )
                for row in cursor:
                    out.setdefault((row[0], row[1]), []).append(row)
        return out

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

        - REF == ALT means this build's reference carries the variant
          allele (F5 c.1601G>A is T>T on GRCh37, whose reference base is
          the Leiden allele). The loader keeps no other REF == ALT rows:
          it drops normal-allele records (HGVS ``c.N=``) written REF ==
          ALT and REF == ALT rows that are not SNVs (cache ``iv:5``). Such
          a row annotates a user who carries its ALT, counting copies of
          it, like any other row.
        - The row's REF agrees with ``variant.ref`` when the input supplies
          one (directly, or as the complement for a reverse-stranded
          input). Together with the locus this identifies the variant by
          ``(chrom, pos, ref, alt)`` rather than by rsID.
        - ADR-0011: ClinVar indels never match single-base array calls.
        - ADR-0035 / ADR-0007: the user carries the ALT, directly or
          via a confirmed strand flip.
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

        The Annotation's rsID is the input's when it has one, otherwise
        the matching row's. identical per-SCV rows collapse.
        """
        if variant.is_no_call:
            return []
        if variant.ref is None and not rsid_loci:
            return []
        if rsid_loci and variant_locus(variant) not in rsid_loci:
            return []
        user_is_multibase = len(variant.allele1) > 1 or len(variant.allele2) > 1
        # ADR-0023: report the user's actual diploid call consistently
        # across annotators, not the matched ALT base alone.
        user_diploid = _user_diploid(variant)
        annotations: list[Annotation] = []
        for row in rows:
            (
                _chrom,
                _pos,
                ref,
                alt,
                clnsig,
                condition,
                gene,
                review_status,
                allele_id,
                row_rsid,
            ) = row
            if variant.ref is not None and variant.ref not in (ref, complement(ref)):
                continue
            clinvar_is_indel = len(ref) > 1 or len(alt) > 1
            if clinvar_is_indel and not user_is_multibase:
                continue
            if not strand_aware_carrier_match(
                variant.ref, variant.allele1, variant.allele2, ref, alt
            ):
                continue
            sig_label = _normalize_clnsig(clnsig) if clnsig else "unknown"
            if not self._include_benign and sig_label in _BENIGN_CLNSIGS:
                continue
            description = (
                f"ClinVar classifies this allele as "
                f"{clnsig.replace('_', ' ') if clnsig else 'unknown significance'}"
            )
            references = [f"clinvar:allele/{allele_id}"] if allele_id else []
            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid if variant.rsid.startswith("rs") else row_rsid,
                    significance=f"clinvar_{sig_label}",
                    category="clinical",
                    magnitude=_magnitude(clnsig),
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition="" if not condition or condition == "." else condition,
                    gene=gene or "",
                    review_status=review_status or "",
                    alt=alt,
                    locus_ref=ref,
                    carried_alleles=carried_in_source_orientation(
                        variant.allele1, variant.allele2, alt
                    ),
                    locus=variant_locus(variant),
                    reference_carries_alt=ref == alt,
                )
            )
        return _collapse_identical_scvs(annotations)


# Column positions in ``_rows_at_loci`` rows used outside the unpacking.
_R_REF = 2
_R_ALT = 3
_R_RSID = 9


def _user_diploid(variant: Variant) -> str:
    """Render the user's diploid call as a sorted two-letter string.

    Used by ClinVar and ClinPGx so the report's "Genotype" column shows
    the same shape for every annotation regardless of source (ADR-0023).
    SNV: `("G", "A") -> "AG"`. Indel passthrough is verbatim.
    """
    a1, a2 = variant.allele1, variant.allele2
    if len(a1) == 1 and len(a2) == 1:
        return "".join(sorted((a1, a2)))
    return f"{a1}/{a2}"

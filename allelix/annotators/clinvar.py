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

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases import manager as _manager_module
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import (
    download,
    fetch_remote_text,
    get_database_info,
    load_clinvar_vcf,
    stamp_existing_clinvar_cache,
    verify_file_hash,
)
from allelix.models import Annotation

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

CLINVAR_SUPPORTED_BUILDS: tuple[str, ...] = ("GRCh37", "GRCh38")

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999

# GH #21: a remote .md5 endpoint can return an HTML error page on a
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


def _vcf_filename_for_url(url: str) -> str:
    """Pick the right local filename suffix based on the URL."""
    return "clinvar.vcf.gz" if url.endswith(".gz") else "clinvar.vcf"


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
        """Download each managed build's ClinVar VCF and ingest atomically."""
        for build in self._builds:
            self._setup_one(build)

    def _setup_one(self, build: str) -> None:
        url = _manager_module.CLINVAR_URL_BY_BUILD[build]
        signal = self._fetch_remote_signal_for(build)
        if signal is None:
            msg = (
                f"clinvar ({build}): cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)
        vcf_path = self.data_dir / _vcf_filename_for_url(url)
        download(url, vcf_path)
        try:
            verify_file_hash(vcf_path, "md5", signal.removeprefix("md5:"))
            load_clinvar_vcf(
                vcf_path,
                self._db_paths[build],
                source_url=url,
                remote_signal=signal,
                record_name=clinvar_record_name(build),
            )
        finally:
            try:
                vcf_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                logger.warning("Could not remove staged VCF at %s", vcf_path)

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

    def reference_for(self, rsid: str, build: str) -> str | None:
        """Return ClinVar's single-base REF allele for `rsid` in `build`, or None.

        ADR-0023: ClinPGx's primary non-finding filter calls this. If the
        return value matches both of the user's alleles, the user is
        homozygous reference and the ClinPGx annotation is a non-finding.

        Lazily builds an in-memory `(rsid -> REF)` map per build on first
        call so subsequent lookups are O(1). Multi-base REFs (indels) are
        skipped — array-based parsers can't call indels, so a multi-base
        REF can't validly suppress a single-base genotype.

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
        rows = conn.execute(
            "SELECT DISTINCT rsid, ref FROM clinvar_variants WHERE length(ref) = 1"
        ).fetchall()
        out: dict[str, str] = {}
        for rsid, ref in rows:
            # If a rsid has multiple single-base REFs (shouldn't happen at
            # one position but defending against future data shapes), keep
            # the first.
            if rsid not in out:
                out[rsid] = ref
        return out

    def fetch_remote_signal(self) -> str | None:
        r"""Composite freshness signal across managed builds.

        Format: `"GRCh37:md5:<hex>|GRCh38:md5:<hex>"`. Returns None if
        ANY managed build's signal probe fails — the CLI then prints
        "can't verify" and skips refresh per ADR-0012's policy.
        """
        parts: list[str] = []
        for build in self._builds:
            sig = self._fetch_remote_signal_for(build)
            if sig is None:
                return None
            parts.append(f"{build}:{sig}")
        return "|".join(parts) if parts else None

    @staticmethod
    def _fetch_remote_signal_for(build: str) -> str | None:
        body = fetch_remote_text(_manager_module.CLINVAR_URL_BY_BUILD[build] + ".md5")
        if not body:
            return None
        first_token = body.strip().split(None, 1)[0] if body.strip() else ""
        if not _MD5_HEX_RE.fullmatch(first_token):
            # CDN error page, redirect interstitial, or empty body. Treat
            # as a transient signal failure rather than poisoning the
            # cache: callers handle `None` as "freshness unknown, skip"
            # in `db update`, and `setup()` raises rather than passing
            # garbage to `verify_file_hash` (which would delete the VCF).
            logger.warning(
                "clinvar(%s): .md5 endpoint returned a body whose first token "
                "is not a 32-char hex digest (got %r); treating as no signal",
                build,
                first_token[:32],
            )
            return None
        return f"md5:{first_token}"

    def cached_remote_signal(self) -> str | None:
        """Composite cached signal across managed builds. None if any missing."""
        parts: list[str] = []
        for build in self._builds:
            info = get_database_info(self._db_paths[build], clinvar_record_name(build))
            if info is None or info["remote_signal"] is None:
                return None
            sig = info["remote_signal"]
            if not sig:
                return None
            parts.append(f"{build}:{sig}")
        return "|".join(parts) if parts else None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return ClinVar annotations whose REF/ALT matches the user's genotype.

        ADR-0007 carrier rule: an entry triggers only if `variant.allele1`
        or `variant.allele2` equals the entry's ALT allele.
        ADR-0011 indel-anchor protection: array-based parsers report
        single-base genotypes; ClinVar's anchor-base indel encoding
        does not match those by string equality.
        ADR-0021: dispatch by `variant.build`. If the matching cache is
        absent, the variant is skipped silently — the user already saw
        the analyze-time build warning.
        """
        if variant.is_no_call:
            return []
        conn = self._connection(variant.build)
        if conn is None:
            return []
        rows = conn.execute(
            "SELECT chromosome, position, ref, alt, clinical_significance, "
            "condition, gene, review_status, allele_id "
            "FROM clinvar_variants WHERE rsid = ?",
            (variant.rsid,),
        ).fetchall()
        annotations: list[Annotation] = []
        carrier_alleles = {variant.allele1, variant.allele2}
        user_is_multibase = len(variant.allele1) > 1 or len(variant.allele2) > 1
        # ADR-0023: report the user's actual diploid call consistently
        # across annotators, not the matched ALT base alone.
        user_diploid = _user_diploid(variant)
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
            ) = row
            clinvar_is_indel = len(ref) > 1 or len(alt) > 1
            if clinvar_is_indel and not user_is_multibase:
                continue
            if alt not in carrier_alleles:
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
                    rsid=variant.rsid,
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
                )
            )
        return annotations

    def bulk_resolve_rsids(self, variants: list[Variant]) -> dict[tuple[str, int, str, str], str]:
        """Resolve rsIDs by (chromosome, position) for variants with no ID.

        Real-world VCFs from variant callers (GATK HaplotypeCaller,
        DeepVariant) emit ``ID=.`` — no rsID. The annotation pipeline is
        rsID-keyed, so without resolution every annotator returns zero
        hits. This method queries ClinVar by ``(chromosome, position)``,
        disambiguates multi-allelic rows by matching ClinVar's ALT
        against the user's carrier alleles (ADR-0007 carrier rule),
        and mutates each input variant's ``rsid`` field in-place when
        a match is found.

        Groups by ``variant.build`` and dispatches per-build (ADR-0021).
        Variants not in ClinVar are left with empty ``rsid`` —
        downstream rsID-keyed annotators (ClinPGx, GWAS, SNPedia)
        correctly skip them.

        The returned ``{(chrom, pos, ref, alt): rsid}`` mapping carries
        the resolved coordinates so callers (e.g. the pipeline's
        enrichment phase) can do position-based gnomAD / AlphaMissense
        lookups for variants whose rsID-keyed enrichment misses.
        """
        by_build: dict[str, list[Variant]] = {}
        for v in variants:
            if v.is_no_call or not v.chromosome or v.position <= 0:
                continue
            by_build.setdefault(v.build, []).append(v)

        resolved: dict[tuple[str, int, str, str], str] = {}
        for build, build_variants in by_build.items():
            conn = self._connection(build)
            if conn is None:
                continue
            by_chrom: dict[str, set[int]] = {}
            for v in build_variants:
                by_chrom.setdefault(v.chromosome, set()).add(v.position)
            clinvar_rows: dict[tuple[str, int], list[tuple[str, str, str]]] = {}
            for chrom, positions in by_chrom.items():
                position_list = list(positions)
                for start in range(0, len(position_list), _BATCH_CHUNK):
                    chunk = position_list[start : start + _BATCH_CHUNK]
                    placeholders = ",".join("?" * len(chunk))
                    cursor = conn.execute(
                        f"SELECT chromosome, position, ref, alt, rsid "
                        f"FROM clinvar_variants "
                        f"WHERE chromosome = ? AND position IN ({placeholders})",
                        (chrom, *chunk),
                    )
                    for c, p, ref, alt, rsid in cursor:
                        if not rsid:
                            continue
                        clinvar_rows.setdefault((c, p), []).append((ref, alt, rsid))
            # Sort each position's rows for deterministic carrier-match
            # selection. Without this, the row iteration order depends on
            # SQLite's physical row order, which can change after VACUUM,
            # version upgrades, or even per-run cache loads. Tuple sort is
            # lexicographic over (ref, alt, rsid) — stable across runs.
            for rows in clinvar_rows.values():
                rows.sort()
            for v in build_variants:
                rows = clinvar_rows.get((v.chromosome, v.position))
                if not rows:
                    continue
                user_alleles = {v.allele1, v.allele2}
                # Stamp only on a UNIQUE consistent match. Multiple matches
                # mean the Variant model can't distinguish them (the 1/1
                # hom-alt case at a shared-anchor position is the canonical
                # example — both rows pass the subset check and there's no
                # REF info in the Variant to break the tie). Abstaining is
                # safer than picking by sort order: the variant just
                # doesn't get rsid-keyed annotations rather than getting
                # the wrong ones. Zero matches → abstain. Single match →
                # stamp. Two+ → abstain. A future Variant.ref field (v2.1+)
                # would let the resolver disambiguate the ambiguous cases.
                matches = [
                    (ref, alt, rsid)
                    for ref, alt, rsid in rows
                    if ref != alt and user_alleles.issubset({ref, alt})
                ]
                if len(matches) == 1:
                    ref, alt, rsid = matches[0]
                    v.rsid = rsid
                    resolved[(v.chromosome, v.position, ref, alt)] = rsid
        return resolved

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Bulk-annotate via chunked ``WHERE rsid IN (...)`` queries per build.

        ClinVar holds one SQLite cache per build (ADR-0021). Variants
        are grouped by ``variant.build`` and each group queried against
        its own connection. Per-variant filters (no-call rejection,
        carrier-allele match, indel-anchor protection per ADR-0011,
        benign suppression) match ``annotate`` exactly.

        Variants whose build has no cache loaded are skipped silently —
        the build warning was already emitted at analyze time.
        """
        variants_list = [v for v in variants if not v.is_no_call]
        if not variants_list:
            return

        # Group by build; ClinVar dispatches by build per ADR-0021
        by_build: dict[str, list[Variant]] = {}
        for v in variants_list:
            by_build.setdefault(v.build, []).append(v)

        # Per-build chunked SQL — collect rsid → rows per build
        rows_by_build_rsid: dict[tuple[str, str], list[tuple]] = {}
        for build, build_variants in by_build.items():
            conn = self._connection(build)
            if conn is None:
                continue
            rsids = list({v.rsid for v in build_variants})
            for start in range(0, len(rsids), _BATCH_CHUNK):
                chunk = rsids[start : start + _BATCH_CHUNK]
                placeholders = ",".join("?" * len(chunk))
                cursor = conn.execute(
                    f"SELECT rsid, chromosome, position, ref, alt, clinical_significance, "
                    f"condition, gene, review_status, allele_id "
                    f"FROM clinvar_variants WHERE rsid IN ({placeholders})",
                    chunk,
                )
                for row in cursor:
                    rows_by_build_rsid.setdefault((build, row[0]), []).append(row[1:])

        for variant in variants_list:
            rows = rows_by_build_rsid.get((variant.build, variant.rsid))
            if not rows:
                continue
            carrier_alleles = {variant.allele1, variant.allele2}
            user_is_multibase = len(variant.allele1) > 1 or len(variant.allele2) > 1
            user_diploid = _user_diploid(variant)
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
                ) = row
                clinvar_is_indel = len(ref) > 1 or len(alt) > 1
                if clinvar_is_indel and not user_is_multibase:
                    continue
                if alt not in carrier_alleles:
                    continue
                sig_label = _normalize_clnsig(clnsig) if clnsig else "unknown"
                if not self._include_benign and sig_label in _BENIGN_CLNSIGS:
                    continue
                description = (
                    f"ClinVar classifies this allele as "
                    f"{clnsig.replace('_', ' ') if clnsig else 'unknown significance'}"
                )
                references = [f"clinvar:allele/{allele_id}"] if allele_id else []
                yield Annotation(
                    source=self.name,
                    rsid=variant.rsid,
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
                )


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

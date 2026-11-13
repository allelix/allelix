# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinPGx annotator. Source-attributed pharmacogenomic annotations (ADR-0003)."""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import (
    Annotator,
    LicenseDescriptor,
    carried_alt,
    forward_alleles_at_call,
    is_clinvar_homref,
    variant_locus,
)
from allelix.databases._versions import PHARMGKB_INTERPRETER_VERSION
from allelix.databases.cache_state import CacheState, CacheStatus, cache_state, connect_readonly
from allelix.databases.cpic_loader import (
    fetch_cpic_allele_functions,
    fetch_cpic_remote_signal,
)
from allelix.databases.manager import (
    PROBE_ERRORS,
    download,
    get_database_info,
    head_request_headers,
    http_signal,
    update_remote_signal,
)
from allelix.databases.pharmgkb_loader import (
    _REQUIRED_PHARMGKB_COLUMNS,
    _REQUIRED_PHARMGKB_TABLES,
    PHARMGKB_CLINICAL_URL,
    PHARMGKB_DB_FILENAME,
    PHARMGKB_MIN_ROWS,
    _normalize_genotype,
    load_pharmgkb_tsv,
    pharmgkb_content_hash,
    schema_is_current,
    stored_content_hash,
)
from allelix.databases.staging import owned_directory
from allelix.models import Annotation, diploid_label
from allelix.utils.allele import (
    STRAND_COMPLEMENT,
    STRAND_DIRECT,
    forward_strand_reading,
    is_single_base,
    strand_aware_genotype_match,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

    from allelix.annotators.base import ForwardSiteProvider, NormalAlleleProvider
    from allelix.annotators.rsid_identity import ForwardSite
    from allelix.models import Variant

logger = logging.getLogger(__name__)

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999

# Allelix-derived magnitude scoring from ClinPGx Level of Evidence. See ADR-0008.
# 1A is the strongest evidence (CPIC guideline-backed); 4 is the weakest.
_LOE_MAGNITUDE: dict[str, float] = {
    "1a": 9.0,
    "1b": 8.0,
    "2a": 7.0,
    "2b": 6.0,
    "3": 4.0,
    "4": 2.0,
}


def _magnitude(level_of_evidence: str) -> float:
    return _LOE_MAGNITUDE.get(level_of_evidence.strip().lower(), 5.0)


def _without_cpic(signal: str) -> str:
    """The composite signal with its CPIC part marked unavailable; the ClinPGx part kept.

    The CPIC part is the last one: a quoted ETag before it may itself hold
    ``|cpic:`` (ADR-0012).
    """
    head, sep, _ = signal.rpartition("|cpic:")
    return f"{head if sep else signal}|cpic:unavailable"


def signal_lacks_cpic(signal: str | None) -> bool:
    """True when a ClinPGx signal's CPIC part says CPIC was unavailable.

    The CPIC part is the last ``|cpic:`` field, or the whole signal when it
    is only the mark; a quoted ETag before it may hold the same text.
    """
    if not signal:
        return False
    _, sep, cpic = signal.rpartition("|cpic:")
    return cpic == "unavailable" if sep else signal == "cpic:unavailable"


def _has_cpic_rows(db_path: Path) -> bool:
    """True when the cache holds at least one CPIC allele-function row (read-only)."""
    try:
        with contextlib.closing(connect_readonly(db_path)) as conn:
            return (
                conn.execute("SELECT 1 FROM pharmgkb_allele_function LIMIT 1").fetchone()
                is not None
            )
    except sqlite3.Error:
        return True  # an unreadable table is cache_status's to report, not a CPIC repair


class PharmGKBAnnotator(Annotator):
    """Annotates variants with ClinPGx's curated drug-gene-variant associations."""

    # `name` is the internal stable identifier — appears in Annotation.source,
    # database_versions.name, and the cache filename. Kept as "pharmgkb" so
    # existing user caches keep working after the ClinPGx rebrand. Display
    # surfaces (display_name, attribution, description text, source_url) are
    # the rebranded ones.
    name: ClassVar[str] = "pharmgkb"
    display_name: ClassVar[str] = "ClinPGx"
    attribution: ClassVar[str] = "ClinPGx"
    requires_download: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="CC-BY-SA-4.0",
        license_url="https://creativecommons.org/licenses/by-sa/4.0/",
        attribution_text=(
            "Pharmacogenomic annotations sourced from ClinPGx (formerly PharmGKB),"
            " used under CC BY-SA 4.0."
        ),
        source_url="https://www.clinpgx.org",
        commercial_ok=True,
    )

    def __init__(
        self,
        data_dir: Path,
        clinvar_ref_provider: Callable[[str, str, tuple[str, int]], str | None] | None = None,
        palindromic_provider: Callable[[set[str]], set[str]] | None = None,
        normal_allele_provider: NormalAlleleProvider | None = None,
        reference_prefetch: Callable[[list[Variant]], None] | None = None,
        forward_site_provider: ForwardSiteProvider | None = None,
    ) -> None:
        """Resolve the ClinPGx SQLite cache path within `data_dir`.

        `clinvar_ref_provider` is a `(rsid, build, locus) -> ref_base | None`
        callable used by the primary non-finding filter (ADR-0023). In
        production it's wired to `ClinVarAnnotator.reference_at`, which
        answers only where ClinVar places the rsID at the call's locus.
        None disables the REF check and falls back to the cache's CPIC-based
        `is_nonfinding` flag for all suppression — the v0.7.1 behavior.

        `palindromic_provider` returns the subset of given rsIDs whose site
        the identity table flags palindromic; complement (strand-flip)
        matching is skipped there. In production it's wired to
        `RsidIdentityAnnotator.palindromic_rsids`. None keeps the per-genotype
        palindrome guard only.

        `normal_allele_provider` gives the normal allele where the build's
        reference carries a ClinVar variant allele, for deriving `alt`
        (`base.carried_alt`). In production it's wired to
        `ClinVarAnnotator.normal_allele_for`. None derives `alt` against
        `Variant.ref` everywhere.

        `reference_prefetch` is given each batch before the two providers
        above are asked about it, so they can answer from one bulk lookup
        (`ReferenceCarriedFallback.prefetch`). None asks them per call.

        `forward_site_provider` gives each rsID's locus and forward alleles
        on a build (`RsidIdentityAnnotator.forward_sites`). Where all of an
        rsID's source genotypes, non-findings included, read only as their
        complement at the call's locus, any genotype there matches by its
        complement, homozygotes included (ADR-0035). None leaves the
        complement match to the call's REF alone.
        """
        super().__init__(data_dir)
        self._db_path = data_dir / PHARMGKB_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._clinvar_ref_provider = clinvar_ref_provider
        self._palindromic_provider = palindromic_provider
        self._normal_allele_provider = normal_allele_provider
        self._reference_prefetch = reference_prefetch
        self._forward_site_provider = forward_site_provider
        # Calls one row met as written and another as its complement, with
        # no site evidence to choose: withheld, counted once each (ADR-0035).
        self._ambiguous_strand_calls = 0
        # Set by setup(): the download held the content of the last build,
        # so the cache was kept and only its signal stored (ADR-0068).
        self.content_unchanged = False

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = connect_readonly(self._db_path)
        return self._conn

    def _site_genotypes(self, rsids: list[str]) -> dict[str, list[str]]:
        """Every source genotype of each rsID, non-findings included.

        The site's orientation evidence: a non-finding row is not reported,
        but its alleles are still the source's alleles at the site.
        """
        conn = self._connection()
        out: dict[str, list[str]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            for rsid, genotype in conn.execute(
                f"SELECT rsid, genotype FROM pharmgkb_annotations WHERE rsid IN ({placeholders})",
                chunk,
            ):
                out.setdefault(rsid, []).append(genotype)
        return out

    def run_counts(self) -> dict[str, int]:
        """``ambiguous_strand_calls``: the calls withheld from two readings (ADR-0035)."""
        return (
            {"ambiguous_strand_calls": self._ambiguous_strand_calls}
            if self._ambiguous_strand_calls
            else {}
        )

    def _site_readings(self, variants: list[Variant]) -> list[str | None]:
        """For each call, in order: how its site reads the rsID's genotypes (``_site_reading``).

        Decided per call, never per rsID: the identity table's site on the
        call's build must be at the call's own locus (a call without a
        position takes the rsID's one site). All None when no identity
        table is wired (ADR-0035).
        """
        if self._forward_site_provider is None or not variants:
            return [None] * len(variants)
        genotypes = self._site_genotypes(sorted({v.rsid for v in variants}))
        sites_by_build: dict[str, dict[str, ForwardSite]] = {}
        for build in sorted({v.build for v in variants}):
            wanted = {v.rsid for v in variants if v.build == build and v.rsid in genotypes}
            sites_by_build[build] = self._forward_site_provider(wanted, build) if wanted else {}
        return [
            _site_reading(v, genotypes.get(v.rsid, []), sites_by_build[v.build]) for v in variants
        ]

    def setup(self, *, force: bool = False) -> None:
        """Download ClinPGx clinical annotations + CPIC allele functions, ingest atomically.

        Two sources are fetched: ClinPGx's `clinicalAnnotations.zip`
        (the annotation rows + per-genotype rows) and CPIC's API (the
        structured per-allele function table per ADR-0020). The primary
        non-finding filter is the ClinVar REF check (ADR-0023); CPIC's
        per-allele function table is the secondary fallback for rsids
        ClinVar doesn't catalog.

        This is also how a cache built by an older interpreter version is
        rebuilt: ``is_ready()`` reports it stale and ``db update`` calls
        this, so the CPIC table is fetched again with the annotations.

        ``force`` (``db update --force``) refreshes when the freshness signal
        is unavailable. ClinPGx publishes no checksum, so the checks are the
        same either way; no signal is stored, so the next ordinary update
        refreshes again.

        When the download and CPIC's data hash the same as the last build's
        and the cache is current, the cache is kept: only the new signal is
        stored, and ``content_unchanged`` is set (ADR-0068). The signal is
        stored only on the cache holding the content it checked: a cache
        another refresh installed meanwhile keeps its own signal.
        """
        url = PHARMGKB_CLINICAL_URL
        signal = self.fetch_remote_signal()
        if signal is None and not force:
            msg = (
                "cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)
        # This refresh's own download directory: the archive it hashes and
        # loads is the one it fetched; another refresh can neither replace
        # nor remove it. It is removed when the refresh ends.
        with owned_directory(self.data_dir, "pharmgkb") as staging:
            zip_path = staging / "clinicalAnnotations.zip"
            # No content-hash verification: ClinPGx publishes no checksum and
            # the content is mutable, so there is nothing to pin or fetch.
            # TLS + Content-Length truncation guard only. See ADR-0029.
            download(url, zip_path)
            try:
                cpic_lookup = fetch_cpic_allele_functions()
            except PROBE_ERRORS as exc:
                # A transport failure, a cut reply, or a partial table
                # (CpicIncompleteError): the cache is built without CPIC rows.
                logger.warning(
                    "CPIC API unavailable (%s) -- proceeding without "
                    "allele function data. Non-finding filter degraded.",
                    exc,
                )
                cpic_lookup = {}
            if not cpic_lookup:
                # The stamp says what was built: without CPIC rows the cache is
                # marked "cpic:unavailable", so status says so and the next
                # update retries even when CPIC's probe answered this time. A
                # forced build without a ClinPGx signal keeps only that mark.
                signal = _without_cpic(signal) if signal is not None else "cpic:unavailable"
            content_hash = pharmgkb_content_hash(zip_path, cpic_lookup)
            self.content_unchanged = False
            if (
                content_hash is not None
                and self.cache_status().ready
                and stored_content_hash(self._db_path) == content_hash
            ):
                if update_remote_signal(
                    self._db_path, "pharmgkb", signal, content_hash=content_hash
                ):
                    self.content_unchanged = True
                else:
                    # Another refresh installed different content after the
                    # check; that cache and its own signal stay as they are.
                    logger.warning(
                        "ClinPGx: another update installed a different cache during this "
                        "one; that cache is kept with its own freshness signal."
                    )
                return
            load_pharmgkb_tsv(
                zip_path,
                self._db_path,
                source_url=url,
                remote_signal=signal,
                allele_function_lookup=cpic_lookup,
                min_rows=PHARMGKB_MIN_ROWS,
                content_hash=content_hash,
            )

    def reference_call_rsids(self) -> Iterable[str]:
        """Every rsID with a row here, streamed (ADR-0064); none when the cache is not ready."""
        if not self.is_ready():
            return ()
        rows = self._connection().execute("SELECT DISTINCT rsid FROM pharmgkb_annotations")
        return (rsid for (rsid,) in rows if rsid)

    def is_ready(self) -> bool:
        """True iff a ClinPGx SQLite cache exists with current schema and interpreter stamp.

        Never touches the network or the file, so ``analyze`` stays
        offline and read-only (ADR-0048). When the interpreter version has
        bumped, the cache is not ready and ``allelix db update`` rebuilds
        it, fetching the ClinPGx annotations and the CPIC allele-function
        table again (``setup``).
        """
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        """State of the ClinPGx cache: tables, columns and interpreter tag (ADR-0048)."""
        return cache_state(
            self._db_path,
            "pharmgkb",
            tag=f"iv:{PHARMGKB_INTERPRETER_VERSION}",
            # The annotation table first: it must hold rows.
            tables=("pharmgkb_annotations", *sorted(_REQUIRED_PHARMGKB_TABLES)),
            columns={"pharmgkb_annotations": _REQUIRED_PHARMGKB_COLUMNS},
        )

    def prepare_cache(self) -> bool:
        """Self-heal a pre-mechanism cache whose baked stamp is the current version.

        Pre-mechanism caches (tag missing, version baked into
        ``remote_signal``) get a one-shot stamp update when the baked
        version is the current one and the schema is current; any other
        outdated cache is left for ``db update`` to rebuild.
        """
        if self.cache_status().state is not CacheState.OUTDATED:
            return False
        info = get_database_info(self._db_path, "pharmgkb")
        if info is None or info.get("local_version_tag"):
            return False
        if not schema_is_current(self._db_path):
            return False
        return _stamp_existing_pharmgkb_cache(self._db_path)

    def version(self) -> str | None:
        """Cached database version (download date, or version supplied to load)."""
        info = get_database_info(self._db_path, "pharmgkb")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Number of (rsid, genotype) annotation rows in the cache."""
        info = get_database_info(self._db_path, "pharmgkb")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open. Safe to call multiple times."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Composite freshness signal for ClinPGx + CPIC (ADR-0020).

        The signal format is `pgkb:<pgkb-signal>|cpic:<cpic-signal>`.

          - ClinPGx portion: its ETag and Last-Modified, each when present
            (ADR-0012).
          - CPIC portion: latest `change_log` date from CPIC's API,
            or ``unavailable`` if the CPIC probe fails.

        Returns None only when ClinPGx itself is unreachable. CPIC
        failure is non-fatal: the signal carries ``cpic:unavailable``
        so the cache is still refreshable (and the mismatch when CPIC
        recovers triggers a re-download automatically).
        """
        headers = head_request_headers(PHARMGKB_CLINICAL_URL)
        if headers is None:
            return None
        pgkb_signal = http_signal(headers)
        if pgkb_signal is None:
            return None

        cpic_signal = fetch_cpic_remote_signal()
        if cpic_signal is None:
            return f"pgkb:{pgkb_signal}|cpic:unavailable"
        return f"pgkb:{pgkb_signal}|cpic:{cpic_signal}"

    def cached_remote_signal(self) -> str | None:
        """Return the remote signal stored at last successful download.

        A cache whose stamp names CPIC data but whose CPIC table is empty
        (2.5.0 stamped the probe's answer even when the data fetch then
        failed) reports ``cpic:unavailable``, so ``db update`` repairs it
        and ``db status`` shows it. One read-only row check; nothing else.
        """
        info = get_database_info(self._db_path, "pharmgkb")
        if not info or not info["remote_signal"]:
            return None
        signal = info["remote_signal"]
        _, sep, cpic = signal.rpartition("|cpic:")
        claims_cpic = bool(sep) and cpic != "unavailable"
        if claims_cpic and not _has_cpic_rows(self._db_path):
            return _without_cpic(signal)
        return signal

    def cpic_unavailable(self) -> bool:
        """True when the cache holds no CPIC allele-function data (ADR-0065).

        The stamp says CPIC was unavailable, or the CPIC table is empty
        (a 2.5.0 stamp claiming CPIC, or no stamp at all). Read-only.
        """
        signal = self.cached_remote_signal()
        if signal is not None and signal.endswith("cpic:unavailable"):
            return True
        return not _has_cpic_rows(self._db_path)

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return ClinPGx annotations for variants the user actually carries.

        Non-finding suppression has two independent signals; either one
        is sufficient to suppress a row:

          1. **ClinVar REF carrier rule (ADR-0023).** If ClinVar has a
             single-base REF for this rsid and the user is homozygous
             for it → suppress before hitting the database.

          2. **CPIC per-allele function (ADR-0020).** The pre-computed
             `is_nonfinding` flag in the cache — set at load time when
             CPIC classifies every user-carried base as Normal function.
             Applied via `AND is_nonfinding = 0` on every query.

        The two checks are additive: ClinVar REF catches genes CPIC
        doesn't cover; CPIC catches rows where both alleles are Normal
        function even when the user isn't homozygous reference per ClinVar
        (e.g. rs1801265 GG in DPYD).

        No-calls and indels are filtered out by `_normalize_genotype()`
        returning None — array-based parsers don't call indels (ADR-0011).
        """
        if variant.is_no_call:
            return []
        user_geno = _normalize_genotype(variant.allele1 + variant.allele2)
        if user_geno is None:
            return []

        if is_clinvar_homref(variant, self._clinvar_ref_provider):
            return []

        # ADR-0035: fetch all rows for the rsid (typically <5) and apply
        # strand-aware genotype match in Python. Forward-stranded users hit
        # the same set of rows as the v0.7-era ``WHERE genotype = ?`` SQL;
        # reverse-stranded users also match when the rsID's source alleles
        # read only as their complement at the call's identity site (any
        # genotype, homozygotes included), or else when Variant.ref confirms
        # orientation AND the source genotype is non-palindromic.
        sql = (
            "SELECT genotype, gene, drugs, phenotype, phenotype_category, "
            "annotation_text, level_of_evidence, score, pgkb_annotation_id "
            "FROM pharmgkb_annotations "
            "WHERE rsid = ? AND is_nonfinding = 0"
        )
        params = (variant.rsid,)

        rows = self._connection().execute(sql, params).fetchall()
        annotations: list[Annotation] = []
        user_diploid = diploid_label(variant)
        palindromic = (
            self._palindromic_provider is not None
            and variant.rsid in self._palindromic_provider({variant.rsid})
        )
        reading = self._site_readings([variant])[0] if rows else None
        met = (
            _genotypes_met(
                variant,
                self._site_genotypes([variant.rsid]).get(variant.rsid, []),
                site_palindromic=palindromic,
                site_reading=reading,
            )
            if rows
            else set()
        )
        if met is None:
            self._ambiguous_strand_calls += 1
            return []
        for row in rows:
            (
                geno,
                gene,
                drugs,
                phenotype,
                _phenotype_category,
                annotation_text,
                level_of_evidence,
                _score,
                pgkb_annotation_id,
            ) = row
            if geno.upper() not in met:
                continue
            sig_label = level_of_evidence.strip().lower() or "unknown"
            description_parts = [f"ClinPGx: {drugs}"] if drugs else ["ClinPGx"]
            if phenotype:
                description_parts.append(phenotype)
            if annotation_text:
                description_parts.append(annotation_text)
            description = " — ".join(description_parts)
            references = (
                [f"pharmgkb:annotation/{pgkb_annotation_id}"] if pgkb_annotation_id else []
            )
            alt, reference_carries_alt = carried_alt(
                variant, variant.allele1, variant.allele2, self._normal_allele_provider
            )
            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"pharmgkb_loe_{sig_label}",
                    category="pharma",
                    magnitude=_magnitude(level_of_evidence),
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=phenotype or "",
                    gene=gene or "",
                    alt=alt,
                    locus_ref=variant.ref,
                    reference_carries_alt=reference_carries_alt,
                    carried_alleles=(variant.allele1, variant.allele2),
                    locus=variant_locus(variant),
                )
            )
        return annotations

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Bulk-annotate via chunked ``WHERE rsid IN (...)`` queries.

        Pre-filters (no-call, indel/non-normalizable genotype, ClinVar
        REF homozygote per ADR-0023) are applied per-variant before
        the SQL batch. The two-key ``(rsid, genotype) AND is_nonfinding=0``
        match from ``annotate`` is preserved by fetching all
        ``is_nonfinding=0`` rows for the rsid and filtering the
        genotype per variant.
        """
        variants = list(variants)
        if self._reference_prefetch is not None:
            self._reference_prefetch(variants)
        # Pre-filter per-variant: only candidates whose genotype is normalizable
        # and that pass the ClinVar REF homref check.
        candidates: list[tuple[Variant, str]] = []
        for variant in variants:
            if variant.is_no_call:
                continue
            user_geno = _normalize_genotype(variant.allele1 + variant.allele2)
            if user_geno is None:
                continue
            if is_clinvar_homref(variant, self._clinvar_ref_provider):
                continue
            candidates.append((variant, user_geno))

        if not candidates:
            return
        rsids = list({v.rsid for v, _ in candidates})
        palindromic_set: set[str] = set()
        if self._palindromic_provider is not None:
            palindromic_set = self._palindromic_provider(set(rsids))
        readings = self._site_readings([v for v, _ in candidates])

        conn = self._connection()
        rows_by_rsid: dict[str, list[tuple]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            cursor = conn.execute(
                f"SELECT rsid, genotype, gene, drugs, phenotype, phenotype_category, "
                f"annotation_text, level_of_evidence, score, pgkb_annotation_id "
                f"FROM pharmgkb_annotations "
                f"WHERE rsid IN ({placeholders}) AND is_nonfinding = 0",
                chunk,
            )
            for row in cursor:
                rows_by_rsid.setdefault(row[0], []).append(row[1:])

        site_genotypes = self._site_genotypes(sorted(rows_by_rsid))
        for (variant, _user_geno), reading in zip(candidates, readings, strict=True):
            rows = rows_by_rsid.get(variant.rsid)
            if not rows:
                continue
            met = _genotypes_met(
                variant,
                site_genotypes.get(variant.rsid, []),
                site_palindromic=variant.rsid in palindromic_set,
                site_reading=reading,
            )
            if met is None:
                self._ambiguous_strand_calls += 1
                continue
            user_diploid = diploid_label(variant)
            for row in rows:
                (
                    geno,
                    gene,
                    drugs,
                    phenotype,
                    _phenotype_category,
                    annotation_text,
                    level_of_evidence,
                    _score,
                    pgkb_annotation_id,
                ) = row
                # ADR-0035: the genotypes this call meets by its one reading
                # (``_genotypes_met``; see annotate()).
                if geno.upper() not in met:
                    continue
                sig_label = level_of_evidence.strip().lower() or "unknown"
                description_parts = [f"ClinPGx: {drugs}"] if drugs else ["ClinPGx"]
                if phenotype:
                    description_parts.append(phenotype)
                if annotation_text:
                    description_parts.append(annotation_text)
                description = " — ".join(description_parts)
                references = (
                    [f"pharmgkb:annotation/{pgkb_annotation_id}"] if pgkb_annotation_id else []
                )
                alt, reference_carries_alt = carried_alt(
                    variant, variant.allele1, variant.allele2, self._normal_allele_provider
                )
                yield Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"pharmgkb_loe_{sig_label}",
                    category="pharma",
                    magnitude=_magnitude(level_of_evidence),
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=phenotype or "",
                    gene=gene or "",
                    alt=alt,
                    locus_ref=variant.ref,
                    reference_carries_alt=reference_carries_alt,
                    carried_alleles=(variant.allele1, variant.allele2),
                    locus=variant_locus(variant),
                )


def _stamp_existing_pharmgkb_cache(db_path: Path) -> bool:
    """One-shot migration: stamp ``local_version_tag`` on a ClinPGx cache.

    Handles legacy caches with ``|iv:N`` baked into ``remote_signal``
    by moving the tag and cleaning the signal, but only when ``N`` is the
    current interpreter version. A cache without the marker, or with
    another ``N``, has unknown or old data: it is left unstamped and the
    caller re-ingests or re-downloads it. Returns True if the current
    interpreter version is now stamped.
    """
    import contextlib

    from allelix.databases.manager import _ensure_local_version_tag_column, parse_baked_version

    if not db_path.exists():
        return False
    tag = f"iv:{PHARMGKB_INTERPRETER_VERSION}"
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            _ensure_local_version_tag_column(conn)
            row = conn.execute(
                "SELECT remote_signal, local_version_tag "
                "FROM database_versions WHERE name='pharmgkb'"
            ).fetchone()
            if not row:
                return False
            sig, existing_tag = row
            if existing_tag == tag:
                return True
            if existing_tag is not None:
                return False
            if parse_baked_version(sig, "iv") != PHARMGKB_INTERPRETER_VERSION:
                return False
            clean_signal = (sig or "").split("|iv:")[0]
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = ? "
                "WHERE name = 'pharmgkb'",
                (clean_signal, tag),
            )
            conn.commit()
        return True
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        return False


def _genotypes_met(
    variant: Variant,
    genotypes: list[str],
    *,
    site_palindromic: bool,
    site_reading: str | None,
) -> set[str] | None:
    """The rsID's genotypes the call meets, by one reading; None when two disagree.

    Every genotype of the rsID counts, non-finding rows included: they are
    evidence of the rows' strand even where they are never reported. Where
    the site decides the reading (``_site_reading``), that is the one
    reading: as written for ``direct``, the complement for ``complement``.
    With no such evidence, a call one genotype meets as written and another
    as its complement (by its REF, ``strand_aware_genotype_match``) cannot be
    read either way, and meets neither (ADR-0035).
    """
    exact = "".join(sorted((variant.allele1.upper(), variant.allele2.upper())))
    if site_reading == STRAND_DIRECT:
        return {g.upper() for g in genotypes if g.upper() == exact}
    site_reads_complement = site_reading == STRAND_COMPLEMENT
    met = {
        g.upper()
        for g in genotypes
        if strand_aware_genotype_match(
            variant.allele1,
            variant.allele2,
            g,
            variant.ref,
            site_palindromic=site_palindromic,
            site_reads_complement=site_reads_complement,
        )
    }
    if site_reads_complement:
        return met
    if exact in met and met - {exact}:
        return None
    return met


def _site_reading(
    variant: Variant, genotypes: list[str], sites: dict[str, ForwardSite]
) -> str | None:
    """How the rsID's source alleles read at the call: ``direct``, ``complement`` or None.

    All the rsID's genotypes together are read against the identity table's
    forward alleles at the call's locus (ADR-0035, ADR-0066's reading):
    ``direct`` where only the alleles as written fit, ``complement`` where
    only their complements do. None (undecided) at a site the table does not
    place at the call, where both readings or neither fit, and where any
    genotype is not two single bases.

    Only for a single-base call: its alleles single bases and its REF one
    base (or not known). A trimmed or padded indel's remaining letters are
    not a sequence allele the source genotype can name, even where they
    spell one of the site's identity alleles (another indel's ALT).
    """
    if not (is_single_base(variant.allele1) and is_single_base(variant.allele2)):
        return None
    if variant.ref is not None and not is_single_base(variant.ref.upper()):
        return None
    forward = forward_alleles_at_call(variant, sites)
    if forward is None or not genotypes:
        return None
    alleles: set[str] = set()
    for genotype in genotypes:
        bases = genotype.upper()
        if len(bases) != 2 or not all(is_single_base(b) for b in bases):
            return None
        alleles.update(bases)
    reading = forward_strand_reading(alleles, forward)
    return reading if reading in (STRAND_DIRECT, STRAND_COMPLEMENT) else None

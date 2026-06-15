# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GWAS Catalog annotator. Source-attributed trait associations."""

from __future__ import annotations

import contextlib
import logging
import sqlite3
import zipfile
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases.gwas_loader import (
    _CATEGORIZER_VERSION,
    _REQUIRED_GWAS_COLUMNS,
    GWAS_CATALOG_URL,
    GWAS_DB_FILENAME,
    GWAS_MIN_ROWS,
    load_gwas_tsv,
    schema_is_current,
)
from allelix.databases.manager import (
    _ensure_local_version_tag_column,
    download,
    get_database_info,
    head_request_headers,
)
from allelix.models import Annotation

_EXCLUDED_TRAIT_CATEGORIES = frozenset(
    {
        "body_measurement",
        "lipid_measurement",
        "hematological_measurement",
        "other_measurement",
        "behavioral",
    }
)

_MUST_INCLUDE_RSIDS = frozenset(
    {
        "rs10737680",  # CFH — age-related macular degeneration
        "rs11209026",  # IL23R — inflammatory bowel disease
        "rs9271366",  # HLA-DRB1 — multiple sclerosis
    }
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999


def _magnitude(p_value: float | None, or_beta: float | None) -> float:
    """Derive magnitude from p-value and optional effect size.

    GH #17: boundary comparisons are inclusive (``<=``) so the canonical
    genome-wide-significance threshold ``p = 5e-8`` lands inside the
    significant bucket rather than the suggestive bucket below it.
    Strict ``<`` made the exact threshold value fall a full magnitude
    below a barely-significant hit.
    """
    if p_value is None:
        base = 2.0
    elif p_value <= 5e-100:
        base = 8.0
    elif p_value <= 5e-20:
        base = 7.0
    elif p_value <= 5e-8:
        base = 6.0
    elif p_value <= 5e-6:
        base = 4.0
    elif p_value <= 5e-4:
        base = 3.0
    else:
        base = 2.0

    if or_beta is not None and or_beta > 0:
        if or_beta >= 3.0 or or_beta <= 0.33:
            base = min(base + 1.0, 9.0)
        elif or_beta >= 2.0 or or_beta <= 0.5:
            base = min(base + 0.5, 9.0)

    return base


_UNKNOWN_RISK_ALLELE_MAG_CAP = 3.0


class GWASCatalogAnnotator(Annotator):
    """Annotates variants with GWAS Catalog trait associations."""

    name: ClassVar[str] = "gwas"
    display_name: ClassVar[str] = "GWAS Catalog"
    attribution: ClassVar[str] = "GWAS Catalog"
    requires_download: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="custom-embl-ebi",
        license_url="https://www.ebi.ac.uk/gwas/docs/about",
        attribution_text=(
            "GWAS Catalog data sourced from NHGRI-EBI GWAS Catalog,"
            " available under EMBL-EBI Terms of Use."
        ),
        source_url="https://www.ebi.ac.uk/gwas/",
        commercial_ok=True,
    )

    def __init__(self, data_dir: Path, *, filter_traits: bool = True) -> None:
        """Initialize with path to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / GWAS_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._filter_traits = filter_traits

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download GWAS Catalog associations ZIP, extract TSV, and ingest."""
        url = GWAS_CATALOG_URL
        signal = self.fetch_remote_signal()
        if signal is None:
            msg = (
                "gwas: cannot verify remote freshness signal. "
                "Refresh aborted to avoid persisting an incomplete cache stamp. "
                "Retry, or pass --force if you accept that next `db update` "
                "will re-download to re-establish the signal."
            )
            raise RuntimeError(msg)
        zip_path = self.data_dir / "gwas_catalog_associations.zip"
        tsv_path = self.data_dir / "gwas_catalog_associations.tsv"
        # No content-hash verification: EBI publishes no checksum for this
        # file and the content is mutable, so there is nothing to pin or
        # fetch. TLS + Content-Length truncation guard only. See ADR-0029.
        download(url, zip_path)
        try:
            with zipfile.ZipFile(zip_path) as zf:
                tsv_names = [n for n in zf.namelist() if n.endswith(".tsv")]
                if not tsv_names:
                    msg = f"No .tsv file found in {zip_path}"
                    raise RuntimeError(msg)
                zf.extract(tsv_names[0], self.data_dir)
                extracted = self.data_dir / tsv_names[0]
                if extracted != tsv_path:
                    extracted.rename(tsv_path)
            load_gwas_tsv(
                tsv_path,
                self._db_path,
                source_url=url,
                remote_signal=signal,
                min_rows=GWAS_MIN_ROWS,
            )
        finally:
            try:
                zip_path.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                logger.warning("Could not remove staged file at %s", zip_path)

    def is_ready(self) -> bool:
        """Return True when the local GWAS cache exists and has current schema.

        Handles three states:

        1. Current tag in ``local_version_tag`` — ready.
        2. No tag (legacy cache or ``|cv:`` still in ``remote_signal``) —
           one-shot migration via ``_stamp_existing_gwas_cache``.
        3. Stale tag (categorizer bumped) — auto-reingest from cached TSV
           if still on disk.
        """
        info = get_database_info(self._db_path, "gwas")
        if info is None:
            return False
        tag = info.get("local_version_tag") or ""
        if tag == f"cv:{_CATEGORIZER_VERSION}":
            return _has_current_gwas_columns(self._db_path)
        if not tag and _stamp_existing_gwas_cache(self._db_path):
            return _has_current_gwas_columns(self._db_path)
        tsv_path = self.data_dir / "gwas_catalog_associations.tsv"
        if tsv_path.exists():
            print(
                "GWAS categorizer changed — re-ingesting from cached TSV...",
                flush=True,
            )
            try:
                load_gwas_tsv(
                    tsv_path,
                    self._db_path,
                    source_url=GWAS_CATALOG_URL,
                    remote_signal=self.cached_remote_signal(),
                    min_rows=GWAS_MIN_ROWS,
                )
            except Exception:
                logger.warning("Auto-reingest from cached TSV failed", exc_info=True)
                return False
            return schema_is_current(self._db_path)
        return False

    def version(self) -> str | None:
        """Return the cached database version string, or None."""
        info = get_database_info(self._db_path, "gwas")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of cached GWAS association records, or None."""
        info = get_database_info(self._db_path, "gwas")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Probe the GWAS Catalog URL for ETag or Last-Modified."""
        headers = head_request_headers(GWAS_CATALOG_URL)
        if headers is None:
            return None
        etag = headers.get("ETag") or headers.get("Etag")
        last_modified = headers.get("Last-Modified") or headers.get("Last-modified")
        if etag:
            return f"etag:{etag.strip()}"
        if last_modified:
            return f"lm:{last_modified.strip()}"
        return None

    def cached_remote_signal(self) -> str | None:
        """Return the remote signal stored during the last successful ingest."""
        info = get_database_info(self._db_path, "gwas")
        if not info or not info["remote_signal"]:
            return None
        return info["remote_signal"] or None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return GWAS Catalog annotations for variants the user carries.

        Carrier matching uses the risk allele when specified. When the
        GWAS entry doesn't specify a risk allele, the annotation fires
        on rsid match alone with a magnitude penalty.
        """
        if variant.is_no_call:
            return []

        sql = (
            "SELECT risk_allele, trait, p_value, or_beta, gene, "
            "study_accession, pubmed_id, trait_category "
            "FROM gwas_associations WHERE rsid = ?"
        )
        rows = self._connection().execute(sql, (variant.rsid,)).fetchall()
        annotations: list[Annotation] = []
        user_diploid = _user_diploid(variant)

        for row in rows:
            (
                risk_allele,
                trait,
                p_value,
                or_beta,
                gene,
                study_accession,
                pubmed_id,
                trait_category,
            ) = row

            if self._filter_traits and trait_category in _EXCLUDED_TRAIT_CATEGORIES:
                continue

            if risk_allele is not None:
                if variant.allele1 != risk_allele and variant.allele2 != risk_allele:
                    continue
                mag = _magnitude(p_value, or_beta)
                risk_note = ""
            else:
                # ADR-0024: unknown risk allele fires on rsID match alone
                # but capped at 3.0 so it doesn't pass typical --min-magnitude
                # thresholds. Without knowing which allele is the risk allele,
                # we can't apply the carrier rule (ADR-0007).
                mag = min(_magnitude(p_value, or_beta), _UNKNOWN_RISK_ALLELE_MAG_CAP)
                risk_note = " (risk allele not specified in study)"

            p_str = f"p={p_value:.1e}" if p_value is not None else "p=N/A"
            gene_str = gene or "—"
            description = f"GWAS Catalog: {trait} ({p_str}, gene: {gene_str}){risk_note}"

            references: list[str] = []
            if pubmed_id:
                references.append(f"pubmed:{pubmed_id}")
            if study_accession:
                references.append(f"gwas:{study_accession}")

            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance="gwas_association",
                    category="trait",
                    magnitude=mag,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=trait,
                    gene=gene or "",
                    alt="",
                    is_must_include=variant.rsid in _MUST_INCLUDE_RSIDS,
                )
            )
        return annotations

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Bulk-annotate via chunked ``WHERE rsid IN (...)`` queries.

        Replaces N per-variant queries with ceil(N/500) batched queries.
        At WGS scale (4-6M variants) this is the difference between
        annotation completing in minutes vs. hours.

        Per-variant filters (no-call rejection, risk-allele carrier
        check, trait-category exclusion, unknown-risk-allele magnitude
        cap) are applied after the SQL — same logic as ``annotate``,
        same outputs, same order.
        """
        variants_list = [v for v in variants if not v.is_no_call]
        if not variants_list:
            return
        rsids = list({v.rsid for v in variants_list})

        conn = self._connection()
        rows_by_rsid: dict[str, list[tuple]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            cursor = conn.execute(
                f"SELECT rsid, risk_allele, trait, p_value, or_beta, gene, "
                f"study_accession, pubmed_id, trait_category "
                f"FROM gwas_associations WHERE rsid IN ({placeholders})",
                chunk,
            )
            for row in cursor:
                rows_by_rsid.setdefault(row[0], []).append(row[1:])

        for variant in variants_list:
            rows = rows_by_rsid.get(variant.rsid)
            if not rows:
                continue
            user_diploid = _user_diploid(variant)
            for row in rows:
                (
                    risk_allele,
                    trait,
                    p_value,
                    or_beta,
                    gene,
                    study_accession,
                    pubmed_id,
                    trait_category,
                ) = row

                if self._filter_traits and trait_category in _EXCLUDED_TRAIT_CATEGORIES:
                    continue

                if risk_allele is not None:
                    if variant.allele1 != risk_allele and variant.allele2 != risk_allele:
                        continue
                    mag = _magnitude(p_value, or_beta)
                    risk_note = ""
                else:
                    mag = min(_magnitude(p_value, or_beta), _UNKNOWN_RISK_ALLELE_MAG_CAP)
                    risk_note = " (risk allele not specified in study)"

                p_str = f"p={p_value:.1e}" if p_value is not None else "p=N/A"
                gene_str = gene or "—"
                description = f"GWAS Catalog: {trait} ({p_str}, gene: {gene_str}){risk_note}"

                references: list[str] = []
                if pubmed_id:
                    references.append(f"pubmed:{pubmed_id}")
                if study_accession:
                    references.append(f"gwas:{study_accession}")

                yield Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance="gwas_association",
                    category="trait",
                    magnitude=mag,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=user_diploid,
                    references=references,
                    condition=trait,
                    gene=gene or "",
                    alt="",
                    is_must_include=variant.rsid in _MUST_INCLUDE_RSIDS,
                )


def _user_diploid(variant: Variant) -> str:
    """Sorted two-letter diploid for SNVs; indel passthrough verbatim."""
    a1, a2 = variant.allele1, variant.allele2
    if len(a1) == 1 and len(a2) == 1:
        return "".join(sorted((a1, a2)))
    return f"{a1}/{a2}"


def _has_current_gwas_columns(db_path: Path) -> bool:
    """True iff the gwas_associations table has the required columns."""
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            cols = {row[1] for row in conn.execute("PRAGMA table_info(gwas_associations)")}
            return _REQUIRED_GWAS_COLUMNS.issubset(cols)
    except sqlite3.DatabaseError:
        return False


def _stamp_existing_gwas_cache(db_path: Path) -> bool:
    """One-shot migration: stamp ``local_version_tag`` on a GWAS cache.

    Handles legacy caches with ``|cv:N`` baked into ``remote_signal``
    by moving the tag and cleaning the signal. Returns True if the
    current categorizer version is now stamped.
    """
    if not db_path.exists():
        return False
    tag = f"cv:{_CATEGORIZER_VERSION}"
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            _ensure_local_version_tag_column(conn)
            row = conn.execute(
                "SELECT remote_signal, local_version_tag FROM database_versions WHERE name='gwas'"
            ).fetchone()
            if not row:
                return False
            sig, existing_tag = row
            if existing_tag == tag:
                return True
            if existing_tag is not None:
                return False
            clean_signal = (sig or "").split("|cv:")[0]
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = ? "
                "WHERE name = 'gwas'",
                (clean_signal, tag),
            )
            conn.commit()
        return True
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        return False

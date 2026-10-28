# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""gnomAD population frequency enrichment.

gnomAD is not a clinical annotator — it does not produce Annotation
objects. It enriches existing annotations with population allele
frequency context. The pipeline calls ``bulk_lookup()`` after all
annotators have run, and stamps each annotation's ``allele_frequency``
field.

License: ODbL v1.0 (Open Database License). We extract only rsID +
allele frequencies (no SpliceAI or other restrictively licensed fields).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.cache_state import (
    CacheStatus,
    cache_state,
    connect_readonly,
)
from allelix.databases.gnomad_loader import (
    GNOMAD_CACHE_URL,
    GNOMAD_DB_FILENAME,
    GNOMAD_EXPECTED_SHA256,
    install_prebuilt_cache,
)
from allelix.databases.manager import (
    download,
    get_database_info,
    verify_file_hash,
)

if TYPE_CHECKING:
    import sqlite3
    from pathlib import Path

    from allelix.models import Annotation, Variant

logger = logging.getLogger(__name__)

_BULK_BATCH_SIZE = 900
# Positions bound per position-lookup statement, plus one chromosome:
# 999 parameters, the SQLite default limit every query here respects.
_POSITIONS_PER_STATEMENT = 998


class GnomadAnnotator(Annotator):
    """Population frequency enrichment from gnomAD.

    Subclasses Annotator for ``db update`` / ``db status`` / ``is_ready()``
    integration. ``annotate()`` always returns ``[]`` — gnomAD does not
    participate in the per-variant annotation loop.
    """

    name: ClassVar[str] = "gnomad"
    display_name: ClassVar[str] = "gnomAD"
    attribution: ClassVar[str] = "gnomAD"
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = False
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="ODbL-1.0",
        license_url="https://opendatacommons.org/licenses/odbl/1-0/",
        attribution_text=("Population frequencies sourced from gnomAD, used under ODbL v1.0."),
        source_url="https://gnomad.broadinstitute.org",
        commercial_ok=True,
    )

    def __init__(self, data_dir: Path) -> None:
        """Bind to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / GNOMAD_DB_FILENAME
        self._conn: sqlite3.Connection | None = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            if not self._db_path.exists():
                raise FileNotFoundError(
                    f"gnomAD cache not found at {self._db_path}. Run `allelix db update` first."
                )
            self._conn = connect_readonly(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built gnomAD exome frequency cache from HuggingFace."""
        gz_path = self.data_dir / "gnomad.sqlite.gz"
        download(GNOMAD_CACHE_URL, gz_path)
        verify_file_hash(gz_path, "sha256", GNOMAD_EXPECTED_SHA256)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=GNOMAD_CACHE_URL,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the gnomAD SQLite cache exists with current schema version.

        A cache with no ``local_version_tag`` is not ready. Accepting it
        would defeat ``GNOMAD_SCHEMA_VERSION``: after a bump, every
        untagged legacy cache would silently pass as the new version.
        Rejecting it tells the user to re-run ``db update``.
        """
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        """State of the gnomAD cache (ADR-0048)."""
        return cache_state(
            self._db_path,
            "gnomad",
            tag=f"sv:{GNOMAD_SCHEMA_VERSION}",
            tables=("gnomad_frequencies",),
        )

    def version(self) -> str | None:
        """Return the cached database version, or None."""
        info = get_database_info(self._db_path, "gnomad")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of rsIDs in the cache, or None."""
        info = get_database_info(self._db_path, "gnomad")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Code-driven source — no runtime freshness probe (ADR-0030)."""
        return None

    def cached_remote_signal(self) -> str | None:
        """Code-driven source — no cached signal to compare (ADR-0030)."""
        return None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Not used — gnomAD enriches, does not annotate. Always returns []."""
        return []

    def lookup(self, rsid: str) -> float | None:
        """Return global allele frequency for a single rsID, or None."""
        conn = self._connection()
        row = conn.execute(
            "SELECT MAX(af) FROM gnomad_frequencies WHERE rsid = ?", (rsid,)
        ).fetchone()
        return row[0] if row else None

    def bulk_lookup(self, rsids: set[str]) -> dict[str, float]:
        """Return ``{rsid: af}`` for all rsIDs found in the cache.

        Fallback for annotations without a known alt allele. Uses MAX to
        resolve multi-allelic sites. Prefer ``bulk_lookup_by_alt`` when alt
        is available.

        Batches into chunks of 900 to stay within SQLite's variable limit.
        """
        if not rsids:
            return {}
        conn = self._connection()
        result: dict[str, float] = {}
        rsid_list = list(rsids)
        for i in range(0, len(rsid_list), _BULK_BATCH_SIZE):
            batch = rsid_list[i : i + _BULK_BATCH_SIZE]
            placeholders = ",".join("?" * len(batch))
            rows = conn.execute(
                f"SELECT rsid, MAX(af) FROM gnomad_frequencies"
                f" WHERE rsid IN ({placeholders}) GROUP BY rsid",
                batch,
            ).fetchall()
            for rsid, af in rows:
                if af is not None:
                    result[rsid] = af
        return result

    def bulk_resolve_coordinates(
        self, rsids: set[str]
    ) -> dict[str, list[tuple[str, int, str, str]]]:
        """Return ``{rsid: [(chrom, pos, ref, alt), ...]}`` from the gnomAD cache.

        Maps rsIDs to genomic coordinates for coordinate-based lookups
        (CADD, future VCF-keyed sources). Multi-allelic sites return
        multiple tuples per rsid.
        """
        if not rsids:
            return {}
        conn = self._connection()
        result: dict[str, list[tuple[str, int, str, str]]] = {}
        rsid_list = list(rsids)
        for i in range(0, len(rsid_list), _BULK_BATCH_SIZE):
            batch = rsid_list[i : i + _BULK_BATCH_SIZE]
            placeholders = ",".join("?" * len(batch))
            rows = conn.execute(
                f"SELECT rsid, chrom, pos, ref, alt FROM gnomad_frequencies"
                f" WHERE rsid IN ({placeholders})",
                batch,
            ).fetchall()
            for rsid, chrom, pos, ref, alt in rows:
                result.setdefault(rsid, []).append((chrom, pos, ref, alt))
        return result

    def bulk_resolve_rsids_from_positions(
        self, positions: set[tuple[str, int]]
    ) -> dict[tuple[str, int], list[tuple[str, str, str]]]:
        """Reverse-lookup rsIDs by ``(chromosome, position)``.

        Counterpart to :meth:`bulk_resolve_coordinates` for the
        flow: a VCF variant arriving with an empty ID column (DeepVariant
        gVCFs, GIAB GRCh38 benchmarks, anything that doesn't stamp dbSNP
        IDs) is parsed as a positional pseudo-ID. The ClinVar resolver
        (:meth:`allelix.annotators.clinvar.ClinVarAnnotator.bulk_resolve_rsids`)
        recovers an rsID only for variants in ClinVar's curated subset —
        which excludes most pharmacogenomic, GWAS-only, and SNPedia-only
        rsIDs commonly found in wellness panels. This method fills the
        gap by querying ``gnomad_frequencies`` (keyed on the full dbSNP
        rsID universe with published allele frequencies).

        Returns ``{(chrom, pos): [(ref, alt, rsid), ...]}``. Multi-allelic
        positions return multiple rows; the caller is responsible for
        applying carrier-rule disambiguation (matching the user's allele
        pair as a subset of ``{ref, alt}`` and abstaining on ties — same
        pattern as the ClinVar resolver).

        Bare chromosomes are expected on both the input set and the
        ``gnomad_frequencies.chrom`` column — Allelix normalizes all
        internal chromosome identifiers to bare form
        (:func:`allelix.parsers._helpers.normalize_chromosome`).
        """
        if not positions:
            return {}
        conn = self._connection()
        result: dict[tuple[str, int], list[tuple[str, str, str]]] = {}
        # One statement per chromosome per chunk of up to
        # ``_POSITIONS_PER_STATEMENT`` positions: ``chrom = ? AND pos IN
        # (...)``. SQLite plans this as a search on both columns of the
        # primary key, ``(chrom, pos, ref, alt)``. The row-value form
        # ``(chrom, pos) IN ((?, ?), ...)`` searched on ``chrom`` alone on
        # SQLite 3.45.1 (Ubuntu 24.04), reading a whole chromosome per
        # statement. The input is grouped by chromosome first, so a batch
        # in file order costs one statement per ~1,000 positions.
        by_chrom: dict[str, list[int]] = {}
        for chrom, pos in positions:
            by_chrom.setdefault(chrom, []).append(pos)
        for chrom, chrom_positions in by_chrom.items():
            chrom_positions.sort()
            for i in range(0, len(chrom_positions), _POSITIONS_PER_STATEMENT):
                chunk = chrom_positions[i : i + _POSITIONS_PER_STATEMENT]
                placeholders = ",".join("?" * len(chunk))
                rows = conn.execute(
                    f"SELECT chrom, pos, ref, alt, rsid FROM gnomad_frequencies"
                    f" WHERE chrom = ? AND pos IN ({placeholders})",
                    (chrom, *chunk),
                ).fetchall()
                for c, p, ref, alt, rsid in rows:
                    if not rsid:
                        continue
                    result.setdefault((c, p), []).append((ref, alt, rsid))
        # Deterministic ordering for the carrier-rule pass.
        for rows in result.values():
            rows.sort()
        return result

    def bulk_lookup_by_alt(self, keys: set[tuple[str, str]]) -> dict[tuple[str, str], float]:
        """Return ``{(rsid, alt): af}`` for exact allele matches."""
        if not keys:
            return {}
        conn = self._connection()
        result: dict[tuple[str, str], float] = {}
        key_list = list(keys)
        batch_size = _BULK_BATCH_SIZE // 2
        for i in range(0, len(key_list), batch_size):
            batch = key_list[i : i + batch_size]
            clauses = " OR ".join(["(rsid = ? AND alt = ?)"] * len(batch))
            params = [v for rsid, alt in batch for v in (rsid, alt)]
            rows = conn.execute(
                f"SELECT rsid, alt, af FROM gnomad_frequencies WHERE {clauses}",
                params,
            ).fetchall()
            for rsid, alt, af in rows:
                if af is not None:
                    result[(rsid, alt)] = af
        return result

    def sites_frequencies(
        self, sites: set[tuple[str, int, str]]
    ) -> dict[tuple[str, int, str], dict[str, float | None]]:
        """Every ALT gnomAD lists at each ``(chrom, pos, ref)``, with its frequency.

        An ALT whose frequency gnomAD does not give is kept, as None; a site
        gnomAD does not list is absent. Each site is a prefix of the
        ``(chrom, pos, ref, alt)`` primary key.
        """
        result: dict[tuple[str, int, str], dict[str, float | None]] = {}
        if not sites:
            return result
        conn = self._connection()
        site_list = sorted(sites)
        batch_size = _BULK_BATCH_SIZE // 3
        for i in range(0, len(site_list), batch_size):
            batch = site_list[i : i + batch_size]
            clauses = " OR ".join(["(chrom = ? AND pos = ? AND ref = ?)"] * len(batch))
            params = [v for site in batch for v in site]
            for chrom, pos, ref, alt, af in conn.execute(
                f"SELECT chrom, pos, ref, alt, af FROM gnomad_frequencies WHERE {clauses}",
                params,
            ):
                result.setdefault((chrom, pos, ref), {})[alt] = af
        return result

    def bulk_lookup_by_position(
        self, keys: set[tuple[str, int, str, str]]
    ) -> dict[tuple[str, int, str, str], float]:
        """Return ``{(chrom, pos, ref, alt): af}`` via primary-key lookup.

        Position-keyed fallback for rsID-less VCFs (ClinVar resolved their
        rsIDs from coordinates; the gnomAD rsID index may not include those
        rsIDs even when the position is in gnomAD). Hits the
        ``(chrom, pos, ref, alt)`` primary key directly, so each row is
        an O(log n) lookup.
        """
        if not keys:
            return {}
        conn = self._connection()
        result: dict[tuple[str, int, str, str], float] = {}
        key_list = list(keys)
        batch_size = _BULK_BATCH_SIZE // 4
        for i in range(0, len(key_list), batch_size):
            batch = key_list[i : i + batch_size]
            clauses = " OR ".join(["(chrom = ? AND pos = ? AND ref = ? AND alt = ?)"] * len(batch))
            params = [v for k in batch for v in k]
            rows = conn.execute(
                f"SELECT chrom, pos, ref, alt, af FROM gnomad_frequencies WHERE {clauses}",
                params,
            ).fetchall()
            for chrom, pos, ref, alt, af in rows:
                if af is not None:
                    result[(chrom, pos, ref, alt)] = af
        return result

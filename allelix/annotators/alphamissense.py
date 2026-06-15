# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""AlphaMissense variant pathogenicity enrichment.

AlphaMissense is not a clinical annotator — it does not produce
Annotation objects. It enriches existing annotations with missense
variant pathogenicity predictions. The pipeline calls
``bulk_lookup()`` after all annotators have run, and stamps each
annotation's ``am_pathogenicity`` and ``am_class`` fields.

License: CC BY 4.0. Attribution: Cheng et al., Science 2023
(doi:10.1126/science.adg7492).
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases._versions import ALPHAMISSENSE_SCHEMA_VERSION
from allelix.databases.alphamissense_loader import (
    ALPHAMISSENSE_CACHE_URL,
    ALPHAMISSENSE_DB_FILENAME,
    ALPHAMISSENSE_EXPECTED_SHA256,
    install_prebuilt_cache,
)
from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
from allelix.databases.manager import (
    download,
    get_database_info,
    verify_file_hash,
)

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation, Variant

logger = logging.getLogger(__name__)

_BULK_BATCH_SIZE = 900


class AlphaMissenseAnnotator(Annotator):
    """Missense variant pathogenicity enrichment from AlphaMissense.

    Subclasses Annotator for ``db update`` / ``db status`` / ``is_ready()``
    integration. ``annotate()`` always returns ``[]`` — AlphaMissense does
    not participate in the per-variant annotation loop.
    """

    name: ClassVar[str] = "alphamissense"
    display_name: ClassVar[str] = "AlphaMissense"
    attribution: ClassVar[str] = "AlphaMissense"
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = False
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="CC-BY-4.0",
        license_url="https://creativecommons.org/licenses/by/4.0/",
        attribution_text=(
            "AlphaMissense predictions from Cheng et al., Science 2023"
            " (doi:10.1126/science.adg7492). Licensed under CC BY 4.0."
        ),
        source_url="https://zenodo.org/records/10813168",
        citation="Cheng et al., Science 2023 (doi:10.1126/science.adg7492)",
        commercial_ok=True,
    )

    def __init__(self, data_dir: Path) -> None:
        """Bind to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / ALPHAMISSENSE_DB_FILENAME
        self._conn: sqlite3.Connection | None = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            if not self._db_path.exists():
                raise FileNotFoundError(
                    f"AlphaMissense cache not found at {self._db_path}. "
                    "Run `allelix db update` first."
                )
            self._conn = sqlite3.connect(self._db_path)
            self._check_gnomad_version()
        return self._conn

    def _check_gnomad_version(self) -> None:
        """Warn if the gnomAD version used to build the AM cache differs from installed."""
        assert self._conn is not None
        row = self._conn.execute(
            "SELECT version FROM database_versions WHERE name = 'alphamissense_gnomad_source'"
        ).fetchone()
        if row is None:
            return
        stamped = row[0]
        if stamped == "no_gnomad":
            logger.warning(
                "AlphaMissense cache was built without gnomAD (--no-gnomad). "
                "rsID lookups will return no results."
            )
            return
        gnomad_info = get_database_info(self.data_dir / GNOMAD_DB_FILENAME, "gnomad")
        if gnomad_info is None:
            return
        installed = gnomad_info["version"]
        if installed and stamped != installed:
            logger.warning(
                "AlphaMissense cache was built against gnomAD %s but installed "
                "gnomAD is %s. rsID mappings may be stale. Rebuild with: "
                "python scripts/build_alphamissense_cache.py",
                stamped,
                installed,
            )

    def setup(self) -> None:
        """Download the pre-built AlphaMissense cache from HuggingFace."""
        gz_path = self.data_dir / "alphamissense.sqlite.gz"
        download(ALPHAMISSENSE_CACHE_URL, gz_path)
        verify_file_hash(gz_path, "sha256", ALPHAMISSENSE_EXPECTED_SHA256)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=ALPHAMISSENSE_CACHE_URL,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the AlphaMissense SQLite cache exists with current schema version."""
        info = get_database_info(self._db_path, "alphamissense")
        if info is None:
            return False
        tag = info.get("local_version_tag") or ""
        return tag == f"sv:{ALPHAMISSENSE_SCHEMA_VERSION}" or not tag

    def version(self) -> str | None:
        """Return the cached database version, or None."""
        info = get_database_info(self._db_path, "alphamissense")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of variants in the cache, or None."""
        info = get_database_info(self._db_path, "alphamissense")
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
        """Not used — AlphaMissense enriches, does not annotate. Always returns []."""
        return []

    def lookup(self, rsid: str) -> tuple[float, str] | None:
        """Return (am_pathogenicity, am_class) for a single rsID, or None."""
        conn = self._connection()
        row = conn.execute(
            "SELECT MAX(am_pathogenicity), am_class FROM alphamissense_scores WHERE rsid = ?",
            (rsid,),
        ).fetchone()
        if row is None or row[0] is None:
            return None
        return (row[0], row[1])

    def bulk_lookup(self, rsids: set[str]) -> dict[str, tuple[float, str]]:
        """Return ``{rsid: (am_pathogenicity, am_class)}`` for found rsIDs.

        Fallback for annotations without a known alt allele. Uses MAX to
        resolve multi-allelic sites. Prefer ``bulk_lookup_by_alt`` when alt
        is available.

        Batches into chunks of 900 to stay within SQLite's variable limit.
        """
        if not rsids:
            return {}
        conn = self._connection()
        result: dict[str, tuple[float, str]] = {}
        rsid_list = list(rsids)
        for i in range(0, len(rsid_list), _BULK_BATCH_SIZE):
            batch = rsid_list[i : i + _BULK_BATCH_SIZE]
            placeholders = ",".join("?" * len(batch))
            rows = conn.execute(
                f"SELECT rsid, MAX(am_pathogenicity), am_class"
                f" FROM alphamissense_scores"
                f" WHERE rsid IN ({placeholders}) GROUP BY rsid",
                batch,
            ).fetchall()
            for rsid, score, cls in rows:
                if score is not None:
                    result[rsid] = (score, cls)
        return result

    def bulk_lookup_by_alt(
        self, keys: set[tuple[str, str]]
    ) -> dict[tuple[str, str], tuple[float, str]]:
        """Return ``{(rsid, alt): (am_pathogenicity, am_class)}`` for exact matches."""
        if not keys:
            return {}
        conn = self._connection()
        result: dict[tuple[str, str], tuple[float, str]] = {}
        key_list = list(keys)
        batch_size = _BULK_BATCH_SIZE // 2
        for i in range(0, len(key_list), batch_size):
            batch = key_list[i : i + batch_size]
            clauses = " OR ".join(["(rsid = ? AND alt = ?)"] * len(batch))
            params = [v for rsid, alt in batch for v in (rsid, alt)]
            rows = conn.execute(
                f"SELECT rsid, alt, am_pathogenicity, am_class"
                f" FROM alphamissense_scores WHERE {clauses}",
                params,
            ).fetchall()
            for rsid, alt, score, cls in rows:
                if score is not None:
                    result[(rsid, alt)] = (score, cls)
        return result

    def bulk_lookup_by_position(
        self, keys: set[tuple[str, int, str, str]]
    ) -> dict[tuple[str, int, str, str], tuple[float, str]]:
        """Return ``{(chrom, pos, ref, alt): (score, class)}`` via PK lookup.

        Position-keyed fallback for rsID-less VCFs whose ClinVar-resolved
        rsIDs aren't indexed in the AlphaMissense cache. Hits the
        ``(chrom, pos, ref, alt)`` primary key directly.
        """
        if not keys:
            return {}
        conn = self._connection()
        result: dict[tuple[str, int, str, str], tuple[float, str]] = {}
        key_list = list(keys)
        batch_size = _BULK_BATCH_SIZE // 4
        for i in range(0, len(key_list), batch_size):
            batch = key_list[i : i + batch_size]
            clauses = " OR ".join(["(chrom = ? AND pos = ? AND ref = ? AND alt = ?)"] * len(batch))
            params = [v for k in batch for v in k]
            rows = conn.execute(
                f"SELECT chrom, pos, ref, alt, am_pathogenicity, am_class"
                f" FROM alphamissense_scores WHERE {clauses}",
                params,
            ).fetchall()
            for chrom, pos, ref, alt, score, cls in rows:
                if score is not None:
                    result[(chrom, pos, ref, alt)] = (score, cls)
        return result

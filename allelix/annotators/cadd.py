# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""CADD variant deleteriousness enrichment.

CADD is not a clinical annotator — it does not produce Annotation
objects. It enriches existing annotations with PHRED-scaled
deleteriousness scores. The pipeline calls ``bulk_lookup()`` after all
annotators have run, and stamps each annotation's ``cadd_phred`` field.

Two modes:

* **Cache mode** (default): pre-built SQLite database from HuggingFace
  containing exome-region CADD scores. Fast, compact (~1 GB).
* **Full mode** (``options.cadd_full = true``): queries the complete
  CADD v1.7 tabix file (``whole_genome_SNVs.tsv.gz``, ~81 GB). Covers
  every scored position in the genome. Requires ``pysam`` and a local
  copy of the tabix file + index. **GRCh38 only.**

License: LicenseRef-CADD — free for non-commercial use only. Commercial
use requires a separate license from University of Washington (CoMotion).
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, Any, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases._versions import CADD_SCHEMA_VERSION
from allelix.databases.cadd_loader import (
    CADD_CACHE_URL,
    CADD_DB_FILENAME,
    CADD_EXPECTED_SHA256,
    install_prebuilt_cache,
)
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

CADD_FULL_FILENAME = "whole_genome_SNVs.tsv.gz"
CADD_INDEL_FILENAME = "gnomad.genomes.r4.0.indel.tsv.gz"


class CaddAnnotator(Annotator):
    """PHRED-scaled deleteriousness enrichment from CADD.

    Subclasses Annotator for ``db update`` / ``db status`` / ``is_ready()``
    integration. ``annotate()`` always returns ``[]`` — CADD does not
    participate in the per-variant annotation loop.
    """

    name: ClassVar[str] = "cadd"
    display_name: ClassVar[str] = "CADD"
    attribution: ClassVar[str] = "CADD"
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = False
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="LicenseRef-CADD",
        license_url="https://cadd.gs.washington.edu/license",
        attribution_text="CADD scores provided by the University of Washington.",
        source_url="https://cadd.gs.washington.edu/",
        citation="Schubach et al., Nucleic Acids Research 2024",
        commercial_ok=False,
        licensable=True,
        purchase_url="https://els2.comotion.uw.edu/product/cadd-scores",
    )

    def __init__(self, data_dir: Path, *, full_mode: bool = False) -> None:
        """Bind to the data directory.

        When ``full_mode`` is True the annotator queries a local tabix
        file instead of the pre-built SQLite cache. The tabix file must
        be placed at ``<data_dir>/whole_genome_SNVs.tsv.gz`` with its
        ``.tbi`` index alongside it.
        """
        super().__init__(data_dir)
        self._db_path = data_dir / CADD_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._full_mode = full_mode
        self._tabix_path = data_dir / CADD_FULL_FILENAME
        # pysam ships no type stubs (see [tool.mypy.overrides]); the tabix
        # handle is held as Any so call sites can use `.fetch()` / `.close()`.
        self._tabix: Any = None
        self._indel_tabix_path = data_dir / CADD_INDEL_FILENAME
        self._indel_tabix: Any = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            if not self._db_path.exists():
                raise FileNotFoundError(
                    f"CADD cache not found at {self._db_path}. "
                    "Run `allelix db update --cadd` first."
                )
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def _open_tabix(self) -> Any:  # noqa: ANN401 — pysam ships no type stubs
        """Open the tabix file for full-mode queries."""
        if self._tabix is None:
            try:
                import pysam
            except ImportError:
                raise ImportError(
                    "Full CADD mode requires pysam. Install with: pip install 'allelix[cadd]'"
                ) from None
            if not self._tabix_path.exists():
                raise FileNotFoundError(
                    f"CADD tabix file not found at {self._tabix_path}. "
                    "Download whole_genome_SNVs.tsv.gz and its .tbi index from "
                    "https://cadd.gs.washington.edu/download"
                )
            self._tabix = pysam.TabixFile(str(self._tabix_path))
        return self._tabix

    def _open_indel_tabix(self) -> Any:  # noqa: ANN401 — pysam ships no type stubs
        """Open the indel tabix file. Returns None if file doesn't exist."""
        if self._indel_tabix is None:
            if not self._indel_tabix_path.exists():
                return None
            tbi = self._indel_tabix_path.parent / (CADD_INDEL_FILENAME + ".tbi")
            if not tbi.exists():
                return None
            try:
                import pysam
            except ImportError:
                return None
            self._indel_tabix = pysam.TabixFile(str(self._indel_tabix_path))
        return self._indel_tabix

    def setup(self) -> None:
        """Download the pre-built CADD cache from HuggingFace."""
        gz_path = self.data_dir / "cadd.sqlite.gz"
        download(CADD_CACHE_URL, gz_path)
        verify_file_hash(gz_path, "sha256", CADD_EXPECTED_SHA256)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=CADD_CACHE_URL,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the active backend is available.

        In full mode, checks for the tabix file. In cache mode, checks
        the SQLite database.
        """
        if self._full_mode:
            return (
                self._tabix_path.exists()
                and (self._tabix_path.parent / (CADD_FULL_FILENAME + ".tbi")).exists()
            )
        info = get_database_info(self._db_path, "cadd")
        if info is None:
            return False
        tag = info.get("local_version_tag") or ""
        return tag == f"sv:{CADD_SCHEMA_VERSION}" or not tag

    def version(self) -> str | None:
        """Return the cached database version, or None."""
        if self._full_mode:
            return "v1.7 (full)" if self.is_ready() else None
        info = get_database_info(self._db_path, "cadd")
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of variants in the cache, or None."""
        if self._full_mode:
            return None
        info = get_database_info(self._db_path, "cadd")
        return info["record_count"] if info else None

    def close(self) -> None:
        """Close the SQLite connection or tabix file if open."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None
        if self._tabix is not None:
            self._tabix.close()
            self._tabix = None
        if self._indel_tabix is not None:
            self._indel_tabix.close()
            self._indel_tabix = None

    def fetch_remote_signal(self) -> str | None:
        """Code-driven source — no runtime freshness probe (ADR-0030)."""
        return None

    def cached_remote_signal(self) -> str | None:
        """Code-driven source — no cached signal to compare (ADR-0030)."""
        return None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Not used — CADD enriches, does not annotate. Always returns []."""
        return []

    def _tabix_lookup(self, chrom: str, pos: int, ref: str, alt: str) -> float | None:
        """Query the tabix file for a single variant.

        SNVs (single-base ref and alt) query the SNV tabix file.
        Indels route to the indel tabix file if available.
        """
        query_chrom = chrom if not chrom.startswith("chr") else chrom[3:]
        is_snv = len(ref) == 1 and len(alt) == 1

        tbx = self._open_tabix() if is_snv else self._open_indel_tabix()

        if tbx is None:
            return None

        try:
            for row in tbx.fetch(query_chrom, pos - 1, pos):
                fields = row.split("\t")
                if len(fields) >= 6 and fields[2] == ref and fields[3] == alt:
                    return float(fields[5])
        except (ValueError, KeyError):
            pass
        return None

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> float | None:
        """Return CADD PHRED score for a single variant, or None."""
        if self._full_mode:
            return self._tabix_lookup(chrom, pos, ref, alt)
        conn = self._connection()
        row = conn.execute(
            "SELECT phred FROM cadd_scores WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?",
            (chrom, pos, ref, alt),
        ).fetchone()
        return row[0] if row else None

    def bulk_lookup(
        self, keys: set[tuple[str, int, str, str]]
    ) -> dict[tuple[str, int, str, str], float]:
        """Return ``{(chrom, pos, ref, alt): phred}`` for all matches.

        In cache mode, batches SQL queries. In full mode, iterates tabix
        lookups (I/O bound on the tabix index, not CPU).
        """
        if not keys:
            return {}
        if self._full_mode:
            return self._tabix_bulk_lookup(keys)
        conn = self._connection()
        result: dict[tuple[str, int, str, str], float] = {}
        key_list = list(keys)
        batch_size = _BULK_BATCH_SIZE // 4
        for i in range(0, len(key_list), batch_size):
            batch = key_list[i : i + batch_size]
            clauses = " OR ".join(["(chrom = ? AND pos = ? AND ref = ? AND alt = ?)"] * len(batch))
            params: list[str | int] = []
            for chrom, pos, ref, alt in batch:
                params.extend([chrom, pos, ref, alt])
            rows = conn.execute(
                f"SELECT chrom, pos, ref, alt, phred FROM cadd_scores WHERE {clauses}",
                params,
            ).fetchall()
            for chrom, pos, ref, alt, phred in rows:
                result[(chrom, pos, ref, alt)] = phred
        return result

    def _tabix_bulk_lookup(
        self, keys: set[tuple[str, int, str, str]]
    ) -> dict[tuple[str, int, str, str], float]:
        """Batch tabix lookups for full mode."""
        result: dict[tuple[str, int, str, str], float] = {}
        for chrom, pos, ref, alt in keys:
            score = self._tabix_lookup(chrom, pos, ref, alt)
            if score is not None:
                result[(chrom, pos, ref, alt)] = score
        return result

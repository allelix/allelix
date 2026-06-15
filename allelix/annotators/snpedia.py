# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""SNPedia annotator. Structured SQL lookups against pre-parsed genotype data.

Reads from the ``snpedia_genotypes`` table in the SNPedia SQLite archive.
The pre-built cache is downloaded from HuggingFace during ``db update``.
It can also be built locally via ``scripts/scrape_snpedia.py`` followed
by ``scripts/parse_snpedia.py``.

SNPedia content is CC-BY-NC-SA 3.0 US. Attribution is required in all
reports.
"""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor, is_clinvar_homref
from allelix.databases.manager import (
    download,
    verify_file_hash,
)
from allelix.databases.snpedia_loader import (
    SNPEDIA_CACHE_URL,
    SNPEDIA_EXPECTED_SHA256,
    install_prebuilt_cache,
)
from allelix.models import Annotation

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

    from allelix.models import Variant

logger = logging.getLogger(__name__)

_BATCH_CHUNK = 500  # SQLite default SQLITE_MAX_VARIABLE_NUMBER is 999

SNPEDIA_DB_FILENAME = "snpedia.sqlite"
SNPEDIA_RECORD_NAME = "snpedia"

_REPUTE_CATEGORY: dict[str, str] = {
    "good": "trait",
    "bad": "clinical",
    "not set": "trait",
    "": "trait",
}

_SUMMARY_SUPPRESS_SUBSTRINGS: tuple[str, ...] = (
    "mis-oriented",
    "mis-orientation",
    "wrong strand",
    "orientation uncertain",
)


class SNPediaAnnotator(Annotator):
    """Annotates variants with SNPedia genotype data via structured SQL lookups."""

    name: ClassVar[str] = "snpedia"
    display_name: ClassVar[str] = "SNPedia"
    attribution: ClassVar[str] = "SNPedia"
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = False
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="CC-BY-NC-SA-3.0-US",
        license_url="https://creativecommons.org/licenses/by-nc-sa/3.0/us/",
        attribution_text=(
            "SNPedia annotations sourced from SNPedia, used under CC BY-NC-SA 3.0 US."
        ),
        source_url="https://www.snpedia.com",
        commercial_ok=False,
    )

    def __init__(
        self,
        data_dir: Path,
        clinvar_ref_provider: Callable[[str, str], str | None] | None = None,
    ) -> None:
        """Initialize with path to the data directory.

        ``clinvar_ref_provider`` is a ``(rsid, build) -> ref_base | None``
        callable used by the ADR-0023 hom-ref check. In production it is
        wired to ``ClinVarAnnotator.reference_for``. ``None`` disables the
        check (tests, standalone use).
        """
        super().__init__(data_dir)
        self._db_path = data_dir / SNPEDIA_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._clinvar_ref_provider = clinvar_ref_provider

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built SNPedia cache from HuggingFace.

        The HuggingFace asset contains raw wiki markup (third-party,
        CC-BY-NC-SA). After download, ``is_ready()`` triggers the
        client-side parse into structured genotype rows and stamps
        ``database_versions`` with proper version metadata.
        """
        gz_path = self.data_dir / "snpedia.sqlite.gz"
        download(SNPEDIA_CACHE_URL, gz_path)
        verify_file_hash(gz_path, "sha256", SNPEDIA_EXPECTED_SHA256)
        install_prebuilt_cache(
            gz_path,
            self._db_path,
            source_url=SNPEDIA_CACHE_URL,
        )
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)
        self.is_ready()

    def is_ready(self) -> bool:
        """Return True when the parsed SNPedia genotype table exists and has data.

        If raw pages exist but the structured table does not, automatically
        parses the raw markup (one-time operation, ~2 minutes).
        """
        if not self._db_path.exists():
            return False
        try:
            from allelix.databases.snpedia_parser import (
                detect_raw_table,
                parse_raw_pages,
                parser_is_current,
            )

            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                has_rows = False
                with contextlib.suppress(sqlite3.OperationalError):
                    has_rows = (
                        conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0] > 0
                    )

                needs_reparse = has_rows and not parser_is_current(conn)
                if has_rows and not needs_reparse:
                    return True

                raw_table = detect_raw_table(conn)
                if raw_table is None:
                    return False

                snp_count = conn.execute(
                    f"SELECT COUNT(*) FROM {raw_table} WHERE category='snp'"
                ).fetchone()[0]
                genotype_count = conn.execute(
                    f"SELECT COUNT(*) FROM {raw_table} WHERE category='genotype'"
                ).fetchone()[0]

            reason = "parser version changed" if needs_reparse else "one-time"
            print(
                f"Parsing {snp_count} SNP pages + {genotype_count} genotype pages"
                f" into structured table ({reason}, ~5 min)...",
                flush=True,
            )
            parsed = parse_raw_pages(str(self._db_path))
            print(f"Parsed {parsed} SNPedia genotype rows.", flush=True)
            return parsed > 0
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return False

    def version(self) -> str | None:
        """Return a version string from the database_versions table."""
        if not self._db_path.exists():
            return None
        try:
            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                row = conn.execute(
                    "SELECT version FROM database_versions WHERE name = ?",
                    (SNPEDIA_RECORD_NAME,),
                ).fetchone()
                if row and row[0]:
                    return row[0]
                return None
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return None

    def record_count(self) -> int | None:
        """Return the number of genotype rows in the structured table."""
        if not self._db_path.exists():
            return None
        try:
            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                count = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
                return count
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            return None

    def close(self) -> None:
        """Close the SQLite connection."""
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
        """Return SNPedia annotations matching the user's genotype."""
        if variant.is_no_call:
            return []

        snp_id = variant.rsid.lower()
        if snp_id.startswith("rs"):
            snp_num = snp_id[2:]
            snp_url_path = f"Rs{snp_num}"
        elif snp_id.startswith("i"):
            snp_num = snp_id[1:]
            snp_url_path = f"I{snp_num}"
        else:
            return []

        if not snp_num or not snp_num.isdigit():
            return []

        if snp_id.startswith("rs") and is_clinvar_homref(variant, self._clinvar_ref_provider):
            return []

        a1, a2 = variant.allele1.upper(), variant.allele2.upper()
        sorted_alleles = (a1, a2) if a1 <= a2 else (a2, a1)

        conn = self._connection()
        rows = conn.execute(
            "SELECT allele1, allele2, magnitude, repute, summary, gene "
            "FROM snpedia_genotypes "
            "WHERE rsid = ? AND allele1 = ? AND allele2 = ?",
            (snp_id, sorted_alleles[0], sorted_alleles[1]),
        ).fetchall()

        annotations: list[Annotation] = []
        for allele1, allele2, magnitude, repute, summary, gene in rows:
            if not summary:
                continue

            summary_lower = summary.lower()
            if any(p in summary_lower for p in _SUMMARY_SUPPRESS_SUBSTRINGS):
                continue

            if magnitude is None:
                magnitude = 0.0

            repute_lower = (repute or "").strip().lower()
            category = _REPUTE_CATEGORY.get(repute_lower, "trait")

            description = f"SNPedia: {summary}"
            genotype_match = f"{allele1}{allele2}"

            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"snpedia_{repute_lower}" if repute_lower else "snpedia_genotype",
                    category=category,
                    magnitude=magnitude,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=genotype_match,
                    references=[f"https://www.snpedia.com/index.php/{snp_url_path}"],
                    gene=gene or "",
                )
            )

        return annotations

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Bulk-annotate via chunked ``WHERE rsid IN (...)`` queries.

        Pre-filters (no-call, non-rs/i prefix, ClinVar-homref for rs
        variants) are applied per-variant before the SQL batch. The
        ``(rsid, allele1, allele2)`` three-key match from ``annotate``
        is preserved by fetching all rows for the rsid and filtering
        the allele tuple per variant — SNPedia typically has <5 rows
        per rsid, so the over-fetch is negligible vs. the saved
        per-variant round-trips.
        """
        # Pre-filter per-variant — only candidates that can possibly match
        candidates: list[tuple[Variant, str, str, tuple[str, str]]] = []
        for variant in variants:
            if variant.is_no_call:
                continue
            snp_id = variant.rsid.lower()
            if snp_id.startswith("rs"):
                snp_num = snp_id[2:]
                snp_url_path = f"Rs{snp_num}"
            elif snp_id.startswith("i"):
                snp_num = snp_id[1:]
                snp_url_path = f"I{snp_num}"
            else:
                continue
            if not snp_num or not snp_num.isdigit():
                continue
            if snp_id.startswith("rs") and is_clinvar_homref(variant, self._clinvar_ref_provider):
                continue
            a1, a2 = variant.allele1.upper(), variant.allele2.upper()
            sorted_alleles = (a1, a2) if a1 <= a2 else (a2, a1)
            candidates.append((variant, snp_id, snp_url_path, sorted_alleles))

        if not candidates:
            return
        rsids = list({snp_id for _, snp_id, _, _ in candidates})

        conn = self._connection()
        rows_by_rsid: dict[str, list[tuple]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            cursor = conn.execute(
                f"SELECT rsid, allele1, allele2, magnitude, repute, summary, gene "
                f"FROM snpedia_genotypes WHERE rsid IN ({placeholders})",
                chunk,
            )
            for row in cursor:
                rows_by_rsid.setdefault(row[0], []).append(row[1:])

        for variant, snp_id, snp_url_path, sorted_alleles in candidates:
            rows = rows_by_rsid.get(snp_id)
            if not rows:
                continue
            for allele1, allele2, magnitude, repute, summary, gene in rows:
                if (allele1, allele2) != sorted_alleles:
                    continue
                if not summary:
                    continue
                summary_lower = summary.lower()
                if any(p in summary_lower for p in _SUMMARY_SUPPRESS_SUBSTRINGS):
                    continue
                if magnitude is None:
                    magnitude = 0.0
                repute_lower = (repute or "").strip().lower()
                category = _REPUTE_CATEGORY.get(repute_lower, "trait")
                description = f"SNPedia: {summary}"
                genotype_match = f"{allele1}{allele2}"
                yield Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"snpedia_{repute_lower}" if repute_lower else "snpedia_genotype",
                    category=category,
                    magnitude=magnitude,
                    description=description,
                    attribution=self.attribution,
                    genotype_match=genotype_match,
                    references=[f"https://www.snpedia.com/index.php/{snp_url_path}"],
                    gene=gene or "",
                )

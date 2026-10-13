# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Cross-build rsID identity table (ADR-0040).

Not a clinical source: it never produces Annotations. The pipeline uses it
for three things, all keyed by locus on the input's own build:

1. ``Variant.ref`` for array data on GRCh37 and GRCh38, from the rows of
   the call's own rsID at its locus (``alleles_at_loci``).
2. rsID recovery for rsID-less VCF records on both builds
   (``alleles_at_loci`` feeding ``resolve_rsid_at_locus``).
3. Translating a GRCh37 ``(chrom, pos, ref, alt)`` to its GRCh38 key so the
   GRCh38-only enrichment sources (gnomAD, AlphaMissense, CADD) can be
   queried for GRCh37 input (``translate_to_grch38``). Palindromic sites
   and alleles that are the GRCh38 reference abstain.

dbSNP is a US government work in the public domain.
"""

from __future__ import annotations

import logging
import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases._versions import RSID_IDENTITY_SCHEMA_VERSION
from allelix.databases.manager import download, get_database_info, verify_file_hash
from allelix.databases.rsid_identity_loader import (
    RSID_IDENTITY_CACHE_URL,
    RSID_IDENTITY_DB_FILENAME,
    RSID_IDENTITY_EXPECTED_SHA256,
    RSID_IDENTITY_RECORD_NAME,
    install_prebuilt_cache,
)

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation, Variant

logger = logging.getLogger(__name__)

_CHUNK = 450  # two parameters per locus; stays under SQLite's 999 limit

# Build -> (position, ref, alt) column names. Fixed mapping, never user input.
_COLUMNS: dict[str, tuple[str, str, str]] = {
    "GRCh37": ("pos37", "ref37", "alt37"),
    "GRCh38": ("pos38", "ref38", "alt38"),
}

Locus = tuple[str, int]
Key = tuple[str, int, str, str]


class RsidIdentityAnnotator(Annotator):
    """Per-build placements and allele pairing for annotator rsIDs."""

    name: ClassVar[str] = "rsid_identity"
    display_name: ClassVar[str] = "dbSNP rsID identity"
    attribution: ClassVar[str] = "dbSNP"
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = False
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="public-domain",
        license_url="https://www.ncbi.nlm.nih.gov/home/about/policies/",
        attribution_text="Cross-build rsID placements derived from NCBI dbSNP.",
        source_url="https://www.ncbi.nlm.nih.gov/snp/",
        commercial_ok=True,
    )

    def __init__(self, data_dir: Path) -> None:
        """Bind to the data directory."""
        super().__init__(data_dir)
        self._db_path = data_dir / RSID_IDENTITY_DB_FILENAME
        self._conn: sqlite3.Connection | None = None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            if not self._db_path.exists():
                raise FileNotFoundError(
                    f"rsID identity cache not found at {self._db_path}. "
                    "Run `allelix db update` first."
                )
            self._conn = sqlite3.connect(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built identity cache and verify its SHA-256."""
        if not RSID_IDENTITY_EXPECTED_SHA256:
            raise RuntimeError(
                "rsid_identity: no published cache is pinned in this build of allelix"
            )
        gz_path = self.data_dir / f"{RSID_IDENTITY_DB_FILENAME}.gz"
        download(RSID_IDENTITY_CACHE_URL, gz_path)
        verify_file_hash(gz_path, "sha256", RSID_IDENTITY_EXPECTED_SHA256)
        install_prebuilt_cache(gz_path, self._db_path, source_url=RSID_IDENTITY_CACHE_URL)
        try:
            gz_path.unlink()
        except OSError:
            logger.warning("Could not remove staged file at %s", gz_path)

    def is_ready(self) -> bool:
        """True when the cache exists with the current schema version."""
        info = get_database_info(self._db_path, RSID_IDENTITY_RECORD_NAME)
        if info is None:
            return False
        return (info.get("local_version_tag") or "") == f"sv:{RSID_IDENTITY_SCHEMA_VERSION}"

    def version(self) -> str | None:
        """Return the dbSNP build the cache was made from, or None."""
        info = get_database_info(self._db_path, RSID_IDENTITY_RECORD_NAME)
        return info["version"] if info else None

    def record_count(self) -> int | None:
        """Return the number of identity rows, or None."""
        info = get_database_info(self._db_path, RSID_IDENTITY_RECORD_NAME)
        return info["record_count"] if info else None

    def close(self) -> None:
        """Release the SQLite connection. Safe to call repeatedly."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        """Pinned prebuilt cache: no remote freshness signal."""
        return None

    def cached_remote_signal(self) -> str | None:
        """Pinned prebuilt cache: no remote freshness signal."""
        return None

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Not used — the identity table never annotates. Always returns []."""
        return []

    def _rows_at_loci(self, build: str, loci: set[Locus], columns: str) -> list[tuple]:
        if build not in _COLUMNS or not loci:
            return []
        pos_col = _COLUMNS[build][0]
        conn = self._connection()
        locus_list = sorted(loci)
        out: list[tuple] = []
        for start in range(0, len(locus_list), _CHUNK):
            chunk = locus_list[start : start + _CHUNK]
            clauses = " OR ".join(f"(chrom = ? AND {pos_col} = ?)" for _ in chunk)
            params = [value for locus in chunk for value in locus]
            out.extend(
                conn.execute(
                    f"SELECT chrom, {pos_col}, {columns} FROM rsid_identity WHERE {clauses}",
                    params,
                )
            )
        return out

    def refs_at_loci(self, build: str, loci: set[Locus]) -> dict[Locus, set[str]]:
        """Distinct reference alleles at each locus on ``build``."""
        ref_col = _COLUMNS.get(build, ("", "", ""))[1]
        out: dict[Locus, set[str]] = {}
        if not ref_col:
            return out
        for chrom, pos, ref in self._rows_at_loci(build, loci, ref_col):
            out.setdefault((chrom, pos), set()).add(ref)
        return out

    def alleles_at_loci(
        self, build: str, loci: set[Locus]
    ) -> dict[Locus, list[tuple[str, str, str]]]:
        """``(ref, alt, rsid)`` rows at each locus on ``build`` (resolver input)."""
        if build not in _COLUMNS:
            return {}
        _, ref_col, alt_col = _COLUMNS[build]
        out: dict[Locus, list[tuple[str, str, str]]] = {}
        for chrom, pos, ref, alt, rsid in self._rows_at_loci(
            build, loci, f"{ref_col}, {alt_col}, rsid"
        ):
            out.setdefault((chrom, pos), []).append((ref, alt, rsid))
        return out

    def palindromic_rsids(self, rsids: set[str]) -> set[str]:
        """Return the subset of rsids flagged palindromic in the identity table."""
        if not rsids or not self.is_ready():
            return set()
        conn = self._connection()
        wanted = sorted(rsids)
        out: set[str] = set()
        for start in range(0, len(wanted), _CHUNK):
            chunk = wanted[start : start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for (rsid,) in conn.execute(
                f"SELECT DISTINCT rsid FROM rsid_identity "
                f"WHERE rsid IN ({marks}) AND palindromic = 1",
                chunk,
            ):
                out.add(rsid)
        return out

    def loci_of_rsids(self, rsids: set[str]) -> dict[str, dict[str, Locus]]:
        """Each rsID's ``{build: (chrom, pos)}`` placement on GRCh37 and GRCh38.

        Used to watch high-value positions in gVCF reference blocks, which
        carry no rsIDs. rsIDs the table does not hold are absent.
        """
        if not rsids:
            return {}
        conn = self._connection()
        wanted = sorted(rsids)
        out: dict[str, dict[str, Locus]] = {}
        for start in range(0, len(wanted), _CHUNK):
            chunk = wanted[start : start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for rsid, chrom, pos37, pos38 in conn.execute(
                f"SELECT DISTINCT rsid, chrom, pos37, pos38 FROM rsid_identity "
                f"WHERE rsid IN ({marks})",
                chunk,
            ):
                out[rsid] = {"GRCh37": (chrom, pos37), "GRCh38": (chrom, pos38)}
        return out

    def translate_to_grch38(self, keys37: set[Key]) -> dict[Key, Key]:
        """Map GRCh37 ``(chrom, pos, ref, alt)`` keys to their GRCh38 keys.

        A key translates only when exactly one GRCh38 key matches it, the
        site is not palindromic (the strand relation between builds would
        be a guess) and the allele is not the GRCh38 reference (there is no
        GRCh38 variant to look up). Everything else abstains.
        """
        loci = {(chrom, pos) for chrom, pos, _ref, _alt in keys37}
        candidates: dict[Key, set[Key]] = {}
        for chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic in self._rows_at_loci(
            "GRCh37", loci, "ref37, alt37, pos38, ref38, alt38, palindromic"
        ):
            key37 = (chrom, pos37, ref37, alt37)
            if key37 not in keys37 or palindromic or alt38 == ref38:
                continue
            candidates.setdefault(key37, set()).add((chrom, pos38, ref38, alt38))
        return {k: next(iter(v)) for k, v in candidates.items() if len(v) == 1}

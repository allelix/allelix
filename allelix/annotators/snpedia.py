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

from allelix.annotators.base import (
    Annotator,
    LicenseDescriptor,
    carried_alt,
    forward_alleles_at_call,
    is_clinvar_homref,
    variant_locus,
)
from allelix.databases.cache_state import (
    READY,
    CacheState,
    CacheStatus,
    cache_state,
    connect_readonly,
)
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
from allelix.utils.allele import (
    STRAND_ABSTAIN_REASONS,
    STRAND_COMPLEMENT,
    STRAND_DIRECT,
    STRAND_NO_FIT,
    STRAND_NOT_IN_IDENTITY,
    STRAND_PALINDROMIC,
    complement,
    forward_strand_reading,
    is_single_base,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

    from allelix.annotators.base import ForwardSiteProvider
    from allelix.annotators.rsid_identity import ForwardSite
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

_PALINDROMIC_SITE_ALLELES: frozenset[frozenset[str]] = frozenset(
    {frozenset({"A", "T"}), frozenset({"C", "G"})}
)


def _site_is_palindromic(rows: list[tuple]) -> bool:
    """Return whether all genotype rows identify an A/T or C/G SNP."""
    alleles = {
        allele
        for row in rows
        for allele in row[:2]
        if len(allele) == 1 and allele in {"A", "C", "G", "T"}
    }
    return frozenset(alleles) in _PALINDROMIC_SITE_ALLELES


def _forward_genotype(
    allele1: str,
    allele2: str,
    orientation: str | None,
) -> tuple[str, str]:
    """Return a SNPedia genotype normalized to the reference-forward strand."""
    if orientation == "minus":
        allele1, allele2 = complement(allele1), complement(allele2)
    return (allele1, allele2) if allele1 <= allele2 else (allele2, allele1)


def _known_orientation(orientation: str | None) -> str | None:
    """``plus`` or ``minus`` when the page states it, else None."""
    normalized = (orientation or "").strip().lower()
    return normalized if normalized in {"plus", "minus"} else None


def _orientation_without_identity(
    orientation: str | None, site_is_palindromic: bool
) -> str | None:
    """A row's orientation from the page alone: unknown reads as plus unless palindromic."""
    known = _known_orientation(orientation)
    if known is not None:
        return known
    return None if site_is_palindromic else "plus"


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
        clinvar_ref_provider: Callable[[str, str, tuple[str, int]], str | None] | None = None,
        normal_allele_provider: Callable[[str, str], str | None] | None = None,
        forward_site_provider: ForwardSiteProvider | None = None,
    ) -> None:
        """Initialize with path to the data directory.

        ``clinvar_ref_provider`` is a ``(rsid, build, locus) -> ref_base | None``
        callable used by the ADR-0023 hom-ref check. In production it is
        wired to ``ClinVarAnnotator.reference_at``, which answers only where
        ClinVar places the rsID at the call's locus. ``None`` disables the
        check (tests, standalone use).

        ``normal_allele_provider`` gives the normal allele where the build's
        reference carries a ClinVar variant allele, for deriving ``alt``
        (``base.carried_alt``); wired to ``ClinVarAnnotator.normal_allele_for``.

        ``forward_site_provider`` gives each rsID's forward alleles on a
        build (``RsidIdentityAnnotator.forward_sites``), the strand check
        for genotype rows. ``None`` keeps the page's orientation
        alone (tests, standalone use).
        """
        super().__init__(data_dir)
        self._db_path = data_dir / SNPEDIA_DB_FILENAME
        self._conn: sqlite3.Connection | None = None
        self._clinvar_ref_provider = clinvar_ref_provider
        self._normal_allele_provider = normal_allele_provider
        self._forward_site_provider = forward_site_provider
        self._strand_abstained = dict.fromkeys(STRAND_ABSTAIN_REASONS, 0)

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = connect_readonly(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built SNPedia cache from HuggingFace.

        The HuggingFace asset contains raw wiki markup (third-party,
        CC-BY-NC-SA). After download, ``prepare_cache()`` runs the
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
        self.prepare_cache()

    def is_ready(self) -> bool:
        """Return True when the parsed SNPedia genotype table exists and has data."""
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        """State of the SNPedia cache (ADR-0048).

        Raw pages not yet parsed, or parsed by an older parser, are
        ``outdated``: the update step parses them (``prepare_cache``).
        """
        status = cache_state(self._db_path, None)
        if not status.ready:
            return status
        from allelix.databases.snpedia_parser import (
            detect_raw_table,
            has_structured_table,
            parser_tag_is_current,
        )

        try:
            with contextlib.closing(connect_readonly(self._db_path)) as conn:
                if has_structured_table(conn):
                    if parser_tag_is_current(conn):
                        return READY
                    return CacheStatus(
                        CacheState.OUTDATED, self._db_path, "parser version changed"
                    )
                if detect_raw_table(conn) is not None:
                    return CacheStatus(CacheState.OUTDATED, self._db_path, "pages not parsed yet")
        except sqlite3.DatabaseError as exc:
            return CacheStatus(CacheState.CORRUPT, self._db_path, str(exc))
        return CacheStatus(
            CacheState.CORRUPT, self._db_path, "no rows in snpedia_genotypes and no raw pages"
        )

    def prepare_cache(self) -> bool:
        """Parse the raw pages into structured rows (one-time, ~2 minutes).

        Also runs when the parser version changed. A legacy cache whose
        ``remote_signal`` carries the current ``|pv:`` marker is stamped
        instead of re-parsed.
        """
        if self.cache_status().state is not CacheState.OUTDATED:
            return False
        try:
            from allelix.databases.snpedia_parser import (
                detect_raw_table,
                has_structured_table,
                parse_raw_pages,
                parser_is_current,
            )

            with contextlib.closing(sqlite3.connect(self._db_path)) as conn:
                has_rows = has_structured_table(conn)
                if has_rows and parser_is_current(conn):
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

            reason = "parser version changed" if has_rows else "one-time"
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
            with contextlib.closing(connect_readonly(self._db_path)) as conn:
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
            with contextlib.closing(connect_readonly(self._db_path)) as conn:
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

    def run_counts(self) -> dict[str, int]:
        """SNPedia sites not assessed for strand, by reason."""
        return dict(self._strand_abstained)

    def _row_orientations(
        self,
        variant: Variant,
        snp_id: str,
        rows: list[tuple],
        sites: dict[str, ForwardSite] | None,
    ) -> list[str | None]:
        """Each row's orientation for matching (``plus``/``minus``), or None to skip it.

        Rows of single bases at an rsID are read against the identity
        table's forward alleles: a stated orientation must put the
        genotype among them; an unknown one takes the direct reading
        of the site's alleles, else the complemented one. A/T and C/G
        sites, rsIDs the table does not place at the call, and genotypes
        that fit neither reading abstain, counted once per site and reason.
        Other rows (indels, i-IDs), and every row when no identity table is
        wired, keep the page's orientation with the palindrome guard.
        """
        site_is_palindromic = _site_is_palindromic(rows)
        out = [_orientation_without_identity(row[6], site_is_palindromic) for row in rows]
        if sites is None or not snp_id.startswith("rs"):
            return out
        snv = [is_single_base(row[0]) and is_single_base(row[1]) for row in rows]
        if not any(snv):
            return out
        forward = forward_alleles_at_call(variant, sites)
        unknown = {
            allele
            for row, is_snv in zip(rows, snv, strict=True)
            if is_snv and _known_orientation(row[6]) is None
            for allele in row[:2]
        }
        unknown_reading = (
            forward_strand_reading(unknown, forward) if unknown and forward is not None else None
        )
        reasons: set[str] = set()
        for i, (row, is_snv) in enumerate(zip(rows, snv, strict=True)):
            if not is_snv:
                continue
            out[i] = None
            if forward is None:
                reasons.add(STRAND_NOT_IN_IDENTITY)
                continue
            known = _known_orientation(row[6])
            if known is not None:
                reading = forward_strand_reading(_forward_genotype(row[0], row[1], known), forward)
                if reading == STRAND_DIRECT:
                    out[i] = known
                else:
                    reasons.add(
                        STRAND_PALINDROMIC if reading == STRAND_PALINDROMIC else STRAND_NO_FIT
                    )
            elif unknown_reading == STRAND_DIRECT:
                out[i] = "plus"
            elif unknown_reading == STRAND_COMPLEMENT:
                out[i] = "minus"
            elif unknown_reading is not None:
                reasons.add(unknown_reading)
        for reason in reasons:
            self._strand_abstained[reason] += 1
        return out

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return SNPedia annotations matching the user's genotype."""
        candidate = self._candidate(variant)
        if candidate is None:
            return []
        snp_id, snp_url_path, sorted_alleles = candidate
        rows = (
            self._connection()
            .execute(
                "SELECT allele1, allele2, magnitude, repute, summary, gene, orientation "
                "FROM snpedia_genotypes "
                "WHERE rsid = ?",
                (snp_id,),
            )
            .fetchall()
        )
        if not rows:
            return []
        sites = self._forward_sites({snp_id}, variant.build)
        return list(
            self._annotations_for(variant, snp_id, snp_url_path, sorted_alleles, rows, sites)
        )

    def _candidate(self, variant: Variant) -> tuple[str, str, tuple[str, str]] | None:
        """``(snp_id, page path, sorted alleles)`` for a call SNPedia can match, else None."""
        if variant.is_no_call:
            return None
        snp_id = variant.rsid.lower()
        if snp_id.startswith("rs"):
            snp_num = snp_id[2:]
            snp_url_path = f"Rs{snp_num}"
        elif snp_id.startswith("i"):
            snp_num = snp_id[1:]
            snp_url_path = f"I{snp_num}"
        else:
            return None
        if not snp_num or not snp_num.isdigit():
            return None
        if snp_id.startswith("rs") and is_clinvar_homref(variant, self._clinvar_ref_provider):
            return None
        a1, a2 = variant.allele1.upper(), variant.allele2.upper()
        sorted_alleles = (a1, a2) if a1 <= a2 else (a2, a1)
        return snp_id, snp_url_path, sorted_alleles

    def _forward_sites(self, rsids: set[str], build: str) -> dict[str, ForwardSite] | None:
        if self._forward_site_provider is None:
            return None
        return self._forward_site_provider({r for r in rsids if r.startswith("rs")}, build)

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
        candidates: list[tuple[Variant, str, str, tuple[str, str]]] = []
        for variant in variants:
            candidate = self._candidate(variant)
            if candidate is not None:
                candidates.append((variant, *candidate))

        if not candidates:
            return
        rsids = list({snp_id for _, snp_id, _, _ in candidates})

        conn = self._connection()
        rows_by_rsid: dict[str, list[tuple]] = {}
        for start in range(0, len(rsids), _BATCH_CHUNK):
            chunk = rsids[start : start + _BATCH_CHUNK]
            placeholders = ",".join("?" * len(chunk))
            cursor = conn.execute(
                f"SELECT rsid, allele1, allele2, magnitude, repute, summary, gene, "
                f"orientation "
                f"FROM snpedia_genotypes WHERE rsid IN ({placeholders})",
                chunk,
            )
            for row in cursor:
                rows_by_rsid.setdefault(row[0], []).append(row[1:])

        sites_by_build: dict[str, dict[str, ForwardSite] | None] = {}
        for build in {variant.build for variant, *_ in candidates}:
            wanted = {
                snp_id
                for variant, snp_id, _, _ in candidates
                if variant.build == build and snp_id in rows_by_rsid
            }
            sites_by_build[build] = self._forward_sites(wanted, build)

        for variant, snp_id, snp_url_path, sorted_alleles in candidates:
            rows = rows_by_rsid.get(snp_id)
            if not rows:
                continue
            yield from self._annotations_for(
                variant, snp_id, snp_url_path, sorted_alleles, rows, sites_by_build[variant.build]
            )

    def _annotations_for(
        self,
        variant: Variant,
        snp_id: str,
        snp_url_path: str,
        sorted_alleles: tuple[str, str],
        rows: list[tuple],
        sites: dict[str, ForwardSite] | None,
    ) -> Iterator[Annotation]:
        """The annotations of one call from its rsID's genotype rows."""
        orientations = self._row_orientations(variant, snp_id, rows, sites)
        for row, orientation in zip(rows, orientations, strict=True):
            if orientation is None:
                continue
            allele1, allele2, magnitude, repute, summary, gene, _ = row
            forward_genotype = _forward_genotype(allele1, allele2, orientation)
            if forward_genotype != sorted_alleles:
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
            forward_allele1, forward_allele2 = forward_genotype
            alt, reference_carries_alt = carried_alt(
                variant, forward_allele1, forward_allele2, self._normal_allele_provider
            )
            genotype_match = f"{forward_allele1}{forward_allele2}"
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
                alt=alt,
                locus_ref=variant.ref,
                reference_carries_alt=reference_carries_alt,
                carried_alleles=(forward_allele1, forward_allele2),
                locus=variant_locus(variant),
            )

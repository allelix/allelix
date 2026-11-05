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
   queried for GRCh37 input (``translate_to_grch38``). A row flagged
   palindromic whose GRCh37 and GRCh38 reference bases differ abstains, and
   so does an allele that is the GRCh38 reference.

dbSNP is a US government work in the public domain.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases._versions import RSID_IDENTITY_SCHEMA_VERSION
from allelix.databases.cache_state import (
    CacheStatus,
    cache_state,
    connect_readonly,
)
from allelix.databases.manager import download, get_database_info, verify_file_hash
from allelix.databases.rsid_identity_loader import (
    RSID_IDENTITY_CACHE_URL,
    RSID_IDENTITY_DB_FILENAME,
    RSID_IDENTITY_EXPECTED_SHA256,
    RSID_IDENTITY_RECORD_NAME,
    install_prebuilt_cache,
)
from allelix.databases.staging import owned_directory
from allelix.utils.allele import complement

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterable, Iterator
    from pathlib import Path

    from allelix.models import Annotation, Variant


_CHUNK = 450  # two parameters per locus; stays under SQLite's 999 limit

# Build -> (position, ref, alt) column names. Fixed mapping, never user input.
_COLUMNS: dict[str, tuple[str, str, str]] = {
    "GRCh37": ("pos37", "ref37", "alt37"),
    "GRCh38": ("pos38", "ref38", "alt38"),
}

Locus = tuple[str, int]
Key = tuple[str, int, str, str]
ForwardSite = tuple[str, int, frozenset[str]]
# rsid -> (chrom, pos37, ref37, pos38, ref38, same_strand); same_strand is
# None where the strand relation between the builds cannot be told.
ReferenceSite = tuple[str, int, str, int, str, bool | None]


def strand_relation(alleles37: set[str], alleles38: set[str], palindromic: bool) -> bool | None:
    """True for the same strand on both builds, False for inverted, None when it cannot be told.

    ``alleles37`` and ``alleles38`` are an rsID's REF and ALTs on each
    build, as ``pair_alleles`` pairs them. A site flagged palindromic, or
    whose allele sets fit both readings or neither, decides nothing.
    """
    if palindromic:
        return None
    same = alleles37 == alleles38
    inverted = all(len(a) == 1 for a in alleles37 | alleles38) and alleles37 == {
        complement(a) for a in alleles38
    }
    if same == inverted:
        return None
    return same


def on_other_strand(allele: str, same_strand: bool) -> str:
    """``allele`` as the other build writes it, given the strand relation."""
    return allele if same_strand else complement(allele)


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
            self._conn = connect_readonly(self._db_path)
        return self._conn

    def setup(self) -> None:
        """Download the pre-built identity cache and verify its SHA-256."""
        if not RSID_IDENTITY_EXPECTED_SHA256:
            raise RuntimeError("no published cache is pinned in this build of allelix")
        # This setup's own download directory: another setup of the same
        # cache can neither replace nor remove the file it verifies.
        with owned_directory(self.data_dir, "rsid_identity") as staging:
            gz_path = staging / f"{RSID_IDENTITY_DB_FILENAME}.gz"
            download(RSID_IDENTITY_CACHE_URL, gz_path)
            verify_file_hash(gz_path, "sha256", RSID_IDENTITY_EXPECTED_SHA256)
            install_prebuilt_cache(gz_path, self._db_path, source_url=RSID_IDENTITY_CACHE_URL)

    def is_ready(self) -> bool:
        """True when the cache exists with the current schema version."""
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        """State of the identity cache (ADR-0048)."""
        return cache_state(
            self._db_path,
            RSID_IDENTITY_RECORD_NAME,
            tag=f"sv:{RSID_IDENTITY_SCHEMA_VERSION}",
            tables=("rsid_identity",),
        )

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

    def reference_differing_sites(self) -> dict[str, ReferenceSite]:
        """The rsIDs whose GRCh37 and GRCh38 REFs are different single bases.

        ``rsid -> (chrom, pos37, ref37, pos38, ref38, same_strand)``. Either
        one build's reference carries the other's variant there (F5 rs6025),
        or the reference strand was inverted between the assemblies and the
        two REFs are one allele written on opposite strands. ``same_strand``
        (``strand_relation`` over all the rsID's rows) tells which spelling
        an allele from one build has on the other; None where it cannot be
        told, and nothing may then be carried across. An rsID whose rows
        disagree on the placement or the REFs is left out. Read once, for
        ``ReferenceCarriedFallback``.
        """
        if not self.is_ready():
            return {}
        found: dict[str, set[tuple[str, int, str, int, str]]] = {}
        alleles: dict[str, tuple[set[str], set[str]]] = {}
        palindromic: set[str] = set()
        rows = self._connection().execute(
            "SELECT rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic "
            "FROM rsid_identity WHERE rsid IN (SELECT rsid FROM rsid_identity "
            "WHERE length(ref37) = 1 AND length(ref38) = 1 AND ref37 <> ref38)"
        )
        for rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, flag in rows:
            found.setdefault(rsid, set()).add((chrom, pos37, ref37, pos38, ref38))
            sets = alleles.setdefault(rsid, (set(), set()))
            sets[0].update((ref37, alt37))
            sets[1].update((ref38, alt38))
            if flag:
                palindromic.add(rsid)
        out: dict[str, ReferenceSite] = {}
        for rsid, sites in found.items():
            if len(sites) != 1:
                continue
            chrom, pos37, ref37, pos38, ref38 = next(iter(sites))
            if len(ref37) != 1 or len(ref38) != 1 or ref37 == ref38:
                continue
            same = strand_relation(*alleles[rsid], rsid in palindromic)
            out[rsid] = (chrom, pos37, ref37, pos38, ref38, same)
        return out

    def normal_alleles_on_other_build(self, carried: dict[str, str], build: str) -> dict[str, str]:
        """The normal allele at rsIDs whose ``build`` reference carries a variant allele.

        ``carried`` maps each rsID to the allele ``build``'s reference
        carries (ClinVar's REF == ALT row; T at F5 rs6025 on GRCh37). The
        normal allele is the other build's REF, written on ``build``'s
        strand: every identity row of the rsID must agree on both REFs,
        ``build``'s REF must be the carried allele, the other REF one base,
        and the strand relation between the builds decided
        (``strand_relation``). A row flagged palindromic decides no strand,
        so such an rsID abstains, rs6025 included (dbSNP lists every base
        there). So does an rsID whose other REF, on ``build``'s strand, is
        the carried allele itself: the reference did not change. Used by
        ClinVar when only one ClinVar build is installed.
        """
        if not carried or build not in _COLUMNS or not self.is_ready():
            return {}
        conn = self._connection()
        wanted = sorted(carried)
        refs: dict[str, set[tuple[str, str]]] = {}
        alleles: dict[str, tuple[set[str], set[str]]] = {}
        palindromic: set[str] = set()
        for start in range(0, len(wanted), _CHUNK):
            chunk = wanted[start : start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for rsid, ref37, alt37, ref38, alt38, flag in conn.execute(
                f"SELECT rsid, ref37, alt37, ref38, alt38, palindromic "
                f"FROM rsid_identity WHERE rsid IN ({marks})",
                chunk,
            ):
                refs.setdefault(rsid, set()).add((ref37, ref38))
                sets = alleles.setdefault(rsid, (set(), set()))
                sets[0].update((ref37, alt37))
                sets[1].update((ref38, alt38))
                if flag:
                    palindromic.add(rsid)
        out: dict[str, str] = {}
        for rsid, pairs in refs.items():
            if len(pairs) != 1:
                continue
            ref37, ref38 = next(iter(pairs))
            ref, other_ref = (ref37, ref38) if build == "GRCh37" else (ref38, ref37)
            if ref != carried[rsid] or len(other_ref) != 1:
                continue
            same = strand_relation(*alleles[rsid], rsid in palindromic)
            if same is None:
                continue
            normal = on_other_strand(other_ref, same)
            if normal != ref:
                out[rsid] = normal
        return out

    def forward_sites(self, rsids: set[str], build: str) -> dict[str, ForwardSite]:
        """Each rsID's ``(chrom, pos, alleles)`` on ``build``: REF plus every ALT.

        The strand reference for GWAS and SNPedia alleles. An rsID
        whose rows disagree on the locus is left out, as are rsIDs the
        table does not hold. Nothing is returned when the cache is not ready.
        """
        if not rsids or build not in _COLUMNS or not self.is_ready():
            return {}
        pos_col, ref_col, alt_col = _COLUMNS[build]
        conn = self._connection()
        wanted = sorted(rsids)
        loci: dict[str, set[Locus]] = {}
        alleles: dict[str, set[str]] = {}
        for start in range(0, len(wanted), _CHUNK):
            chunk = wanted[start : start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for rsid, chrom, pos, ref, alt in conn.execute(
                f"SELECT rsid, chrom, {pos_col}, {ref_col}, {alt_col} "
                f"FROM rsid_identity WHERE rsid IN ({marks})",
                chunk,
            ):
                loci.setdefault(rsid, set()).add((chrom, pos))
                alleles.setdefault(rsid, set()).update((ref, alt))
        return {
            rsid: (*next(iter(places)), frozenset(alleles[rsid]))
            for rsid, places in loci.items()
            if len(places) == 1
        }

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

    def placements(
        self, rsids: set[str]
    ) -> dict[str, list[tuple[str, int, str, str, int, str, str]]]:
        """Every identity row of each rsID: ``(chrom, pos37, ref37, alt37, pos38, ref38, alt38)``.

        rsIDs the table does not hold are absent; nothing is returned when
        the cache is not ready.
        """
        if not rsids or not self.is_ready():
            return {}
        conn = self._connection()
        wanted = sorted(rsids)
        out: dict[str, list[tuple[str, int, str, str, int, str, str]]] = {}
        for start in range(0, len(wanted), _CHUNK):
            chunk = wanted[start : start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for rsid, *row in conn.execute(
                f"SELECT rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38 "
                f"FROM rsid_identity WHERE rsid IN ({marks})",
                chunk,
            ):
                out.setdefault(rsid, []).append(tuple(row))
        return out

    def iter_placements(
        self, rsids: Iterable[str]
    ) -> Iterator[tuple[str, str, int, str, str, int, str, str]]:
        """``placements`` as a stream: ``(rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38)``.

        The rsIDs are read and queried in chunks, so neither they nor the
        rows are held at once. Nothing when the cache is not ready.
        """
        if not self.is_ready():
            return
        conn = self._connection()
        chunk: list[str] = []
        for rsid in rsids:
            chunk.append(rsid)
            if len(chunk) == _CHUNK:
                yield from self._placement_rows(conn, chunk)
                chunk = []
        if chunk:
            yield from self._placement_rows(conn, chunk)

    @staticmethod
    def _placement_rows(
        conn: sqlite3.Connection, chunk: list[str]
    ) -> Iterator[tuple[str, str, int, str, str, int, str, str]]:
        marks = ",".join("?" * len(chunk))
        yield from conn.execute(
            f"SELECT rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38 "
            f"FROM rsid_identity WHERE rsid IN ({marks})",
            chunk,
        )

    def translate_to_grch38(self, keys37: set[Key]) -> dict[Key, Key]:
        """Map GRCh37 ``(chrom, pos, ref, alt)`` keys to their GRCh38 keys.

        A key translates only when exactly one GRCh38 key matches it, the
        allele is not the GRCh38 reference (there is no GRCh38 variant to
        look up), and the row's strand relation is not a guess. A row the
        table flags palindromic is a guess only where its GRCh37 and GRCh38
        reference bases differ; with equal references it translates as
        written. The flag marks a whole site, and the cache sets it on many
        three-allele sites whose pairs are not in doubt. Everything else
        abstains.
        """
        loci = {(chrom, pos) for chrom, pos, _ref, _alt in keys37}
        candidates: dict[Key, set[Key]] = {}
        for chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic in self._rows_at_loci(
            "GRCh37", loci, "ref37, alt37, pos38, ref38, alt38, palindromic"
        ):
            key37 = (chrom, pos37, ref37, alt37)
            if key37 not in keys37 or alt38 == ref38 or (palindromic and ref37 != ref38):
                continue
            candidates.setdefault(key37, set()).add((chrom, pos38, ref38, alt38))
        return {k: next(iter(v)) for k, v in candidates.items() if len(v) == 1}

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Cross-build rsID identity table: builder logic and prebuilt-cache install.

The table (ADR-0040) maps each allele of an rsID between its GRCh37 and
GRCh38 placements. It is built once from the two dbSNP VCFs (GRCh37
``GCF_000001405.25`` and GRCh38 ``GCF_000001405.40``), restricted to the
rsIDs the annotators hold, by ``scripts/build_rsid_identity_cache.py``,
and shipped as a gzipped SQLite cache pinned by SHA-256 like the other
prebuilt caches. dbSNP is a US government work in the public domain.

The pure functions here (accession mapping, VCF row parsing, allele
pairing) carry all of the build logic so they can be tested without
dbSNP.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING, TextIO

from allelix.databases._versions import RSID_IDENTITY_SCHEMA_VERSION
from allelix.databases.loader_utils import install_prebuilt_gz_cache
from allelix.databases.schema import RSID_IDENTITY_SCHEMA
from allelix.utils.allele import complement

if TYPE_CHECKING:
    from collections.abc import Container, Iterable, Iterator
    from pathlib import Path

RSID_IDENTITY_DB_FILENAME = "rsid_identity.sqlite"
RSID_IDENTITY_RECORD_NAME = "rsid_identity"

RSID_IDENTITY_CACHE_URL = (
    "https://huggingface.co/datasets/allelix/allelix-rsid-identity"
    "/resolve/main/rsid_identity.sqlite.gz"
)

RSID_IDENTITY_EXPECTED_SHA256 = "2a00787cd0be4e7572efd21686013c49f3d0180137cd4403730524d50a5db285"

DBSNP_GRCH37_VCF_URL = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
DBSNP_GRCH38_VCF_URL = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"

# (rsid, chrom, pos, ref, alts)
Placement = tuple[str, str, int, str, tuple[str, ...]]
# (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic)
IdentityRow = tuple[str, str, int, str, str, int, str, str, int]


def refseq_to_chrom(accession: str) -> str | None:
    """Map a primary-assembly RefSeq accession to Allelix's chromosome name.

    ``NC_000001.10`` / ``NC_000001.11`` -> ``"1"`` ... ``NC_000023`` -> ``X``,
    ``NC_000024`` -> ``Y``, ``NC_012920`` -> ``MT``. Patches, alternate loci
    and unplaced scaffolds (``NT_`` / ``NW_``) return None: their placements
    are not comparable across assemblies.
    """
    stem = accession.split(".", 1)[0]
    if stem == "NC_012920":
        return "MT"
    if not stem.startswith("NC_0000"):
        return None
    try:
        number = int(stem[len("NC_") :])
    except ValueError:
        return None
    if 1 <= number <= 22:
        return str(number)
    if number == 23:
        return "X"
    if number == 24:
        return "Y"
    return None


def parse_dbsnp_line(line: str, wanted: Container[str]) -> Placement | None:
    """Parse one dbSNP VCF data line; None unless it is a wanted primary placement."""
    if line.startswith("#"):
        return None
    cols = line.split("\t", 5)
    if len(cols) < 5:
        return None
    rsid = cols[2]
    if rsid not in wanted:
        return None
    chrom = refseq_to_chrom(cols[0])
    if chrom is None:
        return None
    try:
        pos = int(cols[1])
    except ValueError:
        return None
    ref = cols[3].upper()
    alts = tuple(a.upper() for a in cols[4].split(",") if a and a != ".")
    if not alts:
        return None
    return rsid, chrom, pos, ref, alts


def pair_alleles(
    ref37: str, alts37: Iterable[str], ref38: str, alts38: Iterable[str]
) -> tuple[list[tuple[str, str, str, str]], bool] | None:
    """Pair an rsID's alleles across builds.

    The allele sets (REF plus ALTs) must be equal on the same strand, or,
    for SNVs, equal after complementing (the reference strand was inverted
    between assemblies). Anything else returns None: the placements cannot
    be paired from sequence alone, and the rsID is left out of the table.

    Returns ``(rows, palindromic)`` where each row is
    ``(ref37, allele37, ref38, allele38)`` for every allele that is
    non-reference on at least one build. ``palindromic`` is True when both
    the same-strand and the complemented pairing fit (A/T, C/G, or all four
    bases): the strand relation is undecidable and translation must abstain.
    """
    set37 = {ref37, *alts37}
    set38 = {ref38, *alts38}
    is_snv = all(len(a) == 1 for a in set37 | set38)
    same = set37 == set38
    flipped = is_snv and set37 == {complement(a) for a in set38}
    if not same and not flipped:
        return None
    rows = []
    for allele37 in sorted(set37):
        allele38 = allele37 if same else complement(allele37)
        if allele37 == ref37 and allele38 == ref38:
            continue
        rows.append((ref37, allele37, ref38, allele38))
    return rows, same and flipped


def pair_grouped_placements(
    groups: Iterable[tuple[str, list[Placement], list[Placement]]],
    stats: dict[str, int] | None = None,
) -> Iterator[IdentityRow]:
    """Turn ``(rsid, grch37 placements, grch38 placements)`` groups into rows.

    An rsID contributes only when it has exactly one primary-assembly
    placement on each build, on the same chromosome, and its alleles pair
    (``pair_alleles``). Every exclusion is counted in ``stats``. Groups are
    consumed one at a time, so memory does not grow with the input.
    """
    counts = stats if stats is not None else {}
    counts.setdefault("only_grch37", 0)
    counts.setdefault("only_grch38", 0)
    for rsid, p37, p38 in groups:
        if not p38:
            counts["only_grch37"] += 1
            continue
        if not p37:
            counts["only_grch38"] += 1
            continue
        if len(p37) != 1 or len(p38) != 1:
            counts["multi_placement"] = counts.get("multi_placement", 0) + 1
            continue
        _, chrom37, pos37, ref37, alts37 = p37[0]
        _, chrom38, pos38, ref38, alts38 = p38[0]
        if chrom37 != chrom38:
            counts["chrom_changed"] = counts.get("chrom_changed", 0) + 1
            continue
        paired = pair_alleles(ref37, alts37, ref38, alts38)
        if paired is None:
            counts["unpairable"] = counts.get("unpairable", 0) + 1
            continue
        rows, palindromic = paired
        counts["rsids"] = counts.get("rsids", 0) + 1
        for r37, a37, r38, a38 in rows:
            counts["rows"] = counts.get("rows", 0) + 1
            yield rsid, chrom37, pos37, r37, a37, pos38, r38, a38, int(palindromic)


def build_identity_rows(
    placements37: dict[str, list[Placement]],
    placements38: dict[str, list[Placement]],
    stats: dict[str, int] | None = None,
) -> Iterator[IdentityRow]:
    """In-memory form of ``pair_grouped_placements`` for small inputs."""
    rsids = sorted(placements37.keys() | placements38.keys())
    yield from pair_grouped_placements(
        ((r, placements37.get(r, []), placements38.get(r, [])) for r in rsids), stats
    )


# ── Low-memory build pipeline (disk + SQLite) ──────────────────

_U64 = (1 << 64) - 1
_K1 = 0x9E3779B97F4A7C15
_K2 = 0xC2B2AE3D27D4EB4F


class RsidBloom:
    """Fixed-size Bloom filter over rsID numbers.

    Lets the dbSNP stream filter test membership against millions of
    rsIDs in a few MB instead of a set. False positives are possible and
    are removed exactly by the SQLite join against the universe table;
    false negatives are not.
    """

    def __init__(self, bits: int = 26, hashes: int = 4) -> None:
        """Allocate ``2**bits`` bits and use ``hashes`` probes per key."""
        self._bits = bits
        self._mask = (1 << bits) - 1
        self._shift = 64 - bits
        self._hashes = hashes
        self._table = bytearray((1 << bits) >> 3)

    def _probes(self, number: int) -> Iterator[int]:
        a = ((number * _K1) & _U64) >> self._shift
        b = (((number * _K2) & _U64) >> self._shift) | 1
        for i in range(self._hashes):
            yield (a + i * b) & self._mask

    def add(self, rsid: str) -> None:
        """Add an ``rs<digits>`` identifier."""
        for x in self._probes(int(rsid[2:])):
            self._table[x >> 3] |= 1 << (x & 7)

    def __contains__(self, rsid: object) -> bool:
        """True when ``rsid`` may be in the set (False means certainly not)."""
        if not isinstance(rsid, str) or not rsid.startswith("rs") or not rsid[2:].isdigit():
            return False
        table = self._table
        return all(table[x >> 3] >> (x & 7) & 1 for x in self._probes(int(rsid[2:])))


def filter_dbsnp_stream(lines: Iterable[str], bloom: RsidBloom, out: TextIO) -> tuple[int, int]:
    """Write ``rsid chrom pos ref alts`` for candidate placements; return (read, written).

    Candidates are primary-assembly placements whose rsID passes ``bloom``.
    The output is a small TSV that ``load_placements`` takes exactly.
    """
    read = written = 0
    for line in lines:
        read += 1
        placement = parse_dbsnp_line(line, bloom)
        if placement is None:
            continue
        rsid, chrom, pos, ref, alts = placement
        out.write(f"{rsid}\t{chrom}\t{pos}\t{ref}\t{','.join(alts)}\n")
        written += 1
    return read, written


def create_work_db(conn: sqlite3.Connection) -> None:
    """Create the build's scratch tables (universe and per-build placements)."""
    conn.execute("PRAGMA temp_store = FILE")
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS universe (rsid TEXT PRIMARY KEY) WITHOUT ROWID;
        CREATE TABLE IF NOT EXISTS placements (
            build TEXT NOT NULL, rsid TEXT NOT NULL, chrom TEXT NOT NULL,
            pos INTEGER NOT NULL, ref TEXT NOT NULL, alts TEXT NOT NULL
        );
        """
    )


def add_universe(conn: sqlite3.Connection, rsids: Iterable[str]) -> None:
    """Insert rsIDs (``rs<digits>`` only) into the universe table, de-duplicated."""
    conn.executemany(
        "INSERT OR IGNORE INTO universe VALUES (?)",
        ((r,) for r in rsids if r and r.startswith("rs") and r[2:].isdigit()),
    )
    conn.commit()


def bloom_from_universe(conn: sqlite3.Connection, bits: int = 26) -> RsidBloom:
    """Build the stream filter from the universe table."""
    bloom = RsidBloom(bits=bits)
    for (rsid,) in conn.execute("SELECT rsid FROM universe"):
        bloom.add(rsid)
    return bloom


def load_placements(conn: sqlite3.Connection, build: str, lines: Iterable[str]) -> int:
    """Load a filtered TSV for ``build``, keeping only exact universe members."""
    conn.execute("CREATE TEMP TABLE IF NOT EXISTS staged (rsid, chrom, pos, ref, alts)")
    conn.execute("DELETE FROM staged")
    conn.executemany(
        "INSERT INTO staged VALUES (?, ?, ?, ?, ?)",
        (
            (rsid, chrom, int(pos), ref, alts)
            for rsid, chrom, pos, ref, alts in (ln.rstrip("\n").split("\t") for ln in lines)
        ),
    )
    cur = conn.execute(
        "INSERT INTO placements SELECT ?, s.rsid, s.chrom, s.pos, s.ref, s.alts "
        "FROM staged s JOIN universe u ON u.rsid = s.rsid",
        (build,),
    )
    conn.execute("DELETE FROM staged")
    conn.commit()
    return cur.rowcount


def grouped_placements(
    conn: sqlite3.Connection,
) -> Iterator[tuple[str, list[Placement], list[Placement]]]:
    """Stream ``(rsid, grch37, grch38)`` placement groups in rsID order."""
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_placements_rsid ON placements(rsid, build, chrom, pos)"
    )
    current: str | None = None
    p37: list[Placement] = []
    p38: list[Placement] = []
    for build, rsid, chrom, pos, ref, alts in conn.execute(
        "SELECT build, rsid, chrom, pos, ref, alts FROM placements "
        "ORDER BY rsid, build, chrom, pos"
    ):
        if rsid != current:
            if current is not None:
                yield current, p37, p38
            current, p37, p38 = rsid, [], []
        placement = (rsid, chrom, pos, ref, tuple(alts.split(",")))
        (p37 if build == "GRCh37" else p38).append(placement)
    if current is not None:
        yield current, p37, p38


def write_identity_db(
    rows: Iterable[IdentityRow],
    db_path: Path,
    *,
    version: str,
    source_url: str,
    downloaded_at: str,
) -> int:
    """Write identity rows into a fresh SQLite cache; return the row count."""
    if db_path.exists():
        db_path.unlink()
    count = 0
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.executescript(RSID_IDENTITY_SCHEMA)
        batch: list[IdentityRow] = []
        for row in rows:
            batch.append(row)
            if len(batch) >= 50_000:
                conn.executemany("INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)", batch)
                count += len(batch)
                batch.clear()
        if batch:
            conn.executemany("INSERT INTO rsid_identity VALUES (?,?,?,?,?,?,?,?,?)", batch)
            count += len(batch)
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES (?, ?, ?, ?, ?, ?)",
            (
                RSID_IDENTITY_RECORD_NAME,
                source_url,
                version,
                downloaded_at,
                count,
                f"sv:{RSID_IDENTITY_SCHEMA_VERSION}",
            ),
        )
        conn.commit()
    return count


def install_prebuilt_cache(
    gz_path: Path,
    db_path: Path,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
) -> None:
    """Decompress a gzipped pre-built identity cache into place."""
    install_prebuilt_gz_cache(
        gz_path,
        db_path,
        RSID_IDENTITY_RECORD_NAME,
        source_url=source_url,
        remote_signal=remote_signal,
        schema_version_tag=f"sv:{RSID_IDENTITY_SCHEMA_VERSION}",
    )

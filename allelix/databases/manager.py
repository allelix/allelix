# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar VCF download, parse, and load into SQLite."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import logging
import os
import sqlite3
import urllib.request
from datetime import UTC, datetime
from typing import TYPE_CHECKING, TypedDict

from allelix import __version__
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

CLINVAR_URL_GRCH37 = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
CLINVAR_URL_GRCH38 = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

# ADR-0021: per-build cache filenames. The annotator holds at most one
# connection per build; each is independent.
CLINVAR_URL_BY_BUILD: dict[str, str] = {
    "GRCh37": CLINVAR_URL_GRCH37,
    "GRCh38": CLINVAR_URL_GRCH38,
}
INSERT_BATCH_SIZE = 5_000
DOWNLOAD_TIMEOUT_SECONDS = 60
SIGNAL_TIMEOUT_SECONDS = 15
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
USER_AGENT = f"allelix/{__version__} (+https://github.com/allelix/allelix)"


class DatabaseInfo(TypedDict):
    """Cached database version metadata."""

    source_url: str
    version: str
    downloaded_at: str
    record_count: int
    remote_signal: str | None
    local_version_tag: str | None


def fetch_remote_text(url: str, timeout: float = SIGNAL_TIMEOUT_SECONDS) -> str | None:
    """Fetch a small text resource (e.g., a `.md5` file) and return its body.

    Returns None on any failure — `db update`'s freshness check treats
    `None` as "can't verify" and falls through to a "skip with notice".
    Never raises.
    """
    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except (OSError, ValueError) as exc:
        if hasattr(exc, "close"):
            exc.close()
        return None


def head_request_headers(
    url: str, timeout: float = SIGNAL_TIMEOUT_SECONDS
) -> dict[str, str] | None:
    """Issue an HTTP HEAD and return the response headers as a plain dict.

    Returns None on any failure. Never raises.
    """
    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD")
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return dict(response.headers.items())
    except (OSError, ValueError) as exc:
        if hasattr(exc, "close"):
            exc.close()
        return None


def download(url: str, dest: Path) -> None:
    """Download `url` to `dest`. Streaming, atomic (.part rename), with timeout.

    - Streams chunks directly to a `.part` sibling file (no in-memory copy).
    - Sets a real User-Agent so CDNs don't reject the default python-urllib UA.
    - `os.replace`s the .part onto `dest` only after a full successful write,
      so a killed mid-download never leaves a half-file at the target name.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part_path = dest.parent / f"{dest.name}.part"
    if part_path.exists():
        part_path.unlink()

    logger.info("Downloading %s -> %s", url, dest)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with (
            urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response,
            part_path.open("wb") as out,
        ):
            expected_size = response.headers.get("Content-Length")
            while chunk := response.read(DOWNLOAD_CHUNK_SIZE):
                out.write(chunk)
            out.flush()
            try:
                os.fsync(out.fileno())
            except OSError:
                logger.debug("fsync unsupported on this filesystem; continuing")

        actual_size = part_path.stat().st_size
        if expected_size is not None:
            expected = int(expected_size)
            if actual_size != expected:
                part_path.unlink(missing_ok=True)
                raise OSError(
                    f"Download truncated: expected {expected:,} bytes, "
                    f"got {actual_size:,} bytes from {url}"
                )
        os.replace(part_path, dest)
    except Exception as exc:
        if hasattr(exc, "close"):
            exc.close()
        if part_path.exists():
            try:
                part_path.unlink()
            except OSError:
                logger.warning("Could not remove stale partial download %s", part_path)
        raise


def verify_file_hash(path: Path, algorithm: str, expected_hex: str) -> None:
    """Verify a file's cryptographic hash and delete it on mismatch.

    Reads the file in streaming chunks to avoid loading multi-GB files
    into memory. On mismatch, deletes the file and raises ``OSError``.
    """
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        while chunk := f.read(DOWNLOAD_CHUNK_SIZE):
            h.update(chunk)
    actual = h.hexdigest()
    if actual != expected_hex:
        path.unlink(missing_ok=True)
        raise OSError(
            f"Integrity check failed for {path.name}: "
            f"expected {algorithm}:{expected_hex}, "
            f"got {algorithm}:{actual}"
        )


def _parse_info(info: str) -> dict[str, str]:
    """Parse a VCF INFO field (`KEY=VALUE;FLAG;...`) into a dict."""
    out: dict[str, str] = {}
    for entry in info.split(";"):
        if "=" in entry:
            key, _, value = entry.partition("=")
            out[key] = value
        else:
            out[entry] = ""
    return out


def _pick(values: list[str], index: int) -> str:
    """Index into a list of parallel-indexed VCF INFO values, padding with last.

    Pads instead of zipping strictly because real ClinVar reliably parallel-
    indexes CLNSIG/CLNDN/ALLELEID with ALT, but third-party VCFs sometimes
    publish a single CLNSIG that applies to all ALTs. Falling through to the
    last value is the more permissive interpretation; a strict zip would drop
    annotations on those rows.
    """
    if not values:
        return ""
    if index < len(values):
        return values[index]
    return values[-1]


def parse_clinvar_version(vcf_path: Path) -> str | None:
    """Extract `##fileDate=YYYYMMDD` from a ClinVar VCF header, or None."""
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            if not raw.startswith("##"):
                return None
            if raw.startswith("##fileDate="):
                return raw.removeprefix("##fileDate=").strip()
    return None


def iter_clinvar_records(vcf_path: Path) -> Iterator[dict[str, object]]:
    """Stream parse a ClinVar VCF (.vcf or .vcf.gz). Skip entries without an RS id.

    Multi-allelic rows (ALT="A,T") are split into one record per ALT. Parallel
    INFO fields (CLNSIG, CLNDN, ALLELEID) are separated by `|` per ClinVar's
    convention and index-paired with the ALTs.
    """
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                logger.warning("Skipping ClinVar line with %d columns", len(parts))
                continue
            chrom, pos_str, _vid, ref, alt_field, _qual, _filter, info = parts[:8]
            try:
                pos = int(pos_str)
            except ValueError:
                logger.warning("Skipping ClinVar entry with non-integer position %r", pos_str)
                continue
            info_dict = _parse_info(info)
            rs = info_dict.get("RS")
            if not rs:
                continue

            alts = alt_field.split(",")
            clnsigs = info_dict.get("CLNSIG", "").split("|") if info_dict.get("CLNSIG") else [""]
            clndns = info_dict.get("CLNDN", "").split("|") if info_dict.get("CLNDN") else [""]
            allele_ids = (
                info_dict.get("ALLELEID", "").split("|") if info_dict.get("ALLELEID") else [""]
            )
            review_status = info_dict.get("CLNREVSTAT", "")
            gene = _extract_gene(info_dict.get("GENEINFO", ""))

            for i, alt in enumerate(alts):
                yield {
                    "rsid": f"rs{rs}",
                    "chromosome": chrom,
                    "position": pos,
                    "ref": ref,
                    "alt": alt,
                    "clinical_significance": _pick(clnsigs, i),
                    "condition": _pick(clndns, i).replace("_", " "),
                    "gene": gene,
                    "review_status": review_status,
                    "allele_id": _safe_int(_pick(allele_ids, i)),
                }


def _extract_gene(geneinfo: str) -> str:
    """`GENEINFO=BRCA1:672|...` → `"BRCA1"`."""
    if not geneinfo:
        return ""
    return geneinfo.split(":", 1)[0].split("|", 1)[0]


def _safe_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def load_clinvar_vcf(
    vcf_path: Path,
    db_path: Path,
    source_url: str = "",
    remote_signal: str | None = None,
    record_name: str = "clinvar",
) -> int:
    """Parse a ClinVar VCF into a fresh SQLite cache atomically.

    Writes to a `.tmp` sibling and `os.replace`s onto `db_path` only after a
    successful commit. If parsing fails mid-load, the previous cache (if any)
    is left intact.

    `remote_signal` is the value `fetch_remote_signal` returned at the time
    of this download; stored alongside version metadata so the next
    `db update` can detect remote changes without re-downloading.

    `record_name` is the key under which the version row is stored. ADR-0021:
    per-build ClinVar caches use `"clinvar.GRCh37"` / `"clinvar.GRCh38"` so
    the same data_dir can hold both. Default `"clinvar"` is the legacy
    single-cache identifier.

    Returns the number of records loaded.
    """
    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    file_date = parse_clinvar_version(vcf_path)
    version = file_date or datetime.now(UTC).strftime("%Y-%m-%d")

    try:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            insert_sql = (
                "INSERT INTO clinvar_variants "
                "(rsid, chromosome, position, ref, alt, clinical_significance, "
                "condition, gene, review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            )
            batch: list[tuple] = []
            count = 0
            for record in iter_clinvar_records(vcf_path):
                batch.append(
                    (
                        record["rsid"],
                        record["chromosome"],
                        record["position"],
                        record["ref"],
                        record["alt"],
                        record["clinical_significance"],
                        record["condition"],
                        record["gene"],
                        record["review_status"],
                        record["allele_id"],
                    )
                )
                if len(batch) >= INSERT_BATCH_SIZE:
                    conn.executemany(insert_sql, batch)
                    count += len(batch)
                    batch.clear()
            if batch:
                conn.executemany(insert_sql, batch)
                count += len(batch)
            from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    record_name,
                    source_url,
                    version,
                    datetime.now(UTC).isoformat(),
                    count,
                    remote_signal or "",
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()
        os.replace(tmp_path, db_path)
        return count
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                logger.warning("Could not remove failed temp DB %s", tmp_path)
        raise


def get_database_info(db_path: Path, name: str) -> DatabaseInfo | None:
    """Return version metadata for a cached database, or None if not present.

    Tolerates older caches that lack ``remote_signal`` or
    ``local_version_tag`` columns by falling back to progressively
    simpler SELECTs.  Missing columns report as None; the next
    ``db update`` self-heals via the annotator's migration path.
    """
    if not db_path.exists():
        return None
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            remote_signal: str | None = None
            local_version_tag: str | None = None
            try:
                row = conn.execute(
                    "SELECT source_url, version, downloaded_at, record_count, "
                    "remote_signal, local_version_tag "
                    "FROM database_versions WHERE name = ?",
                    (name,),
                ).fetchone()
                if row is None:
                    return None
                (
                    source_url,
                    version,
                    downloaded_at,
                    record_count,
                    remote_signal,
                    local_version_tag,
                ) = row
            except sqlite3.OperationalError:
                try:
                    row = conn.execute(
                        "SELECT source_url, version, downloaded_at, record_count, "
                        "remote_signal FROM database_versions WHERE name = ?",
                        (name,),
                    ).fetchone()
                except sqlite3.OperationalError:
                    try:
                        row = conn.execute(
                            "SELECT source_url, version, downloaded_at, record_count "
                            "FROM database_versions WHERE name = ?",
                            (name,),
                        ).fetchone()
                    except sqlite3.DatabaseError:
                        return None
                    if row is None:
                        return None
                    source_url, version, downloaded_at, record_count = row
                except sqlite3.DatabaseError:
                    return None
                else:
                    if row is None:
                        return None
                    source_url, version, downloaded_at, record_count, remote_signal = row
                _ensure_local_version_tag_column(conn)
            except sqlite3.DatabaseError:
                return None
            return DatabaseInfo(
                source_url=source_url,
                version=version,
                downloaded_at=downloaded_at,
                record_count=record_count,
                remote_signal=remote_signal,
                local_version_tag=local_version_tag,
            )
    except sqlite3.DatabaseError:
        return None


def _ensure_local_version_tag_column(conn: sqlite3.Connection) -> None:
    """Add ``local_version_tag`` column if absent (idempotent soft migration)."""
    with contextlib.suppress(sqlite3.OperationalError):
        conn.execute("ALTER TABLE database_versions ADD COLUMN local_version_tag TEXT")


def stamp_remote_signal(
    conn: sqlite3.Connection,
    name: str,
    remote_signal: str,
    source_url: str = "",
) -> None:
    """Ensure ``database_versions`` exists and upsert the remote signal.

    Existing rows keep their version / downloaded_at / record_count
    metadata; only ``remote_signal`` is overwritten.  If the table or row
    is missing (pre-built caches shipped without version metadata), both
    are created with placeholder values that ``db status`` degrades
    gracefully on.
    """
    conn.execute(
        "CREATE TABLE IF NOT EXISTS database_versions ("
        "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, version TEXT, "
        "downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
        "remote_signal TEXT, local_version_tag TEXT)"
    )
    _ensure_local_version_tag_column(conn)
    conn.execute(
        "INSERT INTO database_versions"
        " (name, source_url, version, downloaded_at, record_count, remote_signal)"
        " VALUES (?, ?, NULL, datetime('now'), 0, ?)"
        " ON CONFLICT(name) DO UPDATE SET remote_signal = excluded.remote_signal",
        (name, source_url, remote_signal),
    )


def stamp_existing_clinvar_cache(db_path: Path) -> bool:
    """One-shot migration: stamp ``local_version_tag`` on a ClinVar cache.

    Handles two legacy states:

    1. Pre-mechanism caches (no interpreter stamp at all) — writes the tag.
    2. Old-format caches with ``|iv:N`` baked into ``remote_signal`` —
       moves the tag to ``local_version_tag`` and cleans ``remote_signal``.

    Only stamps when the tag is absent (NULL).  A stale tag (wrong
    version number) means the interpreter changed and the cache needs
    re-downloading — that is NOT self-healable.

    Returns True if the current interpreter version is now stamped.
    Called from ``ClinVarAnnotator.is_ready()`` so existing caches
    self-heal on first run without re-downloading 400 MB.
    """
    if not db_path.exists():
        return False
    import contextlib

    from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

    tag = f"iv:{CLINVAR_INTERPRETER_VERSION}"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        _ensure_local_version_tag_column(conn)
        try:
            rows = conn.execute(
                "SELECT name, remote_signal, local_version_tag "
                "FROM database_versions WHERE name LIKE 'clinvar%'"
            ).fetchall()
        except sqlite3.OperationalError:
            return False
        if not rows:
            return False
        stamped = False
        for name, sig, existing_tag in rows:
            if existing_tag == tag:
                continue
            if existing_tag is not None:
                return False
            clean_signal = (sig or "").split("|iv:")[0]
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = ? WHERE name = ?",
                (clean_signal, tag, name),
            )
            stamped = True
        if stamped:
            conn.commit()
        return True

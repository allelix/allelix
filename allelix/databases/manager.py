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
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

from allelix import __version__
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

# GH #42: ADR-0021 per-build dispatch is preserved at the SQLite cache
# layer (clinvar.GRCh37.sqlite + clinvar.GRCh38.sqlite), but the source
# is a single TSV pair that carries rows for every supported Assembly.
# variant_summary.txt.gz holds per-(VariationID, Assembly) rows with
# position / ref / alt / rsID / gene; submission_summary.txt.gz holds
# per-SCV (VariationID, ClinicalSignificance, ReportedPhenotypeInfo,
# ReviewStatus, SCV ID) rows that join on VariationID. The loader
# emits one cache record per (variant, SCV) — the per-SCV pairing
# that the prior VCF loader's CLNSIG|CLNDN parse could not.
CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
CLINVAR_SUBMISSION_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
)
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


# --- GH #42: per-SCV TSV loader ------------------------------------
#
# Column indices (zero-based) inside the TSV files, taken from the
# verified column-header rows fetched 2026-06-17 from
#   https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
# Both files use `\t` as field separator and `#` as comment-line prefix.

# variant_summary.txt.gz — 43 columns. We only need a handful.
_VS_RS = 9  # "RS# (dbSNP)"     — number only, no `rs` prefix
_VS_GENE_SYMBOL = 4  # "GeneSymbol"
_VS_REVIEW_STATUS = 24  # "ReviewStatus" — aggregate
_VS_ASSEMBLY = 16  # "Assembly"        — "GRCh37" / "GRCh38" / "NCBI36"
_VS_CHROMOSOME = 18  # "Chromosome"
_VS_VARIATION_ID = 30  # "VariationID"   — join key against submission_summary
_VS_POSITION_VCF = 31  # "PositionVCF"
_VS_REF_ALLELE_VCF = 32  # "ReferenceAlleleVCF"
_VS_ALT_ALLELE_VCF = 33  # "AlternateAlleleVCF"

# submission_summary.txt.gz — 16 columns.
_SS_VARIATION_ID = 0  # "VariationID"   — join key
_SS_CLIN_SIG = 1  # "ClinicalSignificance"  — per-SCV (single token)
_SS_REPORTED_PHENO = 5  # "ReportedPhenotypeInfo" — "MedGenID:Name" tuples
_SS_REVIEW_STATUS = 6  # "ReviewStatus" — per-SCV
_SS_SCV = 10  # "SCV"                   — submission accession
_SS_CONTRIBUTES = 15  # "ContributesToAggregateClassification"

# GH #42 follow-up (evaluator defect 5, PR #101): exact set of
# ClinicalSignificance values the TSV loader filters at ingest.
# These are values that ClinVar uses to mark SCVs with no actionable
# classification — placeholders ("-", "", "not specified") AND
# non-classification curatorial categories ("other", "association",
# "association not found") — AND that aren't in
# `allelix.annotators.clinvar._CLNSIG_MAGNITUDE` (so without this
# filter they'd land at the 5.0 default — equal to the analyze
# display floor — and surface to users as bogus annotations).
#
# Values intentionally NOT in this set:
#   - "not provided" / "not_provided" — already maps to 2.0 in
#     _CLNSIG_MAGNITUDE (safe below the floor); dropping at ingest
#     would lose real submitter records.
#   - "no classification for the single variant" — same.
#
# GH #116 added the three non-classification terms ("other",
# "association", "association not found"). The per-SCV TSV
# switch (#42) surfaced these from the underlying submission_summary
# rows where the old summarized data hid them; they're meaningfully
# distinct from placeholders (a curator chose them) but equally
# unactionable as report content.
#
# The protocol's §7b significance-sentinel ship-gate scans the
# same set to verify the loader's commitment against the live cache.
# Keep both in sync.
_CLINVAR_PLACEHOLDER_CLNSIGS: frozenset[str] = frozenset(
    {
        "",
        "-",
        "not specified",
        "no classification provided",
        # GH #116: non-classification curatorial terms — kept out of
        # _CLNSIG_MAGNITUDE deliberately (they aren't classifications)
        # and out of _BENIGN_CLNSIGS (they aren't benign; #56's repute
        # work must not see them as "good").
        "other",
        "association",
        "association not found",
    }
)


def _open_tsv(path: Path) -> Iterator[list[str]]:
    """Stream a (possibly-gzipped) tab-delimited TSV, skipping comment lines."""
    opener: object = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield line.split("\t")


def _decode_reported_phenotype(field: str) -> str:
    """Strip MedGen ID prefixes from a `ReportedPhenotypeInfo` cell.

    submission_summary encodes conditions as ``MedGenID:Name`` tuples
    separated by ``|``. Sometimes ``MedGenID`` is ``na`` when no MedGen
    identifier exists. We want the human-readable name without the
    prefix; multiple conditions stay ``;``-joined to match the VCF
    loader's output convention so the downstream display path is
    unchanged.

    Examples:
        ``"C3150901:Hereditary spastic paraplegia 48"``
            → ``"Hereditary spastic paraplegia 48"``
        ``"na:not provided"``
            → ``"not provided"``
        ``"C1|C2:Cond A|C3:Cond B"``
            → ``"Cond A; Cond B"``  (lone token "C1" dropped — no Name)
    """
    if not field or field == "-":
        return ""
    parts = []
    for token in field.split("|"):
        token = token.strip()
        if not token or token == "-":
            continue
        if ":" not in token:
            # Token without "ID:Name" structure — skip rather than emit
            # a bare MedGen ID that looks like noise.
            continue
        _, _, name = token.partition(":")
        name = name.strip()
        if name and name != "-":
            parts.append(name)
    return "; ".join(parts)


def iter_clinvar_tsv_records(
    variant_summary_path: Path,
    submission_summary_path: Path,
    build: str,
    *,
    aggregate_only: bool = True,
) -> Iterator[dict[str, object]]:
    """Yield one cache record per (variant, SCV submission).

    Joins ``submission_summary.txt.gz`` (per-SCV rows) against
    ``variant_summary.txt.gz`` (per-(VariationID, Assembly) rows) for
    the requested ``build`` and emits one dict per matched submission
    keyed for the downstream cache writer.

    Why two files: ``submission_summary`` carries the true
    (ClinicalSignificance, ReportedPhenotypeInfo) pair for each SCV
    submission. ``variant_summary`` carries the position / ref / alt /
    rsID / gene per (variant, build). The legacy VCF loader's
    CLNSIG|CLNDN Frankenstein pairing (GH #42, removed in stage C) is
    structurally absent here because each row is one submission.

    Memory: the variant_summary index uses a temp on-disk SQLite, NOT
    a Python dict. The full dataset is ~50M variants x ~150 bytes per
    dict entry = ~7.5 GB resident, which would OOM small CI runners
    and is wasteful on any host. SQLite's b-tree pages and OS-level
    caching keep peak memory in the low hundreds of MB regardless of
    input size.

    Args:
        variant_summary_path: Path to ``variant_summary.txt.gz`` (or .txt).
        submission_summary_path: Path to ``submission_summary.txt.gz``.
        build: Target genome build — ``"GRCh37"`` or ``"GRCh38"``. Rows
            from other Assembly values are skipped during the
            variant_summary pass.
        aggregate_only: When True (default), skip SCV submissions whose
            ``ContributesToAggregateClassification`` is not ``"yes"``.
            ClinVar uses that flag to mark submissions that were rolled
            into the aggregate variant-level classification. Submissions
            flagged ``"no"`` (e.g. older submissions superseded by a
            newer one) would otherwise inflate per-variant row counts
            without adding signal. Set False to materialize EVERY
            historical SCV — used by forensic dives, not production
            ingest.

    Yields:
        Dicts shaped for the cache INSERT in ``load_clinvar_tsv``. The
        ``allele_id`` field is repurposed to carry the VariationID
        (per-variant integer, no per-SCV granularity).
    """
    import tempfile

    if build not in ("GRCh37", "GRCh38"):
        msg = f"unsupported build {build!r} — expected 'GRCh37' or 'GRCh38'"
        raise ValueError(msg)

    # Stream variant_summary into a temp SQLite indexed by VariationID
    # for the requested build only. Roughly halves the row count vs.
    # keeping both builds.
    with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        with contextlib.closing(sqlite3.connect(tmp_path)) as tmp_conn:
            tmp_conn.executescript(
                "CREATE TABLE vs ("
                "  variation_id INTEGER PRIMARY KEY, "
                "  rsid TEXT, "
                "  chromosome TEXT, "
                "  position INTEGER, "
                "  ref TEXT, "
                "  alt TEXT, "
                "  gene TEXT, "
                "  review_status TEXT"
                ");"
            )
            batch: list[tuple[object, ...]] = []
            for cols in _open_tsv(variant_summary_path):
                if len(cols) <= _VS_ALT_ALLELE_VCF:
                    continue
                if cols[_VS_ASSEMBLY] != build:
                    continue
                rs = cols[_VS_RS].strip()
                if not rs or rs == "-1":
                    continue
                pos_str = cols[_VS_POSITION_VCF].strip()
                if not pos_str or pos_str == "-1":
                    continue
                try:
                    position = int(pos_str)
                except ValueError:
                    continue
                ref = cols[_VS_REF_ALLELE_VCF].strip()
                alt = cols[_VS_ALT_ALLELE_VCF].strip()
                if not ref or ref in ("-", "na") or not alt or alt in ("-", "na"):
                    # Skip rows without VCF-style ref/alt — these are
                    # complex or copy-number variants the SNV-shaped
                    # cache can't represent. Logged at debug, not warn,
                    # because they're expected (typical clinvar load
                    # drops ~5% of rows here).
                    continue
                try:
                    variation_id = int(cols[_VS_VARIATION_ID].strip())
                except (ValueError, IndexError):
                    continue
                batch.append(
                    (
                        variation_id,
                        f"rs{rs}",
                        cols[_VS_CHROMOSOME].strip(),
                        position,
                        ref,
                        alt,
                        cols[_VS_GENE_SYMBOL].strip(),
                        cols[_VS_REVIEW_STATUS].strip(),
                    )
                )
                if len(batch) >= INSERT_BATCH_SIZE:
                    tmp_conn.executemany(
                        "INSERT OR IGNORE INTO vs VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    batch.clear()
            if batch:
                tmp_conn.executemany(
                    "INSERT OR IGNORE INTO vs VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
            tmp_conn.commit()

            # Second pass: stream submission_summary, look each SCV's
            # VariationID up against the temp index, emit on hit.
            cursor = tmp_conn.cursor()
            for cols in _open_tsv(submission_summary_path):
                if len(cols) <= _SS_CONTRIBUTES:
                    continue
                if aggregate_only and cols[_SS_CONTRIBUTES].strip().lower() != "yes":
                    continue
                try:
                    variation_id = int(cols[_SS_VARIATION_ID].strip())
                except (ValueError, IndexError):
                    continue
                row = cursor.execute(
                    "SELECT rsid, chromosome, position, ref, alt, gene, review_status "
                    "FROM vs WHERE variation_id = ?",
                    (variation_id,),
                ).fetchone()
                if row is None:
                    continue
                rsid, chromosome, position, ref, alt, gene, review_status_agg = row
                clinical_significance = cols[_SS_CLIN_SIG].strip()
                # GH #42 follow-up (evaluator defect 5 + cross-PR value-
                # domain review on PR #101): ClinVar's submission_summary
                # carries several placeholder values that mean "no clinical
                # significance recorded." A placeholder that ISN'T in
                # `_CLNSIG_MAGNITUDE` falls through to that dict's 5.0
                # default — equal to the analyze display floor — and
                # surfaces as a bogus annotation on real rsIDs.
                #
                # The set below is exactly the placeholders that are
                # NOT in the magnitude dict (those that ARE — "not
                # provided", "no_classification_for_the_single_variant"
                # — map to 2.0, safe below the floor; dropping them at
                # ingest would lose real submitter records).
                #
                # Case-insensitive: `_magnitude()` normalizes via
                # `_normalize_clnsig()` (`.lower().replace(" ", "_")`),
                # so a case-sensitive filter here would let
                # "Not Specified" / "NO CLASSIFICATION PROVIDED" slip
                # past the skip and still land at the 5.0 default.
                # `-` (the confirmed real defect) is unaffected — the
                # extra robustness covers the prose placeholders that
                # were added by domain reasoning, not observed
                # in-the-wild casing.
                if clinical_significance.lower() in _CLINVAR_PLACEHOLDER_CLNSIGS:
                    continue
                condition = _decode_reported_phenotype(cols[_SS_REPORTED_PHENO])
                # Per-SCV review_status from submission_summary is more
                # specific than the aggregate; prefer it when present.
                per_scv_review = cols[_SS_REVIEW_STATUS].strip()
                review_status = per_scv_review or review_status_agg
                yield {
                    "rsid": rsid,
                    "chromosome": chromosome,
                    "position": position,
                    "ref": ref,
                    "alt": alt,
                    "clinical_significance": clinical_significance,
                    "condition": condition,
                    "gene": gene,
                    "review_status": review_status,
                    "allele_id": variation_id,
                }
    finally:
        tmp_path.unlink(missing_ok=True)


def load_clinvar_tsv(
    variant_summary_path: Path,
    submission_summary_path: Path,
    db_path: Path,
    build: str,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
    aggregate_only: bool = True,
) -> None:
    """Build a ClinVar SQLite cache from the per-SCV TSV sources.

    Production ingest path for ClinVar (#42). The annotator returns
    MULTIPLE Annotation objects per variant when multiple SCVs target
    the same (chrom, pos, ref, alt); the annotator-side reconciliation
    landed in stage B.
    """
    if db_path.exists():
        db_path.unlink()
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        batch: list[tuple[object, ...]] = []
        for record in iter_clinvar_tsv_records(
            variant_summary_path,
            submission_summary_path,
            build,
            aggregate_only=aggregate_only,
        ):
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
                conn.executemany(
                    "INSERT INTO clinvar_variants "
                    "(rsid, chromosome, position, ref, alt, "
                    "clinical_significance, condition, gene, "
                    "review_status, allele_id) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
                batch.clear()
        if batch:
            conn.executemany(
                "INSERT INTO clinvar_variants "
                "(rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, "
                "review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                batch,
            )
        record_count = conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0]
        if remote_signal:
            from datetime import UTC, datetime

            from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

            stamp_remote_signal(conn, "clinvar", remote_signal, source_url)
            _ensure_local_version_tag_column(conn)
            # GH #42 follow-up (evaluator defect 3): stamp_remote_signal
            # inserts version=NULL by design — it's a freshness-only
            # stamper meant to upsert remote_signal onto existing rows
            # (e.g. the baked-in metadata of the HF .sqlite.gz caches).
            # Caches built from scratch by load_clinvar_tsv have no
            # pre-existing row, so version stays NULL and `db status`
            # shows "version: None". Stamp the build date here so the
            # cache identifies itself. The old VCF loader effectively
            # stamped from ClinVar's VCF ##fileDate header; the per-SCV
            # TSVs don't carry one, so build-date is the practical
            # equivalent. Freshness (remote_signal md5) is unchanged.
            build_date = datetime.now(UTC).strftime("%Y-%m-%d")
            conn.execute(
                "UPDATE database_versions SET version = ?, "
                "local_version_tag = ?, "
                "record_count = ? WHERE name = ?",
                (
                    build_date,
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                    record_count,
                    "clinvar",
                ),
            )
        conn.commit()


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
    """One-shot migration: self-heal a ClinVar cache up to the current iv tag.

    Self-heal is allowed ONLY when there is positive evidence the cache
    was built by the current interpreter version. Absent that evidence
    the function returns False and ``ClinVarAnnotator.is_ready()`` lets
    ``db update`` reingest from upstream. Silent promotion of an
    unknown-version cache is the failure mode this exists to prevent —
    serving pre-format-change data labeled as fresh would be worse than
    paying the redownload cost.

    Decision matrix (for each ``database_versions`` row matching
    ``clinvar%``):

    +-----------------------+---------------------------+-------------------+
    | local_version_tag     | remote_signal `|iv:N`     | Action            |
    +=======================+===========================+===================+
    | == CURRENT            | (any)                     | noop, return True |
    +-----------------------+---------------------------+-------------------+
    | != CURRENT (some N)   | (any)                     | return False      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | baked N == CURRENT        | self-heal: move   |
    |                       |                           | tag to column,    |
    |                       |                           | strip `|iv:`      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | baked N != CURRENT        | return False      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | no `|iv:` marker          | return False      |
    +-----------------------+---------------------------+-------------------+

    The bottom-right case is the v2.2 #42 stage-B safety hole: a
    pre-v2.0.1 cache (iv:1-era, NULL tag, no baked marker) skipping
    straight to v2.2 used to get stamped as iv:CURRENT unconditionally —
    serving old single-row VCF data labeled as fresh per-SCV TSV data.
    Bumping to iv:3 made that promotion newly catastrophic because the
    interpreter changed the data shape, not just the emit rules. NULL
    tag with no marker is now treated as unknown legacy and reingests.

    Returns True only when every clinvar row in the cache now carries
    the current tag.
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
                # iv:N != iv:CURRENT — interpreter changed, reingest.
                return False
            # NULL tag: self-heal only if remote_signal carries a baked
            # `|iv:N` marker matching CURRENT. No marker or mismatched
            # N → unknown / wrong-version legacy → reingest.
            baked = _parse_baked_iv(sig)
            if baked != CLINVAR_INTERPRETER_VERSION:
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


def _parse_baked_iv(remote_signal: str | None) -> int | None:
    """Extract the ``|iv:N`` suffix from a legacy ``remote_signal``, or None.

    Pre-v2.0.1 caches baked the interpreter version into the
    ``remote_signal`` column as ``"<sig>|iv:N"`` before
    ``local_version_tag`` existed. This helper recovers ``N`` so the
    caller can decide between self-heal and reingest. Returns None for
    missing or malformed markers (caller treats as unknown).
    """
    if not remote_signal or "|iv:" not in remote_signal:
        return None
    _, _, after = remote_signal.partition("|iv:")
    head = after.split("|", 1)[0].strip()
    if not head.isdigit():
        return None
    return int(head)

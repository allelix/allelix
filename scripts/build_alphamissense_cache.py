#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Build the AlphaMissense pathogenicity SQLite cache.

Reads the AlphaMissense hg38 TSV (from Zenodo, 71M missense variant
predictions) and joins against the local gnomAD SQLite cache to populate
rsIDs. Output is ~8GB SQLite. The gnomAD dependency is build-time only —
the output cache is self-contained.

IMPORTANT: The gnomAD cache MUST be built first. AlphaMissense source
data uses genomic coordinates (chrom/pos/ref/alt), not rsIDs. The gnomAD
cache provides the coordinate-to-rsID mapping. Without it, all rsID
fields in the output will be NULL and the AlphaMissense annotator (which
queries by rsID) will return zero results.

Build order:
  1. python scripts/build_gnomad_cache.py   (or: allelix db update)
  2. python scripts/build_alphamissense_cache.py

Two data sources:

  (default)     Stream the TSV directly from Zenodo over HTTPS.
                Never saves the source TSV to disk (~613 MB compressed,
                ~3.6 GB uncompressed).

  --tsv PATH    Read a pre-downloaded TSV from a local path.
                Accepts both .tsv and .tsv.gz files.

Source data:
  AlphaMissense_hg38.tsv.gz — 71M missense variant predictions
  https://zenodo.org/records/10813168

Columns in the TSV:
  #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id
  protein_variant  am_pathogenicity  am_class

Usage:
    python scripts/build_alphamissense_cache.py
    python scripts/build_alphamissense_cache.py --tsv AlphaMissense_hg38.tsv.gz
    python scripts/build_alphamissense_cache.py --output am.sqlite
    python scripts/build_alphamissense_cache.py --gnomad-db gnomad.sqlite

The gnomAD cache must already exist (built by build_gnomad_cache.py or
downloaded via `allelix db update`). Default location:
~/.local/share/allelix/gnomad.sqlite

License: AlphaMissense data is CC BY 4.0.
         Cheng et al., Science 2023 (doi:10.1126/science.adg7492)
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import io
import logging
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from allelix.databases.schema import ALPHAMISSENSE_SCHEMA

logger = logging.getLogger(__name__)

ZENODO_TSV_URL = "https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz"

_BATCH_SIZE = 50_000

_AMRow = tuple[
    str,  # chrom
    int,  # pos
    str,  # ref
    str,  # alt
    str | None,  # rsid (from gnomAD join)
    str | None,  # uniprot_id
    str | None,  # transcript_id
    str | None,  # protein_variant
    float,  # am_pathogenicity
    str,  # am_class
]


def _get_gnomad_version(gnomad_db: Path) -> str | None:
    """Read the gnomAD version string from its database_versions table."""
    with contextlib.closing(sqlite3.connect(gnomad_db)) as conn:
        row = conn.execute(
            "SELECT version FROM database_versions WHERE name = 'gnomad'"
        ).fetchone()
        return row[0] if row else None


def _load_gnomad_rsid_map(gnomad_db: Path) -> dict[tuple[str, int, str, str], str]:
    """Load chrom/pos/ref/alt → rsid mapping from gnomAD cache.

    gnomAD stores chrom without 'chr' prefix (e.g. '1', 'X').
    AlphaMissense uses 'chr1', 'chrX' — the caller normalizes.
    """
    if not gnomad_db.exists():
        raise FileNotFoundError(
            f"gnomAD cache not found at {gnomad_db}.\n"
            "The gnomAD cache provides coordinate-to-rsID mapping — without it,\n"
            "all rsIDs will be NULL and the AlphaMissense annotator won't match\n"
            "any variants.\n\n"
            "Build gnomAD first:\n"
            "  allelix db update            (downloads as gnomad.sqlite)\n"
            "  python scripts/build_gnomad_cache.py\n\n"
            "If you already have the gnomAD cache at a different path, pass:\n"
            "  --gnomad-db /path/to/your/gnomad.sqlite\n\n"
            "To build anyway without rsIDs (not useful for allelix), pass --no-gnomad."
        )

    logger.info("Loading rsID map from %s...", gnomad_db)
    t0 = time.monotonic()
    rsid_map: dict[tuple[str, int, str, str], str] = {}
    with contextlib.closing(sqlite3.connect(gnomad_db)) as conn:
        cursor = conn.execute(
            "SELECT chrom, pos, ref, alt, rsid FROM gnomad_frequencies WHERE rsid IS NOT NULL"
        )
        for chrom, pos, ref, alt, rsid in cursor:
            rsid_map[(chrom, pos, ref, alt)] = rsid

    elapsed = time.monotonic() - t0
    logger.info("Loaded %d rsID mappings in %.1fs", len(rsid_map), elapsed)
    return rsid_map


def _open_tsv(tsv_path: Path | None) -> contextlib.AbstractContextManager[io.IOBase]:
    """Open a local TSV or stream from Zenodo. Returns a text-mode context manager."""
    if tsv_path is not None:
        open_fn = gzip.open if tsv_path.name.endswith(".gz") else open
        return open_fn(tsv_path, "rt", encoding="utf-8")
    return _stream_zenodo()


@contextlib.contextmanager
def _stream_zenodo() -> contextlib.AbstractContextManager[io.IOBase]:
    """Stream AlphaMissense TSV from Zenodo over HTTPS."""
    logger.info("Streaming %s", ZENODO_TSV_URL)
    req = urllib.request.Request(ZENODO_TSV_URL)
    with urllib.request.urlopen(req, timeout=600) as resp:
        decompressor = gzip.GzipFile(fileobj=resp)
        text_stream = io.TextIOWrapper(decompressor, encoding="utf-8")
        yield text_stream


def build_cache(
    output_path: Path,
    gnomad_db: Path,
    *,
    tsv_path: Path | None = None,
    skip_gnomad: bool = False,
) -> None:
    """Build the AlphaMissense SQLite cache.

    When *tsv_path* is None (the default), streams the source TSV from
    Zenodo over HTTPS — no local download required.
    """
    if skip_gnomad:
        logger.warning("--no-gnomad: skipping rsID mapping. All rsIDs will be NULL.")
        rsid_map: dict[tuple[str, int, str, str], str] = {}
    else:
        rsid_map = _load_gnomad_rsid_map(gnomad_db)

    tmp_path = output_path.parent / f"{output_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        for stmt in ALPHAMISSENSE_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.commit()

        total_records = 0
        matched_rsids = 0
        t0 = time.monotonic()

        source_label = str(tsv_path) if tsv_path else "Zenodo (streaming)"
        print(f"Processing {source_label}...")
        with _open_tsv(tsv_path) as fh:
            batch: list[_AMRow] = []
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue

                chrom_raw = parts[0]
                chrom = chrom_raw.removeprefix("chr")
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                ref = parts[2]
                alt = parts[3]
                # parts[4] is 'genome' (always 'hg38'), skip
                uniprot_id = parts[5] or None
                transcript_id = parts[6] or None
                protein_variant = parts[7] or None
                try:
                    am_pathogenicity = float(parts[8])
                except ValueError:
                    continue
                am_class = parts[9]

                rsid = rsid_map.get((chrom, pos, ref, alt))
                if rsid is not None:
                    matched_rsids += 1

                batch.append(
                    (
                        chrom,
                        pos,
                        ref,
                        alt,
                        rsid,
                        uniprot_id,
                        transcript_id,
                        protein_variant,
                        am_pathogenicity,
                        am_class,
                    )
                )

                if len(batch) >= _BATCH_SIZE:
                    conn.executemany(
                        "INSERT OR REPLACE INTO alphamissense_scores "
                        "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
                        "protein_variant, am_pathogenicity, am_class) "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    total_records += len(batch)
                    if total_records % 1_000_000 == 0:
                        print(f"  {total_records:,} records...")
                    batch.clear()

            if batch:
                conn.executemany(
                    "INSERT OR REPLACE INTO alphamissense_scores "
                    "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
                    "protein_variant, am_pathogenicity, am_class) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
                total_records += len(batch)
            conn.commit()

        conn.execute(
            "INSERT OR REPLACE INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES (?, ?, ?, datetime('now'), ?)",
            (
                "alphamissense",
                "https://zenodo.org/records/10813168",
                "2023.2",
                total_records,
            ),
        )

        if skip_gnomad:
            gnomad_version = "no_gnomad"
        else:
            gnomad_version = _get_gnomad_version(gnomad_db) or "unknown"
        conn.execute(
            "INSERT OR REPLACE INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES (?, ?, ?, datetime('now'), ?)",
            (
                "alphamissense_gnomad_source",
                str(gnomad_db),
                gnomad_version,
                matched_rsids,
            ),
        )
        conn.commit()

        print("Compacting database...")
        conn.execute("VACUUM")

    os.replace(tmp_path, output_path)

    elapsed = time.monotonic() - t0
    rsid_pct = (matched_rsids / total_records * 100) if total_records else 0
    source = "local" if tsv_path else "streamed"
    print(
        f"\nDone. {total_records:,} records written to {output_path} "
        f"({source}, {elapsed:.0f}s)\n"
        f"rsID coverage: {matched_rsids:,} / {total_records:,} ({rsid_pct:.1f}%)"
    )


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build AlphaMissense pathogenicity cache.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=None,
        help=(
            "Path to a pre-downloaded AlphaMissense_hg38.tsv.gz. "
            "When omitted, streams the TSV directly from Zenodo over HTTPS."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output SQLite path. Default: ~/.local/share/allelix/alphamissense.sqlite",
    )
    parser.add_argument(
        "--gnomad-db",
        type=Path,
        default=None,
        help=(
            "Path to gnomAD SQLite cache for rsID joining. "
            "Default: ~/.local/share/allelix/gnomad.sqlite"
        ),
    )
    parser.add_argument(
        "--no-gnomad",
        action="store_true",
        help="Build without gnomAD rsID mapping. The output cache will have NULL rsIDs.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    data_dir = Path.home() / ".local" / "share" / "allelix"
    data_dir.mkdir(parents=True, exist_ok=True)

    output = args.output or (data_dir / "alphamissense.sqlite")
    gnomad_db = args.gnomad_db or (data_dir / "gnomad.sqlite")

    build_cache(output, gnomad_db, tsv_path=args.tsv, skip_gnomad=args.no_gnomad)


if __name__ == "__main__":
    main()

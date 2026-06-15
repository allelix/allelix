#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Build the CADD PHRED-score SQLite cache.

Reads the full CADD v1.7 SNV and indel prescored files and filters
them to positions present in gnomAD, AlphaMissense, and ClinVar
(GRCh38). SNV keys are packed into int64 for compact storage (~2.5GB
RAM for ~80M positions); indel keys use a tuple set.

Two-pass approach:
  1. Load the position set from all three databases, split into
     packed SNV keys and indel tuple keys.
  2. Stream the CADD SNV file; insert matching rows into SQLite.
  3. (Optional) Stream the CADD indel file; insert matching rows.

Input:
  - whole_genome_SNVs.tsv.gz (~81 GB compressed) — required
  - gnomad.genomes.r4.0.indel.tsv.gz (~1.2 GB) — optional
  Download from: https://cadd.gs.washington.edu/download
  File format (tab-delimited, header lines start with #):
    #Chrom  Pos  Ref  Alt  RawScore  PHRED

Output: cadd.sqlite with the cadd_scores table (see allelix/databases/schema.py).

Usage:
    python scripts/build_cadd_cache.py \
        --snv-input /data/whole_genome_SNVs.tsv.gz
    python scripts/build_cadd_cache.py \
        --snv-input /data/whole_genome_SNVs.tsv.gz \
        --indel-input /data/gnomad.genomes.r4.0.indel.tsv.gz \
        --data-dir ~/.local/share/allelix \
        --output cadd.sqlite

License: CADD data is LicenseRef-CADD (non-commercial use only).
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import logging
import os
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from allelix.databases._versions import CADD_SCHEMA_VERSION
from allelix.databases.schema import CADD_SCHEMA

logger = logging.getLogger(__name__)

_BATCH_SIZE = 100_000

_PROGRESS_INTERVAL = 5_000_000

_CHROM_MAP = {str(i): i for i in range(1, 23)}
_CHROM_MAP.update({"X": 23, "Y": 24, "MT": 25, "M": 25})

_NUC_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}


def _pack(chrom: str, pos: int, ref: str, alt: str) -> int | None:
    """Pack SNV key into int64.

    Returns None for unmapped contigs or non-ACGT alleles.

    Bit layout: chrom in bits 34+ (5 bits, max 31), pos in bits 4-33
    (30 bits, max ~1.07B), ref in bits 2-3, alt in bits 0-1.
    """
    c = _CHROM_MAP.get(chrom)
    if c is None:
        return None
    r = _NUC_MAP.get(ref)
    a = _NUC_MAP.get(alt)
    if r is None or a is None:
        return None
    if pos >= (1 << 30):
        raise ValueError(f"pos {pos} on {chrom} exceeds 30-bit budget")
    return (c << 34) | (pos << 4) | (r << 2) | a


def _normalize_chrom(chrom: str) -> str:
    """Strip chr prefix, normalize M -> MT."""
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


_SOURCE_TABLES = [
    ("gnomad.sqlite", "gnomad_frequencies", "chrom", "pos"),
    ("alphamissense.sqlite", "alphamissense_scores", "chrom", "pos"),
    ("clinvar.GRCh38.sqlite", "clinvar_variants", "chromosome", "position"),
]


def _load_position_set(
    data_dir: Path,
) -> tuple[set[int], set[tuple[str, int, str, str]]]:
    """Load all (chrom, pos, ref, alt) from gnomAD, AlphaMissense, ClinVar.

    Returns (snv_keys, indel_keys) where snv_keys is a set of packed
    int64s and indel_keys is a set of (chrom, pos, ref, alt) tuples.
    """
    snv_keys: set[int] = set()
    indel_keys: set[tuple[str, int, str, str]] = set()
    skipped_contigs = 0

    gnomad_path = data_dir / "gnomad.sqlite"
    if not gnomad_path.exists():
        print(f"Error: gnomAD database not found at {gnomad_path}", file=sys.stderr)
        sys.exit(1)

    for db_name, table, chrom_col, pos_col in _SOURCE_TABLES:
        db_path = data_dir / db_name
        if not db_path.exists():
            print(f"  Warning: {db_name} not found, skipping")
            continue
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            cursor = conn.execute(f"SELECT {chrom_col}, {pos_col}, ref, alt FROM {table}")
            while True:
                rows = cursor.fetchmany(500_000)
                if not rows:
                    break
                for chrom, pos, ref, alt in rows:
                    chrom = _normalize_chrom(str(chrom))
                    pos = int(pos)
                    ref = str(ref)
                    alt = str(alt)
                    if len(ref) == 1 and len(alt) == 1:
                        packed = _pack(chrom, pos, ref, alt)
                        if packed is None:
                            skipped_contigs += 1
                            continue
                        snv_keys.add(packed)
                    else:
                        indel_keys.add((chrom, pos, ref, alt))
        print(f"  {db_name}: done")

    print(f"  Skipped {skipped_contigs:,} unmapped contig/allele rows")
    return snv_keys, indel_keys


def _stream_cadd_file(
    cadd_path: Path,
    conn: sqlite3.Connection,
    snv_keys: set[int],
    indel_keys: set[tuple[str, int, str, str]],
    *,
    is_indel_pass: bool = False,
) -> tuple[int, int, int]:
    """Stream a CADD file and insert matching rows.

    Returns (scanned, inserted, skipped) counts.
    """
    t0 = time.monotonic()
    total_scanned = 0
    total_inserted = 0
    total_skipped = 0
    batch: list[tuple[str, int, str, str, float]] = []

    with gzip.open(cadd_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            total_scanned += 1

            fields = line.split("\t")
            if len(fields) < 6:
                continue

            chrom = _normalize_chrom(fields[0])
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ref = fields[2]
            alt = fields[3]

            if is_indel_pass:
                if (chrom, pos, ref, alt) not in indel_keys:
                    continue
            else:
                packed = _pack(chrom, pos, ref, alt)
                if packed is None:
                    total_skipped += 1
                    continue
                if packed not in snv_keys:
                    continue

            try:
                phred = float(fields[5].rstrip())
            except ValueError:
                continue

            batch.append((chrom, pos, ref, alt, phred))
            if len(batch) >= _BATCH_SIZE:
                conn.executemany(
                    "INSERT OR REPLACE INTO cadd_scores "
                    "(chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
                    batch,
                )
                total_inserted += len(batch)
                batch.clear()

            if total_scanned % _PROGRESS_INTERVAL == 0:
                elapsed = time.monotonic() - t0
                rate = total_scanned / elapsed if elapsed > 0 else 0
                pass_name = "Indel" if is_indel_pass else "SNV"
                print(
                    f"  {pass_name}: scanned {total_scanned:,} rows, "
                    f"matched {total_inserted + len(batch):,} "
                    f"({rate:,.0f} rows/s)"
                )

    if batch:
        conn.executemany(
            "INSERT OR REPLACE INTO cadd_scores "
            "(chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
            batch,
        )
        total_inserted += len(batch)
    conn.commit()

    return total_scanned, total_inserted, total_skipped


def build_cache(
    snv_path: Path,
    indel_path: Path | None,
    output_path: Path,
    data_dir: Path,
) -> None:
    """Build the CADD SQLite cache filtered to database positions."""
    print(f"Loading position set from {data_dir}...")
    t0 = time.monotonic()
    snv_keys, indel_keys = _load_position_set(data_dir)
    elapsed = time.monotonic() - t0
    print(f"  Loaded {len(snv_keys):,} SNV keys + {len(indel_keys):,} indel keys ({elapsed:.1f}s)")

    if not snv_keys and not indel_keys:
        print("Error: no positions loaded from databases.", file=sys.stderr)
        sys.exit(1)

    tmp_path = output_path.parent / f"{output_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        for stmt in CADD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.commit()

        print(f"SNV pass: streaming {snv_path}...")
        snv_scanned, snv_inserted, snv_skipped = _stream_cadd_file(
            snv_path, conn, snv_keys, indel_keys, is_indel_pass=False
        )
        print(
            f"  SNV pass: {snv_scanned:,} scanned, {snv_inserted:,} matched, "
            f"{snv_skipped:,} skipped (unmapped)"
        )

        indel_scanned = 0
        indel_inserted = 0
        indel_skipped = 0
        if indel_path is not None:
            print(f"Indel pass: streaming {indel_path}...")
            indel_scanned, indel_inserted, indel_skipped = _stream_cadd_file(
                indel_path, conn, snv_keys, indel_keys, is_indel_pass=True
            )
            print(
                f"  Indel pass: {indel_scanned:,} scanned, "
                f"{indel_inserted:,} matched, {indel_skipped:,} skipped"
            )

        total_inserted = snv_inserted + indel_inserted
        conn.execute(
            "INSERT OR REPLACE INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count, "
            "local_version_tag) "
            "VALUES (?, ?, ?, datetime('now'), ?, ?)",
            (
                "cadd",
                "https://cadd.gs.washington.edu",
                "v1.7",
                total_inserted,
                f"sv:{CADD_SCHEMA_VERSION}",
            ),
        )
        conn.commit()

        print("Compacting database...")
        conn.execute("VACUUM")

    os.replace(tmp_path, output_path)

    elapsed_total = time.monotonic() - t0
    total_scanned = snv_scanned + indel_scanned
    print(
        f"\nDone. Scanned {total_scanned:,} rows, "
        f"kept {total_inserted:,}. "
        f"Written to {output_path} ({elapsed_total:.0f}s)"
    )


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Build CADD PHRED-score cache filtered to gnomAD, "
            "AlphaMissense, and ClinVar (GRCh38) positions."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--snv-input",
        type=Path,
        required=True,
        help="Path to whole_genome_SNVs.tsv.gz (CADD v1.7, ~81 GB).",
    )
    parser.add_argument(
        "--indel-input",
        type=Path,
        default=None,
        help="Path to gnomad.genomes.r4.0.indel.tsv.gz (optional, ~1.2 GB).",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing gnomad.sqlite, alphamissense.sqlite, "
            "clinvar.GRCh38.sqlite. Default: ~/.local/share/allelix"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output SQLite path. Default: ~/.local/share/allelix/cadd.sqlite",
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

    data_dir = args.data_dir
    if data_dir is None:
        data_dir = Path.home() / ".local" / "share" / "allelix"

    output = args.output
    if output is None:
        data_dir.mkdir(parents=True, exist_ok=True)
        output = data_dir / "cadd.sqlite"

    if not args.snv_input.exists():
        print(f"Error: SNV file not found: {args.snv_input}", file=sys.stderr)
        sys.exit(1)

    if args.indel_input is not None and not args.indel_input.exists():
        print(f"Error: indel file not found: {args.indel_input}", file=sys.stderr)
        sys.exit(1)

    build_cache(args.snv_input, args.indel_input, output, data_dir)


if __name__ == "__main__":
    main()

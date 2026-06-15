#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Build the gnomAD population frequency SQLite cache.

Extracts rsID and allele frequencies from gnomAD v4.1 exome VCFs and
inserts into SQLite. Two data sources:

  --local-dir   Read pre-downloaded VCF files from a local directory.
                Expects gnomad.exomes.v4.1.sites.chr{N}.vcf.bgz files.

  (default)     Stream VCFs over HTTPS from Google Cloud Storage.
                Never saves VCFs to disk.

Two filtering modes:

  --full        Process all records from all 24 chromosomes.
                ~16M rsIDs, ~6GB output.

  (default)     Array-only: filter to rsIDs listed in a manifest file
                (--manifest). ~1M rsIDs covering common genotyping
                arrays.

Usage:
    python scripts/build_gnomad_cache.py --full
    python scripts/build_gnomad_cache.py --full --local-dir /data/gnomad/
    python scripts/build_gnomad_cache.py --manifest data/array_rsid_manifest.txt
    python scripts/build_gnomad_cache.py --full --chromosomes 22,X

The output file is written to the default data directory
(~/.local/share/allelix/gnomad.sqlite) or the path given by --output.

License: gnomAD data is ODbL v1.0. We extract only rsID + allele
frequencies (no SpliceAI or other restrictively licensed fields).
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
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from allelix.databases.schema import GNOMAD_SCHEMA

logger = logging.getLogger(__name__)

GNOMAD_BASE_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
)
GNOMAD_FILENAME = "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"

ALL_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]

_BATCH_SIZE = 50_000


def _parse_info_field(info: str) -> dict[str, str]:
    """Parse a VCF INFO field into a dict of key=value pairs."""
    result: dict[str, str] = {}
    for part in info.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = value
        else:
            result[part] = ""
    return result


def _safe_float(value: str | None) -> float | None:
    """Convert a string to float, returning None on failure."""
    if value is None or value == "" or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


_FreqRow = tuple[
    str,  # chrom
    int,  # pos
    str,  # ref
    str,  # alt
    str,  # rsid
    float | None,  # af
    float | None,  # af_popmax
    str | None,  # popmax
    float | None,  # af_afr
    float | None,  # af_amr
    float | None,  # af_asj
    float | None,  # af_eas
    float | None,  # af_fin
    float | None,  # af_nfe
    float | None,  # af_sas
]


def _iter_vcf_records(
    fileobj: io.IOBase,
    rsid_filter: set[str] | None = None,
) -> Iterator[_FreqRow]:
    """Yield frequency records from a decompressed VCF text stream."""
    for line in fileobj:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue
        rsid = fields[2]
        if rsid == "." or not rsid.startswith("rs"):
            continue
        if rsid_filter is not None and rsid not in rsid_filter:
            continue

        # Strip "chr" prefix — gnomAD VCFs use "chr1", we store "1".
        # CADD also uses "1" (matches). AlphaMissense GRCh38 uses "chr1"
        # — future joins will need to normalize both sides.
        chrom = fields[0].removeprefix("chr")
        try:
            pos = int(fields[1])
        except ValueError:
            continue
        ref = fields[3]
        alt = fields[4]
        info = _parse_info_field(fields[7])

        yield (
            chrom,
            pos,
            ref,
            alt,
            rsid,
            _safe_float(info.get("AF")),
            _safe_float(info.get("AF_grpmax")),
            info.get("grpmax") or None,
            _safe_float(info.get("AF_afr")),
            _safe_float(info.get("AF_amr")),
            _safe_float(info.get("AF_asj")),
            _safe_float(info.get("AF_eas")),
            _safe_float(info.get("AF_fin")),
            _safe_float(info.get("AF_nfe")),
            _safe_float(info.get("AF_sas")),
        )


def _stream_vcf_gz_filtered(
    url: str,
    rsid_filter: set[str] | None = None,
) -> Iterator[_FreqRow]:
    """Stream a bgzipped VCF and yield filtered frequency records."""
    logger.info("Streaming %s", url)
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=600) as resp:
        decompressor = gzip.GzipFile(fileobj=resp)
        text_stream = io.TextIOWrapper(decompressor, encoding="utf-8")
        yield from _iter_vcf_records(text_stream, rsid_filter=rsid_filter)


def _read_local_vcf_gz(
    path: Path,
    rsid_filter: set[str] | None = None,
) -> Iterator[_FreqRow]:
    """Read a local bgzipped VCF and yield filtered frequency records."""
    logger.info("Reading %s", path)
    with path.open("rb") as fh:
        decompressor = gzip.GzipFile(fileobj=fh)
        text_stream = io.TextIOWrapper(decompressor, encoding="utf-8")
        yield from _iter_vcf_records(text_stream, rsid_filter=rsid_filter)


def _load_manifest(path: Path) -> set[str]:
    """Load a newline-delimited rsID manifest file."""
    rsids: set[str] = set()
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line and line.startswith("rs"):
                rsids.add(line)
    logger.info("Loaded %d rsIDs from manifest %s", len(rsids), path)
    return rsids


def build_cache(
    output_path: Path,
    *,
    full: bool = False,
    manifest_path: Path | None = None,
    chromosomes: list[str] | None = None,
    local_dir: Path | None = None,
) -> None:
    """Build the gnomAD SQLite cache from local files or HTTPS streams."""
    chroms = chromosomes or ALL_CHROMOSOMES
    rsid_filter: set[str] | None = None

    if not full:
        if manifest_path is None:
            print(
                "Error: either --full or --manifest is required.",
                file=sys.stderr,
            )
            sys.exit(1)
        rsid_filter = _load_manifest(manifest_path)
        if not rsid_filter:
            print(
                f"Error: manifest {manifest_path} is empty.",
                file=sys.stderr,
            )
            sys.exit(1)

    tmp_path = output_path.parent / f"{output_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
        conn.execute("PRAGMA journal_mode = WAL")
        conn.execute("PRAGMA synchronous = NORMAL")
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.commit()

        total_records = 0
        t0 = time.monotonic()

        for chrom in chroms:
            chrom_t0 = time.monotonic()
            if local_dir is not None:
                vcf_path = local_dir / GNOMAD_FILENAME.format(chrom=chrom)
                print(f"Processing {vcf_path}...")
            else:
                print(f"Processing chr{chrom} (streaming)...")

            try:
                if local_dir is not None:
                    records = _read_local_vcf_gz(vcf_path, rsid_filter)
                else:
                    url = GNOMAD_BASE_URL.format(chrom=chrom)
                    records = _stream_vcf_gz_filtered(url, rsid_filter)

                batch: list[_FreqRow] = []
                chrom_count = 0
                for record in records:
                    batch.append(record)
                    if len(batch) >= _BATCH_SIZE:
                        conn.executemany(
                            "INSERT OR REPLACE INTO gnomad_frequencies "
                            "(chrom, pos, ref, alt, rsid, af, af_popmax, popmax, "
                            "af_afr, af_amr, af_asj, af_eas, af_fin, af_nfe, af_sas) "
                            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                            batch,
                        )
                        chrom_count += len(batch)
                        batch.clear()
                if batch:
                    conn.executemany(
                        "INSERT OR REPLACE INTO gnomad_frequencies "
                        "(chrom, pos, ref, alt, rsid, af, af_popmax, popmax, "
                        "af_afr, af_amr, af_asj, af_eas, af_fin, af_nfe, af_sas) "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    chrom_count += len(batch)
                conn.commit()
                total_records += chrom_count
                elapsed = time.monotonic() - chrom_t0
                print(f"  chr{chrom}: {chrom_count:,} records ({elapsed:.1f}s)")
            except (urllib.error.URLError, FileNotFoundError, OSError) as exc:
                print(f"  ERROR processing chr{chrom}: {exc}", file=sys.stderr)
                continue

        if local_dir is not None:
            source_url = str(local_dir)
        else:
            source_url = GNOMAD_BASE_URL.replace("{chrom}", "*")
        conn.execute(
            "INSERT OR REPLACE INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES (?, ?, ?, datetime('now'), ?)",
            ("gnomad", source_url, "4.1", total_records),
        )
        conn.commit()

        print("Compacting database...")
        conn.execute("VACUUM")

    os.replace(tmp_path, output_path)

    elapsed_total = time.monotonic() - t0
    mode = "--full" if full else "--manifest"
    source = "local" if local_dir else "streamed"
    print(
        f"\nDone. {total_records:,} records written to {output_path} "
        f"({mode}, {source}, {elapsed_total:.0f}s)"
    )


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build gnomAD population frequency cache.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Build full cache from all gnomAD exome VCFs (~185GB download).",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Path to a newline-delimited rsID manifest for filtered mode.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output SQLite path. Default: ~/.local/share/allelix/gnomad.sqlite",
    )
    parser.add_argument(
        "--local-dir",
        type=Path,
        default=None,
        help="Read VCFs from a local directory instead of streaming over HTTPS.",
    )
    parser.add_argument(
        "--chromosomes",
        type=str,
        default=None,
        help="Comma-separated list of chromosomes to process (e.g., '22,X').",
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

    output = args.output
    if output is None:
        data_dir = Path.home() / ".local" / "share" / "allelix"
        data_dir.mkdir(parents=True, exist_ok=True)
        output = data_dir / "gnomad.sqlite"

    chroms = None
    if args.chromosomes:
        chroms = [c.strip() for c in args.chromosomes.split(",")]
        for c in chroms:
            if c not in ALL_CHROMOSOMES:
                print(f"Error: unknown chromosome '{c}'", file=sys.stderr)
                sys.exit(1)

    build_cache(
        output,
        full=args.full,
        manifest_path=args.manifest,
        chromosomes=chroms,
        local_dir=args.local_dir,
    )


if __name__ == "__main__":
    main()

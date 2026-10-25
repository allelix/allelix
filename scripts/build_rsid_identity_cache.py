#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Build the cross-build rsID identity cache (ADR-0040) in bounded memory.

Steps, each on disk rather than in Python data structures, so the build
runs on a 1 GB machine:

1. rsID universe: every rsID held by the installed ClinVar, ClinPGx, GWAS
   Catalog and SNPedia caches, plus the bundled high-value and
   build-detection seed lists, inserted into a scratch SQLite table.
2. Filter: each dbSNP VCF is streamed once (``curl | gzip -dc`` for a URL,
   ``gzip -dc`` for a path; never saved) through an 8 MB Bloom filter of
   the universe, writing primary-assembly candidates to a small per-build
   TSV in the work directory.
3. Pair: the TSVs are loaded into SQLite, joined exactly against the
   universe (removing Bloom false positives), streamed back in rsID order
   and paired across builds, and the rows written to
   ``rsid_identity.sqlite`` (optionally gzipped, with its SHA-256).

``bash``, ``curl`` and ``gzip`` must be on PATH. The work directory is deleted at
the end unless ``--keep-work`` is given.

Usage (one command line):
    python scripts/build_rsid_identity_cache.py --data-dir ~/.local/share/allelix
        [--grch37-vcf PATH_OR_URL] [--grch38-vcf PATH_OR_URL]
        [--output rsid_identity.sqlite] [--work-dir DIR] [--gzip] [--keep-work]

Source data (public domain, NCBI dbSNP):
    https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz  (GRCh37)
    https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz  (GRCh38)
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import gzip
import hashlib
import logging
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path
from typing import TYPE_CHECKING

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from allelix.annotators.clinvar import clinvar_db_filename
from allelix.databases.rsid_identity_loader import (
    DBSNP_GRCH37_VCF_URL,
    DBSNP_GRCH38_VCF_URL,
    RsidBloom,
    add_universe,
    bloom_from_universe,
    create_work_db,
    filter_dbsnp_stream,
    grouped_placements,
    load_placements,
    pair_grouped_placements,
    write_identity_db,
)
from allelix.reports.high_value import load_high_value_snps
from allelix.utils.build_detect import KNOWN_SNP_POSITIONS

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger("build_rsid_identity_cache")

# (cache file, table) pairs whose ``rsid`` column defines the universe.
_RSID_SOURCES: tuple[tuple[str, str], ...] = (
    (clinvar_db_filename("GRCh37"), "clinvar_variants"),
    (clinvar_db_filename("GRCh38"), "clinvar_variants"),
    ("pharmgkb.sqlite", "pharmgkb_annotations"),
    ("gwas.sqlite", "gwas_associations"),
    ("snpedia.sqlite", "snpedia_genotypes"),
)


def collect_universe(work: sqlite3.Connection, data_dir: Path) -> int:
    """Fill the work DB's universe table; return its size."""
    for filename, table in _RSID_SOURCES:
        path = data_dir / filename
        if not path.exists():
            raise SystemExit(f"missing {path}; run `allelix db update --data-dir {data_dir}`")
        with contextlib.closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as src:
            add_universe(work, (r for (r,) in src.execute(f"SELECT DISTINCT rsid FROM {table}")))
        logger.info("%s: universe now %d", filename, _universe_size(work))
    add_universe(work, load_high_value_snps())
    add_universe(work, KNOWN_SNP_POSITIONS)
    return _universe_size(work)


def _universe_size(work: sqlite3.Connection) -> int:
    return int(work.execute("SELECT COUNT(*) FROM universe").fetchone()[0])


def _progress(lines: Iterator[str], label: str) -> Iterator[str]:
    start = time.monotonic()
    for n, line in enumerate(lines, 1):
        if n % 50_000_000 == 0:
            logger.info("%s: %dM lines, %.0fs", label, n // 1_000_000, time.monotonic() - start)
        yield line


def filter_source(source: str, bloom: RsidBloom, out_path: Path) -> None:
    """Stream one dbSNP VCF (URL or path) through the Bloom filter into a TSV."""
    bash = shutil.which("bash")
    if bash is None:
        raise SystemExit("bash is required on PATH to stream dbSNP; install Bash and retry.")
    if source.startswith(("http://", "https://")):
        cmd = f"curl -sSfL --retry 5 '{source}' | gzip -dc"
    else:
        cmd = f"gzip -dc '{source}'"
    label = source.rsplit("/", 1)[-1]
    with (
        subprocess.Popen(
            [bash, "-o", "pipefail", "-c", cmd],
            stdout=subprocess.PIPE,
            text=True,
            bufsize=1 << 20,
        ) as proc,
        out_path.open("w") as out,
    ):
        assert proc.stdout is not None
        read, written = filter_dbsnp_stream(_progress(proc.stdout, label), bloom, out)
    if proc.returncode != 0:
        raise SystemExit(f"streaming {source} failed with exit code {proc.returncode}")
    logger.info("%s: %d lines read, %d candidates written to %s", label, read, written, out_path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    """Build the identity cache from the installed annotator caches and dbSNP."""
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--data-dir", type=Path, required=True, help="annotator cache directory")
    ap.add_argument("--grch37-vcf", default=DBSNP_GRCH37_VCF_URL, help="dbSNP GRCh37 VCF")
    ap.add_argument("--grch38-vcf", default=DBSNP_GRCH38_VCF_URL, help="dbSNP GRCh38 VCF")
    ap.add_argument("--dbsnp-build", default="b157", help="dbSNP build label for the stamp")
    ap.add_argument("--output", type=Path, default=Path("rsid_identity.sqlite"))
    ap.add_argument("--work-dir", type=Path, default=None, help="scratch dir")
    ap.add_argument("--gzip", action="store_true", help="also write OUTPUT.gz and its SHA-256")
    ap.add_argument("--keep-work", action="store_true", help="keep the scratch directory")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    work_dir = args.work_dir or args.output.parent / "rsid_identity_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    work_db = work_dir / "work.sqlite"
    try:
        with contextlib.closing(sqlite3.connect(work_db)) as work:
            create_work_db(work)
            logger.info("rsID universe: %d", collect_universe(work, args.data_dir))
            bloom = bloom_from_universe(work)
            for build, source in (("GRCh37", args.grch37_vcf), ("GRCh38", args.grch38_vcf)):
                tsv = work_dir / f"candidates.{build}.tsv"
                if tsv.exists() and args.keep_work:
                    logger.info("%s: reusing %s", build, tsv)
                else:
                    filter_source(source, bloom, tsv)
                with tsv.open() as fh:
                    logger.info("%s: %d exact placements", build, load_placements(work, build, fh))
            del bloom

            stats: dict[str, int] = {}
            count = write_identity_db(
                pair_grouped_placements(grouped_placements(work), stats),
                args.output,
                version=f"dbSNP {args.dbsnp_build}",
                source_url=f"{args.grch37_vcf} + {args.grch38_vcf}",
                downloaded_at=_dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            )
        with contextlib.closing(sqlite3.connect(args.output)) as conn:
            conn.execute("VACUUM")
        logger.info("wrote %s: %d rows; %s", args.output, count, stats)

        if args.gzip:
            gz_path = args.output.with_name(args.output.name + ".gz")
            with (
                args.output.open("rb") as f_in,
                gzip.open(gz_path, "wb", compresslevel=9) as f_out,
            ):
                shutil.copyfileobj(f_in, f_out)
            logger.info("wrote %s sha256=%s", gz_path, _sha256(gz_path))
    finally:
        if not args.keep_work:
            shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()

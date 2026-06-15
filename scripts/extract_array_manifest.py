#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extract a union rsID manifest from genotype files for the gnomAD cache.

Reads one or more genotype files (any format allelix supports), extracts
every rsID, and writes the deduplicated union to a manifest file. The
manifest is a newline-delimited list of rsIDs used by
``build_gnomad_cache.py --manifest`` to build a filtered gnomAD cache.

Usage:
    python scripts/extract_array_manifest.py \
        test_data/edge_cases/23andme_format_from_genes_for_good_service.txt \
        test_data/edge_cases/ftdna_grch36_positions.csv \
        -o data/array_rsid_manifest.txt

The output contains only rs-prefixed identifiers (I-probes, AX-probes,
and positional IDs are excluded since gnomAD indexes by rsID).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from allelix.parsers import ParserNotFoundError, detect_parser


def extract_rsids(file_path: Path) -> set[str]:
    """Extract all rs-prefixed IDs from a genotype file."""
    try:
        parser = detect_parser(file_path)
    except ParserNotFoundError:
        print(f"WARNING: no parser detected for {file_path}, skipping", file=sys.stderr)
        return set()
    rsids: set[str] = set()
    for variant in parser.parse(file_path):
        if variant.rsid.startswith("rs"):
            rsids.add(variant.rsid)
    return rsids


def main() -> None:
    """CLI entry point."""
    ap = argparse.ArgumentParser(
        description="Extract rsID manifest from genotype files.",
    )
    ap.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="Genotype files to extract rsIDs from.",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("data/array_rsid_manifest.txt"),
        help="Output manifest path (default: data/array_rsid_manifest.txt).",
    )
    args = ap.parse_args()

    all_rsids: set[str] = set()
    for f in args.files:
        if not f.exists():
            print(f"ERROR: {f} not found", file=sys.stderr)
            sys.exit(1)
        rsids = extract_rsids(f)
        print(f"  {f.name}: {len(rsids):,} rsIDs")
        all_rsids |= rsids

    if not all_rsids:
        print("ERROR: no rsIDs extracted from any file", file=sys.stderr)
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    sorted_rsids = sorted(all_rsids)
    args.output.write_text("\n".join(sorted_rsids) + "\n")
    print(f"\nWrote {len(sorted_rsids):,} unique rsIDs to {args.output}")


if __name__ == "__main__":
    main()

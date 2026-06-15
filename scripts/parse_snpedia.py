# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
#!/usr/bin/env python3
"""Parse raw SNPedia wiki markup into structured genotype rows.

Reads from the ``pages`` table (raw archive produced by
``scrape_snpedia.py``) and writes structured fields to a
``snpedia_genotypes`` table in the same SQLite file. Also creates a
``database_versions`` row for ``allelix db status`` symmetry.

No network. No API. Reads existing local data only. Runs once.

This is also run automatically by the annotator on first use. This
script exists for manual/standalone use.

Usage:
    python scripts/parse_snpedia.py
    python scripts/parse_snpedia.py --input /path/to/snpedia.sqlite
"""

from __future__ import annotations

import argparse
import functools
import sys
from pathlib import Path

print = functools.partial(print, flush=True)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from allelix.databases import resolve_data_dir  # noqa: E402
from allelix.databases.snpedia_parser import parse_raw_pages  # noqa: E402


def main() -> None:
    """Parse raw wiki markup into structured genotype rows."""
    default_input = str(resolve_data_dir() / "snpedia.sqlite")
    parser = argparse.ArgumentParser(
        description="Parse raw SNPedia markup into structured genotype table."
    )
    parser.add_argument(
        "--input",
        default=default_input,
        help=f"Input SQLite file (default: {default_input})",
    )
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"ERROR: {args.input} not found. Run scrape_snpedia.py first.")
        sys.exit(1)

    print(f"Parsing {args.input}...")
    count = parse_raw_pages(args.input, verbose=True)
    print(f"\nStructured genotype rows: {count}")


if __name__ == "__main__":
    main()

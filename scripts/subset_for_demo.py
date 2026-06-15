# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Subset a full genotype file to a smaller demo file for the examples/ directory.

Strategy:

1. Run ``allelix analyze`` against the full source file with the requested
   databases. Collect every rsID that produced at least one annotation —
   those are the "must-include" variants. The demo report's signal comes
   from these.
2. Pad the remainder with random variants from the source file so the
   output reaches the target variant count. Padding gives the demo report
   realistic per-chromosome statistics and prevents the file from looking
   suspiciously hand-picked.
3. Emit a 23andMe-format file with a provenance header. The output parses
   identically to a full 23andMe file — same format, same delimiters,
   same column order — it simply has fewer lines.

Designed to be re-runnable. With ``--seed`` pinned, the same input plus the
same database set always produces the same subset.

Usage::

    python3 scripts/subset_for_demo.py \\
        --input PATH-TO-FULL-23ANDME-FILE \\
        --databases PATH-TO-PINNED-DB-DIR \\
        --output examples/sample_input/demo_23andme.txt \\
        --target-variants 10000 \\
        --exclude-snpedia
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Subset a 23andMe file for the examples/ demo.")
    p.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Full 23andMe-format genotype file to subset.",
    )
    p.add_argument(
        "--databases",
        required=True,
        type=Path,
        help="Allelix database cache directory (pinned versions).",
    )
    p.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output 23andMe-format file (gets committed under examples/).",
    )
    p.add_argument(
        "--target-variants",
        type=int,
        default=10_000,
        help="Approximate output variant count. Default: 10000.",
    )
    p.add_argument(
        "--exclude-snpedia",
        action="store_true",
        help="Pass --exclude-snpedia to analyze. Required for the public demo "
        "because SNPedia content is CC BY-NC-SA (non-commercial).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=20260606,
        help="Random seed for padding selection. Pinned for reproducibility.",
    )
    p.add_argument(
        "--allelix",
        default="allelix",
        help="Path to the allelix CLI binary. Default: 'allelix' (must be on PATH).",
    )
    return p.parse_args()


def _run_analyze(
    allelix: str,
    input_path: Path,
    databases: Path,
    exclude_snpedia: bool,
) -> dict:
    """Run ``allelix analyze`` and return the parsed JSON output."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        cmd = [
            allelix,
            "analyze",
            str(input_path),
            "--data-dir",
            str(databases),
            "--output",
            str(tmp_path),
        ]
        if exclude_snpedia:
            cmd.append("--exclude-snpedia")
        print(f"  running: {' '.join(cmd)}", file=sys.stderr)
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            sys.exit(f"allelix analyze failed (exit {result.returncode})")
        return json.loads(tmp_path.read_text())
    finally:
        tmp_path.unlink(missing_ok=True)


def _read_full_source(input_path: Path) -> tuple[list[str], list[tuple[str, str, str, str]]]:
    """Return (header lines, data rows). Each data row is (rsid, chrom, pos, genotype)."""
    header: list[str] = []
    rows: list[tuple[str, str, str, str]] = []
    with input_path.open("r", encoding="utf-8") as fh:
        in_header = True
        for raw in fh:
            line = raw.rstrip("\r\n")
            if in_header and (not line or line.startswith("#")):
                header.append(line)
                continue
            in_header = False
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            rows.append((parts[0], parts[1], parts[2], parts[3]))
    return header, rows


def _build_subset(
    rows: list[tuple[str, str, str, str]],
    annotated_rsids: set[str],
    target: int,
    seed: int,
) -> list[tuple[str, str, str, str]]:
    """Pick rows: every annotated rsID + random padding to reach target.

    Order in the output preserves the source-file order so per-chromosome
    statistics look natural in stats output.
    """
    annotated_rows = [r for r in rows if r[0] in annotated_rsids]
    remaining = [r for r in rows if r[0] not in annotated_rsids]
    need_padding = max(0, target - len(annotated_rows))
    rng = random.Random(seed)
    padding_indices = set(rng.sample(range(len(remaining)), min(need_padding, len(remaining))))
    pad_rows_set = {remaining[i][0] for i in padding_indices}
    keep_rsids = set(annotated_rsids) | pad_rows_set
    return [r for r in rows if r[0] in keep_rsids]


def _write_output(
    output_path: Path,
    header: list[str],
    rows: list[tuple[str, str, str, str]],
    source_count: int,
    annotated_count: int,
) -> None:
    """Write subset to a 23andMe-format file with a provenance preamble."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    provenance = [
        "# Subsetted from openSNP user1190 (CC0 public domain donation).",
        "# Generated by scripts/subset_for_demo.py.",
        f"# Source had {source_count:,} variants; this subset has {len(rows):,} "
        "(including every rsID that produced an annotation, plus random padding).",
        f"# Annotated rsIDs preserved: {annotated_count:,}. "
        "The subset file format is identical to a full 23andMe export.",
        "#",
    ]
    with output_path.open("w", encoding="utf-8") as fh:
        for line in provenance:
            fh.write(line + "\n")
        for line in header:
            fh.write(line + "\n")
        for rsid, chrom, pos, gt in rows:
            fh.write(f"{rsid}\t{chrom}\t{pos}\t{gt}\n")


def main() -> int:  # noqa: D103
    args = _parse_args()
    if not args.input.exists():
        sys.exit(f"--input does not exist: {args.input}")
    if not args.databases.is_dir():
        sys.exit(f"--databases is not a directory: {args.databases}")

    print(f"Analyzing {args.input} to find rsIDs with annotations...", file=sys.stderr)
    payload = _run_analyze(args.allelix, args.input, args.databases, args.exclude_snpedia)
    annotated_rsids = {a["rsid"] for a in payload.get("annotations", [])}
    print(f"  found {len(annotated_rsids):,} unique rsIDs with annotations", file=sys.stderr)

    print(f"Reading source {args.input}...", file=sys.stderr)
    header, rows = _read_full_source(args.input)
    print(f"  source has {len(rows):,} variants, {len(header)} header lines", file=sys.stderr)

    print(f"Building subset (target {args.target_variants:,})...", file=sys.stderr)
    subset = _build_subset(rows, annotated_rsids, args.target_variants, args.seed)
    print(f"  subset has {len(subset):,} variants", file=sys.stderr)

    _write_output(args.output, header, subset, len(rows), len(annotated_rsids))
    size = args.output.stat().st_size
    print(f"Wrote {args.output} ({size:,} bytes)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

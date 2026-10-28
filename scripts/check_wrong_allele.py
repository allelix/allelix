#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Wrong-allele check for the enriched annotations of an Allelix JSON report.

For every annotation that carries a gnomAD AF, an AlphaMissense score or a
CADD score, asks one question: is the attached value the database value for
the allele the person carries, at that locus, with no strand complement?

The check follows the product's enrichment key, never the rsID: the GRCh38
``(chrom, pos, ref, allele)`` of the annotation, where the allele is the
annotation's ``alt`` or, for an alt-less row, the single carried
non-reference allele. A GRCh37 key is translated through the rsID identity
table by locus and alleles, and only when exactly one GRCh38 key matches
that is not a row flagged palindromic whose GRCh37 and GRCh38 references differ
(``allelix/reports/_pipeline.py``, ``allelix/annotators/rsid_identity.py``).

An allele equal to the build's reference base has no variant to look up,
except where the build's reference carries a ClinVar variant allele (a
single-base REF == ALT row in a current, ``iv:5`` or later, ClinVar cache:
Factor V Leiden on GRCh37). There the pipeline keys enrichment by that
allele, and so does this check. Anywhere else such a row with values is
a mismatch.

The JSON report does not carry the locus, so the loci of an rsID come from
the identity table, ClinVar and gnomAD, and, for an rsID no cache places,
from the ID column of the analysed file (``--input``, required).

Only the caches of sources in the report's provenance are required. A value
from a source missing from the provenance is a mismatch.

Each enriched annotation is counted as:
  - allele-direct: the attached values are the database values of its key;
  - non-direct: a value belongs to another allele at the same locus;
  - mismatch: a value matches nothing;
  - unresolved: no locus or no GRCh38 key can be found.

Exit status: 0 when every enriched annotation is allele-direct and there is
at least one (or ``--allow-empty``); 1 otherwise; 2 for a usage error or a
missing cache.

Usage:
    python scripts/check_wrong_allele.py REPORT.json --input ANALYSED_FILE \
        [--data-dir DIR]
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import json
import math
import sqlite3
import sys
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterator

FIELDS = ("allele_frequency", "am_pathogenicity", "cadd_phred")
VALUE_SOURCES = ("gnomad", "alphamissense", "cadd")
CACHE_FILES = {
    "rsid_identity": "rsid_identity.sqlite",
    "clinvar": "clinvar.{build}.sqlite",
    "gnomad": "gnomad.sqlite",
    "alphamissense": "alphamissense.sqlite",
    "cadd": "cadd.sqlite",
}

Locus = tuple[str, int]
Key = tuple[str, int, str, str]
Values = tuple[float | None, float | None, float | None]


class MissingCacheError(Exception):
    """A cache that the report's provenance requires is not installed."""


def _same(a: float | None, b: float | None) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def _carried(annotation: dict[str, Any]) -> list[str]:
    """The alleles a report says the user carries.

    ``A/G`` or ``AG``; a half call's missing allele (``T/.``) and a no-call
    (``-``) are not alleles; a hemizygous call is one allele (``T``).
    """
    genotype = annotation.get("genotype_match") or ""
    parts = genotype.split("/") if "/" in genotype else list(genotype)
    return [x for x in parts if x and x not in (".", "-")]


def scan_input_ids(path: Path, wanted: set[str]) -> dict[str, set[Locus]]:
    """Map each wanted rsID to its ``(chrom, pos)`` loci in the input's ID column.

    Reads a VCF (tab-separated, eight or more columns) or a consumer array
    export (rsID, chromosome, position first; tab- or comma-separated).
    """
    found: dict[str, set[Locus]] = {}
    if not wanted:
        return found
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t" if "\t" in line else ",")
            parts = [p.strip('"') for p in parts]
            if len(parts) >= 8 and parts[1].isdigit():
                chrom, pos, ids = parts[0], parts[1], parts[2].split(";")
            elif len(parts) >= 3 and parts[2].isdigit():
                chrom, pos, ids = parts[1], parts[2], [parts[0]]
            else:
                continue
            for rsid in ids:
                if rsid in wanted:
                    found.setdefault(rsid, set()).add((chrom.removeprefix("chr"), int(pos)))
    return found


class Caches:
    """Read-only connections to the caches; a cache not installed reads as empty."""

    def __init__(self, data_dir: Path, build: str, provenance: set[str]) -> None:
        """Open every installed cache; raise when a provenance source has none."""
        self.build = build
        self._stack = contextlib.ExitStack()
        self._conns: dict[str, sqlite3.Connection] = {}
        for source, pattern in CACHE_FILES.items():
            path = data_dir / pattern.format(build=build)
            if not path.exists():
                if source in provenance:
                    self._stack.close()
                    raise MissingCacheError(f"{source} is in the provenance but {path} is missing")
                continue
            conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
            self._stack.callback(conn.close)
            self._conns[source] = conn
        self.pos, self.ref, self.alt = (
            ("pos38", "ref38", "alt38") if build == "GRCh38" else ("pos37", "ref37", "alt37")
        )

    def close(self) -> None:
        """Close every connection."""
        self._stack.close()

    def rows(self, source: str, sql: str, params: tuple[Any, ...]) -> Iterator[tuple[Any, ...]]:
        """Rows of ``sql`` on ``source``'s cache; none when it is not installed."""
        conn = self._conns.get(source)
        if conn is not None:
            yield from conn.execute(sql, params)

    def loci(self, rsid: str) -> dict[Locus, set[str]]:
        """``(chrom, pos)`` on the report's build -> the caches that place ``rsid`` there."""
        out: dict[Locus, set[str]] = {}
        sql = f"SELECT chrom, {self.pos} FROM rsid_identity WHERE rsid=?"
        for c, p in self.rows("rsid_identity", sql, (rsid,)):
            out.setdefault((c, p), set()).add("identity")
        sql = "SELECT chromosome, position FROM clinvar_variants WHERE rsid=?"
        for c, p in self.rows("clinvar", sql, (rsid,)):
            out.setdefault((c, p), set()).add("clinvar")
        if self.build == "GRCh38":
            sql = "SELECT chrom, pos FROM gnomad_frequencies WHERE rsid=?"
            for c, p in self.rows("gnomad", sql, (rsid,)):
                out.setdefault((c, p), set()).add("gnomad")
        return out

    def refs_at(self, locus: Locus) -> set[str]:
        """Reference alleles the caches give at ``locus`` on the report's build."""
        sql = f"SELECT DISTINCT {self.ref} FROM rsid_identity WHERE chrom=? AND {self.pos}=?"
        refs = {r for (r,) in self.rows("rsid_identity", sql, locus)}
        sql = "SELECT DISTINCT ref FROM clinvar_variants WHERE chromosome=? AND position=?"
        refs |= {r for (r,) in self.rows("clinvar", sql, locus)}
        if self.build == "GRCh38":
            sql = "SELECT DISTINCT ref FROM gnomad_frequencies WHERE chrom=? AND pos=?"
            refs |= {r for (r,) in self.rows("gnomad", sql, locus)}
        return refs

    def reference_carried(self, locus: Locus, allele: str) -> bool:
        """True when the build's reference at ``locus`` carries a ClinVar variant ``allele``.

        The pipeline's rule (since ClinVar cache ``iv:5``, 2.4.2): a single-base REF ==
        ALT row means the reference base is the variant allele. An older
        cache still holds normal-allele records written REF == ALT, so it
        never counts.
        """
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        sql = "SELECT local_version_tag FROM database_versions WHERE name LIKE 'clinvar%'"
        tags = {tag for (tag,) in self.rows("clinvar", sql, ())}
        if tags != {f"iv:{CLINVAR_INTERPRETER_VERSION}"} or len(allele) != 1:
            return False
        sql = (
            "SELECT 1 FROM clinvar_variants "
            "WHERE chromosome=? AND position=? AND ref=? AND alt=? LIMIT 1"
        )
        return next(self.rows("clinvar", sql, (*locus, allele, allele)), None) is not None

    def to38(self, key: Key) -> Key | None:
        """The GRCh38 key of ``key``; GRCh37 keys only when exactly one key translates."""
        if self.build == "GRCh38":
            return key
        sql = (
            "SELECT pos38, ref38, alt38, palindromic FROM rsid_identity "
            "WHERE chrom=? AND pos37=? AND ref37=? AND alt37=?"
        )
        hits = {
            (key[0], p38, r38, a38)
            for (p38, r38, a38, pal) in self.rows("rsid_identity", sql, key)
            if a38 != r38 and not (pal and key[2] != r38)
        }
        return next(iter(hits)) if len(hits) == 1 else None

    def values(self, key38: Key) -> Values:
        """Return the gnomAD AF, AlphaMissense and CADD PHRED stored for ``key38``."""
        where = "WHERE chrom=? AND pos=? AND ref=? AND alt=?"
        found: list[float | None] = []
        for source, sql in (
            ("gnomad", f"SELECT af FROM gnomad_frequencies {where}"),
            ("alphamissense", f"SELECT am_pathogenicity FROM alphamissense_scores {where}"),
            ("cadd", f"SELECT phred FROM cadd_scores {where}"),
        ):
            row = next(self.rows(source, sql, key38), None)
            found.append(row[0] if row else None)
        return (found[0], found[1], found[2])

    def other_alleles(self, key38: Key) -> set[Key]:
        """Other GRCh38 keys at the same locus that carry a value."""
        c, p, r, _a = key38
        sql = "SELECT alt FROM gnomad_frequencies WHERE chrom=? AND pos=?"
        alts = {a for (a,) in self.rows("gnomad", sql, (c, p))}
        sql = "SELECT alt FROM cadd_scores WHERE chrom=? AND pos=? AND ref=?"
        alts |= {a for (a,) in self.rows("cadd", sql, (c, p, r))}
        return {(c, p, r, a) for a in alts} - {key38}


def check_report(
    report: dict[str, Any], caches: Caches, input_path: Path
) -> tuple[Counter[str], list[tuple[str, dict[str, Any], str]]]:
    """Classify every enriched annotation of ``report``; return counts and problems."""
    provenance = {x["name"] for x in report["annotators"]}
    active = [source in provenance for source in VALUE_SOURCES]
    enriched = [a for a in report["annotations"] if any(a.get(f) is not None for f in FIELDS)]

    def agrees(attached: Values, stored: Values) -> bool:
        return all(
            _same(x, y) if on else x is None
            for x, y, on in zip(attached, stored, active, strict=True)
        )

    # An annotation without an rsID (ADR-0046) is placed by its own locus.
    placed = {a["rsid"]: caches.loci(a["rsid"]) for a in enriched if a["rsid"]}
    unplaced = {rsid for rsid, loci in placed.items() if not loci}
    for rsid, loci in scan_input_ids(input_path, unplaced).items():
        for locus in loci:
            placed[rsid].setdefault(locus, set()).add("input-id")

    status: Counter[str] = Counter()
    problems: list[tuple[str, dict[str, Any], str]] = []
    for a in enriched:
        attached: Values = (a.get(FIELDS[0]), a.get(FIELDS[1]), a.get(FIELDS[2]))
        stray = [
            s
            for s, x, on in zip(VALUE_SOURCES, attached, active, strict=True)
            if x is not None and not on
        ]
        if stray:
            status["mismatch"] += 1
            problems.append(("mismatch", a, f"value from {stray}, not in the provenance"))
            continue
        loci = (
            placed.get(a["rsid"], {})
            if a["rsid"]
            else {(a["chromosome"], a["position"]): {"report"}}
            if a.get("chromosome")
            else {}
        )
        keys: list[Key] = []
        ref_is_allele = False
        for locus in sorted(loci):
            for ref in sorted(caches.refs_at(locus)):
                if a.get("alt"):
                    allele = a["alt"]
                else:
                    nonref = {x for x in _carried(a) if x != ref}
                    if len(nonref) != 1:
                        continue
                    allele = nonref.pop()
                if allele == ref and not caches.reference_carried(locus, allele):
                    ref_is_allele = True
                    continue
                key38 = caches.to38((locus[0], locus[1], ref, allele))
                if key38 is not None:
                    keys.append(key38)
        if any(agrees(attached, caches.values(k)) for k in keys):
            status["direct"] += 1
            continue
        wrong = sorted(
            {
                other
                for k in keys
                for other in caches.other_alleles(k)
                if any(
                    x is not None and on and _same(x, y)
                    for x, y, on in zip(attached, caches.values(other), active, strict=True)
                )
            }
        )
        if wrong:
            status["non-direct"] += 1
            problems.append(("non-direct", a, f"value belongs to {wrong[:3]}"))
        elif keys:
            status["mismatch"] += 1
            problems.append(("mismatch", a, f"stored at {keys[0]}: {caches.values(keys[0])}"))
        elif ref_is_allele:
            status["mismatch"] += 1
            problems.append(("mismatch", a, "enriched although the allele is the reference"))
        else:
            status["unresolved"] += 1
            why = "no locus known for the rsID" if not loci else "no GRCh38 key at its loci"
            problems.append(("unresolved", a, why))
    status["enriched"] = len(enriched)
    return status, problems


def main(argv: list[str] | None = None) -> int:
    """Run the check; see the module docstring for the exit status."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("report", type=Path, help="Allelix JSON report")
    ap.add_argument("--input", type=Path, required=True, help="the file the report analysed")
    ap.add_argument("--data-dir", type=Path, help="cache directory (default: Allelix's)")
    ap.add_argument("--allow-empty", action="store_true", help="pass with no enriched row")
    ap.add_argument("--show", type=int, default=10, help="problems to print")
    args = ap.parse_args(argv)
    if not args.input.exists():
        print(f"ERROR: --input {args.input} not found", file=sys.stderr)
        return 2
    data_dir = args.data_dir
    if data_dir is None:
        from allelix.databases import resolve_data_dir

        data_dir = resolve_data_dir()
    name = args.report.name

    report = json.loads(args.report.read_text(encoding="utf-8"))
    build = report["input"]["build"]
    if build not in ("GRCh37", "GRCh38"):
        print(f"{name}: build {build} gets no gnomAD, AlphaMissense or CADD values")
        return 0 if args.allow_empty else 1
    try:
        caches = Caches(data_dir, build, {x["name"] for x in report["annotators"]})
    except MissingCacheError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    try:
        status, problems = check_report(report, caches, args.input)
    finally:
        caches.close()

    total, direct = status["enriched"], status["direct"]
    bad = status["non-direct"] + status["mismatch"] + status["unresolved"]
    for kind, a, why in problems[: args.show]:
        print(
            f"  {kind}: {a['source']} {a['rsid']} alt={a.get('alt')!r} "
            f"gt={a.get('genotype_match')!r} :: {why}"
        )
    print(
        f"{name} ({build}): {direct}/{total} allele-direct, {status['non-direct']} non-direct, "
        f"{status['mismatch']} mismatch, {status['unresolved']} unresolved"
    )
    if bad:
        return 1
    if direct == 0 and not args.allow_empty:
        print(f"FAIL: {name}: no enriched annotation to check", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

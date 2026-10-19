#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
r"""Independent gate on reports of a healthy reference sample.

The protocol's ground truth is Allelix's own output, so a false finding the
tool has always produced passes it (the 2.4.1 Factor V Leiden inversion was
expected by the ground-truth file until 2.4.2). This gate does not compare
with past output. It checks two facts that hold for any healthy reference
sample (the HG002 benchmark, the CC0 arrays):

  1. no ClinVar Pathogenic or Likely pathogenic row is called homozygous;
  2. no ClinVar Pathogenic or Likely pathogenic row sits on a variant whose
     ClinVar aggregate classification is not Pathogenic or Likely
     pathogenic (a missing aggregate counts as not P/LP), except a row the
     report marks ``outside_aggregate``: an outranked P/LP submission that
     ADR-0045 keeps on purpose because the aggregate is not a disease
     classification (drug response, ...). Rule 1 still applies to it.

It reads report output only (JSON fields, not wording) and imports nothing
from Allelix, so a change in Allelix cannot change what the gate means.

Run it on reports made with the default filters. An unfiltered report also
holds the P/LP rows Allelix scores down because their aggregate is not P/LP;
the default floor leaves them out, as a user sees the report.

A row is P/LP when its ``significance`` is ``clinvar_`` followed by parts,
separated by ``/``, that are each Pathogenic or Likely pathogenic (low
penetrance included). An aggregate is P/LP by the same rule on its
germline term (the part before any ``;``), as ClinVar writes it (ADR-0050).
A row is marked only when its ``outside_aggregate`` is JSON ``true``.

Exit status: 0 when every report passes; 1 when any report has a row that
breaks a rule; 2 for a usage error or a report that cannot be read (every
report is still checked and printed).

Usage:
    python scripts/check_reference_sample.py REPORT.json [REPORT.json ...] [--counts-only]

``--counts-only`` prints the counts and no row (for a private sample).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_PLP_PARTS = frozenset(
    {
        "pathogenic",
        "likely pathogenic",
        "pathogenic, low penetrance",
        "likely pathogenic, low penetrance",
    }
)


def _is_plp(term: str) -> bool:
    """True when every ``/`` part of the germline term is Pathogenic or Likely pathogenic."""
    germline = term.split(";")[0].strip().lower().replace("_", " ")
    return bool(germline) and all(part.strip() in _PLP_PARTS for part in germline.split("/"))


def plp_row(annotation: dict[str, Any]) -> bool:
    """True for a ClinVar row whose own classification is Pathogenic or Likely pathogenic."""
    significance = str(annotation.get("significance") or "")
    return (
        annotation.get("source") == "clinvar"
        and significance.startswith("clinvar_")
        and (_is_plp(significance[len("clinvar_") :]))
    )


def check(
    report: dict[str, Any],
) -> tuple[int, list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """``(P/LP rows, homozygous, non-P/LP aggregate, non-P/LP aggregate kept by ADR-0045)``.

    A row in the last list passes rule 2; it is not in the third.
    """
    rows = [a for a in report["annotations"] if plp_row(a)]
    homozygous = [a for a in rows if a.get("zygosity") == "homozygous"]
    not_plp = [a for a in rows if not _is_plp(str(a.get("clinvar_aggregate") or ""))]
    kept = [a for a in not_plp if a.get("outside_aggregate") is True]
    off_aggregate = [a for a in not_plp if a.get("outside_aggregate") is not True]
    return len(rows), homozygous, off_aggregate, kept


def _describe(a: dict[str, Any]) -> str:
    where = a.get("rsid") or f"{a.get('chromosome')}:{a.get('position')}"
    return (
        f"{where} {a.get('gene') or '-'} {a.get('significance')} zygosity={a.get('zygosity')} "
        f"aggregate={a.get('clinvar_aggregate')!r}"
    )


def main(argv: list[str] | None = None) -> int:
    """Check each report; print one verdict line per report (and its rows)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("reports", nargs="+", type=Path, help="Allelix JSON reports")
    ap.add_argument("--counts-only", action="store_true", help="print counts, never a row")
    args = ap.parse_args(argv)

    failed = unreadable = False
    for path in args.reports:
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
            n, homozygous, off, kept = check(report)
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            print(f"ERROR {path.name}: cannot read the report ({exc})")
            unreadable = True
            continue
        verdict = "FAIL" if homozygous or off else "PASS"
        failed = failed or verdict == "FAIL"
        print(
            f"{verdict} {path.name}: {n} ClinVar P/LP row(s); {len(homozygous)} homozygous; "
            f"{len(off)} on a variant whose ClinVar aggregate is not P/LP; "
            f"{len(kept)} outside the aggregate (kept by ADR-0045, passes)"
        )
        if not args.counts_only:
            for a in homozygous:
                print(f"  homozygous P/LP: {_describe(a)}")
            for a in off:
                print(f"  aggregate not P/LP: {_describe(a)}")
            for a in kept:
                print(f"  outside the aggregate (passes): {_describe(a)}")
    if unreadable:
        return 2
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

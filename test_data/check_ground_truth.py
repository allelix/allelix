#!/usr/bin/env python3
"""Ground-truth invariant checker for the §19 VCF battery.

Consumes `HG002_GROUND_TRUTH.yaml` (alongside this file) and asserts a set
of unfalsifiable invariants against an allelix `analyze --report-format
json` output. Sized to catch the failure modes that stale-pin counts
historically masked:

  - Floor invariants — written / total / unique-key counts must clear
    floors sized at ~80% of v2.2.1 observed values. A drop below floor
    is an over-filter or vocab-leak regression — NOT a license to ratchet
    down without root-cause analysis (see §19 prose).

  - Spot-check invariants — published HG002 ground truth at specific
    rsIDs (carrier identity, source-set membership, condition presence).
    Unfalsifiable: the GIAB benchmark file IS the source of truth, the
    checker just verifies allelix reflects it faithfully.

  - Vocabulary invariants — every `significance` value must be in the
    known-vocabulary union. A novel value means the loader stopped
    filtering placeholder CLNSIGs (a #116 / #101 regression).

  - Universal invariants — required-non-null fields, allowed sources,
    attribution mapping, magnitude range.

Usage:
    python3 test_data/check_ground_truth.py <json> <key>

where <key> is one of the keys under `floors:` in the YAML
(`giab_grch38_benchmark`, `giab_grch37_benchmark`, `hg00187_gatkhc_gvcf`).

Exit code 0 on full pass; 1 on any failure, with the failing assertions
printed to stderr.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import yaml


def _key(a: dict) -> tuple:
    return (
        a.get("rsid"),
        a.get("chrom"),
        a.get("pos"),
        a.get("source"),
        a.get("clinical_significance") or a.get("significance") or a.get("trait") or "",
    )


def _check_floors(anns: list[dict], floors: dict, failures: list[str]) -> None:
    n_written = len(anns)
    n_unique = len({_key(a) for a in anns})
    if (
        "written_annotations_mag_0_min" in floors
        and n_written < floors["written_annotations_mag_0_min"]
    ):
        floor = floors["written_annotations_mag_0_min"]
        failures.append(f"floor: written {n_written} < {floor} — over-filter regression?")
    if "unique_dedup_keys_min" in floors and n_unique < floors["unique_dedup_keys_min"]:
        floor = floors["unique_dedup_keys_min"]
        failures.append(f"floor: unique keys {n_unique} < {floor} — lost real classifications?")


def _check_spots(
    anns: list[dict], spec_list: list[dict], file_key: str, failures: list[str]
) -> None:
    for spec in spec_list:
        applies = spec.get("applies_to")
        if applies and file_key not in applies:
            continue
        rs = spec["rsid"]
        rows = [a for a in anns if a.get("rsid") == rs]
        if len(rows) < spec["min_rows"]:
            failures.append(f"spot {rs}: row count {len(rows)} < min_rows {spec['min_rows']}")
            continue
        sources = {a.get("source") for a in rows}
        missing = set(spec["expected_sources"]) - sources
        if missing:
            failures.append(
                f"spot {rs}: missing required sources {sorted(missing)} (have {sorted(sources)})"
            )
        patterns = [re.compile(p) for p in spec["allowed_significance_patterns"]]
        bad = [
            a.get("significance")
            for a in rows
            if not any(p.match(a.get("significance") or "") for p in patterns)
        ]
        if bad:
            failures.append(
                f"spot {rs}: significance values match no allowed pattern: {sorted(set(bad))[:3]}"
            )
        if spec["must_have_condition"] and not any(
            (a.get("condition") or "").strip() for a in rows
        ):
            failures.append(f"spot {rs}: no row has non-empty condition")


def _check_vocab(anns: list[dict], vocab: dict, failures: list[str]) -> None:
    allowed_clinvar = set(vocab["clinvar_allowed"])
    allowed_pharma = set(vocab["pharmgkb_allowed"])
    forbidden = {v for v in vocab["forbidden_significance_values"] if v is not None}
    seen_bad: set[tuple] = set()
    for a in anns:
        sig = a.get("significance")
        src = a.get("source")
        if sig in forbidden or sig is None or sig == "":
            seen_bad.add(("forbidden", sig, src))
        elif src == "clinvar" and sig not in allowed_clinvar:
            seen_bad.add(("clinvar_novel", sig, src))
        elif src == "pharmgkb" and sig not in allowed_pharma:
            seen_bad.add(("pharmgkb_novel", sig, src))
    for kind, sig, src in sorted(seen_bad)[:5]:
        failures.append(f"vocab: {kind} significance {sig!r} (source={src})")


def _check_universal(anns: list[dict], universal: dict, failures: list[str]) -> None:
    required = universal["required_non_null_fields"]
    allowed_sources = set(universal["allowed_source_values"])
    attr_map = universal["attribution_map"]
    mag_lo, mag_hi = universal["magnitude_range"]

    null_fields: dict[str, int] = {}
    bad_sources: set[str] = set()
    bad_attr: set[tuple] = set()
    bad_mag: set[tuple] = set()

    for a in anns:
        for f in required:
            if a.get(f) is None:
                null_fields[f] = null_fields.get(f, 0) + 1
        src = a.get("source")
        if src not in allowed_sources:
            bad_sources.add(src)
        expected_attr = attr_map.get(src)
        if expected_attr and a.get("attribution") != expected_attr:
            bad_attr.add((src, a.get("attribution")))
        mag = a.get("magnitude")
        if mag is None or not (mag_lo <= mag <= mag_hi):
            bad_mag.add((a.get("rsid"), mag))

    for f, n in null_fields.items():
        failures.append(f"universal: {n} row(s) have null required field {f!r}")
    for s in bad_sources:
        failures.append(f"universal: unexpected source {s!r}")
    for src, attr in bad_attr:
        failures.append(f"universal: source={src!r} but attribution={attr!r}")
    for rsid, mag in list(bad_mag)[:5]:
        failures.append(f"universal: magnitude {mag} out of range on {rsid!r}")


def main() -> int:
    """Entry point — see module docstring for usage."""
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2

    report_path = Path(sys.argv[1])
    file_key = sys.argv[2]
    truth_path = Path(__file__).resolve().parent / "HG002_GROUND_TRUTH.yaml"

    truth = yaml.safe_load(truth_path.read_text())
    if file_key not in truth["floors"]:
        print(
            f"unknown file key {file_key!r}; expected one of {list(truth['floors'])}",
            file=sys.stderr,
        )
        return 2

    report = json.loads(report_path.read_text())
    anns = report.get("annotations", [])

    failures: list[str] = []
    _check_floors(anns, truth["floors"][file_key], failures)
    _check_spots(anns, truth.get("analyze", []), file_key, failures)
    _check_vocab(anns, truth["vocabulary"], failures)
    _check_universal(anns, truth["universal"], failures)

    if failures:
        print(f"\n✗ {file_key}: {len(failures)} invariant failure(s)", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    n_unique = len({_key(a) for a in anns})
    print(f"✓ {file_key}: {len(anns)} annotations, {n_unique} unique keys, all invariants hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())

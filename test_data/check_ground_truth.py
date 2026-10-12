#!/usr/bin/env python3
"""Ground-truth invariant checker for the §19 VCF battery.

Consumes `HG002_GROUND_TRUTH.yaml` (alongside this file) and asserts a set
of unfalsifiable invariants against an allelix `analyze --report-format
json` output. Sized to catch the failure modes that stale-pin counts
historically masked:

  - Identity — the report is of the floor entry's `file`, on its `build`,
    and the build came from its `build_source` (ADR-0044 vocabulary).

  - Floor invariants — written and unique-key counts must clear floors
    sized at ~80% of v2.2.1 observed values. A drop below floor is an
    over-filter or vocab-leak regression — NOT a license to ratchet down
    without root-cause analysis (see §19 prose).

  - Spot-check invariants — published HG002 ground truth at specific
    rsIDs (carrier identity, source-set membership, condition presence).
    Unfalsifiable: the GIAB benchmark file IS the source of truth, the
    checker just verifies allelix reflects it faithfully.

  - Vocabulary invariants — every `significance` value must be in the
    known-vocabulary union. A novel value means the loader stopped
    filtering placeholder CLNSIGs (a regression).

  - Universal invariants — required-non-null fields, allowed sources,
    attribution mapping, magnitude range.

  - Report assertions (optional, per floors entry) — report-level fields
    outside the annotation list: `input.gvcf`, the VCF quality gate's
    `input.vcf_quality.low_quality_calls`, and `high_value_no_calls`
    (count and reason). Pins the GQ-floor behavior on low-coverage input.

Every key in the YAML is read: a key this checker does not know fails the
run (exit 2), so the YAML cannot define a check that is never made.

Usage:
    python3 test_data/check_ground_truth.py <json> <key>

where <key> is one of the keys under `floors:` in the YAML
(`giab_grch38_benchmark`, `giab_grch37_benchmark`, `hg00187_gatkhc_gvcf`,
`hg00187_gatkhc_gvcf_default`).

Exit code 0 on full pass; 1 on any failure, with the failing assertions
printed to stderr.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import yaml


# The keys this checker reads, per YAML block. `source_of_truth`, `sample`
# and `note` are descriptive text.
_TOP_KEYS = {"source_of_truth", "sample", "analyze", "vocabulary", "floors", "universal"}
_SPOT_KEYS = {
    "rsid",
    "applies_to",
    "min_rows",
    "expected_sources",
    "allowed_significance_patterns",
    "must_have_condition",
    "note",
}
_FLOOR_KEYS = {
    "file",
    "build",
    "build_source",
    "written_annotations_mag_0_min",
    "unique_dedup_keys_min",
    "report_assertions",
    "note",
}
_FLOOR_REQUIRED = _FLOOR_KEYS - {"report_assertions", "note"}
_ASSERTION_KEYS = {
    "input_gvcf",
    "vcf_quality_low_quality_calls_min",
    "high_value_no_calls_min",
    "high_value_no_calls_reason",
}
_VOCAB_KEYS = {"clinvar_allowed", "pharmgkb_allowed", "forbidden_significance_values"}
_UNIVERSAL_KEYS = {
    "required_non_null_fields",
    "allowed_source_values",
    "attribution_map",
    "magnitude_range",
}


def _schema_errors(truth: dict) -> list[str]:
    """Keys the checker would not read, and floor keys that are missing."""
    errors = []

    def unknown(where: str, block: dict, known: set[str]) -> None:
        errors.extend(f"{where}: unknown key {k!r}" for k in sorted(set(block) - known))

    unknown("top level", truth, _TOP_KEYS)
    for spec in truth.get("analyze", []):
        unknown(f"analyze {spec.get('rsid')}", spec, _SPOT_KEYS)
    unknown("vocabulary", truth.get("vocabulary", {}), _VOCAB_KEYS)
    unknown("universal", truth.get("universal", {}), _UNIVERSAL_KEYS)
    for name, floors in truth.get("floors", {}).items():
        unknown(f"floors {name}", floors, _FLOOR_KEYS)
        errors.extend(
            f"floors {name}: missing key {k!r}" for k in sorted(_FLOOR_REQUIRED - set(floors))
        )
        unknown(
            f"floors {name} report_assertions",
            floors.get("report_assertions") or {},
            _ASSERTION_KEYS,
        )
    return errors


def _check_identity(report: dict, floors: dict, failures: list[str]) -> None:
    input_block = report.get("input") or {}
    expected = (Path(floors["file"]).name, floors["build"], floors["build_source"])
    actual = (input_block.get("file"), input_block.get("build"), input_block.get("build_source"))
    for label, want, got in zip(("file", "build", "build_source"), expected, actual, strict=True):
        if got != want:
            failures.append(f"identity: input.{label} is {got!r}, expected {want!r}")


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


def _check_report_assertions(report: dict, floors: dict, failures: list[str]) -> None:
    assertions = floors.get("report_assertions") or {}
    if not assertions:
        return
    input_block = report.get("input") or {}
    if "input_gvcf" in assertions and input_block.get("gvcf") != assertions["input_gvcf"]:
        failures.append(
            f"report: input.gvcf is {input_block.get('gvcf')!r}, "
            f"expected {assertions['input_gvcf']!r}"
        )
    if "vcf_quality_low_quality_calls_min" in assertions:
        floor = assertions["vcf_quality_low_quality_calls_min"]
        quality = input_block.get("vcf_quality")
        if not isinstance(quality, dict):
            failures.append(f"report: input.vcf_quality is {quality!r}, expected an object")
        elif quality.get("low_quality_calls", 0) < floor:
            failures.append(
                f"report: vcf_quality.low_quality_calls {quality.get('low_quality_calls')} "
                f"< {floor} — GQ floor not applied?"
            )
    # The JSON report omits high_value_no_calls when the list is empty.
    no_calls = report.get("high_value_no_calls") or []
    if "high_value_no_calls_min" in assertions:
        floor = assertions["high_value_no_calls_min"]
        if len(no_calls) < floor:
            failures.append(f"report: high_value_no_calls has {len(no_calls)} entries < {floor}")
    if "high_value_no_calls_reason" in assertions:
        reason = assertions["high_value_no_calls_reason"]
        if not any(entry.get("reason") == reason for entry in no_calls):
            seen = sorted({entry.get("reason") for entry in no_calls}, key=str)
            failures.append(
                f"report: no high_value_no_calls entry has reason {reason!r} (saw {seen})"
            )


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
    schema = _schema_errors(truth)
    if schema:
        print(f"{truth_path.name} has keys this checker does not read:", file=sys.stderr)
        for error in schema:
            print(f"  - {error}", file=sys.stderr)
        return 2
    if file_key not in truth["floors"]:
        print(
            f"unknown file key {file_key!r}; expected one of {list(truth['floors'])}",
            file=sys.stderr,
        )
        return 2

    report = json.loads(report_path.read_text())
    anns = report.get("annotations", [])

    failures: list[str] = []
    _check_identity(report, truth["floors"][file_key], failures)
    _check_floors(anns, truth["floors"][file_key], failures)
    _check_report_assertions(report, truth["floors"][file_key], failures)
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

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report diff engine for comparing analysis runs.

Compares a current analysis run against a previous JSON report to surface
new, removed, and changed annotations. Primary use cases: regression
detection after code changes, QA after database refreshes, and user
version-to-version comparison.

Diff key: ``(source, rsid, condition, description, alt)``. This groups
annotations so that reclassifications (significance changes) appear as
"changed" rather than "removed + added", while two alleles at one variant
stay apart. Under one key, findings with the same ``genotype_match`` pair
first; one left on each side is a finding whose genotype changed, shown as
changed. Row order never matters. A previous entry without ``alt`` (an
older report) is matched on the key without it, as before.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from allelix.models import annotation_label, annotation_to_public_dict, public_dict_label

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_SUPPORTED_SCHEMA_VERSIONS = {str(version) for version in range(1, 11)}

# magnitude comparison tolerance. Magnitudes are scored on a
# 0-10 scale; representation noise from JSON round-tripping (e.g. 7.5
# vs 7.499999999999999) is well below 1e-9. Using exact `!=` previously
# flagged such noise as "changed", filling diffs with non-events.
_MAGNITUDE_TOLERANCE = 1e-9


@dataclass
class ChangedAnnotation:
    """An annotation whose significance or magnitude changed between runs.

    ``previous_magnitude`` is ``None`` when the baseline entry had no
    magnitude recorded (legacy or partial baseline). Treating absence as
    ``0.0`` would flag every such entry as "changed" against any nonzero
    current magnitude.
    """

    current: Annotation
    previous_significance: str
    previous_magnitude: float | None
    # The baseline's genotype; None when the baseline entry has none.
    previous_genotype_match: str | None = None


@dataclass
class DiffResult:
    """The result of comparing current annotations against a previous report."""

    new: list[Annotation] = field(default_factory=list)
    removed: list[dict] = field(default_factory=list)
    changed: list[ChangedAnnotation] = field(default_factory=list)
    previous_generated_at: str = ""

    @property
    def has_changes(self) -> bool:
        """True if any annotations were added, removed, or changed."""
        return bool(self.new or self.removed or self.changed)


# The variant is named by its rsID, or by its locus when it has none
# (ADR-0046), so findings without an rsID do not collide.
def _diff_key_from_annotation(a: Annotation) -> tuple[str, str, str, str]:
    return (a.source, annotation_label(a), a.condition, a.description)


def _diff_key_from_dict(d: dict) -> tuple[str, str, str, str]:
    return (d["source"], public_dict_label(d), d.get("condition", ""), d.get("description", ""))


_DiffKey = tuple[str, str, str, str] | tuple[str, str, str, str, str]


def _previous_order(d: dict) -> tuple[str, str]:
    """A fixed order for baseline entries under one key: genotype, then content."""
    return (str(d.get("genotype_match", "")), json.dumps(d, sort_keys=True, default=str))


def _current_order(a: Annotation) -> tuple[str, str]:
    """The same order for current findings, from their public form."""
    return _previous_order(annotation_to_public_dict(a))


def load_previous_report(path: Path) -> dict:
    """Load and validate a previous JSON report.

    Raises ValueError on invalid JSON or unsupported schema version.
    """
    text = path.read_text(encoding="utf-8")
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        msg = f"Cannot parse {path.name} as JSON: {exc}"
        raise ValueError(msg) from exc

    version = data.get("schema_version")
    if version not in _SUPPORTED_SCHEMA_VERSIONS:
        msg = (
            f"Cannot diff against schema version {version!r} "
            f"(expected one of {sorted(_SUPPORTED_SCHEMA_VERSIONS)}). "
            "Re-generate the baseline report with the current version of Allelix."
        )
        raise ValueError(msg)

    if "annotations" not in data:
        msg = f"{path.name} has no 'annotations' key."
        raise ValueError(msg)

    # validate per-annotation entries. `compute_diff` indexes
    # `d["source"]` / `d["rsid"]` unguarded; a baseline with the right
    # top-level shape but malformed entries (annotations as a dict,
    # entries missing required keys) would raise KeyError or TypeError
    # rather than the documented ValueError. Catch that here.
    annotations = data["annotations"]
    if not isinstance(annotations, list):
        msg = f"{path.name}: 'annotations' must be a list, got {type(annotations).__name__}."
        raise ValueError(msg)
    for i, entry in enumerate(annotations):
        if not isinstance(entry, dict):
            msg = f"{path.name}: annotation #{i} must be an object, got {type(entry).__name__}."
            raise ValueError(msg)
        for required in ("source", "rsid"):
            if required not in entry:
                msg = f"{path.name}: annotation #{i} is missing required key {required!r}."
                raise ValueError(msg)

    return data


def compute_diff(
    current: list[Annotation],
    previous_annotations: list[dict],
    previous_generated_at: str,
) -> DiffResult:
    """Compare current annotations against a previous report's annotation list."""
    # A previous entry with ``alt`` is keyed with it, and every finding under
    # a key is kept (one report can hold the same allele with two genotypes).
    # An entry without ``alt`` (older schema) is matched the old way.
    prev_groups: dict[_DiffKey, list[dict]] = defaultdict(list)
    prev_legacy: dict[_DiffKey, dict] = {}
    for p in previous_annotations:
        old_key = _diff_key_from_dict(p)
        if "alt" in p:
            prev_groups[(*old_key, p["alt"])].append(p)
        else:
            prev_legacy[old_key] = p

    curr_groups: dict[_DiffKey, list[Annotation]] = defaultdict(list)
    curr_legacy: dict[_DiffKey, Annotation] = {}
    for c in current:
        key: _DiffKey = (*_diff_key_from_annotation(c), c.alt)
        if key not in prev_groups and key[:4] in prev_legacy:
            curr_legacy[key[:4]] = c
        else:
            curr_groups[key].append(c)

    # Same genotype first; a finding left on each side is one finding whose
    # genotype changed. Both sides are sorted first, so the pairing never
    # depends on the order of either report.
    pairs: list[tuple[Annotation, dict]] = []
    new: list[Annotation] = []
    removed: list[dict] = []
    for key in [*curr_groups, *(k for k in prev_groups if k not in curr_groups)]:
        previous_left = sorted(prev_groups.get(key, []), key=_previous_order)
        current_left: list[Annotation] = []
        for c in sorted(curr_groups.get(key, []), key=_current_order):
            match = next(
                (
                    i
                    for i, p in enumerate(previous_left)
                    if p.get("genotype_match") == c.genotype_match
                ),
                None,
            )
            if match is None:
                current_left.append(c)
            else:
                pairs.append((c, previous_left.pop(match)))
        pairs.extend(zip(current_left, previous_left, strict=False))
        new.extend(current_left[len(previous_left) :])
        removed.extend(previous_left[len(current_left) :])
    pairs.extend((c, prev_legacy[key]) for key, c in curr_legacy.items())
    removed.extend(p for key, p in prev_legacy.items() if key not in curr_legacy)

    changed: list[ChangedAnnotation] = []
    for c, p in pairs:
        prev_mag = p.get("magnitude")
        significance_changed = c.significance != p.get("significance")
        # Only a report that also names the allele is compared on the
        # genotype; earlier reports wrote it in other forms.
        genotype_changed = (
            "alt" in p and "genotype_match" in p and c.genotype_match != p["genotype_match"]
        )
        # tolerance comparison against representation noise;
        # treat a missing previous magnitude as "no value to compare"
        # rather than the implicit 0.0 the old code substituted.
        magnitude_changed = prev_mag is not None and not math.isclose(
            c.magnitude, prev_mag, abs_tol=_MAGNITUDE_TOLERANCE
        )
        if significance_changed or magnitude_changed or genotype_changed:
            changed.append(
                ChangedAnnotation(
                    current=c,
                    previous_significance=p.get("significance", ""),
                    previous_magnitude=prev_mag,
                    previous_genotype_match=p.get("genotype_match"),
                )
            )

    new.sort(key=lambda a: (-a.magnitude, a.rsid))
    removed.sort(key=lambda d: (-d.get("magnitude", 0.0), d.get("rsid", "")))

    return DiffResult(
        new=new,
        removed=removed,
        changed=changed,
        previous_generated_at=previous_generated_at,
    )


def summarize_diff(diff: DiffResult) -> str:
    """Human-readable one-line summary of changes."""
    parts: list[str] = []

    if diff.new:
        counts: Counter[str] = Counter()
        for a in diff.new:
            counts[a.attribution] += 1
        breakdown = ", ".join(f"{n} {src}" for src, n in counts.most_common())
        parts.append(f"{len(diff.new)} new ({breakdown})")

    if diff.changed:
        parts.append(f"{len(diff.changed)} changed")

    if diff.removed:
        parts.append(f"{len(diff.removed)} removed")

    if not parts:
        return "No changes since previous report."

    date_str = ""
    if diff.previous_generated_at:
        date_str = diff.previous_generated_at[:10]

    summary = "; ".join(parts)
    if date_str:
        return f"Changes since {date_str}: {summary}."
    return f"Changes: {summary}."


def diff_annotation_to_dict(a: ChangedAnnotation) -> dict:
    """Serialize a ChangedAnnotation for JSON output."""
    d = annotation_to_public_dict(a.current)
    d["previous_significance"] = a.previous_significance
    d["previous_magnitude"] = a.previous_magnitude
    d["previous_genotype_match"] = a.previous_genotype_match
    return d

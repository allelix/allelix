# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report diff engine for comparing analysis runs.

Compares a current analysis run against a previous JSON report to surface
new, removed, and changed annotations. Primary use cases: regression
detection after code changes, QA after database refreshes, and user
version-to-version comparison.

Diff key: ``(source, rsid, condition)``. This groups annotations so that
reclassifications (significance changes) appear as "changed" rather than
"removed + added." ``genotype_match`` is excluded because the typical
diff workflow reruns the same genotype file against updated databases.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from allelix.models import annotation_to_public_dict

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_SUPPORTED_SCHEMA_VERSIONS = {"1", "2", "3", "4", "5"}

# GH #25: magnitude comparison tolerance. Magnitudes are scored on a
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
    current magnitude — see GH #25.
    """

    current: Annotation
    previous_significance: str
    previous_magnitude: float | None


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


def _diff_key_from_annotation(a: Annotation) -> tuple[str, str, str, str]:
    return (a.source, a.rsid, a.condition, a.description)


def _diff_key_from_dict(d: dict) -> tuple[str, str, str, str]:
    return (d["source"], d["rsid"], d.get("condition", ""), d.get("description", ""))


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

    # GH #25: validate per-annotation entries. `compute_diff` indexes
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
    prev_by_key: dict[tuple[str, str, str, str], dict] = {}
    for p in previous_annotations:
        key = _diff_key_from_dict(p)
        prev_by_key[key] = p

    curr_by_key: dict[tuple[str, str, str, str], Annotation] = {}
    for c in current:
        key = _diff_key_from_annotation(c)
        curr_by_key[key] = c

    new = [c for key, c in curr_by_key.items() if key not in prev_by_key]
    removed = [p for key, p in prev_by_key.items() if key not in curr_by_key]

    changed: list[ChangedAnnotation] = []
    for key, c in curr_by_key.items():
        if key in prev_by_key:
            p = prev_by_key[key]
            prev_mag = p.get("magnitude")
            significance_changed = c.significance != p.get("significance")
            # GH #25: tolerance comparison against representation noise;
            # treat a missing previous magnitude as "no value to compare"
            # rather than the implicit 0.0 the old code substituted.
            magnitude_changed = prev_mag is not None and not math.isclose(
                c.magnitude, prev_mag, abs_tol=_MAGNITUDE_TOLERANCE
            )
            if significance_changed or magnitude_changed:
                changed.append(
                    ChangedAnnotation(
                        current=c,
                        previous_significance=p.get("significance", ""),
                        previous_magnitude=prev_mag,
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
    return d

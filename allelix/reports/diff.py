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
from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_SUPPORTED_SCHEMA_VERSIONS = {"1", "2", "3", "4"}


@dataclass
class ChangedAnnotation:
    """An annotation whose significance or magnitude changed between runs."""

    current: Annotation
    previous_significance: str
    previous_magnitude: float


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
            if c.significance != p.get("significance") or c.magnitude != p.get("magnitude"):
                changed.append(
                    ChangedAnnotation(
                        current=c,
                        previous_significance=p.get("significance", ""),
                        previous_magnitude=p.get("magnitude", 0.0),
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
    d = {k: v for k, v in asdict(a.current).items() if k != "is_must_include"}
    d["previous_significance"] = a.previous_significance
    d["previous_magnitude"] = a.previous_magnitude
    return d

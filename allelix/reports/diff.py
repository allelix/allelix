# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report diff engine for comparing analysis runs.

Compares a current analysis run against a previous JSON report to surface
new, removed, and changed annotations. Primary use cases: regression
detection after code changes, QA after database refreshes, and user
version-to-version comparison.

A finding's identity is ``(source, variant label, condition, allele,
site)`` (2.5.2; the site 2.7.0): the source, the rsID or, without one, the
locus (ADR-0046), the condition, the allele, and the chromosome and
position the report publishes. One rsID at two sites is two findings, so a
finding that moved is removed and new, never paired across sites. Its
classification and the prose generated from it (``significance``,
``description``, ``magnitude``) are content, so a
reclassification or a reworded SNPedia, ClinPGx or GWAS text is one changed
finding, not a removal and a new one. Two alleles at one variant, or two
conditions, stay apart. Under one identity, findings pair by identical
content first, then by ``genotype_match``, then in a fixed order; one left
on each side is changed. Row order never matters. A previous entry without
``alt``, or without a position (an older report), is matched on the
identity without it, by a current finding whose exact identity no previous
entry has.

A previous finding is removed only when the current run has no finding with
its identity at all, before any display floor: a finding still there but
now scored below ``--min-magnitude`` is changed, and marked
``below_display_floor``. New findings are the ones the report shows.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from allelix.databases.manager import clinvar_condition_names
from allelix.models import annotation_label, annotation_to_public_dict, public_dict_label

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_SUPPORTED_SCHEMA_VERSIONS = {str(version) for version in range(1, 15)}

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
    # The finding is still there but below the report's display floor, so it
    # is not in the annotation list (2.5.2).
    below_display_floor: bool = False


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
def _diff_key_from_annotation(a: Annotation) -> tuple[str, str, str]:
    return (a.source, annotation_label(a), a.condition)


def _diff_key_from_dict(d: dict) -> tuple[str, str, str]:
    return (d["source"], public_dict_label(d), _baseline_condition(d))


def _baseline_condition(d: dict) -> str:
    """A baseline's condition as this version writes it.

    A report from before 2.5.2 can carry ClinVar MedGen IDs in the
    condition; they are removed the way the loader now removes them, so
    the same finding keeps one identity.
    """
    condition = d.get("condition", "")
    if d["source"] == "clinvar" and isinstance(condition, str):
        return clinvar_condition_names(condition)
    return condition


# The identity, then the allele and the site. None stands for a baseline
# entry that does not record them (an older report), or for a current
# finding whose position is not known.
_Site = tuple[str, int] | None
_DiffKey = tuple[str, str, str, str | None, _Site]


def _site(d: dict) -> _Site:
    """The chromosome and position a report entry publishes, or None."""
    position = d.get("position")
    if position is None:
        return None
    return (d.get("chromosome") or "", position)


def _content(d: dict) -> tuple[str, str, str]:
    """What a finding says, compared between runs rather than used to pair them."""
    return (
        str(d.get("significance", "")),
        str(d.get("description", "")),
        str(d.get("genotype_match", "")),
    )


def _previous_order(d: dict) -> tuple[str, str]:
    """A fixed order for baseline entries under one key: genotype, then content."""
    return (str(d.get("genotype_match", "")), json.dumps(d, sort_keys=True, default=str))


def _current_order(a: Annotation) -> tuple[str, str]:
    """The same order for current findings, from their public form."""
    return _previous_order(annotation_to_public_dict(a))


def _same_content(c: Annotation, p: dict) -> bool:
    return _content(annotation_to_public_dict(c)) == _content(p)


def _same_genotype(c: Annotation, p: dict) -> bool:
    return p.get("genotype_match") == c.genotype_match


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

    if not isinstance(data, dict):
        msg = f"{path.name} is not an Allelix JSON report (a JSON object at the top)."
        raise ValueError(msg)
    version = data.get("schema_version")
    generated_at = data.get("generated_at")
    if generated_at is not None and not isinstance(generated_at, str):
        msg = f"{path.name}: 'generated_at' must be a string, got {type(generated_at).__name__}."
        raise ValueError(msg)
    if not isinstance(version, str) or version not in _SUPPORTED_SCHEMA_VERSIONS:
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
        _check_entry_values(path, i, entry)

    return data


# The values `compute_diff` and the diff's output (the terminal's removed
# table, `summarize_diff`) read from a baseline, by type (2.6.1). A field may
# be missing (older reports wrote fewer); one that is present must have the
# type they use, so a malformed baseline stops before any work.
_STRING_FIELDS = ("source", "rsid")
# Present, these must be strings; null is never written for them.
_NON_NULL_STRING_FIELDS = ("significance",)
# Present, these are a string or null (a locus a report does not know).
_OPTIONAL_STRING_FIELDS = (
    "condition",
    "chromosome",
    "ref",
    "alt",
    "description",
    "genotype_match",
    "gene",
    "attribution",
)


def _check_entry_values(path: Path, i: int, entry: dict) -> None:
    def refuse(name: str, kind: str) -> None:
        value = entry[name]
        msg = f"{path.name}: annotation #{i} {name!r} must be {kind}, got {type(value).__name__}."
        raise ValueError(msg)

    for name in _STRING_FIELDS:
        if not isinstance(entry[name], str):
            refuse(name, "a string")
    for name in _NON_NULL_STRING_FIELDS:
        if name in entry and not isinstance(entry[name], str):
            refuse(name, "a string")
    for name in _OPTIONAL_STRING_FIELDS:
        if name in entry and entry[name] is not None and not isinstance(entry[name], str):
            refuse(name, "a string or null")
    position = entry.get("position")
    if position is not None and (isinstance(position, bool) or not isinstance(position, int)):
        refuse("position", "an integer or null")
    if "magnitude" in entry and not _is_finite_number(entry["magnitude"]):
        refuse("magnitude", "a finite number")


def _is_finite_number(value: object) -> bool:
    """A JSON number the diff can compare and print: not bool, null, NaN, or too large."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        return math.isfinite(float(value))
    except OverflowError:
        return False


def compute_diff(
    current: list[Annotation],
    previous_annotations: list[dict],
    previous_generated_at: str,
    *,
    unfiltered: list[Annotation] | None = None,
    below_floor: set[int] | None = None,
) -> DiffResult:
    """Compare current annotations against a previous report's annotation list.

    ``current`` is what the report shows. ``unfiltered`` is every current
    finding in the report's scope before the display floors; a previous
    finding is removed only when it has no counterpart there. A pair whose
    current finding is not shown is changed; it is ``below_display_floor``
    when its ``id`` is in ``below_floor`` (the display floors exclude it),
    or, without ``below_floor``, whenever it is not shown. A finding the GWAS
    rollup hides is not shown but not below the floor. Without
    ``unfiltered``, ``current`` is both.
    """
    shown = {id(a) for a in current}
    pool = current if unfiltered is None else unfiltered
    # A previous entry is keyed with its allele and site; one without either
    # (an older schema) on the identity without it. Every finding under a key
    # is kept: one report can hold the same identity with two genotypes, or
    # (several GWAS studies) two texts.
    prev_groups: dict[_DiffKey, list[dict]] = defaultdict(list)
    for p in previous_annotations:
        prev_groups[(*_diff_key_from_dict(p), p.get("alt"), _site(p))].append(p)

    # A current finding takes the exact key when the baseline has it, else
    # the key of an older entry without the allele, the site or both. The
    # site is the published (trimmed) locus, as the baseline wrote it.
    curr_groups: dict[_DiffKey, list[Annotation]] = defaultdict(list)
    for c in pool:
        identity = _diff_key_from_annotation(c)
        site = _site(annotation_to_public_dict(c))
        candidates: list[_DiffKey] = [
            (*identity, c.alt, site),
            (*identity, None, site),
            (*identity, c.alt, None),
            (*identity, None, None),
        ]
        key = next((k for k in candidates if k in prev_groups), candidates[0])
        curr_groups[key].append(c)

    # Identical content first, then the same genotype, then the rest in a
    # fixed order; one left on each side is one changed finding. Both sides
    # are sorted first, so the pairing never depends on either report's
    # order. A shown finding pairs before one below the display floor.
    pairs: list[tuple[Annotation, dict]] = []
    new: list[Annotation] = []
    removed: list[dict] = []
    for key in [*curr_groups, *(k for k in prev_groups if k not in curr_groups)]:
        previous_left = sorted(prev_groups.get(key, []), key=_previous_order)
        current_left = sorted(
            curr_groups.get(key, []), key=lambda c: (id(c) not in shown, _current_order(c))
        )
        for same in (_same_content, _same_genotype):
            unmatched: list[Annotation] = []
            for c in current_left:
                match = next((i for i, p in enumerate(previous_left) if same(c, p)), None)
                if match is None:
                    unmatched.append(c)
                else:
                    pairs.append((c, previous_left.pop(match)))
            current_left = unmatched
        pairs.extend(zip(current_left, previous_left, strict=False))
        new.extend(current_left[len(previous_left) :])
        removed.extend(previous_left[len(current_left) :])

    changed: list[ChangedAnnotation] = []
    for c, p in pairs:
        prev_mag = p.get("magnitude")
        significance_changed = c.significance != p.get("significance")
        description_changed = "description" in p and c.description != p["description"]
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
        if significance_changed or magnitude_changed or genotype_changed or description_changed:
            changed.append(
                ChangedAnnotation(
                    current=c,
                    previous_significance=p.get("significance", ""),
                    previous_magnitude=prev_mag,
                    previous_genotype_match=p.get("genotype_match"),
                    below_display_floor=(
                        id(c) in below_floor if below_floor is not None else id(c) not in shown
                    ),
                )
            )
    # A current finding without a baseline counterpart is new only when shown.
    new = [c for c in new if id(c) in shown]

    new.sort(key=lambda a: (-a.magnitude, a.rsid))
    removed.sort(key=lambda d: (-d.get("magnitude", 0.0), d.get("rsid", "")))

    return DiffResult(
        new=new,
        removed=removed,
        changed=changed,
        previous_generated_at=previous_generated_at,
    )


def diff_sites(diff: DiffResult) -> dict[str, set[tuple[str, int]]]:
    """The published sites of each rsID across the diff's new, changed and removed findings.

    The site is what the diff pairs on (``_site``): the chromosome and
    position each report publishes. A finding without one adds nothing.
    """
    sites: dict[str, set[tuple[str, int]]] = defaultdict(set)
    entries = [
        *(annotation_to_public_dict(a) for a in diff.new),
        *(annotation_to_public_dict(c.current) for c in diff.changed),
        *diff.removed,
    ]
    for d in entries:
        site = _site(d)
        rsid = d.get("rsid")
        if rsid and site is not None and site[0]:
            sites[rsid].add(site)
    return sites


def diff_site_label(d: dict, sites: dict[str, set[tuple[str, int]]]) -> str | None:
    """``chrom:pos`` of a diff entry whose rsID is at two or more sites in the diff, else None.

    Only then can two rows of one rsID be told apart by their site: a
    finding that moved shows where it left and where it arrived (2.7.0).
    """
    site = _site(d)
    rsid = d.get("rsid")
    if not rsid or site is None or not site[0] or len(sites.get(rsid, ())) < 2:
        return None
    return f"{site[0]}:{site[1]}"


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
    d["below_display_floor"] = a.below_display_floor
    return d

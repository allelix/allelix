# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the report diff engine."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from allelix.models import Annotation
from allelix.reports.diff import (
    ChangedAnnotation,
    DiffResult,
    compute_diff,
    diff_annotation_to_dict,
    load_previous_report,
    summarize_diff,
)


def _ann(
    source: str = "clinvar",
    rsid: str = "rs1801133",
    significance: str = "clinvar_pathogenic",
    magnitude: float = 9.0,
    condition: str = "MTHFR deficiency",
    gene: str = "MTHFR",
    description: str | None = None,
    **kwargs: object,
) -> Annotation:
    return Annotation(
        source=source,
        rsid=rsid,
        significance=significance,
        category="clinical",
        magnitude=magnitude,
        description=description if description is not None else f"{source}: test",
        attribution=source.title(),
        genotype_match="AG",
        condition=condition,
        gene=gene,
        **kwargs,
    )


def _ann_dict(
    source: str = "clinvar",
    rsid: str = "rs1801133",
    significance: str = "clinvar_pathogenic",
    magnitude: float = 9.0,
    condition: str = "MTHFR deficiency",
    gene: str = "MTHFR",
    description: str | None = None,
) -> dict:
    return {
        "source": source,
        "rsid": rsid,
        "significance": significance,
        "category": "clinical",
        "magnitude": magnitude,
        "description": description if description is not None else f"{source}: test",
        "attribution": source.title(),
        "genotype_match": "AG",
        "references": [],
        "condition": condition,
        "gene": gene,
        "is_must_include": False,
    }


class TestComputeDiff:
    """Tests for the core diff algorithm."""

    def test_identical_reports_no_changes(self) -> None:
        current = [_ann()]
        previous = [_ann_dict()]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert not diff.has_changes
        assert diff.new == []
        assert diff.removed == []
        assert diff.changed == []

    def test_new_annotation(self) -> None:
        current = [_ann(), _ann(rsid="rs4680", condition="COMT activity", gene="COMT")]
        previous = [_ann_dict()]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert len(diff.new) == 1
        assert diff.new[0].rsid == "rs4680"
        assert diff.removed == []
        assert diff.changed == []

    def test_removed_annotation(self) -> None:
        current = []
        previous = [_ann_dict()]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert diff.new == []
        assert len(diff.removed) == 1
        assert diff.removed[0]["rsid"] == "rs1801133"
        assert diff.changed == []

    def test_changed_significance(self) -> None:
        current = [_ann(significance="clinvar_likely_benign", magnitude=2.0)]
        previous = [_ann_dict(significance="clinvar_pathogenic", magnitude=9.0)]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert diff.new == []
        assert diff.removed == []
        assert len(diff.changed) == 1
        assert diff.changed[0].current.significance == "clinvar_likely_benign"
        assert diff.changed[0].previous_significance == "clinvar_pathogenic"
        assert diff.changed[0].previous_magnitude == 9.0

    def test_changed_magnitude_only(self) -> None:
        current = [_ann(magnitude=7.0)]
        previous = [_ann_dict(magnitude=9.0)]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert len(diff.changed) == 1
        assert diff.changed[0].current.magnitude == 7.0
        assert diff.changed[0].previous_magnitude == 9.0

    def test_magnitude_representation_noise_not_flagged(self) -> None:
        """GH #25: exact float `!=` flagged JSON round-trip noise as
        changed (7.5 vs 7.499999999999999). Use `math.isclose` so only
        meaningful magnitude differences surface."""
        current = [_ann(magnitude=7.5)]
        previous = [_ann_dict(magnitude=7.5 - 1e-12)]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert diff.changed == []

    def test_missing_previous_magnitude_not_flagged_against_zero(self) -> None:
        """GH #25: a baseline entry without a magnitude key used to be
        treated as previous_magnitude=0.0, so any current magnitude
        looked like a change. Now an absent previous magnitude (and
        unchanged significance) is NOT flagged."""
        prev = _ann_dict(magnitude=5.0)
        del prev["magnitude"]
        current = [_ann(magnitude=5.0)]
        diff = compute_diff(current, [prev], "2026-05-01T00:00:00")
        assert diff.changed == []

    def test_significance_changed_with_missing_previous_magnitude(self) -> None:
        """If significance changes, the entry is still flagged as
        changed even when previous_magnitude is absent.
        ``ChangedAnnotation.previous_magnitude`` carries ``None``."""
        prev = _ann_dict(magnitude=5.0, significance="clinvar_pathogenic")
        del prev["magnitude"]
        current = [_ann(magnitude=5.0, significance="clinvar_benign")]
        diff = compute_diff(current, [prev], "2026-05-01T00:00:00")
        assert len(diff.changed) == 1
        assert diff.changed[0].previous_magnitude is None
        assert diff.changed[0].previous_significance == "clinvar_pathogenic"

    def test_description_change_surfaces_as_removed_plus_added(self) -> None:
        """Description is part of the diff key — wording changes surface honestly."""
        current = [_ann()]
        prev = _ann_dict()
        prev["description"] = "clinvar: updated wording"
        diff = compute_diff(current, [prev], "2026-05-01T00:00:00")
        assert diff.has_changes
        assert len(diff.new) == 1
        assert len(diff.removed) == 1
        assert len(diff.changed) == 0

    def test_mixed_new_removed_changed(self) -> None:
        current = [
            _ann(significance="clinvar_likely_benign", magnitude=2.0),
            _ann(rsid="rs4680", condition="COMT", gene="COMT"),
        ]
        previous = [
            _ann_dict(significance="clinvar_pathogenic", magnitude=9.0),
            _ann_dict(rsid="rs5742904", condition="FH", gene="APOB"),
        ]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert len(diff.new) == 1
        assert len(diff.removed) == 1
        assert len(diff.changed) == 1

    def test_empty_previous_all_new(self) -> None:
        current = [_ann(), _ann(rsid="rs4680", condition="COMT", gene="COMT")]
        diff = compute_diff(current, [], "2026-05-01T00:00:00")
        assert len(diff.new) == 2
        assert diff.removed == []

    def test_empty_current_all_removed(self) -> None:
        previous = [_ann_dict(), _ann_dict(rsid="rs4680", condition="COMT")]
        diff = compute_diff([], previous, "2026-05-01T00:00:00")
        assert diff.new == []
        assert len(diff.removed) == 2

    def test_multi_condition_same_rsid(self) -> None:
        """Same rsid with different conditions should be independent entries."""
        current = [
            _ann(condition="MTHFR deficiency"),
            _ann(condition="Neural tube defects"),
        ]
        previous = [_ann_dict(condition="MTHFR deficiency")]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert len(diff.new) == 1
        assert diff.new[0].condition == "Neural tube defects"

    def test_same_rsid_condition_different_descriptions_no_collision(self) -> None:
        """Multiple annotations sharing (source, rsid, condition) but differing
        in description must not collide — each is tracked independently."""
        current = [
            _ann(condition="HIV", description="efavirenz side-effects"),
            _ann(condition="HIV", description="efavirenz PK"),
            _ann(condition="HIV", description="nevirapine PK"),
        ]
        previous = [
            _ann_dict(condition="HIV", description="efavirenz side-effects"),
            _ann_dict(condition="HIV", description="efavirenz PK"),
            _ann_dict(condition="HIV", description="nevirapine PK"),
        ]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert not diff.has_changes

    def test_one_of_colliding_keys_removed(self) -> None:
        """Removing one of several same-(source, rsid, condition) annotations
        surfaces as exactly 1 removed, not silently dropped."""
        current = [
            _ann(condition="HIV", description="efavirenz PK"),
        ]
        previous = [
            _ann_dict(condition="HIV", description="efavirenz side-effects"),
            _ann_dict(condition="HIV", description="efavirenz PK"),
        ]
        diff = compute_diff(current, previous, "2026-05-01T00:00:00")
        assert len(diff.removed) == 1
        assert diff.removed[0]["description"] == "efavirenz side-effects"
        assert len(diff.new) == 0

    def test_new_sorted_by_magnitude_desc(self) -> None:
        current = [
            _ann(rsid="rs1", magnitude=3.0, condition="low"),
            _ann(rsid="rs2", magnitude=9.0, condition="high"),
        ]
        diff = compute_diff(current, [], "2026-05-01T00:00:00")
        assert diff.new[0].rsid == "rs2"
        assert diff.new[1].rsid == "rs1"

    def test_previous_generated_at_preserved(self) -> None:
        diff = compute_diff([], [], "2026-05-01T12:00:00+00:00")
        assert diff.previous_generated_at == "2026-05-01T12:00:00+00:00"


class TestLoadPreviousReport:
    """Tests for JSON loading and validation."""

    def test_valid_report(self, tmp_path: Path) -> None:
        report = {
            "schema_version": "1",
            "generated_at": "2026-05-01T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        data = load_previous_report(path)
        assert data["schema_version"] == "1"
        assert len(data["annotations"]) == 1

    def test_invalid_json(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.json"
        path.write_text("not json at all {{{")
        with pytest.raises(ValueError, match="Cannot parse"):
            load_previous_report(path)

    def test_wrong_schema_version(self, tmp_path: Path) -> None:
        report = {"schema_version": "99", "annotations": []}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="schema version"):
            load_previous_report(path)

    def test_accepts_v5_baseline(self, tmp_path: Path) -> None:
        """ADR-0035 cluster bump: v5 baselines load without error."""
        report = {
            "schema_version": "5",
            "generated_at": "2026-06-16T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        data = load_previous_report(path)
        assert data["schema_version"] == "5"

    def test_missing_schema_version(self, tmp_path: Path) -> None:
        report = {"annotations": []}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="schema version"):
            load_previous_report(path)

    def test_missing_annotations_key(self, tmp_path: Path) -> None:
        report = {"schema_version": "1"}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="annotations"):
            load_previous_report(path)

    def test_annotations_not_a_list_raises_value_error(self, tmp_path: Path) -> None:
        """GH #25: a baseline with `annotations` as a dict (e.g. someone
        edited the JSON by hand) used to flow through to `compute_diff`
        and raise a less helpful TypeError on iteration. Now caught at
        load with a ValueError."""
        report = {"schema_version": "1", "annotations": {"oops": "dict-not-list"}}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="must be a list"):
            load_previous_report(path)

    def test_annotation_entry_missing_source_raises_value_error(self, tmp_path: Path) -> None:
        """GH #25: previously `compute_diff` indexed `d['source']`
        unguarded and raised KeyError. Now caught at load."""
        ann = _ann_dict()
        del ann["source"]
        report = {"schema_version": "1", "annotations": [ann]}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="missing required key 'source'"):
            load_previous_report(path)

    def test_annotation_entry_missing_rsid_raises_value_error(self, tmp_path: Path) -> None:
        ann = _ann_dict()
        del ann["rsid"]
        report = {"schema_version": "1", "annotations": [ann]}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="missing required key 'rsid'"):
            load_previous_report(path)

    def test_annotation_entry_not_dict_raises_value_error(self, tmp_path: Path) -> None:
        report = {"schema_version": "1", "annotations": ["not-a-dict"]}
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="must be an object"):
            load_previous_report(path)

    def test_v3_report_loads(self, tmp_path: Path) -> None:
        """Schema version 3 reports must load through diff without error."""
        report = {
            "schema_version": "3",
            "generated_at": "2026-06-10T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        data = load_previous_report(path)
        assert data["schema_version"] == "3"
        assert len(data["annotations"]) == 1

    def test_v2_vs_v3_cross_version_diff(self, tmp_path: Path) -> None:
        """Diffing a v2 report against v3 annotations works (diff ignores license_attributions)."""
        previous = {
            "schema_version": "2",
            "generated_at": "2026-05-01T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "prev.json"
        path.write_text(json.dumps(previous))
        data = load_previous_report(path)
        current = [_ann(), _ann(rsid="rs4680", condition="COMT", gene="COMT")]
        diff = compute_diff(current, data["annotations"], data["generated_at"])
        assert len(diff.new) == 1
        assert diff.new[0].rsid == "rs4680"

    def test_v4_report_loads(self, tmp_path: Path) -> None:
        """Schema version 4 reports must load through diff without error."""
        report = {
            "schema_version": "4",
            "generated_at": "2026-06-11T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "report.json"
        path.write_text(json.dumps(report))
        data = load_previous_report(path)
        assert data["schema_version"] == "4"
        assert len(data["annotations"]) == 1

    def test_v3_vs_v4_cross_version_diff(self, tmp_path: Path) -> None:
        """Diffing a v3 report against v4 annotations works."""
        previous = {
            "schema_version": "3",
            "generated_at": "2026-06-01T00:00:00",
            "annotations": [_ann_dict()],
        }
        path = tmp_path / "prev.json"
        path.write_text(json.dumps(previous))
        data = load_previous_report(path)
        current = [_ann(), _ann(rsid="rs4680", condition="COMT", gene="COMT")]
        diff = compute_diff(current, data["annotations"], data["generated_at"])
        assert len(diff.new) == 1
        assert diff.new[0].rsid == "rs4680"


class TestSummarizeDiff:
    """Tests for the human-readable summary."""

    def test_no_changes(self) -> None:
        diff = DiffResult(previous_generated_at="2026-05-01T00:00:00")
        assert summarize_diff(diff) == "No changes since previous report."

    def test_new_only(self) -> None:
        diff = DiffResult(
            new=[_ann(), _ann(rsid="rs4680", condition="COMT", gene="COMT")],
            previous_generated_at="2026-05-01T00:00:00",
        )
        result = summarize_diff(diff)
        assert "2 new" in result
        assert "2026-05-01" in result

    def test_removed_only(self) -> None:
        diff = DiffResult(
            removed=[_ann_dict()],
            previous_generated_at="2026-05-01T00:00:00",
        )
        result = summarize_diff(diff)
        assert "1 removed" in result

    def test_changed_only(self) -> None:
        diff = DiffResult(
            changed=[
                ChangedAnnotation(
                    current=_ann(), previous_significance="old", previous_magnitude=5.0
                )
            ],
            previous_generated_at="2026-05-01T00:00:00",
        )
        result = summarize_diff(diff)
        assert "1 changed" in result

    def test_mixed(self) -> None:
        diff = DiffResult(
            new=[_ann()],
            removed=[_ann_dict(rsid="rs999")],
            changed=[
                ChangedAnnotation(
                    current=_ann(rsid="rs888", condition="X"),
                    previous_significance="old",
                    previous_magnitude=1.0,
                )
            ],
            previous_generated_at="2026-06-01T00:00:00",
        )
        result = summarize_diff(diff)
        assert "1 new" in result
        assert "1 removed" in result
        assert "1 changed" in result
        assert "2026-06-01" in result

    def test_no_date(self) -> None:
        diff = DiffResult(new=[_ann()], previous_generated_at="")
        result = summarize_diff(diff)
        assert result.startswith("Changes:")


class TestDiffAnnotationToDict:
    def test_serializes_changed_annotation(self) -> None:
        changed = ChangedAnnotation(
            current=_ann(significance="clinvar_likely_benign", magnitude=2.0),
            previous_significance="clinvar_pathogenic",
            previous_magnitude=9.0,
        )
        d = diff_annotation_to_dict(changed)
        assert d["significance"] == "clinvar_likely_benign"
        assert d["previous_significance"] == "clinvar_pathogenic"
        assert d["previous_magnitude"] == 9.0
        assert "is_must_include" not in d

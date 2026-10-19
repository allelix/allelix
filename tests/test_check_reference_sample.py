# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The independent reference-sample gate: synthetic report JSON only."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from check_reference_sample import main, plp_row


def _row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "source": "clinvar",
        "rsid": "rs3999990801",
        "gene": "PLACEHOLDER",
        "significance": "clinvar_pathogenic",
        "zygosity": "heterozygous",
        "clinvar_aggregate": "Pathogenic",
    }
    row.update(overrides)
    return row


def _report(tmp_path: Path, rows: list[dict[str, object]], name: str = "r.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps({"schema_version": "9", "annotations": rows}))
    return path


@pytest.mark.parametrize(
    ("significance", "expected"),
    [
        ("clinvar_pathogenic", True),
        ("clinvar_likely_pathogenic", True),
        ("clinvar_pathogenic/likely_pathogenic", True),
        ("clinvar_pathogenic,_low_penetrance", True),
        ("clinvar_likely_pathogenic,_low_penetrance", True),
        ("clinvar_pathogenic/likely_risk_allele", False),
        ("clinvar_conflicting_classifications_of_pathogenicity", False),
        ("clinvar_risk_factor", False),
        ("clinvar_benign", False),
        ("clinvar_uncertain_significance", False),
    ],
)
def test_which_clinvar_rows_are_plp(significance: str, expected: bool) -> None:
    assert plp_row(_row(significance=significance)) is expected


def test_only_clinvar_rows_count() -> None:
    assert not plp_row(_row(source="snpedia", significance="clinvar_pathogenic"))
    assert not plp_row(_row(significance="pathogenic"))


def test_a_healthy_report_passes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rows = [
        _row(),
        _row(rsid="rs3999990802", zygosity="homozygous", significance="clinvar_benign"),
        _row(rsid="rs3999990803", clinvar_aggregate="Pathogenic/Likely pathogenic; risk factor"),
        _row(rsid="rs3999990804", zygosity="hemizygous"),
    ]
    assert main([str(_report(tmp_path, rows))]) == 0
    assert (
        "PASS r.json: 3 ClinVar P/LP row(s); 0 homozygous; 0 on a variant"
        in capsys.readouterr().out
    )


def test_a_homozygous_plp_row_fails(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rows = [_row(), _row(rsid="rs3999990805", zygosity="homozygous")]
    assert main([str(_report(tmp_path, rows))]) == 1
    out = capsys.readouterr().out
    assert "FAIL r.json: 2 ClinVar P/LP row(s); 1 homozygous; 0" in out
    assert "homozygous P/LP: rs3999990805" in out


@pytest.mark.parametrize(
    "aggregate",
    [
        "Conflicting classifications of pathogenicity",
        "Benign",
        "drug response",
        "Pathogenic/Likely risk allele",
        None,
        "",
    ],
)
def test_a_plp_row_whose_aggregate_is_not_plp_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], aggregate: str | None
) -> None:
    rows = [_row(rsid="rs3999990806", clinvar_aggregate=aggregate)]
    assert main([str(_report(tmp_path, rows))]) == 1
    out = capsys.readouterr().out
    assert "1 on a variant whose ClinVar aggregate is not P/LP" in out
    assert "aggregate not P/LP: rs3999990806" in out


def test_a_plp_row_kept_outside_the_aggregate_passes(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rows = [_row(rsid="rs3999990808", clinvar_aggregate="drug response", outside_aggregate=True)]
    assert main([str(_report(tmp_path, rows))]) == 0
    out = capsys.readouterr().out
    assert (
        "PASS r.json: 1 ClinVar P/LP row(s); 0 homozygous; 0 on a variant whose ClinVar "
        "aggregate is not P/LP; 1 outside the aggregate (kept by ADR-0045, passes)" in out
    )
    assert "outside the aggregate (passes): rs3999990808" in out


@pytest.mark.parametrize("mark", [False, None, "true", 1])
def test_the_same_row_without_the_mark_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], mark: object
) -> None:
    rows = [_row(rsid="rs3999990808", clinvar_aggregate="drug response", outside_aggregate=mark)]
    assert main([str(_report(tmp_path, rows))]) == 1
    out = capsys.readouterr().out
    assert "1 on a variant whose ClinVar aggregate is not P/LP; 0 outside the aggregate" in out
    assert "aggregate not P/LP: rs3999990808" in out


def test_a_homozygous_row_kept_outside_the_aggregate_still_fails_rule_1(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rows = [
        _row(
            rsid="rs3999990809",
            zygosity="homozygous",
            clinvar_aggregate="drug response",
            outside_aggregate=True,
        )
    ]
    assert main([str(_report(tmp_path, rows))]) == 1
    out = capsys.readouterr().out
    assert "1 homozygous; 0 on a variant whose ClinVar aggregate is not P/LP; 1 outside" in out
    assert "homozygous P/LP: rs3999990809" in out


def test_counts_only_prints_no_row(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rows = [_row(rsid="rs3999990807", zygosity="homozygous", clinvar_aggregate="Benign")]
    assert main([str(_report(tmp_path, rows)), "--counts-only"]) == 1
    out = capsys.readouterr().out
    assert "1 homozygous; 1 on a variant" in out
    assert "rs3999990807" not in out


def test_every_report_is_checked(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    good = _report(tmp_path, [_row()], "good.json")
    bad = _report(tmp_path, [_row(zygosity="homozygous")], "bad.json")
    assert main([str(good), str(bad)]) == 1
    out = capsys.readouterr().out
    assert "PASS good.json" in out
    assert "FAIL bad.json" in out


def test_an_unreadable_report_is_a_usage_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    broken = tmp_path / "broken.json"
    broken.write_text("{not json")
    assert main([str(broken)]) == 2
    assert main([str(tmp_path / "missing.json")]) == 2
    no_rows = tmp_path / "no_rows.json"
    no_rows.write_text("{}")
    assert main([str(no_rows)]) == 2
    assert "cannot read the report" in capsys.readouterr().out


def test_a_real_location_row_without_rsid_is_named_by_position(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rows = [_row(rsid="", chromosome="1", position=1000101, zygosity="homozygous")]
    assert main([str(_report(tmp_path, rows))]) == 1
    assert "homozygous P/LP: 1:1000101" in capsys.readouterr().out


def test_every_report_is_printed_even_after_an_unreadable_one(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    broken = tmp_path / "broken.json"
    broken.write_text("{not json")
    bad = _report(tmp_path, [_row(zygosity="homozygous")], "bad.json")
    assert main([str(broken), str(bad)]) == 2
    out = capsys.readouterr().out
    assert "ERROR broken.json" in out
    assert "FAIL bad.json" in out


def test_annotations_that_are_not_objects_are_unreadable(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "odd.json"
    path.write_text(json.dumps({"annotations": ["not an object"]}))
    assert main([str(path)]) == 2
    assert "ERROR odd.json: cannot read the report" in capsys.readouterr().out

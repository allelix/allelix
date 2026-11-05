# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Unplaced calls have distinct origins and do not invent a shared position."""

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _identity, _pharm

from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.config import load_config, save_config
from allelix.models import Annotation, Variant
from allelix.parsers import detect_parser
from allelix.reports._pipeline import _mark_haploid, _record_runs, run_analysis

FIRST = "rs3999980501"
SECOND = "rs3999980502"
FILLER = "rs3999980599"


def _fixture(tmp_path, reason="missing", reverse=False):
    data = tmp_path / "data"
    data.mkdir()
    rows = [(FILLER, "2", 1900, "A", "G", 2000, "A", "G", 0)]
    forward, calls = {
        "missing": (None, ("AA", "A")),
        "palindromic": ("AT", ("AA", "A")),
        "both_fit": ("AGT", ("AA", "A")),
        "no_fit": ("AT", ("CC", "C")),
        "di_call": ("AG", ("DD", "D")),
        "placed_same": ("AG", ("AA", "A")),
        "placed_different": ("AG", ("AA", "A")),
    }[reason]
    if forward:
        for i, rsid in enumerate((FIRST, SECOND)):
            pos = 1000 if reason == "placed_same" else 1000 + i
            rows.extend(
                (rsid, "1", pos - 100, forward[0], alt, pos, forward[0], alt, 0)
                for alt in forward[1:]
            )
    _identity(data, rows)
    _pharm(data, FIRST, "MOCK", (calls[0],))
    # The second synthetic source row has a different rsID and the haploid call.
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn, conn:
        conn.execute(
            "INSERT INTO pharmgkb_annotations "
            "(rsid, genotype, gene, drugs, phenotype, phenotype_category, "
            "annotation_text, level_of_evidence, score, pgkb_annotation_id, "
            "function_class, is_nonfinding) "
            "SELECT ?, genotype, gene, drugs, phenotype, phenotype_category, "
            "annotation_text, level_of_evidence, score, 'synthetic-second', "
            "function_class, is_nonfinding FROM pharmgkb_annotations WHERE rsid=?",
            (SECOND, FIRST),
        )
    targets = [f"{FIRST}\t{calls[0]}", f"{SECOND}\t{calls[1]}"]
    if reverse:
        targets.reverse()
    path = tmp_path / "synthetic.txt"
    path.write_text("\n".join([*[f"{FILLER}\tAG"] * 3, *targets]) + "\n")
    return data, path


def _analysis(data, path, annotate=True):
    with contextlib.ExitStack() as stack:
        identity = stack.enter_context(RsidIdentityAnnotator(data))
        pharm = stack.enter_context(PharmGKBAnnotator(data))
        return run_analysis(
            path, detect_parser(path), [pharm] if annotate else [], identity=identity
        )


@pytest.mark.parametrize("reason", ("missing", "palindromic", "both_fit", "no_fit", "di_call"))
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("annotate", (False, True))
def test_unknown_positions_never_mean_two_unrelated_calls_share_a_locus(
    tmp_path, reason, reverse, annotate
):
    data, path = _fixture(tmp_path, reason, reverse)
    result = _analysis(data, path, annotate)
    assert result.placement.unplaced_total == 2
    assert not [w for w in result.report_warnings if w.code == "contradictory_calls"]


@pytest.mark.parametrize("reason,count", (("placed_same", 2), ("placed_different", 0)))
@pytest.mark.parametrize("reverse", (False, True))
def test_known_positions_keep_the_existing_contradiction_rule(tmp_path, reason, count, reverse):
    data, path = _fixture(tmp_path, reason, reverse)
    result = _analysis(data, path)
    assert result.placement.unplaced_total == 0
    warnings = [w for w in result.report_warnings if w.code == "contradictory_calls"]
    assert sum(w.count for w in warnings) == count


@pytest.mark.parametrize("reverse", (False, True))
def test_unknown_locus_origin_uses_each_findings_own_rsid(reverse):
    calls = [
        Variant(FIRST, "", 0, "A", "A", haploid=False, record_index=11),
        Variant(SECOND, "", 0, "C", "C", haploid=True, record_index=12),
    ]
    if reverse:
        calls.reverse()
    findings = [
        Annotation(
            source="pharmgkb",
            rsid=rsid,
            significance="synthetic",
            category="pharma",
            magnitude=6.0,
            description="Synthetic finding",
            attribution="ClinPGx",
            genotype_match=genotype,
            locus=("", 0),
        )
        for rsid, genotype in ((FIRST, "AA"), (SECOND, "CC"))
    ]
    _mark_haploid(findings, calls)
    assert [(a.rsid, a.genotype_match, a.haploid, a.record_index) for a in findings] == [
        (FIRST, "AA", False, 11),
        (SECOND, "C", True, 12),
    ]


@pytest.mark.parametrize("size", (16, 128))
def test_unrelated_unplaced_rows_stay_batched(tmp_path, monkeypatch, size):
    data, path = _fixture(tmp_path)
    path.write_text("".join(f"rs{3999980500 + i}\tA\n" for i in range(size)))
    batches = []
    original = PharmGKBAnnotator.batch_annotate

    def observed(self, variants):
        batch = list(variants)
        batches.append(len(batch))
        yield from original(self, batch)

    monkeypatch.setattr(PharmGKBAnnotator, "batch_annotate", observed)
    result = _analysis(data, path)
    # The accidental shared empty coordinate must not turn one small chunk
    # into a separate database query for every row.
    assert result.total_variants == size
    assert len(batches) <= 2, batches
    assert max(batches) >= size - 1


@pytest.mark.parametrize("placed", (False, True))
def test_true_duplicate_keys_still_split_annotation_runs(placed):
    chrom, pos = ("1", 1000) if placed else ("", 0)
    calls = [
        Variant(FIRST, chrom, pos, "A", "A"),
        Variant(FIRST, chrom, pos, "C", "C", haploid=True),
    ]
    runs = _record_runs(calls)
    assert runs == [[calls[0]], [calls[1]]]


@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("suffix", ("json", "html"))
def test_reports_do_not_publish_a_false_shared_position_warning(tmp_path, reverse, suffix):
    data, path = _fixture(tmp_path, reverse=reverse)
    output = tmp_path / f"report.{suffix}"
    result = CliRunner().invoke(
        main,
        ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(output)],
    )
    assert result.exit_code == 0, result.output
    assert "share a position" not in result.output
    text = output.read_text()
    assert "share a position" not in text
    if suffix == "json":
        report = json.loads(text)
        assert not [w for w in report["report_warnings"] if w["code"] == "contradictory_calls"]
        rows = [a for a in report["annotations"] if a["rsid"] in (FIRST, SECOND)]
        assert len(rows) == 2
        assert all(a["chromosome"] is None and a["position"] is None for a in rows)
        assert {a["rsid"]: a["genotype_match"] for a in rows} == {FIRST: "AA", SECOND: "A"}


@pytest.mark.parametrize("present", (False, True))
@pytest.mark.parametrize("enabled", (False, True))
@pytest.mark.parametrize("command", ("analyze", "pharmacogenomics"))
def test_required_identity_obeys_source_isolation(tmp_path, present, enabled, command):
    data, path = _fixture(tmp_path, "placed_different")
    if not present:
        (data / "rsid_identity.sqlite").unlink()
    config = load_config(data, create=False)
    config.sources["rsid_identity"] = enabled
    save_config(data, config)
    output = tmp_path / "report.json"
    result = CliRunner().invoke(
        main,
        [command, str(path), "--data-dir", str(data), "--no-update", "--output", str(output)],
    )
    if present and enabled:
        assert result.exit_code == 0, result.output
        assert output.is_file()
    else:
        assert result.exit_code == 1, result.output
        assert "identity" in result.output.lower()
        assert not output.exists()

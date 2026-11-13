"""Site evidence belongs to the call whose locus was checked, not its batch."""

import contextlib
import json
from dataclasses import replace

import pytest
from click.testing import CliRunner
from test_clinpgx_reverse_homozygotes import BUILDS, RS, _call, _data, _identity

from allelix.annotators import get_annotators
from allelix.cli import main


@contextlib.contextmanager
def _annotator(data):
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data)
        for item in registry:
            stack.enter_context(item)
        yield next(a for a in registry if a.name == "pharmgkb")


def _elsewhere(call, mismatch):
    if mismatch == "position":
        return replace(call, position=9999)
    if mismatch == "chromosome":
        return replace(call, chromosome="2")
    return replace(call, chromosome="chrUn_synthetic")


def _row_key(row):
    return row.rsid, row.locus, row.genotype_match, row.description


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("mismatch", ("position", "chromosome", "unknown_contig"))
@pytest.mark.parametrize("with_ref", (False, True))
def test_batch_cannot_lend_complement_evidence_to_another_locus(
    tmp_path, build, reverse, mismatch, with_ref
):
    data = _data(tmp_path)
    good = _call(build)
    if not with_ref:
        good = replace(good, ref=None)
    bad = _elsewhere(good, mismatch)
    calls = [bad, good] if reverse else [good, bad]
    with _annotator(data) as annotator:
        good_rows = annotator.annotate(good)
        bad_rows = annotator.annotate(bad)
        assert len(good_rows) == 1, "The placed reverse homozygote is a real positive control"
        assert bad_rows == [], "The other locus supplies no orientation evidence"
        batch = list(annotator.batch_annotate(calls))
    assert [_row_key(row) for row in batch] == [_row_key(row) for row in good_rows], (
        "One placed call must not authorize a complement finding at another locus"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("mismatch", ("position", "chromosome", "unknown_contig"))
def test_unplaced_elsewhere_singleton_does_not_gain_site_evidence(tmp_path, build, mismatch):
    data = _data(tmp_path)
    bad = _elsewhere(_call(build), mismatch)
    with _annotator(data) as annotator:
        assert annotator.annotate(bad) == []
        assert list(annotator.batch_annotate([bad])) == []


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("existing", ("direct", "guarded_heterozygote"))
def test_fix_does_not_remove_existing_rsid_matches_at_other_loci(
    tmp_path, build, reverse, existing
):
    data = _data(tmp_path)
    good = _call(build)
    bad = replace(good, position=9999)
    bad = (
        replace(bad, allele1="T", allele2="T")
        if existing == "direct"
        else replace(bad, allele1="G", allele2="A")
    )
    calls = [bad, good] if reverse else [good, bad]
    with _annotator(data) as annotator:
        expected = [row for call in calls for row in annotator.annotate(call)]
        assert len(expected) == 2
        actual = list(annotator.batch_annotate(calls))
    assert [_row_key(row) for row in actual] == [_row_key(row) for row in expected]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
def test_equivalent_chromosome_spellings_keep_both_calls(tmp_path, build, reverse):
    data = _data(tmp_path)
    first = _call(build)
    second = replace(first, chromosome="chr1")
    calls = [second, first] if reverse else [first, second]
    with _annotator(data) as annotator:
        expected = [row for call in calls for row in annotator.annotate(call)]
        assert len(expected) == 2
        actual = list(annotator.batch_annotate(calls))
    assert [_row_key(row) for row in actual] == [_row_key(row) for row in expected]


@pytest.mark.parametrize("build", BUILDS)
def test_evidence_does_not_survive_into_an_unrelated_later_batch(tmp_path, build):
    data = _data(tmp_path)
    good = _call(build)
    bad = replace(good, position=9999)
    with _annotator(data) as annotator:
        assert len(list(annotator.batch_annotate([good]))) == 1
        assert list(annotator.batch_annotate([bad])) == []


def test_unplaced_list_call_can_use_source_orientation_without_acquiring_a_locus(tmp_path):
    # Placement checks the user's AA; both AA and TT fit G/A/T, so it stays
    # unplaced. ClinPGx checks its full C/T set, which fits only in complement.
    # This is a passing scope control: formats.md must describe the distinction.
    data = _data(tmp_path / "data")
    _identity(
        data,
        [
            (RS, "1", 1000, "G", "A", 2000, "G", "A", 0),
            (RS, "1", 1000, "G", "T", 2000, "G", "T", 0),
        ],
    )
    source = tmp_path / "calls.txt"
    source.write_text(f"{RS}\tAA\n" + "".join(f"rs90050300{i}\tGG\n" for i in range(4)))
    output = tmp_path / "report.json"
    run = CliRunner().invoke(
        main,
        ["analyze", str(source), "--data-dir", str(data), "--no-update", "--output", str(output)],
    )
    assert run.exit_code == 0, (run.output, run.exception)
    report = json.loads(output.read_text())
    assert report["input"]["placement"]["unplaced"]["both_fit"] == 1
    rows = [row for row in report["annotations"] if row["source"] == "pharmgkb"]
    assert len(rows) == 1 and rows[0]["rsid"] == RS and rows[0]["genotype_match"] == "AA"
    assert rows[0]["chromosome"] is None and rows[0]["position"] is None

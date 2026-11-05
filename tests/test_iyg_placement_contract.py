"""Position-free rows use one per-row strand decision and retain unplaced evidence."""

import contextlib
import json

import pytest
from test_annotation_parity_and_named_disclosures import _identity, _pharm

from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.models import Annotation, annotation_to_public_dict
from allelix.parsers import detect_parser
from allelix.reports._pipeline import run_analysis
from allelix.reports.json_report import render_json
from allelix.utils.build_detect import BuildDecisionError

TARGET = "rs3999980101"
FILLER = "rs3999980102"
POS = 1000101


def _input(tmp_path, call, *, filler_count=3, order=False):
    rows = [f"{FILLER}\tAG"] * filler_count
    target = f"{TARGET}\t{call}"
    rows = [target, *rows] if order else [*rows, target]
    path = tmp_path / "synthetic.txt"
    path.write_text("\n".join(rows) + "\n")
    return path


def _data(tmp_path, forward, genotype, *, installed=True):
    data = tmp_path / "data"
    data.mkdir()
    if installed:
        rows = [(FILLER, "2", 2000202, "A", "G", 2000302, "A", "G", 0)]
        if forward:
            ref = forward[0]
            rows.extend(
                (TARGET, "1", POS - 100, ref, alt, POS, ref, alt, int(forward in ("AT", "CG")))
                for alt in forward[1:]
            )
        _identity(data, rows)
    _pharm(data, TARGET, "MOCK", (genotype,))
    return data


def _run(path, data, build=None):
    parser = detect_parser(path)
    with contextlib.ExitStack() as stack:
        identity = stack.enter_context(RsidIdentityAnnotator(data))
        pharm = stack.enter_context(PharmGKBAnnotator(data))
        return run_analysis(path, parser, [pharm], identity=identity, build_override=build)


@pytest.mark.parametrize("forward,call,expected", [("AG", "AG", "AG"), ("AG", "TC", "AG")])
@pytest.mark.parametrize("order", (False, True))
def test_direct_and_complement_rows_are_placed_independently(
    tmp_path, forward, call, expected, order
):
    data = _data(tmp_path, forward, expected)
    result = _run(_input(tmp_path, call, order=order), data)
    rows = [a for a in result.annotations if a.rsid == TARGET]
    assert len(rows) == 1
    assert rows[0].locus == ("1", POS)
    assert rows[0].genotype_match.replace("/", "") == expected
    assert result.build == "GRCh38"
    assert result.placement.placed == 4
    assert result.placement.complemented == int(call != expected)
    assert result.placement.unplaced_total == 0


@pytest.mark.parametrize("forward", ("AG", "CT", "AC", "GT", "AT", "CG", "ACG", "AGT", "ACGT"))
@pytest.mark.parametrize("allele", "ACGT")
@pytest.mark.parametrize("haploid", (False, True))
def test_homozygous_and_haploid_readings_have_the_same_placement_rule(
    tmp_path, forward, allele, haploid
):
    complement = dict(zip("ACGT", "TGCA", strict=True))[allele]
    direct, reverse = allele in forward, complement in forward
    placed = direct != reverse
    expected = complement if placed and reverse else allele
    data = _data(tmp_path, forward, expected * 2)
    result = _run(_input(tmp_path, allele if haploid else allele * 2), data)
    rows = [a for a in result.annotations if a.rsid == TARGET]
    assert rows, "An unplaced row must still reach an rsID-keyed source"
    for row in rows:
        public = annotation_to_public_dict(row)
        assert (public["chromosome"], public["position"]) == (
            ("1", POS) if placed else (None, None)
        )
        assert row.genotype_match.replace("/", "") == expected * (1 if haploid else 2)
        assert row.haploid is haploid
    assert result.placement.placed == 3 + int(placed)
    assert result.placement.complemented == int(placed and reverse)
    assert result.placement.unplaced_total == int(not placed)
    assert sum(result.placement.unplaced.values()) == int(not placed)


@pytest.mark.parametrize("call", ("CG", "GC", "AT", "TA"))
@pytest.mark.parametrize("filler_count", (3, 1000))
def test_both_readings_abstain_even_when_other_rows_unanimously_read_forward(
    tmp_path, call, filler_count
):
    data = _data(tmp_path, "ACGT", call)
    result = _run(_input(tmp_path, call, filler_count=filler_count), data)
    assert result.placement.placed == filler_count
    assert result.placement.unplaced_total == 1
    assert sum(result.placement.unplaced.values()) == 1
    for row in result.annotations:
        public = annotation_to_public_dict(row)
        assert public["chromosome"] is None and public["position"] is None


def test_missing_identity_row_stays_visible_without_invented_coordinates(tmp_path):
    data = _data(tmp_path, "", "AG")
    result = _run(_input(tmp_path, "AG"), data)
    rows = [annotation_to_public_dict(a) for a in result.annotations if a.rsid == TARGET]
    assert rows
    assert all(r["chromosome"] is None and r["position"] is None for r in rows)
    assert result.placement.placed == 3
    assert result.placement.unplaced_total == 1
    assert sum(result.placement.unplaced.values()) == 1


@pytest.mark.parametrize("build", (None, "GRCh38"))
def test_missing_identity_cache_stops_even_with_a_ready_rsid_source(tmp_path, build):
    data = _data(tmp_path, "AG", "AG", installed=False)
    path = _input(tmp_path, "AG")
    # Detection must succeed first: an unknown-format refusal is not this guard.
    assert detect_parser(path).name == "rsid-genotype"
    with pytest.raises(BuildDecisionError, match=r"(?i)identity"):
        _run(path, data, build)


@pytest.mark.parametrize("build", ("GRCh36", "GRCh37"))
def test_position_free_input_never_silently_uses_another_build(tmp_path, build):
    data = _data(tmp_path, "AG", "AG")
    path = _input(tmp_path, "AG")
    assert detect_parser(path).name == "rsid-genotype"
    with pytest.raises(BuildDecisionError, match=r"(?i)grch38"):
        _run(path, data, build)


@pytest.mark.parametrize("placed", (False, True))
def test_json_records_placement_counts_and_only_known_coordinates(tmp_path, placed):
    data = _data(tmp_path, "AG" if placed else "", "AG")
    before = {p.name: p.read_bytes() for p in data.glob("*.sqlite")}
    result = _run(_input(tmp_path, "TC" if placed else "AG"), data)
    output = tmp_path / "report.json"
    render_json(result, output_path=output)
    report = json.loads(output.read_text())
    counts = report["input"]["placement"]
    assert report["input"]["build"] == "GRCh38"
    assert counts["rows"] == 4
    assert counts["placed"] == 3 + int(placed)
    assert counts["complemented"] == int(placed)
    assert sum(counts["unplaced"].values()) == int(not placed)
    rows = [a for a in report["annotations"] if a["rsid"] == TARGET]
    assert rows
    assert all(
        (a["chromosome"], a["position"]) == (("1", POS) if placed else (None, None)) for a in rows
    )
    assert before == {p.name: p.read_bytes() for p in data.glob("*.sqlite")}


@pytest.mark.parametrize("locus,expected", [(("", 0), (None, None)), (("1", POS), ("1", POS))])
def test_public_annotation_coordinates_distinguish_unplaced_from_a_real_site(locus, expected):
    row = Annotation(
        source="pharmgkb",
        rsid=TARGET,
        significance="1a",
        category="drug_response",
        magnitude=6.0,
        description="Synthetic response",
        attribution="ClinPGx",
        genotype_match="AG",
        locus=locus,
    )
    public = annotation_to_public_dict(row)
    assert (public["chromosome"], public["position"]) == expected

"""A rejected sibling cannot outrank the passing split-site half call."""

from dataclasses import replace

import pytest
from rc3_verification_helpers import POS, RS, analyze, caches, extracted, findings, write_vcf

from allelix.cli.utility import _extract_choice
from allelix.models import Variant


@pytest.mark.parametrize("build,indexed", [("GRCh37", False), ("GRCh38", False), ("GRCh38", True)])
@pytest.mark.parametrize("rejection", ["filter", "ft", "gq", "dp", "unknown"])
@pytest.mark.parametrize("rejected_first", [True, False])
def test_extract_agrees_with_analysis_of_the_passing_half(
    tmp_path, monkeypatch, build, indexed, rejection, rejected_first
):
    pairs = [(RS, "A"), (RS, "C")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, [])
    rejected = dict(alt="A", gt="1/0")
    rejected.update(
        {
            "filter": {"filter": "LowQual"},
            "ft": {"ft": "LowQual"},
            "gq": {"gq": 1},
            "dp": {"dp": 0},
            "unknown": {"gt": "./."},
        }[rejection]
    )
    rows = [rejected, dict(alt="C", gt="0/1")]
    if not rejected_first:
        rows.reverse()
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    cv = findings(report, "clinvar")
    assert [(a["alt"], a["genotype_match"], a["zygosity"]) for a in cv] == [
        ("C", "C/.", "half_call")
    ]
    row = extracted(data, source, build, [RS], monkeypatch, indexed=indexed)[RS]
    assert row[3] in ("C/-", "-/C"), "Extract must show the known allele, not -/- or C/G"
    assert row[2].replace(",", "") == str(POS[build])


@pytest.mark.parametrize("reverse", [False, True])
def test_usable_half_beats_full_no_call_at_same_placement(reverse):
    empty = Variant(RS, "1", 100, "-", "-", ref="G", record_index=0)
    half = replace(empty, allele1="C", record_index=1)
    candidates = [empty, half]
    if reverse:
        candidates.reverse()
    assert _extract_choice(candidates, {("1", 100): "G"}, {"1": 0}) is half


def test_placement_still_beats_usable_call_elsewhere_and_ties_use_file_order():
    empty = Variant(RS, "1", 100, "-", "-", ref="G", record_index=0)
    elsewhere = replace(empty, position=200, allele1="C", record_index=1)
    assert _extract_choice([elsewhere, empty], {("1", 100): "G"}, {"1": 0}) is empty
    earlier = replace(empty, allele1="A", record_index=2)
    later = replace(empty, allele1="C", record_index=3)
    assert _extract_choice([later, earlier], {("1", 100): "G"}, {"1": 0}) is earlier

"""Addendum 1: rsID compatibility, allele-only recovery, and record counts."""

import pytest
from split_site_contract_helpers import BUILDS, OTHER, RSID, record, write_vcf

from allelix.models import NO_CALL_MARKER
from allelix.parsers.vcf import VcfParser

# Expected canonical merged ID is independent of the implementation helper.
COMPATIBLE = (
    (RSID, ".", RSID),
    (".", RSID, RSID),
    (RSID.upper(), RSID, RSID),
    ("COSV1;" + RSID.upper(), RSID + ";RCV2", RSID),
    ("22:11924846:G:C", "22:11924847:G:T", ""),
    ("COSV123", "RCV456", ""),
    ("COSV123", "COSV123", ""),
    (".", ".", ""),
    ("COSV123", RSID, RSID),
    (RSID, "RCV456", RSID),
    ("22:11924846:G:C", RSID, RSID),
    (RSID, "22:11924847:G:T", RSID),
)


def _parse(tmp_path, build, rows):
    parser = VcfParser()
    parser.combine_split_sites = True
    calls = list(parser.parse(write_vcf(tmp_path / "input.vcf", build, rows)))
    return parser, calls


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("ids", COMPATIBLE)
@pytest.mark.parametrize("gts", [("1/0", "0/1"), ("1/.", "./1"), ("1|.", ".|1")])
def test_only_distinct_canonical_rsids_can_block_an_exact_pair(tmp_path, build, reverse, ids, gts):
    first, second, expected = ids
    rows = [record("A", gts[0], rsid=first, ps="42"), record("C", gts[1], rsid=second, ps="42")]
    parser, calls = _parse(tmp_path, build, rows[::-1] if reverse else rows)
    assert len(calls) == 1
    call = calls[0]
    assert sorted((call.allele1, call.allele2)) == ["A", "C"]
    assert call.rsid == expected, "Keep the sole canonical rsID; otherwise resolve by position"
    assert not call.is_half_call and parser.split_site_calls == 0


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("gts", [("1/0", "0/1"), ("1/.", "./1"), ("1|0", "0|1")])
def test_two_distinct_rsids_prevent_pair_but_do_not_create_a_genotype_contradiction(
    tmp_path, build, reverse, gts
):
    rows = [record("A", gts[0], rsid=RSID, ps="42"), record("C", gts[1], rsid=OTHER, ps="42")]
    parser, calls = _parse(tmp_path, build, rows[::-1] if reverse else rows)
    assert len(calls) == 2, "Compatible genotype evidence survives as separate half calls"
    assert {(v.rsid, frozenset((v.allele1, v.allele2))) for v in calls} == {
        (RSID, frozenset(("A", NO_CALL_MARKER))),
        (OTHER, frozenset(("C", NO_CALL_MARKER))),
    }
    assert all(v.is_half_call for v in calls)
    assert parser.split_site_calls == 2


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gts", [("1/1", "0/1"), ("1", "0/1"), ("0|1", "0|1")])
def test_distinct_rsids_do_not_rescue_a_true_genotype_contradiction(tmp_path, build, gts):
    parser, calls = _parse(
        tmp_path,
        build,
        [record("A", gts[0], rsid=RSID, ps="42"), record("C", gts[1], rsid=OTHER, ps="42")],
    )
    assert not calls
    assert parser.split_site_calls == 2


# The union of the non-pair guard shapes in test_split_person_genotype,
# test_split_site_phased_pairs and test_split_site_gate_and_vcf_origin.
AGREEING = (
    (record("A", "1/0"), record("C", "1/0")),
    (record("A", "0/1"), record("C", "0/1")),
    (record("A", "1/."), record("C", "0/1")),
    (record("A", "1/0"), record("C", "./1")),
    (record("A", "0|1", ps="300"), record("C", "1|0", ps="301")),
    (record("A", "1|0", ps="300"), record("C", "0|1", ps="301")),
    (record("A", "0|1", ps="300"), record("C", "1/0")),
    (record("A", "1/0"), record("C", "0|1", ps="300")),
)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("rows", AGREEING)
def test_agreeing_guard_shapes_emit_alt_unknown_never_alt_ref(tmp_path, build, rows):
    parser, calls = _parse(tmp_path, build, rows)
    assert len(calls) == 2
    assert {frozenset((v.allele1, v.allele2)) for v in calls} == {
        frozenset(("A", NO_CALL_MARKER)),
        frozenset(("C", NO_CALL_MARKER)),
    }
    assert all(v.is_half_call for v in calls)
    assert parser.split_site_calls == 2


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "rows",
    [
        [record("A", "0/1", filt="q10"), record("C", "0/1", filt="q10")],
        [record("A", "./."), record("C", "./.")],
        [record("A", "0/1"), record("C", "./.")],
        [record("A", "1/0"), record("C", "0/1"), record("T", "0/1", filt="q10")],
    ],
)
def test_non_pair_counts_every_component_including_rejected_records(tmp_path, build, rows):
    parser, _calls = _parse(tmp_path, build, rows)
    assert parser.split_site_calls == len(rows)

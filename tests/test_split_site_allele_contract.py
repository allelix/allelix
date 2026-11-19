"""Allele evidence survives agreeing split representations without inventing a genotype."""

import pytest
from split_site_contract_helpers import (
    BUILDS,
    assert_assessed_halves,
    by_source,
    data_dir,
    record,
    run,
    write_vcf,
)

TABLE = (
    ("merged", (record("A,T", "1/2"),), 2),
    ("complementary", (record("A", "1/0"), record("T", "0/1")), 2),
    ("half_pair", (record("A", "1/."), record("T", "./1")), 2),
    ("same_order", (record("A", "0/1"), record("T", "0/1")), 2),
    ("missing_id", (record("A", "1/0"), record("T", "0/1", rsid=".")), 2),
    ("filtered", (record("A", "0/1"), record("T", "0/1", filt="LowQual")), 1),
    ("low_gq", (record("A", "0/1"), record("T", "0/1", gq=3)), 1),
    ("missing_call", (record("A", "0/1"), record("T", "./.")), 1),
)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("auto", (False, True), ids=("explicit", "auto"))
@pytest.mark.parametrize("name,rows,expected", TABLE, ids=[r[0] for r in TABLE])
def test_eight_row_clinvar_table_keeps_released_allele_findings(
    tmp_path, build, auto, name, rows, expected
):
    data = data_dir(tmp_path / "data", identity=True, alts="AT")
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    report, _ = run(data, source, tmp_path / "out.json", build, auto=auto)
    assert report["input"]["build"] == build
    assert len(by_source(report)["clinvar"]) == expected, name


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "gts,phase",
    [
        (("1/.", "./1"), "."),
        (("./1", "1/."), "."),
        (("1|.", ".|1"), "70"),
        ((".|1", "1|."), "70"),
        (("1|.", ".|1"), "."),
        (("1/0", "0/1"), "."),
    ],
)
@pytest.mark.parametrize("reverse", (False, True))
def test_complete_complementary_halves_are_one_full_call_for_every_source(
    tmp_path, build, gts, phase, reverse
):
    data = data_dir(tmp_path / "data")
    rows = [record("A", gts[0], ps=phase), record("C", gts[1], ps=phase)]
    source = write_vcf(tmp_path / "input.vcf", build, rows[::-1] if reverse else rows)
    report, text = run(data, source, tmp_path / "out.json", build)
    groups = by_source(report)
    assert {s: len(a) for s, a in groups.items()} == {
        "clinvar": 2,
        "gwas": 1,
        "pharmgkb": 1,
        "snpedia": 1,
    }
    assert {a["genotype_match"] for a in report["annotations"]} == {"AC"}
    assert not any(w["code"] == "split_site_calls" for w in report["report_warnings"])
    assert "half call(s)" not in text, "The full call was assessed by genotype sources too"


AGREEING = [
    (record("A", "0/1"), record("T", "0/1")),
    (record("A", "1/0"), record("T", "1/0")),
    (record("A", "0|1", ps="10"), record("T", "1|0", ps="11")),
    (record("A", "0|1", ps="10"), record("T", "0|1", ps="11")),
    (record("A", "1/.", ps="."), record("T", "0/1", ps=".")),
    (record("A", "0|1", ps="10"), record("T", "1/0")),
    (record("A", "1|.", ps="10"), record("T", ".|1", ps="11")),
    (record("A", "1|.", ps="10"), record("T", ".|1")),
    (record("A", "1/0"), record("T", "./1")),
    (record("A", "1/0"), record("T", "0|1", ps="10")),
    (record("A", "1|0", ps="10"), record("T", "0|1", ps="11")),
]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("rows", AGREEING)
def test_agreeing_uncertain_genotypes_give_alleles_only_to_allele_sources(
    tmp_path, build, reverse, rows
):
    rows = [dict(r, alt="C") if r["alt"] == "T" else r for r in rows]
    data = data_dir(tmp_path / "data")
    source = write_vcf(tmp_path / "input.vcf", build, rows[::-1] if reverse else rows)
    report, text = run(data, source, tmp_path / "out.json", build)
    groups = by_source(report)
    assert {a["alt"] for a in groups["clinvar"]} == {"A", "C"}
    assert {a["alt"] for a in groups["gwas"]} == {"A"}
    assert not groups["pharmgkb"] and not groups["snpedia"], (
        "Never give AG and GT for an AT person"
    )
    assert all("." in a["genotype_match"] for a in report["annotations"])
    assert_assessed_halves(text, 2)
    warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warnings) == 1 and warnings[0]["count"] == len(rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "bad", [dict(filt="q10"), dict(gq=3), dict(dp=0), dict(ft="LowAD"), dict(gt="./.")]
)
@pytest.mark.parametrize("reverse", (False, True))
def test_rejected_or_missing_sibling_cancels_no_passing_allele(tmp_path, build, bad, reverse):
    data = data_dir(tmp_path / "data")
    sibling = dict(record("C", "1/1"), **bad)
    rows = [record("A", "0/1"), sibling]
    report, text = run(
        data,
        write_vcf(tmp_path / "input.vcf", build, rows[::-1] if reverse else rows),
        tmp_path / "out.json",
        build,
    )
    groups = by_source(report)
    assert {a["alt"] for a in groups["clinvar"]} == {"A"}
    assert {a["alt"] for a in groups["gwas"]} == {"A"}
    assert not groups["pharmgkb"] and not groups["snpedia"]
    assert_assessed_halves(text, 1)
    warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warnings) == 1 and warnings[0]["count"] == 2


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("bad", [dict(filt="q10"), dict(ft="LowAD"), dict(gq=3), dict(gt="./.")])
def test_rejected_third_record_keeps_alleles_without_inventing_a_full_call(tmp_path, build, bad):
    data = data_dir(tmp_path / "data")
    rows = [record("A", "1/0"), record("C", "0/1"), dict(record("T", "1/1"), **bad)]
    report, text = run(
        data, write_vcf(tmp_path / "input.vcf", build, rows), tmp_path / "out.json", build
    )
    groups = by_source(report)
    assert {s: len(a) for s, a in groups.items()} == {
        "clinvar": 2,
        "gwas": 1,
        "pharmgkb": 0,
        "snpedia": 0,
    }
    assert all("." in a["genotype_match"] for a in report["annotations"])
    assert_assessed_halves(text, 2)
    warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warnings) == 1 and warnings[0]["count"] == 3


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gts", [("1/0", "0/1"), ("1/.", "./1"), ("1|.", ".|1")])
def test_at_person_never_gets_ag_or_gt_genotype_findings(tmp_path, build, gts):
    data = data_dir(tmp_path / "data", alts="AT")
    source = write_vcf(
        tmp_path / "input.vcf", build, [record("A", gts[0], ps="12"), record("T", gts[1], ps="12")]
    )
    report, _ = run(data, source, tmp_path / "out.json", build)
    groups = by_source(report)
    assert [
        (a["source"], a["genotype_match"])
        for a in report["annotations"]
        if a["source"] in ("pharmgkb", "snpedia")
    ] == [("pharmgkb", "AT"), ("snpedia", "AT")]
    assert len(groups["clinvar"]) == 2


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["1/.", "./1"])
def test_standalone_half_call_positive_control(tmp_path, build, gt):
    data = data_dir(tmp_path / "data")
    report, text = run(
        data,
        write_vcf(tmp_path / "input.vcf", build, [record("A", gt)]),
        tmp_path / "out.json",
        build,
    )
    groups = by_source(report)
    assert {s: len(a) for s, a in groups.items()} == {
        "clinvar": 1,
        "gwas": 1,
        "pharmgkb": 0,
        "snpedia": 0,
    }
    assert_assessed_halves(text, 1)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "kind,expected",
    [("full_half_pair", 0), ("contradicting_half", 0), ("agreeing", 2), ("rejected_sibling", 1)],
)
def test_half_call_notice_counts_only_the_calls_actually_assessed(tmp_path, build, kind, expected):
    data = data_dir(tmp_path / "data")
    rows = {
        "full_half_pair": [record("A", "1/."), record("C", "./1")],
        "contradicting_half": [record("A", "1/."), record("C", "1/1")],
        "agreeing": [record("A", "0/1"), record("C", "0/1")],
        "rejected_sibling": [record("A", "0/1"), record("C", "0/1", filt="q10")],
    }[kind]
    _report, text = run(
        data, write_vcf(tmp_path / "input.vcf", build, rows), tmp_path / "out.json", build
    )
    if expected:
        assert_assessed_halves(text, expected)
    else:
        assert "ClinVar and the GWAS Catalog assessed the called allele" not in text

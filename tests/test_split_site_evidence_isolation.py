"""P/LP allele evidence belongs to one withheld site occurrence, not its locus alone."""

import sqlite3
from contextlib import closing

import pytest
from split_site_contract_helpers import BUILDS, OTHER, POS, data_dir, record, run, write_vcf

NONMATCH = "rs999000557"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("order", ["positive_first", "positive_last", "nonmatching_alone"])
def test_revisited_locus_does_not_share_carried_alt_evidence(tmp_path, build, order):
    data = data_dir(tmp_path / "data", identity=False)
    with closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE clinvar_variants SET aggregate='Benign', clinical_significance='Benign' "
            "WHERE alt != 'A'"
        )
    positive = [record("A", "1/1", rsid=OTHER), record("C", "0/1", rsid=OTHER)]
    negative = [record("C", "1/1", rsid=NONMATCH), record("T", "0/1", rsid=NONMATCH)]
    boundary = [record("A", "0/0", rsid=".", pos=POS[build] + 100)]
    rows = (
        negative
        if order == "nonmatching_alone"
        else (
            positive + boundary + negative
            if order == "positive_first"
            else negative + boundary + positive
        )
    )
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    report, _ = run(data, source, tmp_path / "out.json", build, extra=("--show-unchecked",))
    assert not report["annotations"]
    assert next(
        w["count"] for w in report["report_warnings"] if w["code"] == "split_site_calls"
    ) == (2 if order == "nonmatching_alone" else 4)
    named = [s for s in report.get("high_value_unchecked", []) if s["reason"] == "split_site"]
    assert [s["rsid"] for s in named] == ([] if order == "nonmatching_alone" else [OTHER]), (
        "A different withheld group revisiting the locus must not lend its P/LP ALT to this group"
    )

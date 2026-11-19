"""A withheld split site with P/LP evidence for a carried ALT is named inline (2.7.1).

Reviewer A's A-10 and review's extra case: the site is named on the console,
not only counted, whether its records name the P/LP row's rsID, no rsID, or
another rsID. When they name another rsID, the text gives the P/LP row at
the position as evidence and does not say ClinVar classifies the records'
own rsID P/LP. An ALT the records do not carry names nothing.
"""

import sqlite3
from contextlib import closing

import pytest
from split_site_contract_helpers import BUILDS, OTHER, RSID, data_dir, record, run, write_vcf


def _named(report):
    return [s for s in report.get("high_value_unchecked", []) if s["reason"] == "split_site"]


def _contradiction(rsid):
    return [record("A", "1/1", rsid=rsid), record("C", "0/1", rsid=rsid)]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("ids", ["own", "unnamed"])
def test_site_with_its_own_plp_row_is_named_inline(tmp_path, build, ids):
    data = data_dir(tmp_path / "data", identity=False, alts="AC")
    rows = _contradiction(RSID if ids == "own" else ".")
    report, text = run(
        data, write_vcf(tmp_path / "in.vcf", build, rows), tmp_path / "o.json", build
    )
    named = _named(report)
    assert [s["rsid"] for s in named] == [RSID]
    flat = " ".join(text.split())
    assert "0 named here" not in flat, "The site must be named, not only counted"
    assert RSID in flat and "SYNGENE" in flat


@pytest.mark.parametrize("build", BUILDS)
def test_site_naming_another_rsid_is_named_with_the_position_evidence(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=False, alts="AC")
    rows = _contradiction(OTHER)
    report, text = run(
        data, write_vcf(tmp_path / "in.vcf", build, rows), tmp_path / "o.json", build
    )
    named = _named(report)
    assert [s["rsid"] for s in named] == [OTHER]
    # The P/LP row is named as position evidence, not as the records' own rsID's class.
    assert f"{RSID} (SYNGENE)" in named[0]["message"], named[0]["message"]
    assert named[0].get("gene", "") == "", "The P/LP row's gene is not the records' rsID's gene"
    flat = " ".join(text.split())
    assert "0 named here" not in flat
    assert OTHER in flat and RSID in flat and "SYNGENE" in flat
    assert f"SYNGENE ({OTHER})" not in flat, "Do not label the other rsID with the row's gene"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("ids", ["own", "unnamed", "other"])
def test_plp_only_for_an_uncarried_alt_names_nothing(tmp_path, build, ids):
    data = data_dir(tmp_path / "data", identity=False, alts="AC")
    with closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE clinvar_variants SET aggregate='Benign', clinical_significance='Benign' "
            "WHERE alt != 'T'"
        )
    rsid = {"own": RSID, "unnamed": ".", "other": OTHER}[ids]
    report, _ = run(
        data,
        write_vcf(tmp_path / "in.vcf", build, _contradiction(rsid)),
        tmp_path / "o.json",
        build,
    )
    assert not _named(report)

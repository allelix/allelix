"""A contradicting split site at a high-value rsID's placement is named under that rsID (2.7.1).

Reviewer B's B-03: records G>A 1/1 + G>C 0/1 at rs4680's identity placement
contradict each other, so no source reads them. Written rs4680 on both, the
site is named in ``high_value_unchecked``; with dot IDs, positional IDs or
another rsID it was only counted. SUP's ruling: the site is attributed to
high-value rsIDs by placement (chromosome, position, REF, whatever ID or ALT),
as a record is, and named. A wrong REF or another position names nothing.
"""

import sqlite3
from contextlib import closing

import pytest
from split_site_contract_helpers import (
    BUILDS,
    OTHER,
    POS,
    by_source,
    data_dir,
    record,
    run,
    write_vcf,
)

HV = "rs4680"


def _data(path, *, ref="G", shift=0, plp_other=False):
    """rs4680 placed by the identity table at the site (or ``shift`` away), REF ``ref``.

    ClinVar holds no row there, so no P/LP evidence explains the site; with
    ``plp_other`` it holds P/LP rows there under another rsID.
    """
    data = data_dir(path, identity=True, rsid=HV, pathogenic=plp_other, alts="AC")
    for name in ("clinvar.GRCh37.sqlite", "clinvar.GRCh38.sqlite"):
        with closing(sqlite3.connect(data / name)) as conn, conn:
            if plp_other:
                conn.execute("UPDATE clinvar_variants SET rsid=?", (OTHER,))
            else:
                conn.execute("DELETE FROM clinvar_variants")
    with closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE rsid_identity SET ref37=?, ref38=?, pos37=pos37+?, pos38=pos38+?",
            (ref, ref, shift, shift),
        )
    return data


def _ids(style, build):
    return {
        "dot": (".", "."),
        "positional": (f"1:{POS[build]}:G:A", f"1:{POS[build]}:G:C"),
        "other_rsid": (OTHER, OTHER),
        "high_value": (HV, HV),
    }[style]


def _run(tmp_path, build, data, style):
    first, second = _ids(style, build)
    rows = [record("A", "1/1", rsid=first), record("C", "0/1", rsid=second)]
    return run(
        data,
        write_vcf(tmp_path / "input.vcf", build, rows),
        tmp_path / "out.json",
        build,
        extra=("--show-unchecked",),
    )


def _named(report):
    return [r for r in report.get("high_value_unchecked", []) if r.get("reason") == "split_site"]


def _split_count(report):
    return [w["count"] for w in report["report_warnings"] if w["code"] == "split_site_calls"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("style", ["dot", "positional", "other_rsid", "high_value"])
def test_contradicting_site_at_a_high_value_placement_is_named(tmp_path, build, style):
    report, text = _run(tmp_path, build, _data(tmp_path / "data"), style)
    assert all(not found for found in by_source(report).values()), "Still withheld"
    named = _named(report)
    assert [r["rsid"] for r in named] == [HV], "Named once, under the high-value rsID"
    assert named[0]["gene"] == "COMT"
    assert HV in text
    assert _split_count(report) == [2], "The count of records is unchanged"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("style", ["dot", "positional", "other_rsid"])
def test_wrong_ref_names_nothing(tmp_path, build, style):
    report, _text = _run(tmp_path, build, _data(tmp_path / "data", ref="T"), style)
    assert not _named(report)
    assert _split_count(report) == [2]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("style", ["dot", "positional", "other_rsid"])
def test_another_position_names_nothing(tmp_path, build, style):
    report, _text = _run(tmp_path, build, _data(tmp_path / "data", shift=1), style)
    assert not _named(report)
    assert _split_count(report) == [2]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("style", ["dot", "positional"])
def test_plp_naming_keeps_its_own_entry_beside_the_high_value_one(tmp_path, build, style):
    """The P/LP row's rsID is still named for the carried ALT; the high-value rsID is added."""
    report, _text = _run(tmp_path, build, _data(tmp_path / "data", plp_other=True), style)
    named = {r["rsid"]: r for r in _named(report)}
    assert sorted(named) == sorted([HV, OTHER])
    assert named[HV]["gene"] == "COMT"
    assert named[OTHER]["gene"] == "SYNGENE"
    assert _split_count(report) == [2]

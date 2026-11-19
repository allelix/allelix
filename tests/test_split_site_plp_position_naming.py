"""A withheld split site is named by ClinVar's P/LP row at its position for a carried ALT.

SUP's ruling for 2.7.1: named when high-value, or when ClinVar holds a
Pathogenic or Likely pathogenic row at the position for an ALT the person's
records carry, whatever rsIDs the records name. A P/LP row for an ALT the
records do not carry does not count.
"""

import sqlite3
from contextlib import closing

import pytest
from split_site_contract_helpers import BUILDS, OTHER, data_dir, record, run, write_vcf

CONTRADICTION = (record("A", "1/1", rsid=OTHER), record("C", "0/1", rsid=OTHER))


def _add_benign_other(data, build):
    """Give OTHER a Benign row at the same position (ClinVar holds it, not as P/LP)."""
    with closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn, conn:
        info = conn.execute("PRAGMA table_info(clinvar_variants)")
        columns = [r[1] for r in info if r[1] != "id"]
        values = conn.execute(
            "SELECT " + ",".join(columns) + " FROM clinvar_variants LIMIT 1"
        ).fetchone()
        row = dict(zip(columns, values, strict=True))
        row.update(rsid=OTHER, aggregate="Benign", clinical_significance="Benign")
        conn.execute(
            f"INSERT INTO clinvar_variants ({','.join(columns)}) "
            f"VALUES ({','.join('?' for _ in columns)})",
            [row[c] for c in columns],
        )


def _split_named(report):
    return [s for s in report.get("high_value_unchecked", []) if s["reason"] == "split_site"]


def _split_count(report):
    return next(w["count"] for w in report["report_warnings"] if w["code"] == "split_site_calls")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("held_benign", (False, True), ids=("unheld", "benign_only"))
def test_non_plp_rsid_at_a_plp_position_is_named(tmp_path, build, held_benign):
    data = data_dir(tmp_path / "data", identity=False)
    if held_benign:
        _add_benign_other(data, build)
    source = write_vcf(tmp_path / "input.vcf", build, list(CONTRADICTION))
    report, text = run(data, source, tmp_path / "out.json", build, extra=("--show-unchecked",))
    assert not report["annotations"]
    assert _split_count(report) == 2
    named = _split_named(report)
    assert [s["rsid"] for s in named] == [OTHER]
    # The P/LP row is another rsID's: shown as evidence, not as this rsID's gene.
    assert named[0]["gene"] == ""
    assert "SYNGENE" in named[0]["message"]
    assert "SYNGENE" in text


@pytest.mark.parametrize("build", BUILDS)
def test_plp_row_for_an_alt_the_records_do_not_carry_does_not_name_the_site(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=False)
    # ClinVar keeps P/LP only for T; the records carry A and C.
    with closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE clinvar_variants SET aggregate='Benign', clinical_significance='Benign' "
            "WHERE alt IN ('A', 'C')"
        )
    source = write_vcf(tmp_path / "input.vcf", build, list(CONTRADICTION))
    report, _text = run(data, source, tmp_path / "out.json", build)
    assert not report["annotations"]
    assert _split_count(report) == 2
    assert not _split_named(report)

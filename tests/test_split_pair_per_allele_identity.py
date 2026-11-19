"""An exact split pair's identity is decided per ALT (2.7.1 RC3, SUP's per-ALT rule).

Each record of an accepted exact pair keeps its own written rsID; a record
naming none has its ALT looked up alone, by position. Two different rsIDs: the
ALTs are half calls, each under its own rsID (Route S). Exactly one: the full
two-ALT call takes it (D1 point 1), as when one record names it. None: the
unnamed pair is unchanged.

Reviewers' shapes: A-01b and B-02, an unnamed (or positionally named) pair of
which only one ALT is known, lost the GWAS Catalog finding 2.6.3 gave and
extract said "not in file"; B-01, one record naming its rsID hid the other
ALT's rsID, found by position.
"""

import json
import sqlite3
from contextlib import closing

import pysam
import pytest
from click.testing import CliRunner
from split_site_contract_helpers import (
    BUILDS,
    OTHER,
    POS,
    RSID,
    by_source,
    data_dir,
    record,
    run,
    signature,
    write_vcf,
)

from allelix.cli import main, utility

PAIR_SHAPES = {
    "unphased": (("1/0", "."), ("0/1", ".")),
    "half_pair": (("1/.", "."), ("./1", ".")),
    "phased": (("1|0", "7"), ("0|1", "7")),
}


def _data(path, alts, owners):
    """Caches where each ALT in ``owners`` belongs to that rsID, or to none (None)."""
    data = data_dir(path, identity=True, alts=alts)
    for name in ("clinvar.GRCh37.sqlite", "clinvar.GRCh38.sqlite"):
        with closing(sqlite3.connect(data / name)) as conn, conn:
            for alt in "ATC":
                owner = owners.get(alt)
                if owner is None:
                    conn.execute("DELETE FROM clinvar_variants WHERE alt=?", (alt,))
                elif owner != RSID:
                    conn.execute("UPDATE clinvar_variants SET rsid=? WHERE alt=?", (owner, alt))
    with closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn, conn:
        for alt in alts:
            owner = owners.get(alt)
            for column in ("alt37", "alt38"):
                if owner is None:
                    conn.execute(f"DELETE FROM rsid_identity WHERE {column}=?", (alt,))
                elif owner != RSID:
                    conn.execute(f"UPDATE rsid_identity SET rsid=? WHERE {column}=?", (owner, alt))
    with closing(sqlite3.connect(data / "gwas.sqlite")) as conn, conn:
        for alt in alts:
            owner = owners.get(alt)
            if owner is None:
                conn.execute("DELETE FROM gwas_associations WHERE risk_allele=?", (alt,))
            elif owner != RSID:
                conn.execute(
                    "UPDATE gwas_associations SET rsid=? WHERE risk_allele=?", (owner, alt)
                )
    return data


def _rows(alts, ids, shape="unphased"):
    return [
        record(alt, gt, rsid=rsid, ps=ps)
        for alt, (gt, ps), rsid in zip(alts, PAIR_SHAPES[shape], ids, strict=True)
    ]


def _counts(report):
    split = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    return (
        split[0]["count"] if split else 0,
        report["input"]["vcf_quality"]["split_site_half_calls"],
        report["input"]["total_variants"],
    )


def _analyze(tmp_path, data, build, name, rows, *, vcf_out=False):
    extra = ("--vcf-out", str(tmp_path / f"{name}.annotated.vcf")) if vcf_out else ()
    report, _ = run(
        data,
        write_vcf(tmp_path / f"{name}.vcf", build, rows),
        tmp_path / f"{name}.json",
        build,
        extra=extra,
    )
    return report


def _round_trip(tmp_path, data, build, name, report):
    annotated = tmp_path / f"{name}.annotated.vcf"
    ids = [line.split("\t")[2] for line in annotated.read_text().splitlines() if line[0] != "#"]
    again, _ = run(data, annotated, tmp_path / f"{name}.again.json", build)
    assert signature(again) == signature(report)
    assert _counts(again) == _counts(report)
    assert json.dumps(again["report_warnings"]) == json.dumps(report["report_warnings"])
    return ids


def _gwas_rsids(report):
    return {a["rsid"] for a in by_source(report)["gwas"]}


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", sorted(PAIR_SHAPES))
def test_one_known_alt_gives_the_full_call_that_rsid(tmp_path, build, shape):
    """A-01b: only A is known by position; the full A/T call takes its rsID."""
    data = _data(tmp_path / "data", "AT", {"A": RSID, "T": None})
    unnamed = _analyze(
        tmp_path, data, build, "unnamed", _rows("AT", (".", "."), shape), vcf_out=True
    )
    written = _analyze(tmp_path, data, build, "written", _rows("AT", (RSID, "."), shape))
    assert RSID in _gwas_rsids(unnamed), "The GWAS Catalog finding 2.6.3 gave is back"
    assert {a["genotype_match"] for a in by_source(unnamed)["clinvar"]} == {"AT"}
    assert signature(unnamed) == signature(written), "The same as the record naming the rsID"
    assert _counts(unnamed) == _counts(written) == (0, 0, 1), "One full call, no split site"
    assert _round_trip(tmp_path, data, build, "unnamed", unnamed) == [RSID, RSID]


@pytest.mark.parametrize("build", BUILDS)
def test_positional_ids_with_one_known_alt_give_the_full_call_that_rsid(tmp_path, build):
    """B-02: positional IDs are no rsIDs; the one ALT found by position names the call."""
    data = _data(tmp_path / "data", "AC", {"A": RSID, "C": None})
    ids = (f"1:{POS[build]}:G:A", f"1:{POS[build]}:G:C")
    report = _analyze(tmp_path, data, build, "positional", _rows("AC", ids), vcf_out=True)
    assert RSID in _gwas_rsids(report)
    assert {a["genotype_match"] for a in by_source(report)["clinvar"]} == {"AC"}
    assert _counts(report) == (0, 0, 1)
    assert _round_trip(tmp_path, data, build, "positional", report) == [RSID, RSID]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("ids", [(RSID, "."), (".", OTHER)], ids=["first_named", "second_named"])
def test_one_named_record_and_the_other_alts_own_rsid_are_no_pair(tmp_path, build, ids):
    """B-01: the written rsID is that record's only; the other ALT's rsID is found by position."""
    data = _data(tmp_path / "data", "AC", {"A": RSID, "C": OTHER})
    report = _analyze(tmp_path, data, build, "one_named", _rows("AC", ids), vcf_out=True)
    both = _analyze(tmp_path, data, build, "both_named", _rows("AC", (RSID, OTHER)))
    assert {(a["rsid"], a["alt"]) for a in by_source(report)["clinvar"]} == {
        (RSID, "A"),
        (OTHER, "C"),
    }
    assert all("." in a["genotype_match"] for a in report["annotations"]), "Half calls"
    assert signature(report) == signature(both), "The same as both records naming their rsIDs"
    assert _counts(report) == _counts(both) == (2, 2, 2)
    assert _round_trip(tmp_path, data, build, "one_named", report) == [RSID, OTHER]


@pytest.mark.parametrize("build", BUILDS)
def test_both_alts_of_one_rsid_stay_one_full_call(tmp_path, build):
    data = _data(tmp_path / "data", "AC", {"A": RSID, "C": RSID})
    report = _analyze(tmp_path, data, build, "same", _rows("AC", (".", ".")), vcf_out=True)
    assert {a["genotype_match"] for a in report["annotations"]} == {"AC"}
    assert _counts(report) == (0, 0, 1)
    assert _round_trip(tmp_path, data, build, "same", report) == [RSID, RSID]


@pytest.mark.parametrize("build", BUILDS)
def test_no_known_alt_leaves_the_unnamed_pair_unchanged(tmp_path, build):
    data = _data(tmp_path / "data", "AC", {"A": None, "C": None})
    report = _analyze(tmp_path, data, build, "none", _rows("AC", (".", ".")), vcf_out=True)
    assert not report["annotations"]
    assert _counts(report) == (0, 0, 1)
    assert _round_trip(tmp_path, data, build, "none", report) == [".", "."]


def _extract(tmp_path, monkeypatch, data, build, rows, rsid, *, indexed):
    source = write_vcf(tmp_path / "extract.vcf", build, rows)
    if indexed:
        gz = source.with_suffix(".vcf.gz")
        pysam.tabix_compress(str(source), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        source = gz
    captured = []
    monkeypatch.setattr(utility, "_render_extract_table", lambda _p, _w, f: captured.append(f))
    args = ["extract", str(source), "--data-dir", str(data), "--snps", rsid]
    result = CliRunner().invoke(main, [*args, "--build", build.lower()])
    assert result.exit_code == 0, (result.output, result.exception)
    return captured[0].get(rsid)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
def test_extract_finds_the_rsid_of_the_one_known_alt(tmp_path, monkeypatch, build, indexed):
    data = _data(tmp_path / "data", "AT", {"A": RSID, "T": None})
    call = _extract(
        tmp_path, monkeypatch, data, build, _rows("AT", (".", ".")), RSID, indexed=indexed
    )
    assert call is not None, "Not 'not in file'"
    assert {call.allele1, call.allele2} == {"A", "T"}


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
def test_extract_finds_the_other_alts_own_rsid(tmp_path, monkeypatch, build, indexed):
    data = _data(tmp_path / "data", "AC", {"A": RSID, "C": OTHER})
    call = _extract(
        tmp_path, monkeypatch, data, build, _rows("AC", (RSID, ".")), OTHER, indexed=indexed
    )
    assert call is not None, "Not 'not in file'"
    assert call.allele1 == "C" and call.is_half_call

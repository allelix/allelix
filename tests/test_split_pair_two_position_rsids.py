"""A pair whose two ALTs resolve by position to two rsIDs is no pair (2.7.1, Route S).

Reviewer A's A-01: an unnamed exact pair (G>A 1/0 + G>C 0/1, ID '.') whose A
is rs999000555 and whose C is rs999000556 merged into one A/C call with no
rsID, because resolution by position names two rsIDs; the GWAS Catalog, keyed
by rsID, never saw it. SUP ruled Route S: two rsIDs found by position are D1
point 1's two different rsIDs, so the records give their ALTs as half calls,
each under its own rsID, to ClinVar and the GWAS Catalog only; the same as
the records with those rsIDs written.
"""

import json
import sqlite3
from contextlib import closing

import pytest
from split_site_contract_helpers import (
    BUILDS,
    OTHER,
    RSID,
    by_source,
    data_dir,
    record,
    run,
    signature,
    write_vcf,
)


def _two_rsid_data(path):
    """A under RSID and C under OTHER in ClinVar, the identity table and the GWAS Catalog."""
    data = data_dir(path, identity=True, alts="AC")
    for name in ("clinvar.GRCh37.sqlite", "clinvar.GRCh38.sqlite"):
        with closing(sqlite3.connect(data / name)) as conn, conn:
            conn.execute("UPDATE clinvar_variants SET rsid=? WHERE alt='C'", (OTHER,))
    with closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn, conn:
        for column in ("alt37", "alt38"):
            conn.execute(f"UPDATE rsid_identity SET rsid=? WHERE {column}='C'", (OTHER,))
    with closing(sqlite3.connect(data / "gwas.sqlite")) as conn, conn:
        conn.execute("UPDATE gwas_associations SET rsid=? WHERE risk_allele='C'", (OTHER,))
    return data


def _pair(ids):
    return [record("A", "1/0", rsid=ids[0]), record("C", "0/1", rsid=ids[1])]


def _counts(report):
    split = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    quality = report["input"]["vcf_quality"]
    return (split[0]["count"] if split else 0, quality["split_site_half_calls"])


@pytest.mark.parametrize("build", BUILDS)
def test_unnamed_pair_gives_each_allele_its_own_rsid(tmp_path, build):
    data = _two_rsid_data(tmp_path / "data")
    report, _ = run(
        data, write_vcf(tmp_path / "in.vcf", build, _pair((".", "."))), tmp_path / "o.json", build
    )
    groups = by_source(report)
    assert {(a["rsid"], a["alt"]) for a in groups["clinvar"]} == {(RSID, "A"), (OTHER, "C")}
    assert {(a["rsid"], a["alt"]) for a in groups["gwas"]} >= {(RSID, "A")}, (
        "The GWAS Catalog finding 2.6.3 gave for rs999000555 A"
    )
    assert not groups["pharmgkb"] and not groups["snpedia"], "No full genotype for either rsID"
    assert all("." in a["genotype_match"] for a in report["annotations"]), "Half calls, not AC"
    assert _counts(report) == (2, 2), "Not a pair: every record counted, both alleles halves"


@pytest.mark.parametrize("build", BUILDS)
def test_unnamed_and_written_rsids_give_the_same_report(tmp_path, build):
    data = _two_rsid_data(tmp_path / "data")
    reports = {}
    for name, ids in (("unnamed", (".", ".")), ("written", (RSID, OTHER))):
        reports[name], _ = run(
            data,
            write_vcf(tmp_path / f"{name}.vcf", build, _pair(ids)),
            tmp_path / f"{name}.json",
            build,
        )
    assert signature(reports["unnamed"]) == signature(reports["written"])
    assert _counts(reports["unnamed"]) == _counts(reports["written"])


@pytest.mark.parametrize("build", BUILDS)
def test_annotated_vcf_reads_back_the_same(tmp_path, build):
    data = _two_rsid_data(tmp_path / "data")
    annotated = tmp_path / "annotated.vcf"
    first, _ = run(
        data,
        write_vcf(tmp_path / "in.vcf", build, _pair((".", "."))),
        tmp_path / "first.json",
        build,
        extra=("--vcf-out", str(annotated)),
    )
    ids = [line.split("\t")[2] for line in annotated.read_text().splitlines() if line[0] != "#"]
    assert ids == [RSID, OTHER], "Each record carries its own allele's rsID"
    second, _ = run(data, annotated, tmp_path / "second.json", build)
    assert signature(second) == signature(first)
    assert json.dumps(second["report_warnings"]) == json.dumps(first["report_warnings"])
    assert second["input"]["total_variants"] == first["input"]["total_variants"] == 2


@pytest.mark.parametrize("build", BUILDS)
def test_one_rsid_for_both_alleles_stays_one_full_call(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=True, alts="AC")
    report, _ = run(
        data, write_vcf(tmp_path / "in.vcf", build, _pair((".", "."))), tmp_path / "o.json", build
    )
    assert {a["genotype_match"] for a in report["annotations"]} == {"AC"}
    assert _counts(report) == (0, 0)

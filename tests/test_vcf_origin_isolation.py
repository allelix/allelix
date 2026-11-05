"""O1: per-record SNPedia findings must not leak to allele-sharing records.

Outside-review seeds retained: AG-only finding on an AA duplicate, and an
identity-rejected duplicate. All caches and calls below are synthetic/offline.
"""

import contextlib
import json
import sqlite3

import pytest

from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.rsid_identity_loader import write_identity_db
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.parsers.vcf import VcfParser
from allelix.reports import _pipeline
from allelix.reports.json_report import render_json
from allelix.reports.vcf import render_vcf

RS = "rs900000001"


def cache(data, second=False):
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        for a, b, repute in [("A", "G", "Bad"), *([("A", "A", "Good")] if second else [])]:
            conn.execute(
                "INSERT INTO snpedia_genotypes "
                "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
                "VALUES (?,?,?,6,?,'Synthetic finding','TEST','plus')",
                (RS, a, b, repute),
            )
        conn.execute(
            "INSERT INTO database_versions "
            "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
            "VALUES ('snpedia','test://synthetic','synthetic','2026-11-02',?,?)",
            (1 + second, f"pv:{_PARSER_VERSION}"),
        )
    write_identity_db(
        [
            (RS, "1", 1000, "G", "A", 1000, "G", "A", 0),
            ("rs900000002", "1", 2000, "G", "A", 2000, "G", "A", 0),
        ],
        data / "rsid_identity.sqlite",
        version="synthetic",
        source_url="test://synthetic",
        downloaded_at="2026-11-02",
    )


def row(gt, *, rsid=RS, padded=False, gvcf=False, ft="PASS"):
    pos, ref, alt = (999, "CG", "CA") if padded else (1000, "G", "A")
    if gvcf:
        alt += ",<NON_REF>"
    return f"1\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT:FT\t{gt}:{ft}"


def run(tmp_path, data, rows, build, *, keep=False):
    src = tmp_path / "input.vcf"
    src.write_text(
        '##fileformat=VCFv4.2\n##ALT=<ID=NON_REF,Description="Non-reference">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + "\n".join(rows) + "\n"
    )
    result = _pipeline.run_analysis(
        src,
        VcfParser(keep_filtered=keep),
        [SNPediaAnnotator(data)],
        identity=RsidIdentityAnnotator(data),
        build_override=build,
        collect_vcf_write_inputs=True,
    )
    wi = result.vcf_write_inputs
    assert wi is not None
    out = tmp_path / "out.vcf"
    render_vcf(
        input_path=src,
        output_path=out,
        build=result.build,
        annotators_used=result.annotators_used,
        annotations_by_key=wi.annotations_by_key,
        recovered_rsids=wi.recovered_rsids,
        reference_rows=wi.reference_rows,
        record_rows=wi.record_rows,
        record_recovered=wi.record_recovered,
        keep_filtered=keep,
    )
    report = tmp_path / "out.json"
    render_json(result, output_path=report, min_magnitude=0, source_min_magnitudes={"snpedia": 0})
    payload = json.loads(report.read_text())
    assert sorted(a["significance"] for a in payload["annotations"]) == sorted(
        a.significance for a in result.annotations
    )
    actual = []
    for line in out.read_text().splitlines():
        if not line.startswith("#"):
            fields = line.split("\t")
            info = dict(v.split("=", 1) for v in fields[7].split(";") if "=" in v)
            actual.append(info.get("ALLELIX_SNPEDIA", ".").split(",")[0])
    expected = []
    for index in range(len(rows)):
        own = [a for a in result.annotations if a.record_index == index and a.source == "snpedia"]
        assert len(own) <= 1
        expected.append(f"{own[0].rsid}|{own[0].significance}" if own else ".")
    return result, actual, expected


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("other", ["1/1", "1|1", "1", "./.", "0/."])
@pytest.mark.parametrize("shape", ["plain", "padded", "gvcf", "unnamed"])
@pytest.mark.parametrize("batch", [1, 5000])
def test_genotype_specific_info_follows_its_record(
    tmp_path, monkeypatch, build, reverse, other, shape, batch
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch)
    data = tmp_path / "data"
    cache(data)
    kw = dict(
        padded=shape == "padded", gvcf=shape == "gvcf", rsid="." if shape == "unnamed" else RS
    )
    rows = [row("0|1", **kw), row(other, **kw)]
    if reverse:
        rows.reverse()
    result, actual, expected = run(tmp_path, data, rows, build)
    assert len(result.annotations) == 1, "positive AG seed must actually annotate"
    assert actual == expected, "an INFO finding belongs only to its originating input record"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("batch", [1, 5000])
@pytest.mark.parametrize("case", ["two-findings", "identity-rejected", "ft-rejected", "ft-kept"])
def test_distinct_findings_and_exclusions_do_not_pool(
    tmp_path, monkeypatch, build, reverse, batch, case
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch)
    data = tmp_path / "data"
    cache(data, second=case == "two-findings")
    second = row("1/1")
    if case == "identity-rejected":
        second = row("0/1", rsid="rs900000002")
    elif case.startswith("ft-"):
        second = row("1/1", ft="LowAD")
    rows = [row("0/1"), second]
    if reverse:
        rows.reverse()
    result, actual, expected = run(tmp_path, data, rows, build, keep=case == "ft-kept")
    assert len(result.annotations) == (2 if case == "two-findings" else 1)
    if case == "identity-rejected":
        assert any(
            w.code == "identity_locus_conflict" and w.count == 1 for w in result.report_warnings
        )
    assert actual == expected, "source/rank selection must be per record, even across batches"

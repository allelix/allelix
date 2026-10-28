# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Reference-derived ClinVar INFO and record-local stamping eligibility.

All calls and caches are synthetic. The cache layout follows the regular
ClinVar test fixtures. The proposed VCF 0.6.0 contract adds a scalar
ALLELIX_CLINVAR_REF (significance|rsID|gene), for the highest-ranked
reference-derived finding on that record. ALLELIX_CLINVAR remains Number=A
and describes only carried ALT alleles. Both fields may occur together.
"""

from __future__ import annotations

import contextlib
import gzip
import json
import sqlite3
from typing import TYPE_CHECKING

import pysam
import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import (
    ALPHAMISSENSE_SCHEMA_VERSION,
    CADD_SCHEMA_VERSION,
    CLINVAR_INTERPRETER_VERSION,
    GNOMAD_SCHEMA_VERSION,
)
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import (
    ALPHAMISSENSE_SCHEMA,
    CADD_SCHEMA,
    CLINVAR_SCHEMA,
    GNOMAD_SCHEMA,
)

if TYPE_CHECKING:
    from pathlib import Path

BUILDS = ("GRCh37", "GRCh38")
REF_FIELD = "ALLELIX_CLINVAR_REF"
ALT_FIELD = "ALLELIX_CLINVAR"
# Synthetic REF == ALT sites on both builds, spanning all four bases.
SITES = (
    (1000, "rs900252001", "A", "G"),
    (2000, "rs900252002", "C", "T"),
    (3000, "rs900252003", "G", "A"),
    (4000, "rs900252004", "T", "C"),
)
META = """##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=chr1>
##ALT=<ID=NON_REF,Description="Any other allele">
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=KEEP,Number=1,Type=String,Description="Unrelated input information">
##FILTER=<ID=LowQual,Description="Caller rejected">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
"""


def _insert(conn, rsid, pos, ref, alt, *, gene="REFGENE", sig="Pathogenic", allele=1):
    conn.execute(
        "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
        "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
        "submitters_same, submitters_total) VALUES (?, '1', ?, ?, ?, ?, "
        "'Synthetic condition', ?, 'reviewed_by_expert_panel', ?, ?, 1, 1)",
        (rsid, pos, ref, alt, sig, gene, allele, sig.replace("_", " ")),
    )


@pytest.fixture
def data(tmp_path: Path) -> Path:
    directory = tmp_path / "data"
    directory.mkdir()
    for build in BUILDS:
        with contextlib.closing(sqlite3.connect(directory / clinvar_db_filename(build))) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            for pos, rsid, ref, _alt in SITES:
                _insert(conn, rsid, pos, ref, ref, allele=pos)
            # Ordinary alternate-allele finding, with no REF-derived finding.
            _insert(conn, "rs900252005", 5000, "G", "A", gene="ALTGENE", allele=5)
            # Both reference and alternate findings at one multiallelic site.
            _insert(conn, "rs900252006", 6000, "T", "T", allele=6)
            _insert(
                conn,
                "rs900252006",
                6000,
                "T",
                "C",
                gene="ALTGENE",
                sig="Likely_pathogenic",
                allele=7,
            )
            # Same reference allele, lower-ranked SCV: output must pick the top row.
            _insert(
                conn,
                "rs900252006",
                6000,
                "T",
                "T",
                gene="LOWER",
                sig="Uncertain_significance",
                allele=6,
            )
            # The published F5 placement, with invented clinical submission content.
            pos, ref = (169519049, "T") if build == "GRCh37" else (169549811, "C")
            _insert(conn, "rs6025", pos, ref, "T", gene="F5", allele=8)
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
                "record_count, local_version_tag) VALUES (?, 'test://synthetic', "
                "'20261025', '2026-10-25', 9, ?)",
                (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    return directory


def _row(pos, rsid, ref, alt, gt, *, gq=99, dp=20, filt="PASS", second=None, info="KEEP=x"):
    sample = f"{gt}:{gq}:{dp}" + (f"\t{second}" if second is not None else "")
    return f"chr1\t{pos}\t{rsid}\t{ref}\t{alt}\t50\t{filt}\t{info}\tGT:GQ:DP\t{sample}\n"


def _input(tmp_path: Path, rows, *, samples="S", compressed=False, owned=False) -> Path:
    path = tmp_path / ("input.vcf.gz" if compressed else "input.vcf")
    stale = (
        (
            '##INFO=<ID=ALLELIX_CLINVAR_REF,Number=1,Type=String,Description="Old">\n'
            '##INFO=<ID=ALLELIX_CLINVAR,Number=A,Type=String,Description="Old">\n'
            "##ALLELIX_VCF_SchemaVersion=0.5.0\n"
        )
        if owned
        else ""
    )
    meta = (
        META
        if any("<NON_REF>" in row for row in rows)
        else "\n".join(
            line for line in META.split("\n") if not line.startswith("##ALT=<ID=NON_REF,")
        )
    )
    text = meta + stale + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    text += samples + "\n" + "".join(rows)
    if compressed:
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _run(data: Path, src: Path, tmp_path: Path, build: str, *options: str, name="out"):
    out, report = tmp_path / f"{name}.vcf", tmp_path / f"{name}.json"
    args = [
        "analyze",
        str(src),
        "--no-update",
        "--data-dir",
        str(data),
        "--build",
        build.lower(),
        "--output",
        str(report),
        "--vcf-out",
        str(out),
        *options,
    ]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    body = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    with pysam.VariantFile(str(out)) as stream:
        header = stream.header.copy()
        decoded = [dict(record.info) for record in stream]
    return json.loads(report.read_text()), body, decoded, header, out


def _value(rsid, gene="REFGENE", sig="clinvar_pathogenic"):
    return f"{sig}|{rsid}|{gene}"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gvcf", [False, True], ids=["vcf", "gvcf"])
@pytest.mark.parametrize("gt", ["0/1", "0|1", "1|0", "0/.", ".|0", "0/0", "0|0", "0"])
def test_reference_findings_have_a_separate_scalar_field(data, tmp_path, build, gvcf, gt):
    rows = [
        _row(pos, rsid, ref, alt + (",<NON_REF>" if gvcf else ""), gt)
        for pos, rsid, ref, alt in SITES
    ]
    src = _input(tmp_path, rows, compressed=gvcf)
    report, body, info, _header, _out = _run(data, src, tmp_path, build)
    expected = {_value(rsid) for _pos, rsid, _ref, _alt in SITES}
    assert {
        f"{a['significance']}|{a['rsid']}|{a['gene']}" for a in report["annotations"]
    } == expected
    assert all(a["magnitude"] == 9.0 for a in report["annotations"])
    assert report["input"]["gvcf"] is gvcf
    if "." in gt:
        assert all(a["zygosity"] == "half_call" for a in report["annotations"])
    assert [item.get(REF_FIELD) for item in info] == [
        _value(rsid) for _pos, rsid, _ref, _alt in SITES
    ]
    assert all(ALT_FIELD not in item for item in info), "REF findings must never occupy ALT slots"
    assert [row[:7] + row[8:] for row in body] == [
        row.rstrip("\n").split("\t")[:7] + row.rstrip("\n").split("\t")[8:] for row in rows
    ]
    assert all(item["KEEP"] == "x" for item in info)


@pytest.mark.parametrize("gt", ["0/1", "0|1", "0/.", ".|0"])
def test_factor_v_leiden_reference_finding_matches_json(data, tmp_path, gt):
    src = _input(tmp_path, [_row(169519049, "rs6025", "T", "C", gt)])
    report, _body, info, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert [(a["rsid"], a["significance"], a["magnitude"]) for a in report["annotations"]] == [
        ("rs6025", "clinvar_pathogenic", 9.0)
    ]
    assert info[0].get(REF_FIELD) == _value("rs6025", "F5")
    assert ALT_FIELD not in info[0]


def test_vcf_schema_declares_reference_field_and_alt_only_semantics(data, tmp_path):
    src = _input(tmp_path, [_row(6000, "rs900252006", "T", "C,G,<NON_REF>", "0/1")])
    _report, _body, _info, header, out = _run(data, src, tmp_path, "GRCh37")
    assert "##ALLELIX_VCF_SchemaVersion=0.6.0\n" in out.read_text()
    assert REF_FIELD in header.info
    declaration = header.info[REF_FIELD]
    assert declaration.number == 1 and declaration.type == "String"
    description = declaration.description.lower()
    assert "reference" in description
    for component in ("significance", "rsid", "gene"):
        assert component in description
    assert header.info[ALT_FIELD].number == "A"
    assert "homozygous-reference" not in header.info[ALT_FIELD].description.lower()


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt,alt", [("0/1", "C,G,<NON_REF>"), ("0|2", "G,C,<NON_REF>")])
def test_reference_and_alternate_findings_coexist_without_slot_collision(
    data, tmp_path, build, gt, alt
):
    src = _input(tmp_path, [_row(6000, "rs900252006", "T", alt, gt)])
    report, _body, info, _header, _out = _run(data, src, tmp_path, build)
    assert {a["alt"] for a in report["annotations"]} == {"T", "C"}
    expected_alt = _value("rs900252006", "ALTGENE", "clinvar_likely_pathogenic")
    assert info[0].get(ALT_FIELD) == tuple(
        expected_alt if a == "C" else "." for a in alt.split(",")
    )
    assert info[0].get(REF_FIELD) == _value("rs900252006")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True], ids=["eligible-first", "eligible-last"])
@pytest.mark.parametrize(
    "bad_gt,gq,dp,filt",
    [
        ("0/1", 3, 20, "PASS"),
        ("./.", 99, 20, "PASS"),
        (".|.", 99, 20, "PASS"),
        ("0/0", 99, 20, "PASS"),
        ("0|0", 99, 20, "PASS"),
        ("2/2", 99, 20, "PASS"),
        ("0/1", 99, 0, "PASS"),
        ("0/1", 99, 20, "LowQual"),
    ],
    ids=[
        "low-gq",
        "no-call",
        "phased-no-call",
        "hom-ref",
        "phased-hom-ref",
        "other-alt",
        "zero-depth",
        "filtered",
    ],
)
def test_duplicate_record_does_not_borrow_an_alt_finding(
    data, tmp_path, build, reverse, bad_gt, gq, dp, filt
):
    good = _row(5000, "rs900252005", "G", "A,T,<NON_REF>", "0/1")
    bad = _row(5000, "rs900252005", "G", "A,T,<NON_REF>", bad_gt, gq=gq, dp=dp, filt=filt)
    rows = [bad, good] if reverse else [good, bad]
    report, body, info, _header, _out = _run(data, _input(tmp_path, rows), tmp_path, build)
    assert any(a["rsid"] == "rs900252005" for a in report["annotations"])
    assert info[int(reverse)].get(ALT_FIELD) == (_value("rs900252005", "ALTGENE"), ".", ".")
    assert not any(k.startswith("ALLELIX_") for k in info[int(not reverse)])
    assert all(REF_FIELD not in item for item in info)
    assert [r[8:] for r in body] == [r.rstrip("\n").split("\t")[8:] for r in rows]


@pytest.mark.parametrize("bad_gt,gq", [("0/1", 3), ("./.", 99), ("1/1", 99)])
@pytest.mark.parametrize("reverse", [False, True])
def test_duplicate_record_does_not_borrow_a_reference_finding(data, tmp_path, bad_gt, gq, reverse):
    good = _row(4000, "rs900252004", "T", "C", "0|1")
    bad = _row(4000, "rs900252004", "T", "C", bad_gt, gq=gq)
    rows = [bad, good] if reverse else [good, bad]
    _report, _body, info, _header, _out = _run(data, _input(tmp_path, rows), tmp_path, "GRCh37")
    assert info[int(reverse)].get(REF_FIELD) == _value("rs900252004")
    assert not any(k.startswith("ALLELIX_") for k in info[int(not reverse)])
    assert all(ALT_FIELD not in item for item in info)


@pytest.mark.parametrize("sample", ["FIRST", "SECOND"])
@pytest.mark.parametrize("min_gq,expected", [(20, [True, False]), (0, [True, True])])
def test_selected_sample_and_configured_quality_decide_each_record(
    data, tmp_path, sample, min_gq, expected
):
    selected = ["0|1:99:20", "0|1:3:20"]
    rows = []
    for call in selected:
        first, second = (call, "0/0:99:20") if sample == "FIRST" else ("0/0:99:20", call)
        gt, gq, dp = first.split(":")
        rows.append(_row(5000, "rs900252005", "G", "A", gt, gq=int(gq), dp=int(dp), second=second))
    src = _input(tmp_path, rows, samples="FIRST\tSECOND")
    _report, _body, info, _header, _out = _run(
        data, src, tmp_path, "GRCh37", "--sample", sample, "--min-gq", str(min_gq)
    )
    assert [ALT_FIELD in item for item in info] == expected


@pytest.mark.parametrize("keep", [False, True])
def test_keep_filtered_does_not_bypass_sample_quality(data, tmp_path, keep):
    rows = [
        _row(5000, "rs900252005", "G", "A", "0/1"),
        _row(5000, "rs900252005", "G", "A", "0/1", filt="LowQual"),
        _row(5000, "rs900252005", "G", "A", "0/1", filt="LowQual", gq=3),
    ]
    options = ("--keep-filtered",) if keep else ()
    _report, _body, info, _header, _out = _run(
        data, _input(tmp_path, rows), tmp_path, "GRCh37", *options
    )
    assert [ALT_FIELD in item for item in info] == [True, keep, False]
    if keep:
        assert info[1]["ALLELIX_FILTER"] == "LowQual"
    assert not any(k.startswith("ALLELIX_") for k in info[2])


@pytest.mark.parametrize("gt", ["1/1", "1|1", "1/.", "./1", "./."])
def test_no_reference_allele_means_no_reference_finding(data, tmp_path, gt):
    src = _input(tmp_path, [_row(4000, "rs900252004", "T", "C,<NON_REF>", gt)])
    report, _body, info, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert report["annotations"] == []
    assert REF_FIELD not in info[0] and ALT_FIELD not in info[0]


@pytest.mark.parametrize("gt", ["0/1", "0|1", "1/."])
def test_grch38_factor_v_leiden_remains_an_alt_finding(data, tmp_path, gt):
    src = _input(tmp_path, [_row(169549811, "rs6025", "C", "T,<NON_REF>", gt)])
    report, _body, info, _header, _out = _run(data, src, tmp_path, "GRCh38")
    assert [(a["rsid"], a["alt"]) for a in report["annotations"]] == [("rs6025", "T")]
    assert info[0].get(ALT_FIELD) == (_value("rs6025", "F5"), ".")
    assert REF_FIELD not in info[0]


def test_reannotation_replaces_both_fields_without_leaking_stale_stamps(data, tmp_path):
    stale = "KEEP=x;ALLELIX_CLINVAR_REF=STALE_REF;ALLELIX_CLINVAR=STALE_ALT"
    rows = [
        _row(4000, "rs900252004", "T", "C", "0/1", info=stale),
        _row(4000, "rs900252004", "T", "C", "./.", info=stale),
    ]
    first = _run(data, _input(tmp_path, rows, owned=True), tmp_path, "GRCh37", name="first")
    second = _run(data, first[4], tmp_path, "GRCh37", name="second")
    assert first[1:3] == second[1:3]
    assert second[2][0].get(REF_FIELD) == _value("rs900252004")
    assert second[2][1] == {"KEEP": "x"}
    text = second[4].read_text()
    assert "STALE_" not in text
    assert text.count(f"##INFO=<ID={REF_FIELD},") == 1


@pytest.mark.parametrize("pos,ref,alt", [(4000, "TC", "T"), (3999, "AT", "A")])
def test_multibase_reference_finding_is_explicit_not_an_alt_finding(data, tmp_path, pos, ref, alt):
    src = _input(tmp_path, [_row(pos, ".", ref, alt, "0/0")])
    report, body, info, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert [(a["rsid"], a["genotype_match"]) for a in report["annotations"]] == [
        ("rs900252004", "TT")
    ]
    assert body[0][3:5] == [ref, alt]
    assert info[0].get(REF_FIELD) == _value("rs900252004")
    assert ALT_FIELD not in info[0]


@pytest.mark.parametrize("gt", ["0/1", "0/.", "0/0"])
def test_wrong_reference_base_cannot_create_a_reference_finding(data, tmp_path, gt):
    src = _input(tmp_path, [_row(4000, "rs900252004", "A", "C", gt)])
    report, _body, info, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert report["annotations"] == []
    assert REF_FIELD not in info[0] and ALT_FIELD not in info[0]


@pytest.fixture
def enriched_data(data):
    (data / "config.toml").write_text("[sources]\ncadd = true\n", encoding="utf-8")
    write_identity_db(
        [("rs6025", "1", 169519049, "T", "T", 169549811, "C", "T", 0)],
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://identity",
        downloaded_at="2026-10-25",
    )
    sites = [("1", 169549811, "C", "T", "rs6025"), ("1", 5000, "G", "A", "rs900252005")]
    for name, schema, version in (
        ("gnomad", GNOMAD_SCHEMA, GNOMAD_SCHEMA_VERSION),
        ("alphamissense", ALPHAMISSENSE_SCHEMA, ALPHAMISSENSE_SCHEMA_VERSION),
        ("cadd", CADD_SCHEMA, CADD_SCHEMA_VERSION),
    ):
        with contextlib.closing(sqlite3.connect(data / f"{name}.sqlite")) as conn:
            conn.executescript(schema)
            if name == "gnomad":
                conn.executemany(
                    "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
                    "VALUES (?,?,?,?,?,0.01)",
                    sites,
                )
            elif name == "alphamissense":
                conn.executemany(
                    "INSERT INTO alphamissense_scores VALUES "
                    "(?, ?, ?, ?, ?, 'P1', 'T1', 'A1V', 0.91, 'likely_pathogenic')",
                    sites,
                )
            else:
                conn.executemany(
                    "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?,?,?,?,23.5)",
                    [site[:4] for site in sites],
                )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
                "record_count, local_version_tag) "
                "VALUES (?, 'test://synthetic', '1', '2026-10-25', 2, ?)",
                (name, f"sv:{version}"),
            )
            conn.commit()
    return data


@pytest.mark.parametrize("gt", ["0/0", "0|1", "0/."])
def test_reference_enrichment_never_occupies_alt_slots(enriched_data, tmp_path, gt):
    src = _input(tmp_path, [_row(169519049, "rs6025", "T", "C,G,<NON_REF>", gt)])
    report, _body, info, _header, _out = _run(enriched_data, src, tmp_path, "GRCh37")
    assert len(report["annotations"]) == 1
    annotation = report["annotations"][0]
    assert annotation["allele_frequency"] == 0.01
    assert annotation["am_pathogenicity"] == 0.91
    assert annotation["cadd_phred"] == 23.5
    assert not (
        {ALT_FIELD, "ALLELIX_GNOMAD_AF", "ALLELIX_ALPHAMISSENSE", "ALLELIX_CADD"} & info[0].keys()
    )
    assert info[0].get(REF_FIELD) == _value("rs6025", "F5")


@pytest.mark.parametrize("bad_gt,gq", [("0/1", 3), ("./.", 99), ("0/0", 99)])
def test_enrichment_also_respects_source_record_eligibility(enriched_data, tmp_path, bad_gt, gq):
    rows = [
        _row(5000, "rs900252005", "G", "A", "0/1"),
        _row(5000, "rs900252005", "G", "A", bad_gt, gq=gq),
    ]
    _report, _body, info, _header, _out = _run(
        enriched_data, _input(tmp_path, rows), tmp_path, "GRCh38"
    )
    assert info[0]["ALLELIX_GNOMAD_AF"] == pytest.approx((0.01,))
    assert info[0]["ALLELIX_ALPHAMISSENSE"] == ("likely_pathogenic|0.91",)
    assert info[0]["ALLELIX_CADD"] == (23.5,)
    assert not any(k.startswith("ALLELIX_") for k in info[1])


@pytest.mark.parametrize("bad_gt,gq", [("0/1", 3), ("./.", 99), ("0/0", 99)])
def test_recovered_id_and_finding_do_not_spill_to_an_ineligible_record(data, tmp_path, bad_gt, gq):
    rows = [
        _row(5000, ".", "G", "A", "0/1"),
        _row(5000, ".", "G", "A", bad_gt, gq=gq),
    ]
    _report, body, info, _header, _out = _run(data, _input(tmp_path, rows), tmp_path, "GRCh38")
    assert body[0][2] == "rs900252005"
    assert info[0][ALT_FIELD] == (_value("rs900252005", "ALTGENE"),)
    assert body[1][2] == "."
    assert not any(k.startswith("ALLELIX_") for k in info[1])


def test_two_eligible_duplicate_records_both_keep_their_findings(data, tmp_path):
    rows = [
        _row(5000, "rs900252005", "G", "A", "0/1"),
        _row(5000, "rs900252005", "G", "A", "1|1"),
    ]
    _report, body, info, _header, _out = _run(data, _input(tmp_path, rows), tmp_path, "GRCh38")
    assert [item.get(ALT_FIELD) for item in info] == [(_value("rs900252005", "ALTGENE"),)] * 2
    assert [row[-1] for row in body] == ["0/1:99:20", "1|1:99:20"]


@pytest.mark.parametrize(
    "source,field",
    [
        ("clinvar", ALT_FIELD),
        ("gwas", "ALLELIX_GWAS"),
        ("pharmgkb", "ALLELIX_CLINPGX"),
        ("snpedia", "ALLELIX_SNPEDIA"),
    ],
)
def test_each_primary_source_obeys_the_same_record_eligibility(
    data, tmp_path, monkeypatch, source, field
):
    # Inject a payload after real parsing/matching, to isolate the shared
    # writer dispatch from each source's separate biological matching rules.
    from allelix.cli import _helpers

    original = _helpers.run_analysis

    def with_source_payload(*args, **kwargs):
        result = original(*args, **kwargs)
        for annotation in result.annotations:
            annotation.source = source
            annotation.trait = "Synthetic trait"
            annotation.p_value = 1e-9
        return result

    monkeypatch.setattr(_helpers, "run_analysis", with_source_payload)
    rows = [
        _row(5000, "rs900252005", "G", "A", "0/1"),
        _row(5000, "rs900252005", "G", "A", "0/1", gq=3),
    ]
    _report, _body, info, _header, _out = _run(data, _input(tmp_path, rows), tmp_path, "GRCh38")
    assert field in info[0]
    assert not any(k.startswith("ALLELIX_") for k in info[1])


@pytest.mark.parametrize(
    "pos,ref,alt,info",
    [
        (3999, "A", "<NON_REF>", "END=4001;KEEP=x"),
        (4000, "T", "<NON_REF>", "END=4001;KEEP=x"),
    ],
)
def test_end_reference_blocks_keep_their_existing_exclusion(data, tmp_path, pos, ref, alt, info):
    """Span summaries retain their documented JSON-only findings."""
    src = _input(tmp_path, [_row(pos, ".", ref, alt, "0/0", info=info)])
    report, body, decoded, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert [(a["rsid"], a["genotype_match"]) for a in report["annotations"]] == [
        ("rs900252004", "TT")
    ]
    assert body[0][7] == info
    assert not any(k.startswith("ALLELIX_") for k in decoded[0])


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0|1", "0/."])
def test_requesting_vcf_output_does_not_change_json_findings(data, tmp_path, build, gt):
    rows = [_row(pos, rsid, ref, alt, gt) for pos, rsid, ref, alt in SITES]
    src = _input(tmp_path, rows)
    with_vcf, _body, _info, _header, _out = _run(data, src, tmp_path, build)
    plain = tmp_path / "without-vcf.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(src),
            "--no-update",
            "--data-dir",
            str(data),
            "--build",
            build.lower(),
            "--output",
            str(plain),
        ],
    )
    assert result.exit_code == 0, result.output
    without_vcf = json.loads(plain.read_text())
    assert with_vcf["annotations"] == without_vcf["annotations"]
    assert with_vcf["input"]["vcf_quality"] == without_vcf["input"]["vcf_quality"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0/0", "0|0", "0", "0/.", ".|0"])
def test_single_position_altless_record_carries_its_reference_finding(data, tmp_path, build, gt):
    rows = [_row(pos, rsid, ref, ".", gt) for pos, rsid, ref, _alt in SITES]
    report, body, info, header, _out = _run(data, _input(tmp_path, rows), tmp_path, build)
    assert {a["rsid"] for a in report["annotations"]} == {site[1] for site in SITES}
    assert all(a["magnitude"] == 9.0 for a in report["annotations"])
    if "." in gt:
        assert all(a["zygosity"] == "half_call" for a in report["annotations"])
    assert [item.get(REF_FIELD) for item in info] == [
        _value(rsid) for _pos, rsid, _ref, _alt in SITES
    ]
    assert header.info[REF_FIELD].number == 1
    assert all(ALT_FIELD not in item for item in info)
    assert all(row[4] == "." for row in body)
    assert [row[-1] for row in body] == [f"{gt}:99:20"] * len(SITES)


def test_factor_v_leiden_altless_homozygous_reference_record_is_stamped(data, tmp_path):
    src = _input(tmp_path, [_row(169519049, "rs6025", "T", ".", "0/0")])
    report, body, info, _header, _out = _run(data, src, tmp_path, "GRCh37")
    assert [(a["rsid"], a["genotype_match"], a["magnitude"]) for a in report["annotations"]] == [
        ("rs6025", "TT", 9.0)
    ]
    assert body[0][4] == "."
    assert info[0].get(REF_FIELD) == _value("rs6025", "F5")
    assert ALT_FIELD not in info[0]

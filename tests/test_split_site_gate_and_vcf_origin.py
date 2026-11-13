"""Split-site decisions retain rejected records and each ALT's output origin."""

import json
import sqlite3
from contextlib import closing

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _vcf
from test_split_person_genotype import SOURCES, _calls, _data, _row
from validated_call_helpers import BUILDS, RSID, analyze, identity

from allelix.cli import main


def _gate(path, row, kind):
    lines = path.read_text().splitlines()
    cols = lines[2 + row].split("\t")
    if kind == "FILTER":
        cols[6] = "q10"
    else:
        cols[8] = "GT:FT"
        cols[9] += ":q10"
    lines[2 + row] = "\t".join(cols)
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("kind", ("FILTER", "FT"))
@pytest.mark.parametrize("shape", ("pair", "third_alt"))
def test_rejected_split_record_cannot_leave_a_partial_or_combined_finding(
    tmp_path, build, reverse, kind, shape
):
    data = tmp_path / "data"
    pos = _data(data, build)
    identity(data, "ACGT")
    rows = [_row(pos, "A", "1/0"), _row(pos, "T", "0/1")]
    bad = 1
    if shape == "third_alt":
        rows.append(_row(pos, "C", "0/1"))
        bad = 2
    if reverse:
        rows.reverse()
        bad = len(rows) - 1 - bad
    path = _vcf(tmp_path / "input.vcf", rows)
    _gate(path, bad, kind)
    result = analyze(data, path, build, SOURCES)
    assert not result.annotations, "A rejected component cannot leave AG, GT or inferred AT"
    assert any(w.code == "split_site_calls" and w.count for w in result.report_warnings)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("kind", ("FILTER", "FT"))
def test_keep_filtered_restores_the_exact_complementary_pair(tmp_path, build, reverse, kind):
    from allelix.parsers.vcf import VcfParser

    data = tmp_path / "data"
    pos = _data(data, build)
    rows = [_row(pos, "A", "1/0"), _row(pos, "T", "0/1")]
    if reverse:
        rows.reverse()
    path = _vcf(tmp_path / "input.vcf", rows)
    _gate(path, 1, kind)
    parser = VcfParser(keep_filtered=True)
    result = analyze(data, path, build, SOURCES, parser=parser)
    assert _calls(result) == [("pharmgkb", "AT"), ("snpedia", "AT")]
    assert not any(w.code == "split_site_calls" for w in result.report_warnings)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
def test_each_combined_alt_finding_is_written_on_its_own_input_record(tmp_path, build, reverse):
    data = tmp_path / "data"
    pos = _data(data, build)
    with closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn, conn:
        conn.execute(
            "INSERT INTO clinvar_variants "
            "(rsid,chromosome,position,ref,alt,gene,clinical_significance,condition,"
            "review_status,allele_id,aggregate,submitters_same,submitters_total) "
            "SELECT rsid,chromosome,position,ref,'T',gene,clinical_significance,condition,"
            "review_status,allele_id+100,aggregate,submitters_same,submitters_total "
            "FROM clinvar_variants WHERE rsid=? AND alt='A'",
            (RSID,),
        )
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\nsnpedia=true\npharmgkb=true\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=true\n"
    )
    rows = [_row(pos, "A", "1/0"), _row(pos, "T", "0/1")]
    source = _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows)
    report, output = tmp_path / "report.json", tmp_path / "output.vcf"
    run = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            build.lower(),
            "--min-magnitude",
            "0",
            "--output",
            str(report),
            "--vcf-out",
            str(output),
        ],
    )
    assert run.exit_code == 0, (run.output, run.exception)
    findings = json.loads(report.read_text())["annotations"]
    assert sorted(
        (a["alt"], a["genotype_match"]) for a in findings if a["source"] == "clinvar"
    ) == [("A", "AT"), ("T", "AT")]
    written = [
        line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")
    ]
    assert [row[:7] for row in written] == [
        line.split("\t")[:7]
        for line in source.read_text().splitlines()
        if not line.startswith("#")
    ]
    assert {row[4] for row in written if "ALLELIX_CLINVAR=" in row[7]} == {"A", "T"}, (
        "Both carried ALT findings must survive in their own record's INFO, in either row order"
    )


@pytest.mark.parametrize("build", BUILDS)
def test_a_reference_record_at_another_locus_ends_the_consecutive_run(tmp_path, build):
    data = tmp_path / "data"
    pos = _data(data, build)
    path = _vcf(
        tmp_path / "unsorted.vcf",
        [_row(pos, "A", "1/0"), _row(pos + 100, "C", "0/0", "."), _row(pos, "T", "0/1")],
    )
    result = analyze(data, path, build, SOURCES)
    assert _calls(result) == [
        ("pharmgkb", "AG"),
        ("pharmgkb", "GT"),
        ("snpedia", "AG"),
        ("snpedia", "GT"),
    ], "Non-adjacent split records are outside the documented consecutive-run rule"

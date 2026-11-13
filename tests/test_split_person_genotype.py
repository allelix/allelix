"""Only the exact complementary split pair represents one two-ALT genotype."""

import json
import re
import sqlite3
from contextlib import closing

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, identity, setup

from allelix.cli import main

SOURCES = ("snpedia", "pharmgkb")


def _data(path, build):
    pos = setup(path, build)
    identity(path, "AGT")
    with closing(sqlite3.connect(path / "snpedia.sqlite")) as conn, conn:
        conn.execute("DELETE FROM snpedia_genotypes")
        for call in ("AG", "GT", "AT", "AA", "TT", "GG"):
            conn.execute(
                "INSERT INTO snpedia_genotypes "
                "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
                "VALUES (?,?,?,6,'Bad',?,'SYNGENE','plus')",
                (RSID, *call, f"Synthetic {call}"),
            )
    with closing(sqlite3.connect(path / "pharmgkb.sqlite")) as conn, conn:
        conn.execute("DELETE FROM pharmgkb_annotations")
        for n, call in enumerate(("AG", "GT", "AT", "AA", "TT", "GG")):
            conn.execute(
                "INSERT INTO pharmgkb_annotations "
                "(rsid,genotype,gene,drugs,phenotype,phenotype_category,annotation_text,"
                "level_of_evidence,score,pgkb_annotation_id,function_class,is_nonfinding) "
                "VALUES (?,?,'SYNGENE','Synthetic drug','Synthetic response','Efficacy',"
                "?,'1A',10,?,'decreased',0)",
                (RSID, call, f"Synthetic {call}", str(900450000 + n)),
            )
    return pos


def _row(pos, alt, gt, rsid=RSID):
    return rsid, "1", pos, "G", alt, gt, "."


def _calls(result):
    return sorted((a.source, a.genotype_match) for a in result.annotations)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("swap_alts", (False, True))
@pytest.mark.parametrize("named", (False, True))
def test_exact_complementary_pair_matches_one_merged_two_alt_call(
    tmp_path, build, reverse, swap_alts, named
):
    data = tmp_path / "data"
    pos = _data(data, build)
    rid = RSID if named else "."
    merged = analyze(
        data, _vcf(tmp_path / "merged.vcf", [_row(pos, "A,T", "1/2", rid)]), build, SOURCES
    )
    assert _calls(merged) == [("pharmgkb", "AT"), ("snpedia", "AT")]
    alts = ("T", "A") if swap_alts else ("A", "T")
    rows = [_row(pos, alts[0], "1/0", rid), _row(pos, alts[1], "0/1", rid)]
    result = analyze(
        data, _vcf(tmp_path / "split.vcf", rows[::-1] if reverse else rows), build, SOURCES
    )
    assert _calls(result) == _calls(merged), "AG and GT are not this person's AT genotype"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize(
    "first,second",
    (
        ("1/0", "1/0"),
        ("0/1", "0/1"),
        ("1/1", "0/1"),
        ("1/0", "1/1"),
        ("1/.", "0/1"),
        ("1/0", "./1"),
        ("1", "0/1"),
    ),
)
def test_noncomplementary_split_calls_abstain_and_are_counted(
    tmp_path, build, reverse, first, second
):
    data = tmp_path / "data"
    pos = _data(data, build)
    rows = [_row(pos, "A", first), _row(pos, "T", second)]
    result = analyze(
        data, _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows), build, SOURCES
    )
    assert not result.annotations, "Do not infer either partial genotype as the person's call"
    assert any(w.count and w.count > 0 for w in result.report_warnings), (
        "The split-call abstention must be counted in the report, not silently dropped"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
def test_a_third_distinct_alt_cannot_be_folded_into_a_diploid_pair(tmp_path, build, reverse):
    data = tmp_path / "data"
    pos = _data(data, build)
    identity(data, "ACGT")
    rows = [_row(pos, "A", "1/0"), _row(pos, "T", "0/1"), _row(pos, "C", "0/1")]
    result = analyze(
        data, _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows), build, SOURCES
    )
    assert not result.annotations
    assert any(w.count and w.count > 0 for w in result.report_warnings)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ("single", "duplicate", "merged", "different_loci"))
def test_scope_controls_do_not_treat_every_pair_of_records_as_a_split_call(tmp_path, build, shape):
    data = tmp_path / "data"
    pos = _data(data, build)
    rows = [_row(pos, "A", "0/1")]
    expected = [("pharmgkb", "AG"), ("snpedia", "AG")]
    if shape == "duplicate":
        rows *= 2
        expected = sorted(expected * 2)
    elif shape == "merged":
        rows = [_row(pos, "A,T", "1/2")]
        expected = [("pharmgkb", "AT"), ("snpedia", "AT")]
    elif shape == "different_loci":
        # Another, unnamed non-target record must not contaminate the site's call.
        rows += [(".", "2", pos + 100, "G", "T", "0/1", ".")]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, SOURCES)
    assert _calls(result) == expected


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("surface", ("json", "html", "terminal"))
def test_cli_reports_the_person_genotype_in_auto_and_explicit_builds(
    tmp_path, build, explicit, surface
):
    data = tmp_path / "data"
    pos = _data(data, build)
    (data / "config.toml").write_text(
        "[sources]\nclinvar=false\nsnpedia=true\npharmgkb=true\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=true\n"
    )
    source = _vcf(tmp_path / "split.vcf", [_row(pos, "A", "1/0"), _row(pos, "T", "0/1")])
    length = 249250621 if build == "GRCh37" else 248956422
    source.write_text(
        source.read_text().replace(
            "##fileformat=VCFv4.2\n", f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
        )
    )
    args = ["analyze", str(source), "--data-dir", str(data), "--no-update", "--min-magnitude", "0"]
    if explicit:
        args += ["--build", build.lower()]
    output = tmp_path / f"report.{surface}"
    if surface != "terminal":
        args += ["--output", str(output)]
    run = CliRunner(env={"COLUMNS": "160"}).invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    if surface == "json":
        report = json.loads(output.read_text())
        assert report["input"]["build"] == build
        assert sorted((r["source"], r["genotype_match"]) for r in report["annotations"]) == [
            ("pharmgkb", "AT"),
            ("snpedia", "AT"),
        ]
    elif surface == "terminal":
        # Terminal shows GT and condition, not the sources' summary text.
        # Read the rendered cells so an AT elsewhere in the report cannot
        # conceal partial AG/GT findings or duplicated source rows.
        rows = [
            [cell.strip() for cell in re.split("[│┃]", line)[1:-1]]
            for line in run.output.splitlines()
            if line.lstrip().startswith(("│", "┃"))
        ]
        header = next(row for row in rows if "GT" in row and "Source" in row)
        source_column, gt_column = header.index("Source"), header.index("GT")
        findings = [
            (row[source_column], row[gt_column])
            for row in rows
            if len(row) == len(header) and row[source_column] in ("ClinPGx", "SNPedia")
        ]
        assert sorted(findings) == [("ClinPGx", "AT"), ("SNPedia", "AT")], run.output
    else:
        text = output.read_text()
        assert "Synthetic AT" in text
        assert "Synthetic AG" not in text
        assert "Synthetic GT" not in text


@pytest.mark.parametrize("build", BUILDS)
def test_split_pair_is_not_cut_at_the_pipeline_batch_boundary(tmp_path, build, monkeypatch):
    from allelix.reports import _pipeline

    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", 1)
    data = tmp_path / "data"
    pos = _data(data, build)
    path = _vcf(tmp_path / "input.vcf", [_row(pos, "A", "1/0"), _row(pos, "T", "0/1")])
    result = analyze(data, path, build, SOURCES)
    assert _calls(result) == [("pharmgkb", "AT"), ("snpedia", "AT")]


@pytest.mark.parametrize("build", BUILDS)
def test_abstention_is_counted_even_without_a_matching_genotype_page(tmp_path, build):
    data = tmp_path / "data"
    pos = _data(data, build)
    with closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.execute("DELETE FROM snpedia_genotypes WHERE allele1 != allele2")
    path = _vcf(tmp_path / "input.vcf", [_row(pos, "A", "0/1"), _row(pos, "T", "0/1")])
    result = analyze(data, path, build, ("snpedia",))
    assert not result.annotations
    assert any(w.count and w.count > 0 for w in result.report_warnings)

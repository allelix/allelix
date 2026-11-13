"""Every placement of a ClinVar VariationID keeps its own submissions and locus."""

import contextlib
import gzip
import json
import sqlite3
from collections import Counter

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import (
    ClinVarAnnotator,
    clinvar_db_filename,
    clinvar_record_name,
)
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import iter_clinvar_tsv_records, load_clinvar_tsv
from allelix.models import Variant
from tests.databases.test_clinvar_tsv_loader import _SS_HEADER, _VS_HEADER, _ss_row, _vs_row

RS = "rs900535001"
VID = "900535001"
BUILDS = ("GRCh37", "GRCh38")


def _variant(build, chrom, pos, *, ref="G", alt="A", rs=RS, vid=VID, **kwargs):
    return _vs_row(
        assembly=build,
        chrom=chrom,
        position_vcf=str(pos),
        ref=ref,
        alt=alt,
        rs=rs.removeprefix("rs") if rs else "-1",
        variation_id=vid,
        gene="SYNGENE",
        name="NM_SYNTHETIC:c.1G>A",
        **kwargs,
    )


def _files(path, vs_rows, ss_rows=None, compressed=False):
    if ss_rows is None:
        ss_rows = [_ss_row(variation_id=VID, scv="SCV900535001.1")]
    suffix = ".tsv.gz" if compressed else ".tsv"
    paths = []
    for name, header, rows in (
        ("variants", _VS_HEADER, vs_rows),
        ("submissions", _SS_HEADER, ss_rows),
    ):
        target = path / (name + suffix)
        payload = ("\n".join([header, *rows]) + "\n").encode()
        target.write_bytes(gzip.compress(payload) if compressed else payload)
        paths.append(target)
    return paths


def _source_rows(build, reverse=False, *, same_position=False, duplicate=False, rs=RS):
    rows = [
        _variant(build, "X", 100001, rs=rs),
        _variant(build, "Y", 100001 if same_position else 200001, rs=rs),
    ]
    if duplicate:
        rows = [rows[0], rows[0], rows[1], rows[1]]
    return rows[::-1] if reverse else rows


def _load(data, vs, ss, build):
    path = data / clinvar_db_filename(build)
    load_clinvar_tsv(
        vs, ss, path, build, record_name=clinvar_record_name(build), remote_signal="md5:synthetic"
    )
    return path


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("same_position", (False, True))
@pytest.mark.parametrize("duplicate", (False, True))
@pytest.mark.parametrize("compressed", (False, True))
def test_each_submission_is_joined_once_to_each_distinct_par_placement(
    tmp_path, build, reverse, same_position, duplicate, compressed
):
    submissions = [
        _ss_row(
            variation_id=VID,
            scv="SCV900535001.1",
            reported_phenotype="C900001:First synthetic condition",
        ),
        _ss_row(
            variation_id=VID,
            scv="SCV900535002.1",
            significance="Likely pathogenic",
            reported_phenotype="C900002:Second synthetic condition",
        ),
    ]
    # The other assembly must not supply a third placement in this build.
    source_rows = _source_rows(build, reverse, same_position=same_position, duplicate=duplicate)
    source_rows += [_variant("GRCh38" if build == "GRCh37" else "GRCh37", "1", 777777)]
    vs, ss = _files(tmp_path, source_rows, submissions, compressed)
    rows = list(iter_clinvar_tsv_records(vs, ss, build))
    expected = Counter(
        (c, p, sig)
        for c, p in (("X", 100001), ("Y", 100001 if same_position else 200001))
        for sig in ("Pathogenic", "Likely pathogenic")
    )
    actual = Counter((r["chromosome"], r["position"], r["clinical_significance"]) for r in rows)
    assert actual == expected, "VariationID is not a unique genomic placement"
    assert all(r["allele_id"] == int(VID) and r["rsid"] == RS for r in rows)
    assert all((r["submitters_same"], r["submitters_total"]) == (1, 2) for r in rows)
    assert {r["condition"] for r in rows} == {
        "First synthetic condition",
        "Second synthetic condition",
    }


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("rsid", (RS, ""))
@pytest.mark.parametrize("batch", (False, True))
def test_both_par_copies_are_stored_counted_and_match_their_own_loci(
    tmp_path, build, reverse, rsid, batch
):
    vs, ss = _files(tmp_path, _source_rows(build, reverse, rs=rsid))
    db = _load(tmp_path, vs, ss, build)
    with contextlib.closing(sqlite3.connect(db)) as conn:
        stored = conn.execute("SELECT chromosome,position FROM clinvar_variants").fetchall()
        count, tag = conn.execute(
            "SELECT record_count,local_version_tag FROM database_versions"
        ).fetchone()
    assert Counter(stored) == Counter({("X", 100001): 1, ("Y", 200001): 1})
    assert count == 2 and tag == f"iv:{CLINVAR_INTERPRETER_VERSION}"
    calls = [
        Variant(rsid, c, p, "G", "A", build=build, ref="G", ref_from_input=True)
        for c, p in (("X", 100001), ("Y", 200001))
    ]
    with ClinVarAnnotator(tmp_path, builds=[build]) as ann:
        assert ann.is_ready()
        rows = (
            list(ann.batch_annotate(calls))
            if batch
            else [r for call in calls for r in ann.annotate(call)]
        )
    assert Counter(r.locus for r in rows) == Counter({("X", 100001): 1, ("Y", 200001): 1})
    assert all(r.source == "clinvar" and r.alt == "A" for r in rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("prefetch", (False, True))
def test_reference_queries_use_each_stored_placements_own_alleles(
    tmp_path, build, reverse, prefetch
):
    source = [
        _variant(build, "X", 100001, ref="G", alt="A"),
        _variant(build, "Y", 200001, ref="C", alt="T"),
    ]
    vs, ss = _files(tmp_path, source[::-1] if reverse else source)
    _load(tmp_path, vs, ss, build)
    calls = [
        Variant(RS, c, p, alt, alt, build=build, ref=ref)
        for c, p, ref, alt in (("X", 100001, "G", "A"), ("Y", 200001, "C", "T"))
    ]
    with ClinVarAnnotator(tmp_path, builds=[build]) as ann:
        if prefetch:
            ann.prefetch_rows_at(calls)
        assert [ann.reference_at(RS, build, (v.chromosome, v.position)) for v in calls] == [
            "G",
            "C",
        ]
        assert ann.reference_at(RS, build, ("Y", 100001)) is None
        assert ann.reference_at(RS, build, ("X", 200001)) is None


@pytest.mark.parametrize("reverse", (False, True))
def test_reference_carried_translation_reads_all_placements_on_the_other_build(tmp_path, reverse):
    source = [
        _variant("GRCh37", "X", 100001, ref="A", alt="A"),
        _variant("GRCh37", "Y", 200001, ref="T", alt="T"),
        _variant("GRCh38", "X", 300001, ref="G", alt="A"),
        _variant("GRCh38", "Y", 400001, ref="C", alt="T"),
    ]
    vs, ss = _files(tmp_path, source[::-1] if reverse else source)
    for build in BUILDS:
        _load(tmp_path, vs, ss, build)
    with ClinVarAnnotator(tmp_path) as ann:
        assert ann.normal_allele_for(RS, "GRCh37", ("X", 100001)) == "G"
        assert ann.normal_allele_for(RS, "GRCh37", ("Y", 200001)) == "C"
        assert ann.reference_at(RS, "GRCh37", ("X", 100001)) == "G"
        assert ann.reference_at(RS, "GRCh37", ("Y", 200001)) == "C"


@pytest.mark.parametrize("build", BUILDS)
def test_released_iv9_cache_is_outdated_without_rewriting_it(tmp_path, build):
    vs, ss = _files(tmp_path, [_variant(build, "X", 100001)])
    db = _load(tmp_path, vs, ss, build)
    with contextlib.closing(sqlite3.connect(db)) as conn, conn:
        conn.execute("UPDATE database_versions SET local_version_tag='iv:9'")
    before = db.read_bytes()
    with ClinVarAnnotator(tmp_path, builds=[build]) as ann:
        assert not ann.is_ready(), "The release's old first-placement cache must rebuild"
    assert db.read_bytes() == before
    assert CLINVAR_INTERPRETER_VERSION > 9


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("haploid", (False, True))
def test_analyze_finds_both_par_copies_with_auto_and_explicit_builds(
    tmp_path, build, reverse, explicit, haploid
):
    data = tmp_path / "data"
    data.mkdir()
    vs, ss = _files(tmp_path, _source_rows(build, reverse))
    _load(data, vs, ss, build)
    xlen, ylen = (155270560, 59373566) if build == "GRCh37" else (156040895, 57227415)
    source = tmp_path / "input.vcf"
    gt = "1" if haploid else "0/1"
    source.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=X,length={xlen}>\n##contig=<ID=Y,length={ylen}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nX\t100001\t{RS}\tG\tA\t.\tPASS\t.\tGT\t{gt}\nY\t200001\t{RS}\tG\tA\t.\tPASS\t.\tGT\t{gt}\n"
    )
    output = tmp_path / "report.json"
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(output),
    ]
    if explicit:
        args += ["--build", build.lower()]
    run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    report = json.loads(output.read_text())
    assert report["input"]["build"] == build
    rows = [r for r in report["annotations"] if r["source"] == "clinvar"]
    assert Counter((r["chromosome"], r["position"]) for r in rows) == Counter(
        {("X", 100001): 1, ("Y", 200001): 1}
    )


@pytest.mark.parametrize("build", BUILDS)
def test_identical_duplicate_source_rows_do_not_multiply_submissions(tmp_path, build):
    row = _variant(build, "X", 100001)
    vs, ss = _files(tmp_path, [row, row, row])
    records = list(iter_clinvar_tsv_records(vs, ss, build))
    assert len(records) == 1 and records[0]["chromosome"] == "X"


@pytest.mark.parametrize("build", BUILDS)
def test_invalid_second_placement_does_not_erase_valid_first(tmp_path, build):
    source = [_variant(build, "X", 100001), _variant(build, "Y", -1)]
    vs, ss = _files(tmp_path, source)
    records = list(iter_clinvar_tsv_records(vs, ss, build))
    assert len(records) == 1 and records[0]["position"] == 100001


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("alt,vid", (("A", VID), ("C", "900535002")))
@pytest.mark.parametrize("batch", (False, True))
def test_submissions_do_not_cross_between_variation_ids_at_the_same_par_loci(
    tmp_path, build, reverse, alt, vid, batch
):
    source = [
        _variant(build, chrom, pos, alt=allele, vid=number)
        for number, allele in ((VID, "A"), ("900535002", "C"))
        for chrom, pos in (("X", 100001), ("Y", 200001))
    ]
    submissions = [
        _ss_row(variation_id=VID, reported_phenotype="C900001:First synthetic condition"),
        _ss_row(
            variation_id="900535002",
            scv="SCV900535002.1",
            reported_phenotype="C900002:Second synthetic condition",
        ),
    ]
    vs, ss = _files(tmp_path, source[::-1] if reverse else source, submissions)
    _load(tmp_path, vs, ss, build)
    calls = [
        Variant(RS, chrom, pos, alt, alt, ref="G", ref_from_input=True, build=build)
        for chrom, pos in (("X", 100001), ("Y", 200001))
    ]
    with ClinVarAnnotator(tmp_path, builds=[build]) as ann:
        rows = (
            list(ann.batch_annotate(calls))
            if batch
            else [r for call in calls for r in ann.annotate(call)]
        )
    assert Counter(r.locus for r in rows) == Counter({("X", 100001): 1, ("Y", 200001): 1})
    expected = "First synthetic condition" if vid == VID else "Second synthetic condition"
    assert all(r.condition == expected and r.alt == alt for r in rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("forensic", (False, True))
def test_multiple_placements_keep_the_same_submission_inclusion_policy(tmp_path, build, forensic):
    submissions = [
        _ss_row(variation_id=VID),
        _ss_row(variation_id=VID, scv="SCV900535002.1", significance="Benign", contributes="no"),
    ]
    vs, ss = _files(tmp_path, _source_rows(build), submissions)
    rows = list(iter_clinvar_tsv_records(vs, ss, build, aggregate_only=not forensic))
    expected = {("X", "Pathogenic"): 1, ("Y", "Pathogenic"): 1}
    if forensic:
        expected.update({("X", "Benign"): 1, ("Y", "Benign"): 1})
    assert Counter((r["chromosome"], r["clinical_significance"]) for r in rows) == Counter(
        expected
    )
    assert all(
        r["submitters_total"] == 1 for r in rows if r["clinical_significance"] == "Pathogenic"
    )

"""Padding in identity rows must not erase the same single-base source findings."""

import json
import sqlite3
from contextlib import closing

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _identity, _vcf
from validated_call_helpers import BUILDS, POSITIONS, RSID, analyze, setup, signatures

from allelix.cli import main

SOURCES = ("gwas", "snpedia")
SHAPES = {
    "trimmed": (0, "G", "A"),
    "prefix": (-1, "CG", "CA"),
    "suffix": (0, "GC", "AC"),
    "both": (-1, "CGC", "CAC"),
}


def _replace_identity(data, shape, *, kind="snv"):
    delta, ref, alt = SHAPES[shape]
    if kind == "indel":
        alt = "C"
    elif kind == "other_locus":
        delta += 20
    (data / "rsid_identity.sqlite").unlink()
    _identity(
        data,
        [
            (
                RSID,
                "1",
                POSITIONS["GRCh37"] + delta,
                ref,
                alt,
                POSITIONS["GRCh38"] + delta,
                ref,
                alt,
                0,
            )
        ],
    )


def _input(path, build, shape, named):
    delta, ref, alt = SHAPES[shape]
    return _vcf(
        path, [(RSID if named else ".", "1", POSITIONS[build] + delta, ref, alt, "0/1", ".")]
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("cache_shape", ("prefix", "suffix", "both"))
@pytest.mark.parametrize("input_shape", SHAPES)
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("source", SOURCES)
def test_padded_identity_matches_trimmed_identity_for_the_same_call(
    tmp_path, build, cache_shape, input_shape, named, source
):
    data = tmp_path / "data"
    setup(data, build)
    path = _input(tmp_path / "input.vcf", build, input_shape, named)
    before = analyze(data, path, build, (source,))
    assert [a.source for a in before.annotations] == [source], "The unpadded control must match"
    _replace_identity(data, cache_shape)
    after = analyze(data, path, build, (source,))
    assert signatures(after) == signatures(before), (
        "Only identity's storage spelling changed, not the allele or its locus"
    )
    assert [(a.locus, a.record_key, a.record_index) for a in after.annotations] == [
        (a.locus, a.record_key, a.record_index) for a in before.annotations
    ]
    assert after.strand_abstentions == before.strand_abstentions


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("input_shape", ("trimmed", "prefix"))
def test_cli_recovers_and_annotates_with_padded_identity_auto_or_explicit(
    tmp_path, build, explicit, named, input_shape
):
    data = tmp_path / "data"
    setup(data, build)
    _replace_identity(data, "prefix")
    (data / "config.toml").write_text(
        "[sources]\nclinvar=false\nsnpedia=true\ngwas=true\npharmgkb=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=true\n"
    )
    path = _input(tmp_path / "input.vcf", build, input_shape, named)
    length = 249250621 if build == "GRCh37" else 248956422
    path.write_text(
        path.read_text().replace(
            "##fileformat=VCFv4.2\n", f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
        )
    )
    out = tmp_path / "report.json"
    args = [
        "analyze",
        str(path),
        "--data-dir",
        str(data),
        "--no-update",
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--output",
        str(out),
    ]
    if explicit:
        args += ["--build", build.lower()]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    report = json.loads(out.read_text())
    assert report["input"]["build"] == build
    assert sorted((r["source"], r["rsid"]) for r in report["annotations"]) == [
        ("gwas", RSID),
        ("snpedia", RSID),
    ]
    assert not any(any(report["strand_abstentions"][source].values()) for source in SOURCES)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("kind", ("snv", "indel", "other_locus"))
def test_only_equivalent_single_base_identity_rows_are_evidence(tmp_path, build, kind):
    data = tmp_path / "data"
    setup(data, build)
    if kind != "snv":
        _replace_identity(data, "prefix", kind=kind)
    result = analyze(data, _input(tmp_path / "input.vcf", build, "trimmed", True), build, SOURCES)
    assert {a.source for a in result.annotations} == (set(SOURCES) if kind == "snv" else set())


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("orientation", (None, "minus"))
@pytest.mark.parametrize("input_shape", ("trimmed", "prefix"))
def test_trimmed_identity_evidence_serves_unknown_and_stated_reverse_pages(
    tmp_path, build, orientation, input_shape
):
    data = tmp_path / "data"
    setup(data, build)
    with closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.execute("UPDATE snpedia_genotypes SET orientation=?", (orientation,))
        if orientation == "minus":
            conn.execute("UPDATE snpedia_genotypes SET allele1='C', allele2='T'")
    path = _input(tmp_path / "input.vcf", build, input_shape, True)
    before = analyze(data, path, build, ("snpedia",))
    assert len(before.annotations) == 1
    _replace_identity(data, "prefix")
    after = analyze(data, path, build, ("snpedia",))
    assert signatures(after) == signatures(before)


@pytest.mark.parametrize("build", BUILDS)
def test_trimming_does_not_turn_a_both_fit_site_into_decisive_strand_evidence(tmp_path, build):
    data = tmp_path / "data"
    pos = setup(data, build)
    with closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.execute("UPDATE snpedia_genotypes SET allele1='C', allele2='G', orientation=NULL")
    with closing(sqlite3.connect(data / "gwas.sqlite")) as conn, conn:
        conn.execute("UPDATE gwas_associations SET risk_allele='C'")
    (data / "rsid_identity.sqlite").unlink()
    _identity(
        data,
        [(RSID, "1", POSITIONS["GRCh37"] - 1, "AC", "AG", POSITIONS["GRCh38"] - 1, "AC", "AG", 1)],
    )
    source = _vcf(tmp_path / "input.vcf", [(RSID, "1", pos, "C", "G", "0/1", ".")])
    result = analyze(data, source, build, SOURCES)
    assert not result.annotations
    assert all(sum(result.strand_abstentions[source].values()) == 1 for source in SOURCES)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", ("clinvar", "pharmgkb"))
@pytest.mark.parametrize("input_shape", ("trimmed", "prefix"))
def test_existing_named_clinvar_and_clinpgx_findings_are_unchanged(
    tmp_path, build, source, input_shape
):
    data = tmp_path / "data"
    setup(data, build)
    path = _input(tmp_path / "input.vcf", build, input_shape, True)
    before = analyze(data, path, build, (source,))
    assert len(before.annotations) == 1
    _replace_identity(data, "prefix")
    after = analyze(data, path, build, (source,))
    assert signatures(after) == signatures(before)

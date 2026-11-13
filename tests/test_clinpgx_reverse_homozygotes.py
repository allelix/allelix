"""ClinPGx hom-ALT matches need site evidence, not REF in that genotype.

All caches and calls are synthetic. A definitive reverse site selects that
reading; preserve suppression and the existing guarded heterozygote path.
"""

import contextlib
import json
import sqlite3
from dataclasses import replace

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _identity, _pharm

from allelix.annotators import get_annotators
from allelix.cli import main
from allelix.models import Variant
from allelix.utils.allele import complement

RS = "rs900502001"
BUILDS = ("GRCh37", "GRCh38")
PAIRS = [(r, a) for r in "ACGT" for a in "ACGT" if a not in (r, complement(r))]


def _data(path, ref="G", alt="A", *, identity="present", level="1A", reverse=False):
    path.mkdir(exist_ok=True)
    genotypes = [complement(ref) * 2, "".join(sorted(complement(ref + alt))), complement(alt) * 2]
    _pharm(path, RS, "SYNGENE", genotypes[::-1] if reverse else genotypes, level=level)
    if identity != "missing":
        rows = [(RS, "1", 1000, ref, alt, 2000, ref, alt, int(identity == "palindromic"))]
        if identity == "ambiguous":
            rows.append((RS, "1", 1000, ref, complement(alt), 2000, ref, complement(alt), 0))
            # Only the ambiguous target row is evidence in this control.
            with contextlib.closing(sqlite3.connect(path / "pharmgkb.sqlite")) as conn, conn:
                conn.execute(
                    "DELETE FROM pharmgkb_annotations WHERE genotype != ?", (complement(alt) * 2,)
                )
        if identity == "multiple":
            rows.append((RS, "2", 3000, ref, alt, 4000, ref, alt, 0))
        _identity(path, rows)
    return path


def _call(build, ref="G", alt="A"):
    return Variant(RS, "1", 1000 if build == "GRCh37" else 2000, alt, alt, build=build, ref=ref)


def _find(data, call, path, *, switched_off=()):
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, switched_off=switched_off)
        for item in registry:
            stack.enter_context(item)
        annotator = next(a for a in registry if a.name == "pharmgkb")
        assert annotator.is_ready(), "The synthetic cache must be usable before matching"
        return (
            annotator.annotate(call)
            if path == "single"
            else list(annotator.batch_annotate([call]))
        )


def _assert_match(rows, call, level):
    assert len(rows) == 1, "Hom-ALT must match the independently reverse-only site's source row"
    row = rows[0]
    assert row.source == "pharmgkb" and row.rsid == RS
    assert row.significance == f"pharmgkb_loe_{level.lower()}"
    assert row.genotype_match == call.allele1 * 2
    assert row.alt == call.allele1 and row.locus_ref == call.ref
    assert row.carried_alleles == (call.allele1, call.allele2)
    assert row.locus == (call.chromosome, call.position)
    assert row.gene == "SYNGENE" and "Synthetic drug" in row.description


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("ref,alt", PAIRS)
@pytest.mark.parametrize("level", ["1A", "3"])
def test_reverse_only_site_matches_homozygous_alt(tmp_path, build, path, ref, alt, level):
    data = _data(tmp_path, ref, alt, level=level)
    call = _call(build, ref, alt)
    _assert_match(_find(data, call, path), call, level)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("support_suppressed", [False, True])
def test_site_evidence_uses_all_genotypes_not_only_the_homozygote(
    tmp_path, build, path, reverse, support_suppressed
):
    data = _data(tmp_path, reverse=reverse)
    # TT alone is ambiguous at G/A/T; C in the site's CC/CT rows establishes
    # the reverse-only reading for the full source allele set C/T.
    _identity(
        data,
        [
            (RS, "1", 1000, "G", "A", 2000, "G", "A", 0),
            (RS, "1", 1000, "G", "T", 2000, "G", "T", 0),
        ],
    )
    if support_suppressed:
        with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn, conn:
            conn.execute("UPDATE pharmgkb_annotations SET is_nonfinding=1 WHERE genotype!='TT'")
    call = _call(build)
    _assert_match(_find(data, call, path), call, "1A")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize(
    "state", ["missing", "disabled", "palindromic", "ambiguous", "multiple", "wrong_locus"]
)
def test_uncertain_or_unavailable_evidence_does_not_admit_new_hom_alt(
    tmp_path, build, path, state
):
    data = _data(tmp_path, identity=state)
    call = _call(build)
    if state == "wrong_locus":
        call = replace(call, position=9999)
    off = {"rsid_identity"} if state == "disabled" else set()
    assert _find(data, call, path, switched_off=off) == []


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize(
    "control", ["direct", "heterozygote", "cpic_nonfinding", "clinvar_homref", "no_call", "di"]
)
def test_existing_matching_and_suppression_stay_in_scope(tmp_path, build, path, control):
    data = _data(tmp_path)
    call = _call(build)
    expected = 0
    if control == "direct":
        # Identity resolves this source site as reverse-only: TT describes
        # forward AA, not this forward TT call. Exact spelling cannot override it.
        call = replace(call, allele1="T", allele2="T")
        expected = 0
    elif control == "heterozygote":
        call = replace(call, allele1="G", allele2="A")
        expected = 1
    elif control == "cpic_nonfinding":
        with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn, conn:
            conn.execute("UPDATE pharmgkb_annotations SET is_nonfinding=1 WHERE genotype='TT'")
    elif control == "clinvar_homref":
        from test_annotation_parity_and_named_disclosures import _clinvar

        _clinvar(
            data,
            build,
            [(RS, "1", call.position, "A", "G", "SYNGENE", "Pathogenic", "Pathogenic", 1, 1)],
        )
    elif control == "no_call":
        call = replace(call, allele1="-", allele2="-")
    else:
        call = replace(call, allele1="D", allele2="I")
    assert len(_find(data, call, path)) == expected


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_cli_reports_the_reverse_hom_alt_in_input_orientation(tmp_path, build, command, surface):
    data = _data(tmp_path / "data")
    call = _call(build)
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t{call.position}\t{RS}\tG\tA\t.\tPASS\t.\tGT\t1/1\n"
    )
    args = [command, str(src), "--build", build.lower(), "--data-dir", str(data), "--no-update"]
    output = tmp_path / f"report.{surface}"
    if surface != "terminal":
        args += ["--output", str(output)]
    run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    if surface == "json":
        rows = [
            r for r in json.loads(output.read_text())["annotations"] if r["source"] == "pharmgkb"
        ]
        assert len(rows) == 1, "The missing reverse homozygote must reach the report"
        assert rows[0]["genotype_match"] == "AA" and rows[0]["alt"] == "A"
    elif surface == "html":
        text = output.read_text()
        marker = '<script id="variant-data" type="application/json">'
        assert marker in text, "Missing reverse finding: no variant data was rendered"
        start = text.index(marker)
        payload = text[text.index(">", start) + 1 : text.index("</script", start)]
        rows = json.loads(payload)
        assert len(rows) == 1 and rows[0]["rsid"] == RS and rows[0]["genotype"] == "AA"
    else:
        assert RS in run.output and "SYNGENE" in run.output and "AA" in run.output


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
def test_orientation_evidence_is_read_on_the_calls_build(tmp_path, build, path):
    data = _data(tmp_path)
    _identity(data, [(RS, "1", 1000, "G", "A", 2000, "C", "T", 0)])
    call = _call(build) if build == "GRCh37" else _call(build, "C", "T")
    _assert_match(_find(data, call, path), call, "1A")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("state", ["missing", "palindromic", "ambiguous"])
def test_existing_exact_and_legacy_heterozygote_matches_do_not_need_new_evidence(
    tmp_path, build, path, state
):
    data = _data(tmp_path, identity=state)
    direct = replace(_call(build), allele1="T", allele2="T")
    assert len(_find(data, direct, path)) == 1
    if state == "missing":
        heterozygote = replace(_call(build), allele1="G", allele2="A")
        assert len(_find(data, heterozygote, path)) == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
def test_batch_never_borrows_another_rsids_orientation(tmp_path, build, reverse):
    data = _data(tmp_path)
    other = "rs900502002"
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn, conn:
        conn.execute(
            "INSERT INTO pharmgkb_annotations "
            "(rsid, genotype, gene, drugs, phenotype, annotation_text, "
            "level_of_evidence, function_class, is_nonfinding) "
            "VALUES (?, 'TT', 'OTHER', 'Other synthetic drug', "
            "'Other response', 'Other', '1A', 'decreased', 0)",
            (other,),
        )
    _identity(
        data,
        [
            (RS, "1", 1000, "G", "A", 2000, "G", "A", 0),
            (other, "1", 3000, "G", "A", 4000, "G", "A", 0),
            (other, "1", 3000, "G", "T", 4000, "G", "T", 0),
        ],
    )
    first = _call(build)
    second = replace(first, rsid=other, position=3000 if build == "GRCh37" else 4000)
    calls = [second, first] if reverse else [first, second]
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data)
        for item in registry:
            stack.enter_context(item)
        annotator = next(a for a in registry if a.name == "pharmgkb")
        rows = list(annotator.batch_annotate(calls))
    assert [r.rsid for r in rows] == [RS], (
        "Only the independently resolved site's row is a finding"
    )

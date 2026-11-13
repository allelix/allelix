"""One call cannot receive two ClinPGx genotypes from undecided strand readings."""

import contextlib
import json
import re
import sqlite3
from dataclasses import replace
from html import unescape

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity, _pharm

from allelix.annotators import get_annotators
from allelix.cli import main
from allelix.models import Variant
from allelix.reports.json_report import SCHEMA_VERSION

RS = "rs900583011"
OTHER = "rs900583012"
BUILDS = ("GRCh37", "GRCh38")
POS = {"GRCh37": 1000000, "GRCh38": 5000000}
LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}
COUNT = "ambiguous_strand_calls"
WARNING = "clinpgx_ambiguous_strand"
MODES = ("missing", "disabled", "no_fit", "elsewhere", "multi_locus")


def _data(
    path,
    *,
    rows=("AG", "CT"),
    mode="missing",
    reverse=False,
    nonfinding=(),
    clinvar=False,
    site="GA",
    duplicate=False,
):
    path.mkdir(exist_ok=True)
    rows = tuple(rows) * (2 if duplicate else 1)
    _pharm(path, RS, "SYNGENE", rows[::-1] if reverse else rows)
    with contextlib.closing(sqlite3.connect(path / "pharmgkb.sqlite")) as conn, conn:
        conn.execute("UPDATE pharmgkb_annotations SET annotation_text='row ' || genotype")
        for genotype in nonfinding:
            conn.execute(
                "UPDATE pharmgkb_annotations SET is_nonfinding=1 WHERE genotype=?", (genotype,)
            )
    if mode != "missing":
        offset = 100 if mode == "elsewhere" else 0
        identity = [
            (
                RS,
                "1",
                POS["GRCh37"] + offset,
                site[0],
                alt,
                POS["GRCh38"] + offset,
                site[0],
                alt,
                int(mode == "palindromic"),
            )
            for alt in site[1:]
        ]
        if mode == "multi_locus":
            identity.append((RS, "2", 2000000, "G", "A", 6000000, "G", "A", 0))
        _identity(path, identity)
    if clinvar:
        for build in BUILDS:
            _clinvar(
                path,
                build,
                [(RS, "1", POS[build], "G", "A", "SYNGENE", "Pathogenic", "Pathogenic", 1, 1)],
            )
    enabled = {"pharmgkb"} | ({"rsid_identity"} if mode not in ("missing", "disabled") else set())
    if clinvar:
        enabled.add("clinvar")
    names = (
        "clinvar",
        "pharmgkb",
        "rsid_identity",
        "snpedia",
        "gwas",
        "gnomad",
        "alphamissense",
        "cadd",
    )
    (path / "config.toml").write_text(
        "[sources]\n" + "".join(f"{n}={str(n in enabled).lower()}\n" for n in names)
    )
    return path


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    root = tmp_path_factory.mktemp("ambiguous-pgx")
    memo = {}

    def get(**options):
        key = tuple(sorted(options.items()))
        if key not in memo:
            memo[key] = _data(root / str(len(memo)), **options)
        return memo[key]

    return get


def _call(build, genotype="AG", ref="G", **changes):
    return replace(
        Variant(RS, "1", POS[build], genotype[0], genotype[1], build=build, ref=ref), **changes
    )


def _run(data, calls, path):
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data)
        for annotator in registry:
            stack.enter_context(annotator)
        pgx = next(a for a in registry if a.name == "pharmgkb")
        assert pgx.is_ready(), "Synthetic cache is ready before checking the behaviour"
        assert pgx.run_counts().get(COUNT, 0) == 0
        if path == "batch":
            findings = list(pgx.batch_annotate(calls))
        else:
            findings = [a for call in calls for a in pgx.annotate(call)]
        return findings, pgx.run_counts().get(COUNT, 0)


def _labels(rows):
    return sorted(a.description.split(" — ")[-1] for a in rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize(
    "genotype,ref,rows",
    (
        ("AG", "G", ("AG", "CT")),
        ("CT", "C", ("AG", "CT")),
        ("AC", "A", ("AC", "GT")),
        ("GT", "G", ("AC", "GT")),
    ),
)
def test_two_undecided_readings_abstain_independent_of_row_order(
    cache, build, path, mode, reverse, genotype, ref, rows
):
    findings, count = _run(
        cache(mode=mode, reverse=reverse, rows=rows), [_call(build, genotype, ref)], path
    )
    assert not findings, "Do not prefer exact spelling or emit both readings at an undecided site"
    assert count == 1, "Count the call once, not one count per matched source row"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("clinvar", (False, True))
@pytest.mark.parametrize("nonfinding", ((), ("AG",), ("CT",)))
@pytest.mark.parametrize("reverse", (False, True))
def test_nonfinding_filter_cannot_decide_ambiguous_orientation(
    cache, build, path, clinvar, nonfinding, reverse
):
    findings, count = _run(
        cache(nonfinding=nonfinding, clinvar=clinvar, reverse=reverse), [_call(build)], path
    )
    assert not findings, (
        "A nonfinding row is still evidence; reporting filters cannot choose strand"
    )
    assert count == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("reverse", (False, True))
def test_duplicate_source_rows_do_not_inflate_the_call_count(cache, build, path, reverse):
    findings, count = _run(cache(duplicate=True, reverse=reverse), [_call(build)], path)
    assert count == 1, "Four matched source rows still represent one undecidable input call"
    assert not findings


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("reverse", (False, True))
def test_one_calls_ambiguity_does_not_leak_to_a_later_call(cache, build, path, reverse):
    # This other genotype has only one candidate row. Its call must not
    # inherit the AG/CT ambiguity merely because it has the same rsID.
    ambiguous = _call(build)
    exact_only = _call(build, "AA", chromosome="", position=0, ref=None)
    calls = [ambiguous, exact_only]
    if reverse:
        calls.reverse()
    findings, count = _run(cache(rows=("AG", "CT", "AA")), calls, path)
    assert count == 1
    assert _labels(findings) == ["row AA"]
    assert findings[0].locus == ("", 0)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("reading", ("direct", "complement"))
@pytest.mark.parametrize("reverse", (False, True))
def test_nonfinding_row_can_supply_definitive_site_evidence(cache, build, path, reading, reverse):
    # CC/GG alone fit in both readings. The nonfinding AA or TT decides
    # orientation; do not discard it before making the site decision.
    evidence = "AA" if reading == "direct" else "TT"
    data = cache(
        rows=("CC", "GG", evidence),
        nonfinding=(evidence,),
        mode="no_fit",
        site="GAC",
        reverse=reverse,
    )
    findings, count = _run(data, [_call(build, "CC")], path)
    assert _labels(findings) == ["row CC" if reading == "direct" else "row GG"]
    assert count == 0, "A resolved reading is not an ambiguity"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize(
    "case",
    (
        "exact_only",
        "complement_only",
        "same_genotype_rows",
        "palindromic",
        "no_ref",
        "no_call",
        "absent_rsid",
    ),
)
def test_single_reading_and_existing_exclusions_are_not_counted(cache, build, path, case):
    options = {}
    call = _call(build)
    expected = ["row AG"]
    if case == "exact_only":
        options["rows"] = ("AG",)
    elif case == "complement_only":
        options["rows"] = ("CT",)
        expected = ["row CT"]
    elif case == "same_genotype_rows":
        options.update(rows=("AG",), duplicate=True)
        expected = ["row AG", "row AG"]
    elif case == "palindromic":
        options.update(rows=("AT",), mode="palindromic", site="AT")
        call = _call(build, "AT", "A")
        expected = ["row AT"]
    elif case == "no_ref":
        options["rows"] = ("AG",)
        call = _call(build, ref=None)
    elif case == "no_call":
        call = _call(build, "--")
        expected = []
    elif case == "absent_rsid":
        call = _call(build, rsid=OTHER)
        expected = []
    findings, count = _run(cache(**options), [call], path)
    assert _labels(findings) == expected
    assert count == 0


def _cli(
    tmp_path, build, *, command="analyze", explicit=True, surface="json", clinvar=False, floor="0"
):
    data = _data(tmp_path / "data", clinvar=clinvar, duplicate=True)
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=1,length={LENGTH[build]}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t{POS[build]}\t{RS}\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
    )
    args = [command, str(source), "--data-dir", str(data), "--no-update", "--min-magnitude", floor]
    if explicit:
        args += ["--build", build.lower()]
    out = tmp_path / ("report." + surface)
    if surface != "terminal":
        args += ["--output", str(out)]
    result = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    return out.read_text() if surface != "terminal" else result.output


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("command", ("analyze", "pharmacogenomics"))
@pytest.mark.parametrize("clinvar", (False, True))
def test_pipeline_reports_one_ambiguity_in_the_current_schema_and_keeps_other_sources(
    tmp_path, build, explicit, command, clinvar
):
    report = json.loads(_cli(tmp_path, build, explicit=explicit, command=command, clinvar=clinvar))
    assert report["schema_version"] == SCHEMA_VERSION
    assert report["input"]["build"] == build
    warnings = [w for w in report["report_warnings"] if w["code"] == WARNING]
    assert len(warnings) == 1 and warnings[0]["count"] == 1
    text = warnings[0]["message"].lower()
    assert "clinpgx" in text and ("strand" in text or "reading" in text)
    assert not [a for a in report["annotations"] if a["source"] == "pharmgkb"]
    if clinvar and command == "analyze":
        assert [a for a in report["annotations"] if a["source"] == "clinvar"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("surface", ("terminal", "html"))
def test_human_report_explains_the_abstention(tmp_path, build, surface):
    text = _cli(tmp_path, build, surface=surface).lower()
    if surface == "html":
        # Inspect visible warnings, not static JavaScript/CSS/legend prose.
        notices = re.findall(r'<div class="notice-warn">(.*?)</div>', text, flags=re.S)
        text = unescape(" ".join(re.sub(r"<[^>]*>", "", n) for n in notices))
    assert "clinpgx" in text and (
        "ambiguous" in text or "cannot decide" in text or "could not decide" in text
    )
    assert "strand" in text or "reading" in text


@pytest.mark.parametrize("build", BUILDS)
def test_hidden_findings_floor_does_not_hide_the_abstention(tmp_path, build):
    report = json.loads(_cli(tmp_path, build, floor="10"))
    warnings = [w for w in report["report_warnings"] if w["code"] == WARNING]
    assert len(warnings) == 1 and warnings[0]["count"] == 1

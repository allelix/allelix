# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A source choice isolates its cache from every field of a report."""

import contextlib
import sqlite3

import pytest
from test_annotation_parity_and_named_disclosures import _clinvar, _run, _stamp
from validated_call_helpers import BUILDS, POSITIONS, RSID, call, reference_effect, setup, sql

from allelix.annotators import get_annotators
from allelix.config import load_config, save_config
from allelix.databases._versions import ALPHAMISSENSE_SCHEMA_VERSION, CADD_SCHEMA_VERSION
from allelix.databases.schema import ALPHAMISSENSE_SCHEMA, CADD_SCHEMA

SOURCES = (
    "clinvar",
    "pharmgkb",
    "gwas",
    "snpedia",
    "gnomad",
    "alphamissense",
    "cadd",
    "rsid_identity",
)
FLAGS = {
    "gnomad": "--no-gnomad",
    "alphamissense": "--no-alphamissense",
    "cadd": "--no-cadd",
    "snpedia": "--exclude-snpedia",
}
CHOICES = [(source, "config") for source in SOURCES] + [(source, "flag") for source in FLAGS]
FLOORS = ("--min-magnitude", "0", "--gwas-min-magnitude", "0", "--snpedia-min-magnitude", "0")


def _scores(data):
    for source, schema, tag, statement in (
        (
            "alphamissense",
            ALPHAMISSENSE_SCHEMA,
            ALPHAMISSENSE_SCHEMA_VERSION,
            "INSERT INTO alphamissense_scores "
            "(chrom,pos,ref,alt,rsid,am_pathogenicity,am_class) "
            "VALUES ('1',11796321,'G','A','rs1801133',0.9,'likely_pathogenic')",
        ),
        (
            "cadd",
            CADD_SCHEMA,
            CADD_SCHEMA_VERSION,
            "INSERT INTO cadd_scores (chrom,pos,ref,alt,phred) VALUES ('1',11796321,'G','A',35.0)",
        ),
    ):
        with contextlib.closing(sqlite3.connect(data / (source + ".sqlite"))) as conn:
            conn.executescript(schema)
            conn.execute(statement)
            _stamp(conn, source, 1, f"sv:{tag}")
            conn.commit()


def _park(data, source):
    pattern = "clinvar.*.sqlite" if source == "clinvar" else source + ".sqlite"
    paths = list(data.glob(pattern))
    assert paths, source
    for path in paths:
        path.rename(path.with_suffix(".parked"))


def _floors(command):
    return ("--min-magnitude", "0") if command == "pharmacogenomics" else FLOORS


def _configure(data, source=None, method="config", *, command="analyze"):
    # Pharma refuses --exclude-snpedia; its supported source switch is config.
    if command == "pharmacogenomics" and source == "snpedia":
        method = "config"
    config = load_config(data, create=False)
    config.sources.update(dict.fromkeys(SOURCES, True))
    if source is not None and method == "config":
        config.sources[source] = False
    save_config(data, config)
    return (FLAGS[source],) if source is not None and method == "flag" else ()


def _prepare(root, build, scenario="ordinary"):
    data = root / "data"
    reference = scenario.startswith("reference")
    pos = reference_effect(data, build) if reference else setup(data, build)
    _scores(data)
    other = "GRCh38" if build == "GRCh37" else "GRCh37"
    _clinvar(
        data,
        other,
        [
            (
                RSID,
                "1",
                POSITIONS[other],
                "A" if reference else "G",
                "G" if reference else "A",
                "MTHFR",
                "Pathogenic",
                "Pathogenic",
                1,
                1,
            )
        ],
    )
    if reference:
        for db, statement in (
            ("pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='AA'"),
            ("snpedia.sqlite", "UPDATE snpedia_genotypes SET allele1='A',allele2='A'"),
            (
                "gnomad.sqlite",
                "UPDATE gnomad_frequencies SET ref='A',alt='G',af=0.01"
                if build == "GRCh37"
                else "UPDATE gnomad_frequencies SET ref='G',alt='A',af=0.99",
            ),
        ):
            sql(data, db, statement)
        if scenario == "reference-input":
            (data / f"clinvar.{other}.sqlite").unlink()
            _park(data, "gnomad")
        elif scenario == "reference-other":
            (data / f"clinvar.{build}.sqlite").unlink()
            _park(data, "gnomad")
        elif scenario == "reference-frequency":
            _park(data, "clinvar")
    if scenario == "reverse":
        sql(data, "gwas.sqlite", "UPDATE gwas_associations SET risk_allele='T'")
        sql(
            data,
            "snpedia.sqlite",
            "UPDATE snpedia_genotypes SET allele1='C',allele2='T',orientation=NULL",
        )
        sql(data, "pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='CT'")
    if scenario == "palindromic":
        sql(data, "rsid_identity.sqlite", "UPDATE rsid_identity SET palindromic=1")
        sql(data, "pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='CT'")
    path = call(root / "input.vcf", pos, gt="1/1" if reference else "0/1")
    return data, path


def _bytes(root, data, path, build, command, extra, label):
    report, _ = _run(root, data, path, build, "json", command=command, extra=extra, label=label)
    raw = (root / (label + ".json")).read_bytes()
    raw = b"".join(
        line for line in raw.splitlines(keepends=True) if not line.startswith(b'  "generated_at":')
    )
    return report, raw


def _pair(root, data, path, build, command, source, method, extra=()):
    options = _floors(command) + _configure(data, source, method, command=command) + extra
    if command == "methylation":
        options += ("--include-gwas",)
    before, first = _bytes(root, data, path, build, command, options, "present")
    _park(data, source)
    after, second = _bytes(root, data, path, build, command, options, "absent")
    for report in (before, after):
        state = next(row for row in report["sources"] if row["name"] == source)
        assert state["state"] == "skipped" and state["version"] is None
        assert all(row["source"] != source for row in report["annotations"])
    assert first == second


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics", "methylation"])
@pytest.mark.parametrize("source,method", CHOICES)
def test_every_disabled_source_has_cache_independent_complete_reports(
    tmp_path, build, command, source, method
):
    data, path = _prepare(tmp_path, build)
    _pair(tmp_path, data, path, build, command, source, method)


# Each callback edge gets a call that reaches it, including the two ClinVar builds.
EDGES = [
    ("clinvar", "reference-input"),
    ("clinvar", "reference-other"),
    ("gnomad", "reference-frequency"),
    ("rsid_identity", "reference-input"),
    ("rsid_identity", "reference-other"),
    ("rsid_identity", "reference-frequency"),
    ("rsid_identity", "reverse"),
    ("rsid_identity", "palindromic"),
]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source,scenario", EDGES)
def test_each_callback_edge_is_isolated_in_complete_reports(tmp_path, build, source, scenario):
    data, path = _prepare(tmp_path, build, scenario)
    _pair(tmp_path, data, path, build, "analyze", source, "config")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source,scenario", EDGES)
def test_disabled_callback_provider_is_never_queried(
    tmp_path, monkeypatch, build, source, scenario
):
    data, path = _prepare(tmp_path, build, scenario)
    options = FLOORS + _configure(data, source)
    registry = get_annotators(data)
    cls = type(next(a for a in registry if a.name == source))
    queried = []
    original = cls._connection

    def track(self, *args, **kwargs):
        queried.append(source)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(cls, "_connection", track)
    _bytes(tmp_path, data, path, build, "analyze", options, "report")
    assert not queried


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", ["clinvar", "gnomad", "rsid_identity"])
@pytest.mark.parametrize("shape", ["plain", "prefix"])
def test_resolution_and_panel_do_not_consult_a_disabled_provider(tmp_path, build, source, shape):
    data, _ = _prepare(tmp_path, build)
    for other in ("clinvar", "gnomad", "rsid_identity"):
        if other != source:
            _park(data, other)
    path = call(tmp_path / "input.vcf", POSITIONS[build], rsid=".", shape=shape)
    panel = tmp_path / "panel.txt"
    panel.write_text(RSID + "\n")
    _pair(tmp_path, data, path, build, "analyze", source, "config", ("--filter-file", str(panel)))


@pytest.mark.parametrize("build", BUILDS)
def test_gnomad_fallback_flag_has_complete_report_parity(tmp_path, build):
    data, path = _prepare(tmp_path, build, "reference-frequency")
    _pair(tmp_path, data, path, build, "analyze", "gnomad", "flag")


@pytest.mark.parametrize("build", BUILDS)
def test_all_on_control_queries_every_installed_source(tmp_path, build):
    data, path = _prepare(tmp_path, build)
    options = FLOORS + _configure(data)
    report, first = _bytes(tmp_path, data, path, build, "analyze", options, "first")
    again, second = _bytes(tmp_path, data, path, build, "analyze", options, "second")
    assert first == second and report["annotations"] == again["annotations"]
    assert {r["name"] for r in report["sources"] if r["state"] == "used"} == set(SOURCES)
    assert {r["source"] for r in report["annotations"]} == {
        "clinvar",
        "pharmgkb",
        "gwas",
        "snpedia",
    }
    assert all(r["allele_frequency"] is not None for r in report["annotations"])
    assert any(r["cadd_phred"] == 35.0 for r in report["annotations"])
    assert any(r["am_pathogenicity"] == 0.9 for r in report["annotations"])


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", ["clinvar", "gnomad", "rsid_identity"])
def test_common_allele_cap_uses_only_enabled_evidence(tmp_path, build, source):
    data, path = _prepare(tmp_path, build)
    sql(
        data,
        f"clinvar.{build}.sqlite",
        "UPDATE clinvar_variants SET review_status='no assertion criteria provided'",
    )
    sql(data, "gnomad.sqlite", "UPDATE gnomad_frequencies SET af=0.2")
    _configure(data)
    control, _ = _bytes(tmp_path, data, path, build, "analyze", FLOORS, "on")
    assert control["clinvar"]["common_allele_rows_scored_down"] > 0
    _pair(tmp_path, data, path, build, "analyze", source, "config")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", ["pharmacogenomics", "methylation"])
def test_focused_default_gwas_exclusion_has_complete_report_parity(tmp_path, build, command):
    data, path = _prepare(tmp_path, build)
    options = _floors(command) + _configure(data)
    before, first = _bytes(tmp_path, data, path, build, command, options, "present")
    _park(data, "gwas")
    after, second = _bytes(tmp_path, data, path, build, command, options, "absent")
    for report in (before, after):
        state = next(row for row in report["sources"] if row["name"] == "gwas")
        assert state["state"] == "skipped" and state["reason"] == "excluded"
    assert first == second


@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
@pytest.mark.parametrize(
    "source,scenario",
    [
        ("clinvar", "reference-input"),
        ("gnomad", "reference-frequency"),
        ("rsid_identity", "reference-frequency"),
    ],
)
def test_each_unavailable_dependency_is_disclosed(tmp_path, source, scenario, surface):
    data, path = _prepare(tmp_path, "GRCh37", scenario)
    options = FLOORS + _configure(data, source)
    report, console = _run(tmp_path, data, path, "GRCh37", surface, extra=options)
    if surface == "json":
        notices = [
            w for w in report["report_warnings"] if w["code"] == "source_dependency_unavailable"
        ]
        assert notices
        text = " ".join(w["message"] for w in notices)
    else:
        text = report
    word = "identity" if source == "rsid_identity" else source
    assert word in text.lower()
    assert "normal allele" in text.lower()
    assert "normal allele" in console.lower()


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", ["clinvar", "gnomad", "rsid_identity"])
@pytest.mark.parametrize(
    "shape,named,gt",
    [
        ("plain", True, "0/1"),
        ("prefix", False, "0/1"),
        ("prefix", True, "1"),
        ("plain", False, "1/1"),
        ("prefix", True, "0/0"),
        ("plain", False, "0/."),
        ("prefix", False, "./."),
        ("plain", True, "./1"),
    ],
)
def test_source_isolation_survives_call_representation(tmp_path, build, source, shape, named, gt):
    data, _ = _prepare(tmp_path, build)
    path = call(
        tmp_path / "input.vcf", POSITIONS[build], rsid=RSID if named else ".", shape=shape, gt=gt
    )
    _pair(tmp_path, data, path, build, "analyze", source, "config")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("format_name", ["vcf", "23andme"])
def test_disabled_identity_cannot_reject_a_contradictory_placement(tmp_path, build, format_name):
    data, _ = _prepare(tmp_path, build)
    pos = POSITIONS[build] + 10
    if format_name == "vcf":
        path = call(tmp_path / "input.vcf", pos)
    else:
        path = tmp_path / "input.txt"
        path.write_text(
            "# This data file generated by 23andMe\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            f"{RSID}\t1\t{pos}\tAG\n"
        )
    _pair(tmp_path, data, path, build, "analyze", "rsid_identity", "config")

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Diff real synthetic-cache annotations by finding identity, not generated prose."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION, PHARMGKB_INTERPRETER_VERSION
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GWAS_SCHEMA, PHARMGKB_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.models import Variant, annotation_to_public_dict
from allelix.reports.diff import compute_diff, load_previous_report

RSID = "rs900252301"
SOURCES = ("clinvar", "snpedia", "pharmgkb", "gwas")
BUILDS = ("GRCh37", "GRCh38")


def _cache(
    data,
    source,
    build,
    *,
    revised=False,
    alt="G",
    condition="Synthetic condition",
    aggregate=None,
    rsid=RSID,
    position=None,
):
    """Write loader-shaped rows; every finding below comes from an actual annotator."""
    data.mkdir(parents=True, exist_ok=True)
    pos = position if position is not None else (1000 if build == "GRCh37" else 2000)
    if source == "clinvar":
        filename, schema = clinvar_db_filename(build), CLINVAR_SCHEMA
        name, tag = clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"
    elif source == "pharmgkb":
        filename, schema = "pharmgkb.sqlite", PHARMGKB_SCHEMA
        name, tag = source, f"iv:{PHARMGKB_INTERPRETER_VERSION}"
    elif source == "gwas":
        filename, schema = "gwas.sqlite", GWAS_SCHEMA
        name, tag = source, f"cv:{_CATEGORIZER_VERSION}"
    else:
        filename, schema = "snpedia.sqlite", _STRUCTURED_SCHEMA
        name, tag = source, f"pv:{_PARSER_VERSION}"
    with contextlib.closing(sqlite3.connect(data / filename)) as conn:
        conn.executescript(schema)
        if source == "clinvar":
            conn.execute(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id, "
                "aggregate, submitters_same, submitters_total) "
                "VALUES (?, '1', ?, 'A', ?, ?, ?, 'MTHFR', "
                "'criteria provided, single submitter', 900252301, ?, 1, 2)",
                (
                    rsid,
                    pos,
                    alt,
                    "Likely pathogenic" if revised else "Pathogenic",
                    condition,
                    aggregate,
                ),
            )
        elif source == "snpedia":
            conn.execute(
                "INSERT INTO snpedia_genotypes "
                "(rsid, allele1, allele2, magnitude, repute, summary, gene, orientation) "
                "VALUES (?, 'A', ?, 6.0, 'Bad', ?, 'MTHFR', 'plus')",
                (
                    rsid,
                    alt,
                    "Revised synthetic interpretation"
                    if revised
                    else "Original synthetic interpretation",
                ),
            )
        elif source == "pharmgkb":
            conn.execute(
                "INSERT INTO pharmgkb_annotations "
                "(rsid, genotype, gene, drugs, phenotype, phenotype_category, annotation_text, "
                "level_of_evidence, score, pgkb_annotation_id, function_class, is_nonfinding) "
                "VALUES (?, ?, 'MTHFR', 'Synthetic drug', ?, 'Efficacy', ?, '1A', 10, "
                "'900252301', 'decreased', 0)",
                (
                    rsid,
                    "A" + alt,
                    condition,
                    "Revised synthetic interpretation"
                    if revised
                    else "Original synthetic interpretation",
                ),
            )
        else:
            conn.execute(
                "INSERT INTO gwas_associations "
                "(rsid, risk_allele, trait, p_value, or_beta, ci_text, gene, study_accession, "
                "pubmed_id, trait_category) VALUES (?, ?, ?, 1e-30, 1.5, '[1.2-1.8]', ?, "
                "'GCST900252301', '900252301', 'disease')",
                (rsid, alt, condition, "GENE_REVISED" if revised else "GENE_ORIGINAL"),
            )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count, local_version_tag) "
            "VALUES (?, 'test://synthetic', 'synthetic', '2001-01-01', 1, ?)",
            (name, tag),
        )
        conn.commit()
    if source in ("snpedia", "gwas"):
        write_identity_db(
            [(rsid, "1", 1000, "A", alt, 2000, "A", alt, 0)],
            data / RSID_IDENTITY_DB_FILENAME,
            version="synthetic",
            source_url="test://identity",
            downloaded_at="2001-01-01",
        )
    return Variant(rsid, "1", pos, "A", alt, build=build, ref="A", ref_from_input=True)


def _findings(data, source, build, *, batch=True, **kwargs):
    variant = _cache(data, source, build, **kwargs)
    with contextlib.ExitStack() as stack:
        if source in ("snpedia", "gwas"):
            identity = stack.enter_context(RsidIdentityAnnotator(data))
            cls = SNPediaAnnotator if source == "snpedia" else GWASCatalogAnnotator
            annotator = cls(data, forward_site_provider=identity.forward_sites)
        else:
            annotator = (
                ClinVarAnnotator(data, builds=(build,))
                if source == "clinvar"
                else PharmGKBAnnotator(data)
            )
        stack.enter_context(annotator)
        assert annotator.is_ready(), annotator.cache_status()
        findings = (
            list(annotator.batch_annotate([variant])) if batch else annotator.annotate(variant)
        )
    assert len(findings) == 1, (source, build, findings)
    assert findings[0].alt == kwargs.get("alt", "G")
    return findings


def _counts(diff):
    return len(diff.changed), len(diff.new), len(diff.removed)


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
def test_generated_prose_and_reclassification_are_changed_not_replaced(
    tmp_path, source, build, batch
):
    old = _findings(tmp_path / "old", source, build, batch=batch)
    new = _findings(tmp_path / "new", source, build, batch=batch, revised=True)
    assert old[0].description != new[0].description
    assert (old[0].rsid, old[0].condition, old[0].alt) == (
        new[0].rsid,
        new[0].condition,
        new[0].alt,
    )
    if source == "clinvar":
        assert old[0].significance == "clinvar_pathogenic"
        assert new[0].significance == "clinvar_likely_pathogenic"
    else:
        # These changes are prose-only, so dropping description from the key alone is insufficient.
        assert (old[0].significance, old[0].magnitude) == (new[0].significance, new[0].magnitude)
    for previous, current in ((old, new), (new, old)):
        diff = compute_diff(
            current, [annotation_to_public_dict(a) for a in previous], "2001-01-01"
        )
        assert _counts(diff) == (1, 0, 0)
        assert diff.changed[0].current == current[0]


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("build", BUILDS)
def test_genuinely_different_allele_is_new_and_removed(tmp_path, source, build):
    old = _findings(tmp_path / "old", source, build, alt="G")
    new = _findings(tmp_path / "new", source, build, alt="C")
    diff = compute_diff(new, [annotation_to_public_dict(old[0])], "2001-01-01")
    assert _counts(diff) == (0, 1, 1)
    assert diff.new[0].alt == "C" and diff.removed[0]["alt"] == "G"


@pytest.mark.parametrize("source", ["clinvar", "pharmgkb", "gwas"])
@pytest.mark.parametrize("build", BUILDS)
def test_genuinely_different_condition_is_new_and_removed(tmp_path, source, build):
    old = _findings(tmp_path / "old", source, build)
    new = _findings(tmp_path / "new", source, build, condition="Another synthetic condition")
    diff = compute_diff(new, [annotation_to_public_dict(old[0])], "2001-01-01")
    assert _counts(diff) == (0, 1, 1)


@pytest.mark.parametrize("source", SOURCES)
def test_identical_actual_annotator_output_is_quiet(tmp_path, source):
    rows = _findings(tmp_path / "data", source, "GRCh38")
    assert _counts(compute_diff(rows, [annotation_to_public_dict(rows[0])], "2001-01-01")) == (
        0,
        0,
        0,
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "aggregate,expected",
    [
        ("Conflicting classifications of pathogenicity", 4.0),
        ("Benign", 1.0),
        ("Likely benign", 2.0),
    ],
)
def test_schema8_uncapped_baseline_is_changed_when_aggregate_caps_the_same_row(
    tmp_path, build, aggregate, expected
):
    # With no aggregate metadata the annotator produces the uncapped Pathogenic 9.0
    # row written by 2.4.1. Keep its real emitted description, allele and genotype.
    old = _findings(tmp_path / "old", "clinvar", build)
    new = _findings(tmp_path / "new", "clinvar", build, aggregate=aggregate)
    assert old[0].magnitude == 9.0 and not old[0].aggregate_capped
    assert new[0].magnitude == expected and new[0].aggregate_capped
    previous = annotation_to_public_dict(old[0])
    for key in (
        "clinvar_aggregate",
        "submitters_same",
        "submitters_total",
        "aggregate_capped",
        "outside_aggregate",
        "common_allele_capped",
        "clinvar_allele_af",
    ):
        previous.pop(key, None)
    baseline = tmp_path / "baseline.json"
    baseline.write_text(
        json.dumps({"schema_version": "8", "allelix_version": "2.4.1", "annotations": [previous]}),
        encoding="utf-8",
    )
    loaded = load_previous_report(baseline)
    diff = compute_diff(new, loaded["annotations"], "2001-01-01")
    assert _counts(diff) == (1, 0, 0)
    assert diff.changed[0].previous_magnitude == 9.0
    assert diff.changed[0].current.aggregate_capped


def _input(path, build="GRCh37"):
    pos = 1000 if build == "GRCh37" else 2000
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{pos}\t{RSID}\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )


def _args(command, source, data, build="GRCh37"):
    return [
        command,
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build.lower(),
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--snpedia-min-magnitude",
        "0",
    ]


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("output_kind", ["json", "html", "terminal"])
@pytest.mark.parametrize("same_input", [False, True], ids=["different-input", "same-input"])
def test_cli_diff_discloses_a_different_baseline_input(tmp_path, command, output_kind, same_input):
    data = tmp_path / "data"
    _cache(data, "pharmgkb", "GRCh37")
    source = tmp_path / "current.vcf"
    _input(source)
    previous_source = source if same_input else tmp_path / "other.vcf"
    if not same_input:
        _input(previous_source)
    baseline = tmp_path / "baseline.json"
    runner = CliRunner()
    first = runner.invoke(
        main, [*_args(command, previous_source, data), "--output", str(baseline)]
    )
    assert first.exit_code == 0, first.output
    assert json.loads(baseline.read_text())["input"]["file"] == previous_source.name
    args = [*_args(command, source, data), "--diff", str(baseline)]
    if output_kind != "terminal":
        args += ["--output", str(tmp_path / f"report.{output_kind}")]
    result = runner.invoke(main, args)
    assert result.exit_code == 0, result.output
    warnings = [
        line
        for line in result.output.splitlines()
        if "warn" in line.lower() and source.name in line and previous_source.name in line
    ]
    assert bool(warnings) is not same_input, result.output


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("output_kind", ["json", "html", "terminal"])
def test_cli_reclassification_uses_actual_annotations_in_every_diff_renderer(
    tmp_path, build, output_kind
):
    data = tmp_path / "data"
    _cache(data, "clinvar", build)
    source = tmp_path / "sample.vcf"
    _input(source, build)
    runner = CliRunner()
    baseline = tmp_path / "baseline.json"
    first = runner.invoke(
        main, [*_args("analyze", source, data, build), "--output", str(baseline)]
    )
    assert first.exit_code == 0, first.output
    assert len(json.loads(baseline.read_text())["annotations"]) == 1
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.execute("UPDATE clinvar_variants SET clinical_significance = 'Likely pathogenic'")
        conn.commit()
    args = [*_args("analyze", source, data, build), "--diff", str(baseline)]
    output = tmp_path / f"changed.{output_kind}"
    if output_kind != "terminal":
        args += ["--output", str(output)]
    result = runner.invoke(main, args)
    assert result.exit_code == 0, result.output
    if output_kind == "json":
        diff = json.loads(output.read_text())["diff"]
        assert (len(diff["changed"]), len(diff["new"]), len(diff["removed"])) == (1, 0, 0)
    else:
        text = result.output if output_kind == "terminal" else output.read_text()
        assert "1 changed" in text and "1 new" not in text and "1 removed" not in text


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("build", BUILDS)
def test_genuinely_different_variant_label_is_new_and_removed(tmp_path, source, build):
    old = _findings(tmp_path / "old", source, build)
    new = _findings(tmp_path / "new", source, build, rsid="rs900252302")
    assert _counts(compute_diff(new, [annotation_to_public_dict(old[0])], "2001-01-01")) == (
        0,
        1,
        1,
    )


@pytest.mark.parametrize("build", BUILDS)
def test_unnamed_clinvar_reclassification_keeps_locus_identity(tmp_path, build):
    old = _findings(tmp_path / "old", "clinvar", build, rsid="")
    new = _findings(tmp_path / "new", "clinvar", build, rsid="", revised=True)
    assert old[0].rsid == new[0].rsid == ""
    assert _counts(compute_diff(new, [annotation_to_public_dict(old[0])], "2001-01-01")) == (
        1,
        0,
        0,
    )


@pytest.mark.parametrize("build", BUILDS)
def test_different_unnamed_loci_stay_new_and_removed(tmp_path, build):
    old = _findings(tmp_path / "old", "clinvar", build, rsid="", position=3000)
    new = _findings(tmp_path / "new", "clinvar", build, rsid="", position=4000)
    assert _counts(compute_diff(new, [annotation_to_public_dict(old[0])], "2001-01-01")) == (
        0,
        1,
        1,
    )


def _schema8_cli_baseline(tmp_path, build):
    """A real CLI-produced uncapped finding in the older report envelope."""
    data, source, baseline = (tmp_path / name for name in ("data", "sample.vcf", "baseline.json"))
    _cache(data, "clinvar", build)
    _input(source, build)
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build.lower(),
    ]
    result = CliRunner().invoke(main, [*args, "--output", str(baseline)])
    assert result.exit_code == 0, result.output
    payload = json.loads(baseline.read_text())
    assert len(payload["annotations"]) == 1
    assert payload["annotations"][0]["magnitude"] == 9.0
    payload["schema_version"] = "8"
    payload["allelix_version"] = "2.4.1"
    for row in payload["annotations"]:
        for key in (
            "clinvar_aggregate",
            "submitters_same",
            "submitters_total",
            "aggregate_capped",
            "outside_aggregate",
            "common_allele_capped",
            "clinvar_allele_af",
        ):
            row.pop(key, None)
    baseline.write_text(json.dumps(payload), encoding="utf-8")
    return data, baseline, args


def _visible_text(document):
    """Check findings, excluding styles, scripts and the static scoring tutorial."""
    from html.parser import HTMLParser

    class Text(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self.hidden = 0
            self.education_depth = 0
            self.parts = []

        def handle_starttag(self, tag, attrs):
            if tag == "details":
                if self.education_depth:
                    self.education_depth += 1
                elif "education" in (dict(attrs).get("class") or "").split():
                    self.education_depth = 1
            if tag in ("script", "style"):
                self.hidden += 1

        def handle_endtag(self, tag):
            if tag == "details" and self.education_depth:
                self.education_depth -= 1
            if tag in ("script", "style"):
                self.hidden -= 1

        def handle_data(self, data):
            if not self.hidden and not self.education_depth:
                self.parts.append(data)

    parser = Text()
    parser.feed(document)
    parser.close()
    return " ".join(" ".join(parser.parts).split())


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("output_kind", ["json", "html", "terminal"])
@pytest.mark.parametrize("floor", [None, 0], ids=["default-floor", "floor-zero"])
def test_schema8_capped_finding_stays_changed_below_the_display_floor(
    tmp_path, build, output_kind, floor
):
    """A display filter cannot turn an extant, scored-down finding into a removal."""
    import re

    data, baseline, args = _schema8_cli_baseline(tmp_path, build)
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.execute(
            "UPDATE clinvar_variants SET aggregate = "
            "'Conflicting classifications of pathogenicity'"
        )
        conn.commit()
    pos = 1000 if build == "GRCh37" else 2000
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        rows = annotator.annotate(
            Variant(RSID, "1", pos, "A", "G", build=build, ref="A", ref_from_input=True)
        )
    assert len(rows) == 1 and rows[0].magnitude == 4.0 and rows[0].aggregate_capped
    args += ["--diff", str(baseline)]
    if floor is not None:
        args += ["--min-magnitude", str(floor)]
    output = tmp_path / f"report.{output_kind}"
    if output_kind != "terminal":
        args += ["--output", str(output)]
    result = CliRunner().invoke(main, args, env={"COLUMNS": "240", "NO_COLOR": "1"})
    assert result.exit_code == 0, result.output
    assert "scored down" in " ".join(result.output.split()).lower()
    if output_kind == "json":
        payload = json.loads(output.read_text())
        diff = payload["diff"]
        assert (len(diff["changed"]), len(diff["new"]), len(diff["removed"])) == (1, 0, 0)
        changed = diff["changed"][0]
        assert changed["rsid"] == RSID and changed["alt"] == "G"
        assert changed["previous_magnitude"] == 9.0 and changed["magnitude"] == 4.0
        assert changed["aggregate_capped"] is True
        # An explicit per-entry flag makes the display status machine-readable.
        assert changed.get("below_display_floor", False) is (floor is None)
        assert len(payload["annotations"]) == (0 if floor is None else 1)
        assert payload["clinvar"]["pathogenic_rows_scored_down"] == 1
    else:
        text = (
            _visible_text(output.read_text())
            if output_kind == "html"
            else " ".join(result.output.split())
        )
        assert "1 changed" in text and "1 new" not in text and "1 removed" not in text
        # Existing HTML keeps the rsID in its detail payload; either location is valid.
        assert RSID in (output.read_text() if output_kind == "html" else text)
        assert re.search(r"(?<![\d.])9(?:\.0)?(?![\d.])", text), text
        assert re.search(r"(?<![\d.])4(?:\.0)?(?![\d.])", text), text
        below = bool(re.search(r"below\s+(?:the\s+)?display\s+floor", text, re.IGNORECASE))
        assert below is (floor is None), text


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("output_kind", ["json", "html", "terminal"])
def test_schema8_truly_absent_identity_is_removed_at_default_filters(tmp_path, build, output_kind):
    """Keep a ready cache, but move its only finding away from the input call."""
    data, baseline, args = _schema8_cli_baseline(tmp_path, build)
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.execute("UPDATE clinvar_variants SET position = position + 1")
        conn.commit()
    pos = 1000 if build == "GRCh37" else 2000
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        assert annotator.is_ready()
        assert (
            annotator.annotate(
                Variant(RSID, "1", pos, "A", "G", build=build, ref="A", ref_from_input=True)
            )
            == []
        )
    args += ["--diff", str(baseline)]
    output = tmp_path / f"absent.{output_kind}"
    if output_kind != "terminal":
        args += ["--output", str(output)]
    result = CliRunner().invoke(main, args, env={"COLUMNS": "240", "NO_COLOR": "1"})
    assert result.exit_code == 0, result.output
    if output_kind == "json":
        payload = json.loads(output.read_text())
        diff = payload["diff"]
        assert (len(diff["changed"]), len(diff["new"]), len(diff["removed"])) == (0, 0, 1)
        assert diff["removed"][0]["rsid"] == RSID
        assert diff["removed"][0]["magnitude"] == 9.0
        assert payload["annotations"] == []
    else:
        text = (
            _visible_text(output.read_text())
            if output_kind == "html"
            else " ".join(result.output.split())
        )
        assert "1 removed" in text and "1 changed" not in text and "1 new" not in text

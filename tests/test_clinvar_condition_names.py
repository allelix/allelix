# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Decode multi-condition ClinVar submissions without stripping genuine name text."""

from __future__ import annotations

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.databases.manager import (
    _SS_REQUIRED_COLUMNS,
    _VS_REQUIRED_COLUMNS,
    _decode_reported_phenotype,
    iter_clinvar_tsv_records,
    load_clinvar_tsv,
)
from allelix.models import Variant

# The F5 multi-condition shape, with an invented first identifier and a synthetic variant.
NAMES = ("Congenital factor V deficiency", "Budd-Chiari syndrome", "Thrombophilia")
TOKENS = ("C900252801:" + NAMES[0], "C0856761:" + NAMES[1], "C900252802:" + NAMES[2])
EXPECTED = "; ".join(NAMES)
RAW_CASES = (
    ("|".join(TOKENS), EXPECTED),
    (";".join(TOKENS), EXPECTED),
    (TOKENS[0] + ";" + TOKENS[1] + "|" + TOKENS[2], EXPECTED),
    (TOKENS[0] + "|" + TOKENS[1] + ";" + TOKENS[2], EXPECTED),
    (" ; ".join(TOKENS), EXPECTED),
    ("na:Disorder: childhood onset", "Disorder: childhood onset"),
    ("Disorder: childhood onset", "Disorder: childhood onset"),
    ("C3 deficiency: severe", "C3 deficiency: severe"),
    ("na:C3 deficiency", "C3 deficiency"),
    ("C3 deficiency", "C3 deficiency"),
    ("C1234567-like syndrome", "C1234567-like syndrome"),
    ("na:C1234567-like syndrome", "C1234567-like syndrome"),
    ("Syndrome with C1234567 marker", "Syndrome with C1234567 marker"),
    ("C1|C2:Condition A;C3;na:Condition B", "Condition A; Condition B"),
    ("C1:-;C2:Condition A|-", "Condition A"),
    ("", ""),
    ("-", ""),
)


@pytest.mark.parametrize("raw,expected", RAW_CASES)
def test_only_identifier_tokens_are_removed(raw, expected):
    assert _decode_reported_phenotype(raw) == expected


def _sources(tmp_path, raw):
    vs = tmp_path / "variant_summary.tsv"
    ss = tmp_path / "submission_summary.tsv"
    header = [_VS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(43)]
    header[0] = "#AlleleID"
    lines = ["\t".join(header)]
    for build, pos in (("GRCh37", 1000), ("GRCh38", 2000)):
        row = [""] * 43
        values = {
            0: "900252801",
            1: "single nucleotide variant",
            2: "NM_SYNTHETIC(F5):c.1A>G",
            4: "F5",
            6: "Pathogenic",
            9: "900252801",
            16: build,
            18: "1",
            24: "criteria provided, single submitter",
            30: "900252801",
            31: str(pos),
            32: "A",
            33: "G",
        }
        for index, value in values.items():
            row[index] = value
        lines.append("\t".join(row))
    vs.write_text("\n".join(lines) + "\n")
    header = [_SS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(16)]
    header[0] = "#VariationID"
    row = [""] * 16
    for index, value in {
        0: "900252801",
        1: "Pathogenic",
        5: raw,
        6: "criteria provided, single submitter",
        10: "SCV900252801.1",
        15: "yes",
    }.items():
        row[index] = value
    ss.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")
    return vs, ss


@pytest.mark.parametrize("raw,expected", RAW_CASES[:15])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_tsv_cache_and_actual_annotation_preserve_condition_names(tmp_path, raw, expected, build):
    vs, ss = _sources(tmp_path, raw)
    records = list(iter_clinvar_tsv_records(vs, ss, build))
    assert len(records) == 1
    data = tmp_path / "data"
    data.mkdir()
    load_clinvar_tsv(
        vs, ss, data / clinvar_db_filename(build), build, record_name=clinvar_record_name(build)
    )
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        assert annotator.is_ready(), annotator.cache_status()
        variant = Variant(
            "rs900252801",
            "1",
            1000 if build == "GRCh37" else 2000,
            "A",
            "G",
            ref="A",
            build=build,
            ref_from_input=True,
        )
        rows = list(annotator.batch_annotate([variant]))
    assert len(rows) == 1, rows
    assert rows[0].condition == expected
    assert records[0]["condition"] == expected


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_f5_shaped_multi_condition_names_reach_reports_without_ids(tmp_path, build, surface):
    import json

    from click.testing import CliRunner

    from allelix.cli import main

    vs, ss = _sources(tmp_path, ";".join(TOKENS))
    data = tmp_path / "data"
    data.mkdir()
    load_clinvar_tsv(
        vs, ss, data / clinvar_db_filename(build), build, record_name=clinvar_record_name(build)
    )
    source = tmp_path / "input.vcf"
    pos = 1000 if build == "GRCh37" else 2000
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{pos}\trs900252801\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    )
    args = [
        "analyze",
        str(source),
        "--build",
        build.lower(),
        "--data-dir",
        str(data),
        "--no-update",
    ]
    output = tmp_path / f"report.{surface}"
    if surface != "terminal":
        args += ["--output", str(output)]
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
    assert result.exit_code == 0, result.output
    if surface == "json":
        rows = [
            row
            for row in json.loads(output.read_text())["annotations"]
            if row["source"] == "clinvar"
        ]
        assert len(rows) == 1
        assert rows[0]["condition"] == EXPECTED
    else:
        text = output.read_text() if surface == "html" else result.output
        for name in NAMES:
            assert name in text
        assert "C0856761" not in text and "C900252802" not in text


# Literal output of the old loader: it removed only the first identifier.
# Do not derive this from the decoder under test, or a loader fix would hide read-time defects.
LEGACY_CONDITION = NAMES[0] + ";" + TOKENS[1] + ";" + TOKENS[2]


def _existing_cache(tmp_path, build, stored_condition):
    import contextlib
    import sqlite3

    vs, ss = _sources(tmp_path, "na:Synthetic placeholder")
    data = tmp_path / "data"
    data.mkdir()
    cache = data / clinvar_db_filename(build)
    load_clinvar_tsv(vs, ss, cache, build, record_name=clinvar_record_name(build))
    with contextlib.closing(sqlite3.connect(cache)) as conn:
        conn.execute("UPDATE clinvar_variants SET condition = ?", (stored_condition,))
        # The released cache's stamp is deliberate: these reads must need neither rebuild nor bump.
        conn.execute("UPDATE database_versions SET local_version_tag = 'iv:9'")
        conn.commit()
    source = tmp_path / "input.vcf"
    pos = 1000 if build == "GRCh37" else 2000
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{pos}\trs900252801\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    )
    return data, cache, source


def _read_existing_cache(tmp_path, monkeypatch, build, surface, stored_condition):
    from click.testing import CliRunner

    from allelix.cli import main

    data, cache, source = _existing_cache(tmp_path, build, stored_condition)
    before = (cache.read_bytes(), cache.stat().st_mtime_ns)

    def no_rebuild(*args, **kwargs):
        pytest.fail("Reading an existing ClinVar cache must not rebuild it")

    monkeypatch.setattr(ClinVarAnnotator, "setup", no_rebuild)
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        assert annotator.is_ready(), annotator.cache_status()
    args = [
        "analyze",
        str(source),
        "--build",
        build.lower(),
        "--data-dir",
        str(data),
        "--no-update",
    ]
    output = tmp_path / f"report.{surface}"
    if surface != "terminal":
        args += ["--output", str(output)]
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert (cache.read_bytes(), cache.stat().st_mtime_ns) == before
    return result, output, args, cache


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
@pytest.mark.parametrize(
    "stored,expected",
    [
        (LEGACY_CONDITION, EXPECTED),
        (EXPECTED, EXPECTED),
        (
            "C3 deficiency: severe; C1234567-like syndrome",
            "C3 deficiency: severe; C1234567-like syndrome",
        ),
    ],
    ids=["old-id-laden-cache", "already-clean-cache", "legitimate-colon-and-c-word"],
)
def test_existing_cache_names_are_clean_without_rebuild(
    tmp_path, monkeypatch, build, surface, stored, expected
):
    import json

    result, output, _, _ = _read_existing_cache(tmp_path, monkeypatch, build, surface, stored)
    if surface == "json":
        rows = [
            row
            for row in json.loads(output.read_text())["annotations"]
            if row["source"] == "clinvar"
        ]
        assert len(rows) == 1
        assert rows[0]["condition"] == expected
    else:
        text = output.read_text() if surface == "html" else result.output
        for name in expected.split("; "):
            assert name in text
        if stored == LEGACY_CONDITION:
            assert "C0856761" not in text and "C900252802" not in text


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_old_baseline_condition_and_clean_current_row_keep_one_diff_identity(
    tmp_path, monkeypatch, build
):
    import json

    from click.testing import CliRunner

    from allelix.cli import main

    # Generate the actual current annotation, then put only its baseline condition back
    # into the literal old-cache spelling. Classification, allele and input stay identical.
    _, output, args, cache = _read_existing_cache(tmp_path, monkeypatch, build, "json", EXPECTED)
    baseline = json.loads(output.read_text())
    assert len(baseline["annotations"]) == 1
    assert baseline["annotations"][0]["source"] == "clinvar"
    assert baseline["annotations"][0]["condition"] == EXPECTED
    baseline["schema_version"] = "10"
    baseline["allelix_version"] = "2.5.1"
    baseline["annotations"][0]["condition"] = LEGACY_CONDITION
    previous = tmp_path / "baseline-2.5.1.json"
    previous.write_text(json.dumps(baseline))
    previous_bytes = previous.read_bytes()
    cache_before = (cache.read_bytes(), cache.stat().st_mtime_ns)
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, [*args, "--diff", str(previous)])
    assert result.exit_code == 0, (result.output, result.exception)
    assert previous.read_bytes() == previous_bytes
    assert (cache.read_bytes(), cache.stat().st_mtime_ns) == cache_before
    report = json.loads(output.read_text())
    assert len(report["annotations"]) == 1
    assert report["annotations"][0]["condition"] == EXPECTED
    diff = report["diff"]
    assert diff["new"] == [] and diff["removed"] == [], diff
    assert len(diff["changed"]) in (0, 1), diff
    if diff["changed"]:
        assert diff["changed"][0]["condition"] == EXPECTED

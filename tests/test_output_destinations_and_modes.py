# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Output destinations are distinct before work starts; published modes honor the user."""

from __future__ import annotations

import contextlib
import os
import sqlite3
import stat
from itertools import combinations
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import _helpers, main, utility
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.exporters.plink import default_plink_prefix, export_plink, plink_paths
from allelix.models import Variant
from allelix.reports import atomic_write_text
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from allelix.reports.vcf import render_vcf

VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    "1\t1000\trs900252401\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
)
PAIRS = list(combinations((".fam", ".bim", ".bed"), 2))


def _ready_cache(data):
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.execute(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, allele_id) "
            "VALUES ('rs900252401', '1', 1000, 'A', 'G', 'Pathogenic', "
            "'Synthetic condition', 'MTHFR', 900252401)"
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count, local_version_tag) "
            "VALUES (?, 'test://synthetic', 'synthetic', '2001-01-01', 1, ?)",
            (clinvar_record_name("GRCh37"), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _snapshot(root):
    """Catch writes, replacements, mode changes, created files and altered links."""
    result = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            result[str(path.relative_to(root))] = ("symlink", os.readlink(path))
        elif path.is_file():
            info = path.stat()
            result[str(path.relative_to(root))] = (
                path.read_bytes(),
                info.st_ino,
                info.st_mtime_ns,
                stat.S_IMODE(info.st_mode),
            )
        else:
            result[str(path.relative_to(root))] = ("directory",)
    return result


def _forbid_work(monkeypatch):
    calls = []

    def forbidden(*args, **kwargs):
        calls.append("analysis/cache setup started")
        raise AssertionError("Conflicting destinations reached analysis/cache setup")

    monkeypatch.setattr(_helpers, "resolve_data_dir", forbidden)
    monkeypatch.setattr(utility, "resolve_data_dir", forbidden)
    monkeypatch.setattr(_helpers, "run_analysis", forbidden)
    return calls


def _assert_refused(result, calls, root, before):
    assert calls == [], result.output
    assert result.exit_code != 0, result.output
    assert not isinstance(result.exception, AssertionError), result.output
    assert "error" in result.output.lower(), result.output
    assert _snapshot(root) == before


@pytest.mark.parametrize("report_format", ["json", "html"])
@pytest.mark.parametrize("alias", ["same", "relative", "symlink", "parent-symlink", "hardlink"])
@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reverse"])
def test_report_and_vcf_output_aliases_refused_before_analysis(
    tmp_path, monkeypatch, report_format, alias, reverse
):
    source = tmp_path / "input.vcf"
    source.write_text(VCF, encoding="utf-8")
    first = tmp_path / f"report.{report_format}"
    if alias == "same":
        second = first
    elif alias == "relative":
        monkeypatch.chdir(tmp_path)
        second = Path(first.name)
    elif alias == "symlink":
        # A dangling link is an alias even before either destination exists.
        second = tmp_path / "alias.vcf"
        second.symlink_to(first.name)
    elif alias == "parent-symlink":
        link = tmp_path / "directory-link"
        link.symlink_to(tmp_path, target_is_directory=True)
        second = link / first.name
    else:
        first.write_bytes(b"existing output must survive")
        second = tmp_path / "alias.vcf"
        os.link(first, second)
    report, annotated = (second, first) if reverse else (first, second)
    before = _snapshot(tmp_path)
    calls = _forbid_work(monkeypatch)
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--output",
            str(report),
            "--report-format",
            report_format,
            "--vcf-out",
            str(annotated),
        ],
    )
    _assert_refused(result, calls, tmp_path, before)


@pytest.mark.parametrize("pair", PAIRS)
@pytest.mark.parametrize("alias", ["symlink", "dangling-symlink", "hardlink"])
@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize(
    "explicit_prefix", [False, True], ids=["default-prefix", "explicit-prefix"]
)
def test_every_plink_output_pair_is_checked_before_analysis(
    tmp_path, monkeypatch, pair, alias, reverse, explicit_prefix
):
    source = tmp_path / "sample.vcf"
    source.write_text(VCF, encoding="utf-8")
    prefix = tmp_path / "output.v1" if explicit_prefix else default_plink_prefix(source)
    paths = plink_paths(prefix)
    a, b = (paths[ext] for ext in pair)
    if reverse:
        a, b = b, a
    if alias != "dangling-symlink":
        a.write_bytes(b"prior output")
    if alias == "hardlink":
        os.link(a, b)
    else:
        b.symlink_to(a.name)
    before = _snapshot(tmp_path)
    calls = _forbid_work(monkeypatch)
    args = ["export", "plink", str(source)]
    if explicit_prefix:
        args += ["--output-prefix", str(prefix)]
    _assert_refused(CliRunner().invoke(main, args), calls, tmp_path, before)


@pytest.mark.parametrize("pair", PAIRS)
@pytest.mark.parametrize("alias", ["symlink", "hardlink"])
def test_direct_plink_writer_refuses_collisions_before_consuming_input(tmp_path, pair, alias):
    prefix = tmp_path / "output"
    paths = plink_paths(prefix)
    first, second = (paths[ext] for ext in pair)
    first.write_bytes(b"prior output")
    if alias == "hardlink":
        os.link(first, second)
    else:
        second.symlink_to(first.name)
    before = _snapshot(tmp_path)
    consumed = []

    def variants():
        consumed.append(True)
        yield Variant("rs900252401", "1", 1000, "A", "G")

    with pytest.raises(ValueError):
        export_plink(variants(), prefix, "GRCh37")
    assert consumed == []
    assert _snapshot(tmp_path) == before


def _write(writer, source, destination):
    if writer == "plink":
        export_plink(iter([Variant("rs900252401", "1", 1000, "A", "G")]), destination, "GRCh37")
    elif writer == "vcf":
        render_vcf(
            input_path=source,
            output_path=destination,
            build="GRCh37",
            annotators_used=[],
            annotations_by_key={},
        )
    elif writer == "atomic-text":
        atomic_write_text(destination, "synthetic report\n", input_path=source)
    else:
        result = AnalysisResult(
            file_path=source,
            parser_name="vcf",
            parser_display_name="VCF",
            sample_id="SAMPLE",
            build="GRCh37",
            total_variants=1,
            skipped_count=0,
            annotators_used=[],
        )
        renderer = render_json if writer == "json" else render_html
        renderer(result, output_path=destination)


@pytest.mark.parametrize("writer", ["json", "html", "vcf", "plink", "atomic-text"])
@pytest.mark.parametrize("mask", [0o022, 0o077], ids=["umask022", "umask077"])
@pytest.mark.parametrize(
    "prior_mode", [None, 0o640, 0o664], ids=["new", "existing640", "existing664"]
)
def test_published_outputs_follow_umask_or_keep_existing_mode(tmp_path, writer, mask, prior_mode):
    source, destination = tmp_path / "input.vcf", tmp_path / "report"
    source.write_text(VCF, encoding="utf-8")
    outputs = list(plink_paths(destination).values()) if writer == "plink" else [destination]
    if prior_mode is not None:
        for path in outputs:
            path.write_bytes(b"prior report")
            path.chmod(prior_mode)
    previous_mask = os.umask(mask)
    try:
        _write(writer, source, destination)
    finally:
        os.umask(previous_mask)
    expected = prior_mode if prior_mode is not None else 0o666 & ~mask
    assert {path.name: stat.S_IMODE(path.stat().st_mode) for path in outputs} == {
        path.name: expected for path in outputs
    }
    assert all(path.read_bytes() != b"prior report" for path in outputs)
    assert source.read_text(encoding="utf-8") == VCF
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("extension", ["json", "html"])
@pytest.mark.parametrize("mask", [0o022, 0o077], ids=["umask022", "umask077"])
@pytest.mark.parametrize("prior_mode", [None, 0o640], ids=["new", "existing640"])
def test_cli_report_modes(tmp_path, command, extension, mask, prior_mode):
    source, output, data = (
        tmp_path / "input.vcf",
        tmp_path / f"report.{extension}",
        tmp_path / "data",
    )
    source.write_text(VCF, encoding="utf-8")
    _ready_cache(data)
    if prior_mode is not None:
        output.write_bytes(b"prior report")
        output.chmod(prior_mode)
    previous_mask = os.umask(mask)
    try:
        result = CliRunner().invoke(
            main,
            [
                command,
                str(source),
                "--no-update",
                "--build",
                "grch37",
                "--data-dir",
                str(data),
                "--output",
                str(output),
            ],
        )
    finally:
        os.umask(previous_mask)
    assert result.exit_code == 0, result.output
    assert stat.S_IMODE(output.stat().st_mode) == (
        prior_mode if prior_mode is not None else 0o666 & ~mask
    )
    assert source.read_text(encoding="utf-8") == VCF


def test_distinct_report_and_vcf_destinations_remain_supported(tmp_path):
    source, report, annotated, data = (
        tmp_path / name for name in ("input.vcf", "report.json", "annotated.vcf", "data")
    )
    source.write_text(VCF, encoding="utf-8")
    _ready_cache(data)
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--no-update",
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--output",
            str(report),
            "--vcf-out",
            str(annotated),
        ],
    )
    assert result.exit_code == 0, result.output
    assert report.read_text().lstrip().startswith("{")
    assert annotated.read_text().startswith("##fileformat=VCFv4.2")
    assert source.read_text(encoding="utf-8") == VCF

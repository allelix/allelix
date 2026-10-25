# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Input aliases must be refused before cache or output mutation."""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import _helpers, main, utility

_PAYLOAD = (Path(__file__).parent / "fixtures" / "mock_23andme.txt").read_bytes()
_ALIASES = [
    "direct",
    "relative",
    "symlink",
    "parent-symlink",
    "hardlink",
    "legacy",
    "legacy-link",
    "input-symlink",
    "legacy-hardlink",
]


def _paths(tmp_path, kind, destination):
    source = tmp_path / "genotypes.txt"
    source.write_bytes(_PAYLOAD)
    if kind == "direct":
        source.unlink()
        source = destination
        source.write_bytes(_PAYLOAD)
    elif kind == "relative":
        source.unlink()
        source = destination
        source.write_bytes(_PAYLOAD)
        (tmp_path / "spelling").mkdir()
        destination = tmp_path / "spelling" / ".." / destination.name
    elif kind == "symlink":
        destination.symlink_to(source.name)
    elif kind == "parent-symlink":
        source.unlink()
        source = destination
        source.write_bytes(_PAYLOAD)
        parent = tmp_path / "parent-link"
        parent.symlink_to(tmp_path, target_is_directory=True)
        destination = parent / destination.name
    elif kind == "hardlink":
        os.link(source, destination)
    elif kind == "legacy":
        source.unlink()
        source = destination.with_name(destination.name + ".tmp")
        source.write_bytes(_PAYLOAD)
    elif kind == "legacy-link":
        destination.with_name(destination.name + ".tmp").symlink_to(source.name)
    elif kind == "input-symlink":
        source.rename(destination)
        source.symlink_to(destination.name)
    elif kind == "legacy-hardlink":
        os.link(source, destination.with_name(destination.name + ".tmp"))
    return source, destination


def _forbid_cache_mutation(monkeypatch):
    calls = []

    def mutation(*args, **kwargs):
        calls.append("cache mutation")
        raise AssertionError("Input alias reached cache setup")

    monkeypatch.setattr(_helpers, "resolve_data_dir", mutation)
    monkeypatch.setattr(utility, "resolve_data_dir", mutation)
    return calls


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("report_format", ["json", "html"])
@pytest.mark.parametrize("kind", _ALIASES)
def test_report_alias_refused_before_cache_setup(
    tmp_path, monkeypatch, command, report_format, kind
):
    source, output = _paths(tmp_path, kind, tmp_path / f"report.{report_format}")
    before = source.read_bytes()
    calls = _forbid_cache_mutation(monkeypatch)
    result = CliRunner().invoke(
        main,
        [
            command,
            str(source),
            "--format",
            "23andme",
            "--output",
            str(output),
            "--report-format",
            report_format,
        ],
    )
    assert calls == [], result.output
    assert result.exit_code != 0
    assert "input" in result.output.lower() and "alias" in result.output.lower(), result.output
    assert source.read_bytes() == before


@pytest.mark.parametrize("extension", [".fam", ".bim", ".bed"])
@pytest.mark.parametrize("kind", _ALIASES)
def test_every_actual_plink_destination_checked_before_setup(
    tmp_path, monkeypatch, extension, kind
):
    source, output = _paths(tmp_path, kind, tmp_path / ("release" + extension))
    prefix = output.with_suffix("")  # release.fam/.bim/.bed come from the prefix "release".
    before = source.read_bytes()
    calls = _forbid_cache_mutation(monkeypatch)
    result = CliRunner().invoke(
        main,
        ["export", "plink", str(source), "--format", "23andme", "--output-prefix", str(prefix)],
    )
    assert calls == [], result.output
    assert result.exit_code != 0
    assert "input" in result.output.lower() and "alias" in result.output.lower(), result.output
    assert source.read_bytes() == before


def test_complete_report_and_vcf_set_checked_before_setup(tmp_path, monkeypatch):
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    report = tmp_path / "report.json"
    before = source.read_bytes()
    calls = _forbid_cache_mutation(monkeypatch)
    result = CliRunner().invoke(
        main, ["analyze", str(source), "--output", str(report), "--vcf-out", str(source)]
    )
    assert calls == [], result.output
    assert result.exit_code != 0
    assert not report.exists()
    assert source.read_bytes() == before


@pytest.mark.parametrize(
    "command", [["analyze"], ["methylation"], ["pharmacogenomics"], ["export", "plink"]]
)
def test_output_help_describes_input_alias_refusal(command):
    result = CliRunner().invoke(main, [*command, "--help"])
    assert result.exit_code == 0
    assert "input aliases" in " ".join(result.output.split())
    assert "atomic" in result.output

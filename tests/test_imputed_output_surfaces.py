# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Imputed input remains identified when a report or export leaves the run."""

from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import detect_parser

if TYPE_CHECKING:
    from pathlib import Path

    from click.testing import Result


def _input(tmp_path: Path, kind: str, imputed: bool) -> Path:
    if kind == "gfg-text":
        path = tmp_path / "calls.txt"
        qualifier = "imputed " if imputed else ""
        path.write_text(
            f"# Genes for Good v0.0: placeholder, {qualifier}genotype calls\n"
            "# Human genome reference build: GRCh37\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs900260281\t1\t12345\tAG\n"
        )
        return path
    path = tmp_path / "calls.vcf"
    extra = "##Genes_for_Good_version=v0.0 placeholder\n" if kind == "gfg-vcf" else ""
    if kind == "generic-vcf" and imputed:
        extra += '##INFO=<ID=DR2,Number=A,Type=Float,Description="Synthetic imputation quality">\n'
    # Phase alone is not an imputation marker on an ordinary VCF.
    gt = "0/1" if kind == "gfg-vcf" and not imputed else "0|1"
    path.write_text(
        "##fileformat=VCFv4.2\n##reference=GRCh37\n##contig=<ID=1,length=249250621>\n"
        + extra
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMOCK\n"
        + f"1\t12345\trs900260281\tA\tG\t.\tPASS\t.\tGT\t{gt}\n"
    )
    return path


def _method(kind: str, imputed: bool) -> str | None:
    return ("vcf-imputed" if kind == "generic-vcf" else "array-imputed") if imputed else None


def _analyze(path: Path, data: Path, *extra: str) -> Result:
    return CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--no-update",
            *extra,
        ],
    )


@pytest.mark.parametrize("kind", ["gfg-text", "gfg-vcf", "generic-vcf"])
@pytest.mark.parametrize("imputed", [False, True])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_report_imputation_surfaces(tmp_path, clinvar_data_dir, kind, imputed, surface):
    path = _input(tmp_path, kind, imputed)
    expected = _method(kind, imputed)
    assert detect_parser(path).get_metadata(path).get("method") == expected
    out = tmp_path / f"report.{surface}"
    args = () if surface == "terminal" else ("--output", str(out))
    result = _analyze(path, clinvar_data_dir, *args)
    assert result.exit_code == 0, result.output
    if surface == "json":
        assert json.loads(out.read_text())["input"].get("method") == expected
    elif surface == "html":
        page = out.read_text()
        banner = page.split('<details class="tech-details">')[0].lower()
        assert ("imput" in banner) is imputed
        if imputed:
            assert "confirm any finding by direct genotyping" in banner
    else:
        text = " ".join(result.output.lower().split())
        assert ("input method:" in text and "imput" in text) is imputed
        if imputed:
            assert "confirm any finding by direct genotyping" in text


@pytest.mark.parametrize("kind", ["gfg-vcf", "generic-vcf"])
@pytest.mark.parametrize("imputed", [False, True])
def test_annotated_vcf_carries_its_input_method(tmp_path, clinvar_data_dir, kind, imputed):
    import pysam

    path = _input(tmp_path, kind, imputed)
    out = tmp_path / "annotated.vcf"
    result = _analyze(path, clinvar_data_dir, "--vcf-out", str(out))
    assert result.exit_code == 0, result.output
    with pysam.VariantFile(str(out)) as vcf:
        marks = [
            dict(row.items()) for row in vcf.header.records if row.key == "ALLELIX_InputMethod"
        ]
        calls = list(vcf)
    assert len(calls) == 1 and calls[0].samples["MOCK"]["GT"] == (0, 1)
    if imputed:
        assert len(marks) == 1, out.read_text()
        assert marks[0]["ID"] == _method(kind, imputed)
        assert "imput" in marks[0]["Description"].lower()
        schema = re.search(
            r"^##ALLELIX_VCF_SchemaVersion=(\d+)\.(\d+)\.(\d+)$", out.read_text(), re.M
        )
        assert schema and tuple(map(int, schema.groups())) > (0, 6, 0)
    else:
        assert not marks, "A measured input must not acquire an imputed-method mark"


@pytest.mark.parametrize("kind", ["gfg-text", "gfg-vcf", "generic-vcf"])
@pytest.mark.parametrize("imputed", [False, True])
def test_plink_export_discloses_imputation(tmp_path, clinvar_data_dir, kind, imputed):
    path = _input(tmp_path, kind, imputed)
    prefix = tmp_path / "export"
    result = CliRunner().invoke(
        main,
        [
            "export",
            "plink",
            str(path),
            "--build",
            "grch37",
            "--data-dir",
            str(clinvar_data_dir),
            "--output-prefix",
            str(prefix),
        ],
    )
    assert result.exit_code == 0, result.output
    assert all(prefix.with_suffix(suffix).is_file() for suffix in (".bed", ".bim", ".fam"))
    text = " ".join(result.output.lower().split())
    assert bool(re.search(r"\bimputed\b", text)) is imputed, result.output

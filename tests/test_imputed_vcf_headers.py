# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A VCF whose header declares an imputation tool's INFO field is marked imputed.

Minimac declares the INFO fields ``IMPUTED``, ``TYPED`` and ``R2``; Beagle
``DR2`` and ``IMP``; Minimac, Beagle and GLIMPSE the FORMAT dosage ``DS``, and
GLIMPSE the FORMAT ``HS``. Only declarations count, never values. FORMAT
``GP`` alone and GATK's INFO ``DS`` (downsampled) do not. Measured VCFs
(benchmark call sets, GATK and DeepVariant gVCFs, a population call set)
declare none of them and stay unmarked.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from test_sources import data  # noqa: F401  (the data-folder fixture)

from allelix.cli import main
from allelix.parsers.vcf import IMPUTATION_INFO_IDS, METHOD_VCF_IMPUTED, VcfParser
from allelix.reports import input_method_label, is_imputed

if TYPE_CHECKING:
    from pathlib import Path

# The tools' declarations as their VCFs write them.
_MINIMAC = {
    "IMPUTED": (
        '##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">'
    ),
    "TYPED": (
        '##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">'
    ),
    "R2": (
        '##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">'
    ),
}
_BEAGLE = {
    "DR2": (
        '##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated '
        'squared correlation between estimated REF dose and true REF dose">'
    ),
    "IMP": '##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">',
}


def _info(ids: list[str]) -> list[str]:
    return [f'##INFO=<ID={i},Number=.,Type=String,Description="placeholder">' for i in ids]


def _format(ids: list[str]) -> list[str]:
    return [f'##FORMAT=<ID={i},Number=.,Type=String,Description="placeholder">' for i in ids]


# Measured VCFs: the INFO and FORMAT IDs their headers declare (the IDs only;
# descriptions are placeholders).
_MEASURED = {
    "benchmark-call-set": _info(
        [
            "DPSum",
            "platforms",
            "platformnames",
            "platformbias",
            "datasets",
            "datasetnames",
            "datasetsmissingcall",
            "callsets",
            "callsetnames",
            "varType",
            "filt",
            "callable",
            "difficultregion",
            "arbitrated",
            "callsetwiththisuniqgenopassing",
            "callsetwithotheruniqgenopassing",
        ]
    )
    + _format(["DP", "GQ", "ADALL", "AD", "PS"]),
    "gatk-gvcf": _info(
        [
            "BaseQRankSum",
            "DP",
            "DS",
            "END",
            "ExcessHet",
            "InbreedingCoeff",
            "MLEAC",
            "MLEAF",
            "MQRankSum",
            "RAW_MQandDP",
            "ReadPosRankSum",
        ]
    )
    + _format(["AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB"]),
    "deepvariant-gvcf": _info(["END"])
    + _format(["GQ", "DP", "MIN_DP", "AD", "VAF", "PL", "PS", "MED_DP"]),
    "population-call-set": _info(
        [
            "AC",
            "AF",
            "AN",
            "DP",
            "DS",
            "END",
            "FS",
            "HWE",
            "MQ",
            "QD",
            "SOR",
            "VQSLOD",
            "VariantType",
            "culprit",
            "ExcHet",
        ]
    ),
}


def _vcf(path: Path, header: list[str], gts: tuple[str, ...] = ("0/1", "1/1")) -> Path:
    body = "".join(
        f"1\t{3999300 + i}\trs{3999990300 + i}\tA\tG\t.\tPASS\t.\tGT\t{gt}\n"
        for i, gt in enumerate(gts)
    )
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + "".join(line + "\n" for line in header)
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + body
    )
    return path


def test_the_five_fields() -> None:
    assert frozenset({"IMPUTED", "TYPED", "R2", "DR2", "IMP"}) == IMPUTATION_INFO_IDS
    assert is_imputed(METHOD_VCF_IMPUTED)
    assert "imputed" in input_method_label(METHOD_VCF_IMPUTED)


@pytest.mark.parametrize(("field", "line"), [*_MINIMAC.items(), *_BEAGLE.items()])
def test_each_declared_field_marks_the_file(tmp_path: Path, field: str, line: str) -> None:
    path = _vcf(tmp_path / f"{field}.vcf", [line])
    assert VcfParser().get_metadata(path).get("method") == METHOD_VCF_IMPUTED


_GLIMPSE = [
    '##INFO=<ID=RAF,Number=A,Type=Float,Description="ALT frequency in the panel">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT frequency from DS/GP">',
    '##INFO=<ID=INFO,Number=A,Type=Float,Description="IMPUTE info score">',
    '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Genotype dosage">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posteriors">',
    '##FORMAT=<ID=HS,Number=1,Type=Integer,Description="Sampled haplotype pairs">',
]


@pytest.mark.parametrize(
    "header",
    [_GLIMPSE, _format(["DS"]), _format(["HS"]), _format(["GP", "DS"])],
    ids=["glimpse", "format-ds", "format-hs", "format-gp-ds"],
)
def test_format_dosage_or_hs_marks_the_file(tmp_path: Path, header: list[str]) -> None:
    path = _vcf(tmp_path / "x.vcf", header, gts=("0|1", "1|1"))
    assert VcfParser().get_metadata(path).get("method") == METHOD_VCF_IMPUTED


@pytest.mark.parametrize("tool", ["minimac", "beagle"])
def test_a_tools_full_header_marks_the_file(tmp_path: Path, tool: str) -> None:
    lines = list((_MINIMAC if tool == "minimac" else _BEAGLE).values())
    path = _vcf(tmp_path / "x.vcf", lines, gts=("0|1", "1|1"))
    assert VcfParser().get_metadata(path).get("method") == METHOD_VCF_IMPUTED


@pytest.mark.parametrize("shape", list(_MEASURED))
def test_measured_vcfs_stay_unmarked(tmp_path: Path, shape: str) -> None:
    path = _vcf(tmp_path / "x.vcf", _MEASURED[shape])
    assert "method" not in VcfParser().get_metadata(path)


@pytest.mark.parametrize(
    "header",
    [
        _format(["GP"]),  # GP alone: a variant caller can add it
        _info(["DS", "HS"]),  # GATK's INFO DS is "downsampled": INFO DS, HS do not count
        _info(["r2", "imp", "Imputed"]),  # case matters: the tools write upper case
        _info(["MAF_R2", "IMPUTED_FRACTION", "DR2X"]),  # a longer ID is another field
        ['##FILTER=<ID=IMP,Description="placeholder">', '##ALT=<ID=R2,Description="placeholder">'],
        ['##source="imputed with a placeholder tool"', "##contig=<ID=IMP>"],
    ],
    ids=[
        "format-gp-only",
        "info-ds-hs",
        "other-case",
        "longer-ids",
        "other-line-kinds",
        "words-not-ids",
    ],
)
def test_only_an_info_declaration_counts(tmp_path: Path, header: list[str]) -> None:
    path = _vcf(tmp_path / "x.vcf", header)
    assert "method" not in VcfParser().get_metadata(path)


def test_values_are_not_read(tmp_path: Path) -> None:
    """An R2 value in the records without the declaration marks nothing."""
    path = tmp_path / "x.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t3999301\trs3999990301\tA\tG\t.\tPASS\tR2=0.95;IMPUTED\tGT\t0/1\n"
    )
    assert "method" not in VcfParser().get_metadata(path)


def test_genes_for_good_keeps_its_own_method(tmp_path: Path) -> None:
    path = _vcf(
        tmp_path / "x.vcf",
        ["##Genes_for_Good_version=v0.0 placeholder", _MINIMAC["R2"]],
        gts=("0|1", "1|1"),
    )
    assert VcfParser().get_metadata(path).get("method") == "array-imputed"


def test_report_says_imputed(tmp_path: Path, data: Path) -> None:  # noqa: F811
    path = _vcf(tmp_path / "x.vcf", [_BEAGLE["DR2"]])
    out = tmp_path / "r.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(out.read_text())["input"]["method"] == METHOD_VCF_IMPUTED
    flat = " ".join(result.output.split())
    label = "imputed: the VCF header declares imputation fields (Minimac, Beagle or GLIMPSE)"
    assert label in flat
    assert "confirm any finding by direct genotyping" in flat


def test_html_report_reads_as_any_imputed_file(tmp_path: Path, data: Path) -> None:  # noqa: F811
    """The caution in view, and "Imputed genotypes" in place of the pseudogene paragraph."""
    path = _vcf(tmp_path / "x.vcf", [_MINIMAC["R2"]])
    out = tmp_path / "r.html"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    page = out.read_text()
    start = page.index('<details class="tech-details">')
    top = page[:start]
    assert "<strong>Input method.</strong> Imputed: the VCF header declares imputation" in top
    assert "confirm any finding by direct genotyping" in top
    assert "<h3>Imputed genotypes</h3>" in page
    assert "Most genotypes in this file were filled in (imputed)" in page
    assert "<h3>Pseudogene cross-hybridization</h3>" not in page

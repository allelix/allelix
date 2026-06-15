# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Focused report subcommands: `methylation` and `pharmacogenomics`."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.cli import _helpers, main
from allelix.cli._options import (
    _BUILD_OPT,
    _DATA_DIR_OPT,
    _DIFF_OPT,
    _EXCLUDE_SNPEDIA_OPT,
    _FILE_ARG,
    _FORMAT_OPT,
    _GWAS_ALL_OPT,
    _GWAS_MIN_MAG_OPT,
    _INCLUDE_BENIGN_OPT,
    _INCLUDE_GWAS_OPT,
    _MIN_MAG_OPT,
    _NO_ALPHAMISSENSE_OPT,
    _NO_CADD_OPT,
    _NO_GNOMAD_OPT,
    _NO_UPDATE_OPT,
    _OUTPUT_OPT,
    _REPORT_FORMAT_OPT,
    _SAMPLE_OPT,
    _SNPEDIA_MIN_MAG_OPT,
)
from allelix.reports.methylation import METHYLATION_PANEL_GENES

if TYPE_CHECKING:
    from pathlib import Path


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_INCLUDE_GWAS_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
@_NO_CADD_OPT
@_SAMPLE_OPT
def methylation(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    include_gwas: bool,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
    no_cadd: bool,
    sample: str | None,
) -> None:
    """Methylation-pathway-focused report (MTHFR, MTR, MTRR, COMT, CBS, …)."""
    excluded: set[str] = set()
    if not include_gwas:
        excluded.add("gwas")
    if exclude_snpedia:
        excluded.add("snpedia")
    _helpers._run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category=None,
        genes=METHYLATION_PANEL_GENES,
        build=_helpers._normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset(excluded) if excluded else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
        no_cadd=no_cadd,
        sample=sample,
    )


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_INCLUDE_GWAS_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
@_NO_CADD_OPT
@_SAMPLE_OPT
def pharmacogenomics(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    include_gwas: bool,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
    no_cadd: bool,
    sample: str | None,
) -> None:
    """Pharmacogenomics-focused report (annotations from ClinPGx-style sources)."""
    excluded: set[str] = set()
    if not include_gwas:
        excluded.add("gwas")
    if exclude_snpedia:
        excluded.add("snpedia")
    _helpers._run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category="pharma",
        genes=None,
        build=_helpers._normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset(excluded) if excluded else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
        no_cadd=no_cadd,
        sample=sample,
    )

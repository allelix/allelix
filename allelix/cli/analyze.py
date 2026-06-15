# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The `analyze` subcommand — full annotation against all ready databases."""

from __future__ import annotations

from typing import TYPE_CHECKING

import click

from allelix.cli import _helpers, main
from allelix.cli._options import (
    _BUILD_OPT,
    _DATA_DIR_OPT,
    _DIFF_OPT,
    _EXCLUDE_SNPEDIA_OPT,
    _FILE_ARG,
    _FILTER_FILE_OPT,
    _FORMAT_OPT,
    _GWAS_ALL_OPT,
    _GWAS_MIN_MAG_OPT,
    _INCLUDE_BENIGN_OPT,
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

if TYPE_CHECKING:
    from pathlib import Path


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_MIN_MAG_OPT
@click.option(
    "--category",
    type=str,
    default=None,
    help="Filter to a single bucket (clinical, pharma).",
)
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_GWAS_MIN_MAG_OPT
@_SNPEDIA_MIN_MAG_OPT
@_GWAS_ALL_OPT
@_EXCLUDE_SNPEDIA_OPT
@_DIFF_OPT
@_FILTER_FILE_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
@_NO_CADD_OPT
@_SAMPLE_OPT
def analyze(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    min_magnitude: float,
    category: str | None,
    output: Path | None,
    report_format: str | None,
    build: str,
    include_benign: bool,
    gwas_min_magnitude: float,
    snpedia_min_magnitude: float,
    gwas_all: bool,
    exclude_snpedia: bool,
    diff_path: Path | None,
    filter_file: Path | None,
    no_update: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
    no_cadd: bool,
    sample: str | None,
) -> None:
    """Annotate a genotype file against all ready reference databases."""
    filter_genes: frozenset[str] | None = None
    filter_rsids: frozenset[str] | None = None
    if filter_file is not None:
        filter_genes, filter_rsids = _helpers._parse_filter_file(filter_file)
        # Empty sets (file had only comments/blanks) still apply — they
        # mean "match nothing", producing an empty report.

    _helpers._run_analysis_command(
        file_path=file_path,
        fmt=fmt,
        data_dir=data_dir,
        output=output,
        report_format=report_format,
        min_magnitude=min_magnitude,
        category=category,
        genes=filter_genes,
        rsids=filter_rsids,
        build=_helpers._normalize_cli_build(build),
        include_benign=include_benign,
        gwas_min_magnitude=gwas_min_magnitude,
        snpedia_min_magnitude=snpedia_min_magnitude,
        exclude_sources=frozenset({"snpedia"}) if exclude_snpedia else None,
        gwas_all=gwas_all,
        diff_path=diff_path,
        no_update=no_update,
        no_gnomad=no_gnomad,
        no_alphamissense=no_alphamissense,
        no_cadd=no_cadd,
        sample=sample,
    )

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared click option decorator constants.

Each option is defined once here as a decorator that subcommand modules
apply. Centralizing the option definitions keeps help text consistent
across commands and makes it cheap to add a new option to every command
that takes a particular concept (e.g., ``--no-cadd``).
"""

from __future__ import annotations

from pathlib import Path

import click

_FILE_ARG = click.argument(
    "file_path", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
_FORMAT_OPT = click.option(
    "--format", "fmt", default=None, help="Force a specific parser. Default: auto-detect."
)
_DATA_DIR_OPT = click.option(
    "--data-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override database cache location.",
)
_MIN_MAG_OPT = click.option(
    "--min-magnitude",
    type=float,
    default=5.0,
    show_default=True,
    help="Filter annotations below this magnitude. Use 0 for the full unfiltered set.",
)
# The focused subcommands (`allelix methylation`, `allelix pharmacogenomics`)
# default to a lower floor than the general `analyze` report. The intent of
# a focused report is to surface focused-level signal — ClinPGx LoE 3 hits
# on COMT / MTHFR / MTR / CBS for `methylation`, similar for the pharma
# subcommand — which the magnitude-5 floor that protects `analyze` from
# noise would filter out as if they were noise. Floor 3.0 keeps LoE 3
# in, keeps ClinVar Benign (mag 1.0) out.
_FOCUSED_MIN_MAG_OPT = click.option(
    "--min-magnitude",
    type=float,
    default=3.0,
    show_default=True,
    help=(
        "Filter annotations below this magnitude. Focused subcommands default "
        "to 3.0 (vs. analyze's 5.0) so ClinPGx LoE 3 hits on panel genes — "
        "the intended signal of a focused report — aren't filtered out. "
        "Use 0 for the full unfiltered set."
    ),
)
_OUTPUT_OPT = click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write a report file (.html or .json). Omit for terminal output.",
)
_REPORT_FORMAT_OPT = click.option(
    "--report-format",
    type=click.Choice(["html", "json"], case_sensitive=False),
    default=None,
    help="Override report format detection (otherwise inferred from --output extension).",
)
_INCLUDE_BENIGN_OPT = click.option(
    "--include-benign",
    is_flag=True,
    default=False,
    help="Include ClinVar Benign/Likely_benign annotations (suppressed by default).",
)
_GWAS_MIN_MAG_OPT = click.option(
    "--gwas-min-magnitude",
    type=float,
    default=9.0,
    show_default=True,
    help="Magnitude floor for GWAS Catalog annotations (overrides --min-magnitude for GWAS).",
)
_SNPEDIA_MIN_MAG_OPT = click.option(
    "--snpedia-min-magnitude",
    type=float,
    default=2.0,
    show_default=True,
    help="Magnitude floor for SNPedia annotations (overrides --min-magnitude for SNPedia).",
)
_INCLUDE_GWAS_OPT = click.option(
    "--include-gwas",
    is_flag=True,
    default=False,
    help="Include GWAS Catalog annotations (excluded by default in focused reports).",
)
_EXCLUDE_SNPEDIA_OPT = click.option(
    "--exclude-snpedia",
    is_flag=True,
    default=False,
    help="Exclude SNPedia annotations. Required for commercial use (CC BY-NC-SA 3.0).",
)
_GWAS_ALL_OPT = click.option(
    "--gwas-all",
    is_flag=True,
    default=False,
    help="Include all GWAS trait categories (disables default noise filtering).",
)
_DIFF_OPT = click.option(
    "--diff",
    "diff_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Dev/QA tool: compare current output against a previous JSON report "
        "to detect regressions from code changes, database refreshes, or "
        "filter adjustments. Shows new, changed, and removed annotations. "
        "Not a monitoring tool — use for version-to-version validation."
    ),
)
_FILTER_FILE_OPT = click.option(
    "--filter-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Plain text file with rsIDs and/or gene names (one per line) to "
        "filter the report. Lines matching '^rs\\d+$' are rsIDs; everything "
        "else is a gene name. Comments (#) and blank lines are ignored."
    ),
)
_NO_UPDATE_OPT = click.option(
    "--no-update",
    is_flag=True,
    default=False,
    help="Skip the pre-analysis database freshness check.",
)
_NO_GNOMAD_OPT = click.option(
    "--no-gnomad",
    is_flag=True,
    default=False,
    help="Skip gnomAD population frequency enrichment.",
)
_NO_ALPHAMISSENSE_OPT = click.option(
    "--no-alphamissense",
    is_flag=True,
    default=False,
    help="Skip AlphaMissense variant pathogenicity enrichment.",
)
_NO_CADD_OPT = click.option(
    "--no-cadd",
    is_flag=True,
    default=False,
    help="Skip CADD deleteriousness score enrichment.",
)
_BUILD_OPT = click.option(
    "--build",
    type=click.Choice(["grch37", "grch38", "auto"], case_sensitive=False),
    default="auto",
    help=(
        "Genome build of the input file. 'auto' detects from position data "
        "and ignores the file header. 'grch37' / 'grch38' force a "
        "specific build, skipping detection."
    ),
)
_SAMPLE_OPT = click.option(
    "--sample",
    type=str,
    default=None,
    help=(
        "VCF only: which sample column to read from a multi-sample VCF. "
        "Required when the VCF has more than one sample; ignored on "
        "single-sample VCFs and on array formats. Use `allelix stats <file>` "
        "or any plain VCF tool to list available samples."
    ),
)

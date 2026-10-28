# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared click option decorator constants.

Each option is defined once here as a decorator that subcommand modules
apply. Centralizing the option definitions keeps help text consistent
across commands and makes it cheap to add a new option to every command
that takes a particular concept (e.g., ``--no-cadd``).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import click


class _FiniteFloat(click.types.FloatParamType):
    """A float option that refuses nan, inf and -inf as a usage error (2.5.2)."""

    def convert(
        self,
        value: Any,  # noqa: ANN401 — click passes a string or an already-parsed number
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> float:
        """Parse like ``float`` and fail on a value that is not a finite number."""
        number: float = super().convert(value, param, ctx)
        if not math.isfinite(number):
            self.fail(f"{value!r} is not a finite number.", param, ctx)
        return number


_FINITE_FLOAT = _FiniteFloat()

_FILE_ARG = click.argument(
    "file_path", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
_FORMAT_OPT = click.option(
    "--format",
    "fmt",
    default=None,
    help=(
        "Read the file with this reader even if it does not recognise the file, at your "
        "own risk (when it does not, every report carries a warning). Default: auto-detect."
    ),
)
_DATA_DIR_OPT = click.option(
    "--data-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Override database cache location.",
)
_MIN_MAG_OPT = click.option(
    "--min-magnitude",
    type=_FINITE_FLOAT,
    default=5.0,
    show_default=True,
    help=(
        "Filter annotations below this magnitude. 0 removes this floor; the per-source "
        "floors (--gwas-min-magnitude, --snpedia-min-magnitude) and the other filters "
        "still apply."
    ),
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
    type=_FINITE_FLOAT,
    default=3.0,
    show_default=True,
    help=(
        "Filter annotations below this magnitude. Focused subcommands default "
        "to 3.0 (vs. analyze's 5.0) so ClinPGx LoE 3 hits on panel genes — "
        "the intended signal of a focused report — aren't filtered out. "
        "0 removes this floor; the per-source floors and the other filters still apply."
    ),
)
_OUTPUT_OPT = click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Write a report file (.html or .json). Refuses input aliases (and, for "
        "analyze, the --vcf-out file) and publishes atomically; a new file gets "
        "the umask's mode. Omit for terminal output."
    ),
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
    type=_FINITE_FLOAT,
    default=9.0,
    show_default=True,
    help=(
        "Magnitude floor for GWAS Catalog annotations; GWAS rows must meet the"
        " larger of this and --min-magnitude."
    ),
)
_SNPEDIA_MIN_MAG_OPT = click.option(
    "--snpedia-min-magnitude",
    type=_FINITE_FLOAT,
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
        "else is a gene name. A # starts a comment, on its own line or after "
        "an entry; the rest of the line and blank lines are ignored."
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
    type=click.Choice(["grch37", "grch38", "grch36", "auto"], case_sensitive=False),
    default="auto",
    help=(
        "Genome build of the input file. 'auto' decides from the evidence "
        "the file carries (VCF contig lengths, known-SNP positions) and stops "
        "when that evidence contradicts the file's label or is missing. "
        "'grch37' / 'grch38' force a specific build; the run still stops when "
        "the file's contig lengths or confident known-SNP positions say "
        "another build. Where the file's known-SNP positions match no build (0 of "
        "N), a forced build still runs, with a Build warning on the console and "
        "in every report. 'grch36' forces an old GRCh36 / hg18 file, with "
        "position-keyed lookups off and never checked (see docs/grch36-liftover.md)."
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
# VCF quality gate. Rows the variant caller filtered are skipped
# unless --keep-filtered; calls below the GQ floor become no-calls.
_KEEP_FILTERED_OPT = click.option(
    "--keep-filtered",
    is_flag=True,
    default=False,
    help=(
        "VCF only: keep rows whose FILTER is not PASS (for example LowQual or a "
        "VQSR tranche). By default the variant caller's filters are honoured "
        "and those rows are skipped and counted."
    ),
)
_MIN_GQ_OPT = click.option(
    "--min-gq",
    type=click.IntRange(min=0),
    default=20,
    show_default=True,
    help=(
        "VCF only: genotype quality floor. A call or gVCF reference block "
        "with GQ below this value, or with no reads (DP 0), is treated as a "
        "no-call. 0 turns the GQ floor off."
    ),
)
# annotated VCF output. Requires VCF/gVCF input (any other format
# errors out early). Existing output paths are refused (no clobber, no
# --force in v2.3.0). See ADR-0036 for the output contract.
_VCF_OUT_OPT = click.option(
    "--vcf-out",
    "vcf_out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Write an annotated VCF at this path (v2.3.0+; experimental output "
        "contract, ADR-0036). Requires VCF or gVCF input. Refuses to "
        "overwrite an existing file or input aliases; publishes atomically."
    ),
)

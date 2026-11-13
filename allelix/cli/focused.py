# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Focused report subcommands: `methylation` and `pharmacogenomics`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import click
from click.core import ParameterSource

from allelix.cli import _helpers, main
from allelix.cli._options import (
    _BUILD_OPT,
    _DATA_DIR_OPT,
    _DIFF_OPT,
    _EXCLUDE_SNPEDIA_OPT,
    _FILE_ARG,
    _FOCUSED_MIN_MAG_OPT,
    _FORMAT_OPT,
    _GWAS_ALL_OPT,
    _GWAS_MIN_MAG_OPT,
    _INCLUDE_BENIGN_OPT,
    _INCLUDE_GWAS_OPT,
    _NO_ALPHAMISSENSE_OPT,
    _NO_CADD_OPT,
    _NO_GNOMAD_OPT,
    _NO_UPDATE_OPT,
    _OUTPUT_OPT,
    _REPORT_FORMAT_OPT,
    _SAMPLE_OPT,
    _SHOW_UNCHECKED_OPT,
    _SNPEDIA_MIN_MAG_OPT,
)
from allelix.reports.methylation import METHYLATION_PANEL_GENES

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path


def _refuse_explicit(ctx: click.Context, param: click.Parameter, value: object) -> object:
    """Stop with a usage error when an option that never applies was given (2.7.0)."""
    if param.name is not None and ctx.get_parameter_source(param.name) not in (
        ParameterSource.DEFAULT,
        None,
    ):
        option = param.opts[0]
        raise click.UsageError(
            f"{option} never applies to pharmacogenomics: its report shows only ClinPGx "
            "findings, so no GWAS Catalog or SNPedia finding can reach it. Leave the "
            "option out, or use analyze for those sources.",
            ctx,
        )
    return value


# The options pharmacogenomics refuses (``_not_for_pharma``).
_REFUSED_ON_PHARMA: set[str] = set()


class _PharmacogenomicsCommand(click.Command):
    """A command whose typo hints leave out the options it refuses (2.7.0).

    Click suggests the closest of all a command's long options for a
    mistyped one, hidden ones included, so ``--include-gw`` would suggest
    ``--include-gwas``, which pharmacogenomics then refuses.
    """

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        """Parse as Click does; a mistyped option's hints name only options it accepts."""
        try:
            return super().parse_args(ctx, args)
        except click.NoSuchOption as exc:
            # get_params, not self.params: it includes the generated --help.
            offered = [
                opt
                for param in self.get_params(ctx)
                for opt in (*param.opts, *param.secondary_opts)
                if opt.startswith("--") and opt not in _REFUSED_ON_PHARMA
            ]
            raise click.NoSuchOption(
                exc.option_name, message=exc.message, possibilities=offered, ctx=exc.ctx
            ) from None


def _not_for_pharma(name: str, default: object) -> Callable:
    """``name`` on pharmacogenomics: hidden, at its default, refused when given."""
    _REFUSED_ON_PHARMA.add(name)
    if isinstance(default, bool):
        return click.option(
            name, is_flag=True, default=default, hidden=True, callback=_refuse_explicit
        )
    return click.option(name, type=float, default=default, hidden=True, callback=_refuse_explicit)


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_FOCUSED_MIN_MAG_OPT
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
@_SHOW_UNCHECKED_OPT
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
    show_unchecked: bool,
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
        show_unchecked=show_unchecked,
    )


@main.command(cls=_PharmacogenomicsCommand)
@_FILE_ARG
@_FORMAT_OPT
@_DATA_DIR_OPT
@_FOCUSED_MIN_MAG_OPT
@_OUTPUT_OPT
@_REPORT_FORMAT_OPT
@_BUILD_OPT
@_INCLUDE_BENIGN_OPT
@_not_for_pharma("--gwas-min-magnitude", 9.0)
@_not_for_pharma("--snpedia-min-magnitude", 2.0)
@_not_for_pharma("--include-gwas", False)
@_not_for_pharma("--gwas-all", False)
@_not_for_pharma("--exclude-snpedia", False)
@_DIFF_OPT
@_NO_UPDATE_OPT
@_NO_GNOMAD_OPT
@_NO_ALPHAMISSENSE_OPT
@_NO_CADD_OPT
@_SAMPLE_OPT
@_SHOW_UNCHECKED_OPT
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
    show_unchecked: bool,
) -> None:
    """Pharmacogenomics-focused report (annotations from ClinPGx-style sources).

    The GWAS Catalog and SNPedia options of analyze and methylation are
    refused here: the report keeps ClinPGx findings only, so they never
    change it.
    """
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
        show_unchecked=show_unchecked,
    )

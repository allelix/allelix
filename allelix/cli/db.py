# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The `db` subcommand group — local reference database cache management."""

from __future__ import annotations

from typing import TYPE_CHECKING

import click
from rich.table import Table

from allelix.annotators import get_annotators
from allelix.cli import _helpers, main
from allelix.cli._helpers import console
from allelix.cli._options import _DATA_DIR_OPT
from allelix.databases import resolve_data_dir

if TYPE_CHECKING:
    from pathlib import Path


@main.group()
def db() -> None:
    """Manage local reference database cache."""


def _confirm_cadd_license(*, license_held: bool = False) -> bool:
    """Show the CADD license notice and ask for confirmation."""
    if license_held:
        console.print(
            "\n[bold yellow]CADD License Notice[/bold yellow]\n"
            "Commercial license asserted. Proceeding with CADD download.\n"
        )
        return True
    console.print(
        "\n[bold yellow]CADD License Notice[/bold yellow]\n"
        "CADD scores are provided by the University of Washington.\n"
        "Commercial use requires a license from UW CoMotion\n"
        "([link=https://els2.comotion.uw.edu/product/cadd-scores]"
        "https://els2.comotion.uw.edu/product/cadd-scores[/link]).\n"
        "By continuing, you confirm that your use is non-commercial\n"
        "or that you hold a valid commercial license.\n"
    )
    return click.confirm("Continue with CADD download?", default=False)


@db.command("update")
@_DATA_DIR_OPT
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-download even if the local cache appears current.",
)
@click.option(
    "--no-gnomad",
    is_flag=True,
    default=False,
    help="Skip gnomAD population frequency database.",
)
@click.option(
    "--no-alphamissense",
    is_flag=True,
    default=False,
    help="Skip AlphaMissense pathogenicity database.",
)
@click.option(
    "--cadd",
    "include_cadd",
    is_flag=True,
    default=False,
    help="Download CADD deleteriousness scores (non-commercial use only; disabled by default).",
)
@click.option(
    "--build",
    type=click.Choice(["grch37", "grch38", "both"], case_sensitive=False),
    default="both",
    help=(
        "Which ClinVar genome build(s) to download. 'both' (default) keeps "
        "GRCh37 and GRCh38 caches in sync so `analyze` can dispatch by "
        "detected build (ADR-0021). 'grch37' / 'grch38' restrict to one to "
        "save bandwidth."
    ),
)
def db_update(
    data_dir: Path | None,
    force: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
    include_cadd: bool,
    build: str,
) -> None:
    """Download or refresh reference databases.

    For each annotator:
      - no cache → download
      - --force → download
      - cache + remote signal matches cache → skip
      - cache + remote signal differs (or legacy v0.4.1 cache with no
        stored signal) → download
      - cache + remote signal can't be fetched → skip with notice (use
        --force to override)

    `--build` selects which ClinVar build(s) to manage. Default 'both'
    downloads GRCh37 and GRCh38 caches.
    """
    resolved = resolve_data_dir(data_dir)
    console.print(f"Data directory: [cyan]{resolved}[/cyan]")

    from allelix.config import load_config

    cfg = load_config(resolved)

    clinvar_builds = _helpers._resolve_clinvar_builds(build)
    for annotator in get_annotators(
        resolved, clinvar_builds=clinvar_builds, cadd_full=cfg.cadd_full
    ):
        with annotator:
            if no_gnomad and annotator.name == "gnomad":
                console.print(f"  [dim]{annotator.name}: skipped (--no-gnomad)[/dim]")
                continue
            if no_alphamissense and annotator.name == "alphamissense":
                console.print(f"  [dim]{annotator.name}: skipped (--no-alphamissense)[/dim]")
                continue

            if annotator.name == "cadd":
                if not include_cadd and not cfg.is_enabled("cadd"):
                    console.print(
                        f"  [dim]{annotator.name}: disabled "
                        "(enable with `allelix config set sources.cadd true` "
                        "or pass `--cadd`)[/dim]"
                    )
                    continue
                if (not annotator.is_ready() or force) and not _confirm_cadd_license(
                    license_held=cfg.license_held("cadd"),
                ):
                    console.print(f"  [dim]{annotator.name}: skipped (declined)[/dim]")
                    continue

            if not annotator.requires_download:
                if annotator.is_ready():
                    console.print(
                        f"  [dim]{annotator.name}: ready "
                        f"({annotator.version() or 'unknown'})[/dim]"
                    )
                continue

            if not annotator.is_ready():
                console.print(f"  [bold]{annotator.name}[/bold]: downloading…")
                if _helpers._run_setup(annotator):
                    console.print(
                        f"  [green]✓ {annotator.name} ready[/green] "
                        f"(version {annotator.version() or '(unknown)'})"
                    )
                continue

            if force:
                console.print(f"  [bold]{annotator.name}[/bold]: --force; refreshing…")
                if _helpers._run_setup(annotator):
                    console.print(
                        f"  [green]✓ {annotator.name} refreshed[/green] "
                        f"(version {annotator.version() or '(unknown)'})"
                    )
                continue

            # Code-driven sources (commit-pinned HF caches) are updated
            # only via code changes — no runtime freshness probe needed.
            if not annotator.server_driven_freshness:
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {annotator.version() or '(unknown)'})[/dim]"
                )
                continue

            remote = annotator.fetch_remote_signal()
            if remote is None:
                console.print(
                    f"  [yellow]{annotator.name}: cache present, but remote "
                    "freshness can't be verified (network error or no signal). "
                    "Pass --force to refresh anyway.[/yellow]"
                )
                continue

            cached = annotator.cached_remote_signal()
            if cached == remote:
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {annotator.version() or '(unknown)'})[/dim]"
                )
                continue

            if cached is None:
                # GH #20: a cache with no stored freshness signal almost
                # always predates the signal mechanism — i.e., it is old.
                # The previous behavior was to stamp the live remote signal
                # onto the cache and call it current, which permanently
                # marked stale data as fresh (only `--force` would escape).
                # Treat tagless caches as needing a refresh.
                console.print(
                    f"  [bold]{annotator.name}[/bold]: cache predates the "
                    "freshness signal; re-downloading…"
                )
            else:
                console.print(
                    f"  [bold]{annotator.name}[/bold]: remote signal changed; refreshing…"
                )
            if _helpers._run_setup(annotator):
                console.print(
                    f"  [green]✓ {annotator.name} refreshed[/green] "
                    f"(version {annotator.version() or '(unknown)'})"
                )


@db.command("status")
@_DATA_DIR_OPT
def db_status(data_dir: Path | None) -> None:
    """Show installed reference database versions and freshness."""
    from allelix.config import load_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)
    table = Table(title=f"Reference Databases ({resolved})")
    table.add_column("Annotator", style="cyan", no_wrap=True)
    table.add_column("Ready", justify="center")
    table.add_column("Version")
    table.add_column("Records", justify="right")
    for annotator in get_annotators(resolved, cadd_full=cfg.cadd_full):
        with annotator:
            ready = annotator.is_ready()
            ready_marker = "[green]yes[/green]" if ready else "[red]no[/red]"
            version = annotator.version() or "—"
            sig = getattr(annotator, "cached_remote_signal", lambda: None)()
            if sig and "cpic:unavailable" in sig:
                version += " (no CPIC)"
            records = "—"
            count_fn = getattr(annotator, "record_count", None)
            if callable(count_fn):
                count = count_fn()
                if count is not None:
                    records = f"{count:,}"
            table.add_row(annotator.display_name, ready_marker, version, records)
    console.print(table)

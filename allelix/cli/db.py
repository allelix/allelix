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


# Cache items considered part of the CADD opt-in footprint. Matched by
# prefix against each entry in the data dir so a future renaming of the
# sqlite file (e.g. cadd.GRCh38.sqlite) is still caught.
_CADD_PREFIXES: tuple[str, ...] = ("cadd",)


# Names that identify a directory as an allelix data dir for `db clean`'s
# pre-deletion guard. Matched by prefix so future variants (e.g. multi-
# build forms like clinvar.GRCh37.sqlite, sidecar files like
# snpedia.sqlite-bak, the PharmGKB raw zip `clinicalAnnotations.zip`) are
# all caught. `config.toml` is included as a presence-only marker — it
# never appears in the deletion entries (it's preserved), but its
# presence is the strongest signal "yes, this is an allelix data dir."
_KNOWN_ARTIFACT_PREFIXES: frozenset[str] = frozenset(
    {
        "alphamissense",
        "cadd",
        "clinicalAnnotations",
        "clinvar",
        "config.toml",
        "gnomad",
        "gwas",
        "pharmgkb",
        "snpedia",
    }
)


def _is_cadd_cache(name: str) -> bool:
    return any(name.startswith(p) for p in _CADD_PREFIXES)


def _looks_like_data_dir(resolved: Path) -> bool:
    """Return True if the dir contains at least one recognized allelix artifact.

    The guard for `db clean`: a fat-fingered `--data-dir ~/Documents` should
    not result in `rm -rf ~/Documents`. We scan the dir for any name matching
    a known cache-file prefix; absence of any match means the user almost
    certainly mistyped the path.
    """
    return any(
        any(entry.name.startswith(p) for p in _KNOWN_ARTIFACT_PREFIXES)
        for entry in resolved.iterdir()
    )


def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{n:,} {unit}"
        n /= 1024  # type: ignore[assignment]
    return f"{n:,.1f} TB"


def _iter_cache_entries(resolved: Path, *, keep_cadd: bool) -> list[Path]:
    """Return data-dir entries to be deleted by `db clean`.

    Excludes `config.toml` (per ADR-0006: the cache dir doubles as the
    XDG-shaped config home in older installs; never delete config). When
    `keep_cadd` is True, also excludes CADD cache files.
    """
    out: list[Path] = []
    for entry in sorted(resolved.iterdir()):
        if entry.name == "config.toml":
            continue
        if keep_cadd and _is_cadd_cache(entry.name):
            continue
        out.append(entry)
    return out


def _entry_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


@db.command("clean")
@_DATA_DIR_OPT
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show what would be deleted without removing anything.",
)
@click.option(
    "--keep-cadd",
    is_flag=True,
    default=False,
    help="Preserve the CADD cache (the largest opt-in download, ~5.8 GB).",
)
@click.option(
    "--yes",
    "skip_confirm",
    is_flag=True,
    default=False,
    help="Skip the confirmation prompt (for scripted use).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help=(
        "Bypass the looks-like-an-allelix-cache safety guard. Required when "
        "the target directory contains content but no recognized allelix "
        "cache files (e.g. you've staged caches under a non-standard name)."
    ),
)
def db_clean(
    data_dir: Path | None,
    dry_run: bool,
    keep_cadd: bool,
    skip_confirm: bool,
    force: bool,
) -> None:
    """Remove downloaded reference database caches to reclaim disk space.

    The caches are disposable — every file removed here can be re-fetched
    by `allelix db update`. Your config.toml (license assertions, source
    toggles) is preserved; it lives outside the cache.

    Pass `--keep-cadd` to preserve the CADD cache (largest opt-in
    download, ~5.8 GB; expensive to re-download). Pass `--dry-run` to
    preview the deletion without acting. Pass `--yes` to skip the
    confirmation prompt for scripted use.
    """
    import shutil
    import sys

    # resolve_data_dir creates the dir if missing — existence is guaranteed
    # but the dir may be empty (matched as a no-op below).
    resolved = resolve_data_dir(data_dir)
    entries = _iter_cache_entries(resolved, keep_cadd=keep_cadd)

    # Safety guard: refuse to delete content from a directory that doesn't
    # look like an allelix cache. Runs BEFORE --yes / --dry-run so a typo'd
    # --data-dir cannot result in `rm -rf` against an unrelated tree even
    # in the scripted-skip-confirm path. --force is the explicit opt-out.
    # Plain print to stderr (same routing as `db path --check`) — Rich
    # console.print would wrap the message at the terminal width and
    # would route to stdout, both of which are wrong for an error.
    if entries and not force and not _looks_like_data_dir(resolved):
        print(
            f"error: {resolved} doesn't look like an allelix cache directory "
            "(no recognized cache files found). Refusing to delete. "
            "Pass --force if this is intended.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not entries:
        console.print(
            f"[dim]Nothing to clean at {resolved}"
            f"{' (CADD preserved via --keep-cadd)' if keep_cadd else ''}[/dim]"
        )
        return

    table = Table(title=f"{'Would delete' if dry_run else 'Will delete'} ({resolved})")
    table.add_column("Item", style="cyan", no_wrap=True)
    table.add_column("Size", justify="right")
    total = 0
    for entry in entries:
        size = _entry_size(entry)
        total += size
        table.add_row(entry.name, _human_bytes(size))
    table.add_row("[bold]total[/bold]", f"[bold]{_human_bytes(total)}[/bold]")
    console.print(table)

    if keep_cadd:
        console.print(
            "[dim]CADD cache preserved (--keep-cadd). Re-running `db update --cadd` "
            "would otherwise re-download ~5.8 GB.[/dim]"
        )

    if dry_run:
        console.print("[dim]Dry run — no files removed.[/dim]")
        return

    if not skip_confirm and not click.confirm(
        f"Delete {len(entries)} item(s) ({_human_bytes(total)})?",
        default=False,
    ):
        console.print("[yellow]Aborted.[/yellow]")
        return

    removed = 0
    for entry in entries:
        if entry.is_file() or entry.is_symlink():
            entry.unlink()
        else:
            shutil.rmtree(entry)
        removed += 1
    console.print(
        f"[green]✓ Removed {removed} item(s) ({_human_bytes(total)}).[/green] "
        "Re-populate with [bold]allelix db update[/bold]."
    )


@db.command("path")
@_DATA_DIR_OPT
@click.option(
    "--check",
    is_flag=True,
    default=False,
    help="Verify the path exists and is writable (non-zero exit on failure).",
)
def db_path(data_dir: Path | None, check: bool) -> None:
    """Print the resolved reference-database cache directory.

    Useful for scripting / backup integration:

        ALLELIX_DATA=$(allelix db path)
        du -sh "$ALLELIX_DATA"

    With `--check`, additionally verify the path is writable; exit
    non-zero (and print a diagnostic to stderr) if it is not.
    Existence is guaranteed — the path is created on resolution.
    """
    import os
    import sys

    resolved = resolve_data_dir(data_dir)
    # Plain print (not rich console.print) so the output is shell-safe
    # for `$(allelix db path)` capture — no styling escape sequences.
    print(str(resolved))

    if check and not os.access(resolved, os.W_OK):
        print(f"error: path is not writable: {resolved}", file=sys.stderr)
        sys.exit(1)


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

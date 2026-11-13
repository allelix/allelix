# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The `db` subcommand group — local reference database cache management."""

from __future__ import annotations

from typing import TYPE_CHECKING

import click
from rich.table import Table

from allelix.annotators import get_annotators
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator, signal_lacks_cpic
from allelix.cli import _helpers, main
from allelix.cli._helpers import console
from allelix.cli._options import _DATA_DIR_OPT
from allelix.databases import resolve_data_dir
from allelix.databases.cache_state import CacheState
from allelix.databases.manager import compare_signals
from allelix.utils.markup import literal

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.annotators.base import Annotator


@main.group()
def db() -> None:
    """Manage local reference database cache."""


def _print_refreshed(annotator: Annotator) -> None:
    """The line after a successful refresh; a download with unchanged content says so."""
    version = literal(annotator.version() or "(unknown)")
    if getattr(annotator, "content_unchanged", False):
        console.print(
            f"  [green]✓ {annotator.name}: content unchanged; cache kept[/green] "
            f"(version {version})"
        )
    else:
        console.print(f"  [green]✓ {annotator.name} refreshed[/green] (version {version})")
    _print_source_skips(annotator)


def _print_source_skips(annotator: Annotator) -> None:
    """The saved skip counts under a source's line, where it keeps them."""
    skip_line = _source_skip_line(annotator)
    if skip_line:
        console.print(f"    {skip_line}")


def _source_skip_line(annotator: Annotator) -> str | None:
    """The release's skipped source rows, for a source that counts them (GWAS)."""
    counts_fn = getattr(annotator, "source_skip_counts", None)
    counts = counts_fn() if callable(counts_fn) else None
    if counts is None:
        return None
    return (
        f"{annotator.display_name}: source rows skipped because the risk allele "
        f"is not the row's own SNP's: other_snp {counts['other_snp']:,}, "
        f"unattributed {counts['unattributed']:,} (not in the record count)"
    )


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
    help=(
        "Re-download even if the local cache appears current. ClinVar is still "
        "refreshed only when NCBI's checksums can be fetched; CADD full mode "
        "downloads nothing."
    ),
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
    help=(
        "Download CADD deleteriousness scores (non-commercial use, or commercial use "
        "with a CADD license asserted as license.cadd; disabled by default)."
    ),
)
@click.option(
    "--build",
    type=click.Choice(["grch37", "grch38", "both"], case_sensitive=False),
    default="both",
    help=(
        "Which ClinVar genome build(s) to download. 'both' (default) keeps "
        "GRCh37 and GRCh38 caches in sync so `analyze` can dispatch by "
        "the build it decides (ADR-0044). 'grch37' / 'grch38' restrict to one to "
        "save bandwidth."
    ),
)
@click.option(
    "--only",
    "only_source",
    type=click.Choice(["clinvar"], case_sensitive=False),
    default=None,
    help="Update only the named database (currently supported: clinvar).",
)
@click.option(
    "--allow-clinvar-shrink",
    is_flag=True,
    default=False,
    help=(
        "For this run only, accept a ClinVar release with fewer rows than 90% of the "
        "cache it replaces. The fixed row floor and the zero-row check still apply; "
        "old and new counts are printed and logged."
    ),
)
def db_update(
    data_dir: Path | None,
    force: bool,
    no_gnomad: bool,
    no_alphamissense: bool,
    include_cadd: bool,
    build: str,
    only_source: str | None,
    allow_clinvar_shrink: bool,
) -> None:
    """Download or refresh reference databases.

    \b
    For each annotator:
      - no cache → download
      - --force → download
      - CADD full mode → nothing is downloaded, even with --force: ready
        when its tabix files are in place and readable, otherwise the
        update fails naming the file (db update cannot fetch them)
      - cache + remote signal matches cache → skip
      - cache + remote signal differs (or legacy v0.4.1 cache with no
        stored signal) → download
      - cache + remote sends a different freshness header than the one
        stored (ETag or Last-Modified) → keep with notice (use --force)
      - cache + remote signal can't be fetched → skip with notice (use
        --force to override), except ClinVar: it is never refreshed without
        NCBI's checksums, even with --force
      - ClinPGx: a download whose content matches the last build keeps
        the cache ("content unchanged")

    `--build` selects which ClinVar build(s) to manage. Default 'both'
    downloads GRCh37 and GRCh38 caches.
    """  # noqa: D301 - plain string: click's \b marker is the backspace character
    resolved = resolve_data_dir(data_dir)
    console.print(f"Data directory: [cyan]{literal(str(resolved))}[/cyan]")

    from allelix.config import load_config

    cfg = load_config(resolved)

    clinvar_builds = _helpers._resolve_clinvar_builds(build)
    failed: list[str] = []
    annotators = get_annotators(resolved, clinvar_builds=clinvar_builds, cadd_full=cfg.cadd_full)
    # ADR-0048: convert WAL caches and bring outdated ones up to date from
    # local data before deciding what to download.
    _helpers._prepare_caches(
        [a for a in annotators if only_source is None or a.name == only_source], resolved
    )
    for annotator in annotators:
        if isinstance(annotator, ClinVarAnnotator):
            annotator.allow_shrink = allow_clinvar_shrink
        with annotator:
            if only_source is not None and annotator.name != only_source:
                continue
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
                if (
                    annotator.update_supplies_backend()
                    and (not annotator.is_ready() or force)
                    and not _confirm_cadd_license(
                        license_held=cfg.license_held("cadd"),
                    )
                ):
                    console.print(f"  [dim]{annotator.name}: skipped (declined)[/dim]")
                    continue

            # ADR-0048: a cache still in WAL mode after the update step
            # could not be converted. It is not downloaded again.
            status = annotator.cache_status()
            if not force and status.state in (CacheState.WAL, CacheState.READONLY_WAL):
                console.print(
                    f"  [yellow]{annotator.name}: {literal(str(status.path))} is still in WAL "
                    "mode and could not be converted (no write access to its directory, or "
                    "another process has it open). Nothing was downloaded.[/yellow]"
                )
                failed.append(annotator.name)
                continue

            if not annotator.update_supplies_backend():
                # Nothing to download: the files it reads are the user's (CADD
                # full mode). If they cannot be used (missing, unreadable), the
                # reason is printed and the source fails.
                status = annotator.cache_status()
                if status.ready:
                    console.print(
                        f"  [dim]{annotator.name}: ready "
                        f"({literal(annotator.version() or 'unknown')})[/dim]"
                    )
                elif annotator.requires_download:
                    # A downloaded source whose backend update cannot supply
                    # (CADD full mode): it fails. A source that downloads
                    # nothing is not managed here and is skipped, as before.
                    # soft_wrap: a long path is never broken across lines.
                    console.print(
                        f"  [yellow]{literal(_helpers.unready_line(annotator, status))}[/yellow]",
                        soft_wrap=True,
                    )
                    failed.append(annotator.name)
                continue

            if not annotator.is_ready():
                console.print(f"  [bold]{annotator.name}[/bold]: downloading…")
                _check_setup_destinations(annotator, resolved)
                if _helpers._run_setup(annotator, force=force):
                    console.print(
                        f"  [green]✓ {annotator.name} ready[/green] "
                        f"(version {literal(annotator.version() or '(unknown)')})"
                    )
                    _print_source_skips(annotator)
                else:
                    failed.append(annotator.name)
                continue

            if force:
                console.print(f"  [bold]{annotator.name}[/bold]: --force; refreshing…")
                _check_setup_destinations(annotator, resolved)
                if _helpers._run_setup(annotator, force=force):
                    _print_refreshed(annotator)
                else:
                    failed.append(annotator.name)
                continue

            # Code-driven sources (commit-pinned HF caches) are updated
            # only via code changes — no runtime freshness probe needed.
            if not annotator.server_driven_freshness:
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {literal(annotator.version() or '(unknown)')})[/dim]"
                )
                continue

            remote = annotator.fetch_remote_signal()
            if remote is None:
                # ClinVar's signal is the checksum that verifies its download:
                # --force cannot refresh without it either.
                advice = (
                    "; cache kept. ClinVar is never refreshed unverified; try again later."
                    if annotator.name == "clinvar"
                    else ". Pass --force to refresh anyway."
                )
                console.print(
                    f"  [yellow]{annotator.name}: cache present, but remote "
                    f"freshness can't be verified (network error or no signal){advice}[/yellow]"
                )
                _print_source_skips(annotator)
                continue

            cached = annotator.cached_remote_signal()
            verdict = compare_signals(cached, remote) if cached is not None else "changed"
            if verdict == "same":
                console.print(
                    f"  [dim]{annotator.name}: already current "
                    f"(version {literal(annotator.version() or '(unknown)')})[/dim]"
                )
                _print_source_skips(annotator)
                continue
            if verdict == "unknown":
                console.print(
                    f"  [yellow]{annotator.name}: cache kept; the remote sends different "
                    "freshness fields than the cache recorded, so a change cannot be "
                    "verified. Pass --force to refresh anyway.[/yellow]"
                )
                _print_source_skips(annotator)
                continue

            if cached is None:
                # a cache with no stored freshness signal almost
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
            _check_setup_destinations(annotator, resolved)
            if _helpers._run_setup(annotator):
                _print_refreshed(annotator)
            else:
                failed.append(annotator.name)

    if failed:
        import sys

        console.print(f"[red]Failed to install or refresh: {', '.join(failed)}[/red]")
        sys.exit(1)


def _owned_cache_names() -> tuple[str, ...]:
    """The file name of every cache allelix builds or installs, legacy ones included."""
    from allelix.annotators.clinvar import CLINVAR_SUPPORTED_BUILDS, clinvar_db_filename
    from allelix.annotators.snpedia import SNPEDIA_DB_FILENAME
    from allelix.databases.alphamissense_loader import ALPHAMISSENSE_DB_FILENAME
    from allelix.databases.cadd_loader import CADD_DB_FILENAME
    from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
    from allelix.databases.gwas_loader import GWAS_DB_FILENAME
    from allelix.databases.pharmgkb_loader import PHARMGKB_DB_FILENAME
    from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME

    prebuilt = (
        GNOMAD_DB_FILENAME,
        ALPHAMISSENSE_DB_FILENAME,
        CADD_DB_FILENAME,
        SNPEDIA_DB_FILENAME,
        RSID_IDENTITY_DB_FILENAME,
    )
    return (
        *(clinvar_db_filename(build) for build in CLINVAR_SUPPORTED_BUILDS),
        # v0.7.x's single ClinVar cache, before the per-build files (2.6.1).
        "clinvar.sqlite",
        *prebuilt,
        GWAS_DB_FILENAME,
        PHARMGKB_DB_FILENAME,
    )


def _owned_final_names() -> frozenset[str]:
    """Every name allelix writes into the data directory, without staging names.

    The caches; the files kept beside them; and the downloads earlier releases
    wrote into the data directory itself (prebuilt ``.gz`` files, the GWAS and
    ClinPGx zips, ClinVar's VCF before the tab-delimited source).
    """
    from allelix.annotators.gwas import GWAS_TSV_FILENAME
    from allelix.annotators.snpedia import SNPEDIA_DB_FILENAME
    from allelix.databases.alphamissense_loader import ALPHAMISSENSE_DB_FILENAME
    from allelix.databases.cadd_loader import CADD_DB_FILENAME
    from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
    from allelix.databases.manager import CLINVAR_BUILD_LOG
    from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME

    prebuilt = (
        GNOMAD_DB_FILENAME,
        ALPHAMISSENSE_DB_FILENAME,
        CADD_DB_FILENAME,
        SNPEDIA_DB_FILENAME,
        RSID_IDENTITY_DB_FILENAME,
    )
    downloads = (
        *(f"{name}.gz" for name in prebuilt),
        "gwas_catalog_associations.zip",
        "clinicalAnnotations.zip",
        "clinvar.vcf.gz",
        "clinvar.vcf",
    )
    return frozenset((*_owned_cache_names(), GWAS_TSV_FILENAME, CLINVAR_BUILD_LOG, *downloads))


# Prefixes of the staging directories (``owned_directory``) allelix makes.
_STAGING_DIRECTORY_PREFIXES: frozenset[str] = frozenset(
    {
        "alphamissense",
        "cadd",
        "clinvar",
        "gnomad",
        "gwas",
        "pharmgkb",
        "rsid_identity",
        "snpedia",
    }
)

# Suffixes allelix adds to a name above: SQLite's side files, and the staging
# names earlier releases shared between writers (``<name>.tmp``, ``.part``,
# ``.partial``), with their own side files.
_SIDE_SUFFIXES: tuple[str, ...] = ("-wal", "-shm", "-journal")
_SHARED_STAGING_SUFFIXES: tuple[str, ...] = (".tmp", ".part", ".partial")


def _owned_base(name: str, finals: frozenset[str]) -> str | None:
    """The name in ``finals`` that ``name`` belongs to, if allelix created it."""
    from allelix.databases.staging import staging_file_target

    target = staging_file_target(name)
    if target is not None:
        return target if target in finals else None
    for side in ("", *_SIDE_SUFFIXES):
        stem = name.removesuffix(side) if side else name
        if side and stem == name:
            continue
        if stem in finals:
            return stem
        for shared in _SHARED_STAGING_SUFFIXES:
            if stem.endswith(shared) and stem.removesuffix(shared) in finals:
                return stem.removesuffix(shared)
    return None


# What a staging directory holds while a download runs: the downloads
# themselves, the GWAS table once moved out of ``extracted/``, their partial
# (``owned_staging``) files, and a cache's partial file while it builds there
# (``.<cache>.<16 hex>.partial`` and its side files). A directory holding
# anything else (a finished cache built or copied there, an extraction, files
# of the user's own) is not removed, only listed (2.6.1).
_STAGING_DOWNLOAD_NAMES: frozenset[str] = frozenset(
    {
        "variant_summary.txt.gz",
        "submission_summary.txt.gz",
        "gnomad.sqlite.gz",
        "alphamissense.sqlite.gz",
        "cadd.sqlite.gz",
        "snpedia.sqlite.gz",
        "rsid_identity.sqlite.gz",
        "gwas_catalog_associations.zip",
        "gwas_catalog_associations.tsv",
        "clinicalAnnotations.zip",
    }
)


def _is_staging_download(entry: Path) -> bool:
    """An exact name allelix writes in a staging directory.

    A download, a partial file of a download or of a cache, or an empty
    ``extracted/``. A non-empty ``extracted/`` holds a table whose name comes
    from EBI's zip, so it cannot be told from a user's file; a directory
    holding one is kept.
    """
    from allelix.databases.staging import staging_file_target

    if entry.is_symlink():
        return False
    if entry.is_dir():
        return entry.name == "extracted" and not any(entry.iterdir())
    name = entry.name
    target = staging_file_target(name)
    return (
        name in _STAGING_DOWNLOAD_NAMES
        or target in _STAGING_DOWNLOAD_NAMES
        or target in _owned_cache_names()
    )


def _has_staging_shape(entry: Path) -> bool:
    """A real directory named like allelix's staging (``.<source>-<8 characters>``)."""
    from allelix.databases.staging import staging_directory_prefix

    return (
        staging_directory_prefix(entry.name) in _STAGING_DIRECTORY_PREFIXES
        and entry.is_dir()
        and not entry.is_symlink()
    )


def _is_owned_staging_directory(entry: Path) -> bool:
    """A real directory of allelix's staging shape holding only names allelix writes there.

    The shape alone (``.<source>-<8 characters>``) is not enough: a user's
    dated backup such as ``.clinvar-20261101/`` has it too.
    """
    if not _has_staging_shape(entry):
        return False
    try:
        return all(_is_staging_download(child) for child in entry.iterdir())
    except OSError:
        return False


def _is_cadd_name(entry: Path, base: str | None) -> bool:
    from allelix.databases.cadd_loader import CADD_DB_FILENAME
    from allelix.databases.staging import staging_directory_prefix

    if base is not None:
        return base.startswith(CADD_DB_FILENAME)
    return staging_directory_prefix(entry.name) == "cadd"


def _check_setup_destinations(annotator: Annotator, resolved: Path) -> None:
    """Stop before a download whose results could not be written (2.5.2)."""
    import os

    from allelix.reports import UnwritableDestinationError, check_writable

    try:
        for path in annotator.setup_destinations():
            check_writable(path)
        if not os.access(resolved, os.W_OK | os.X_OK):
            msg = f"Cannot write to the data directory {resolved}: it is not writable."
            raise UnwritableDestinationError(msg)
    except OSError as exc:
        raise click.ClickException(str(exc)) from exc


def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{n:,} {unit}"
        n /= 1024  # type: ignore[assignment]
    return f"{n:,.1f} TB"


def _clean_plan(resolved: Path, *, keep_cadd: bool) -> tuple[list[Path], list[Path]]:
    """Split the data directory into the entries `db clean` removes and those it keeps.

    It removes only entries allelix creates: caches, their side files, the
    files kept beside them, staging entries and earlier releases' downloads.
    A name allelix gives a file counts only for a file or a symbolic link (a
    link is removed, never followed); a directory is removed only when it
    has the shape of a staging directory allelix makes and holds nothing but
    names allelix writes there (2.6.1). Every other entry is
    kept, among them CADD full-mode files, which the user places there, a
    real directory under a cache's name, and `config.toml` (ADR-0006). With
    `keep_cadd`, the CADD cache's own entries are kept too. The kept list
    leaves out `config.toml` and those CADD entries, which `db clean` names
    on its own.
    """
    finals = _owned_final_names()
    remove: list[Path] = []
    keep: list[Path] = []
    for entry in sorted(resolved.iterdir()):
        is_file_like = entry.is_symlink() or not entry.is_dir()
        base = _owned_base(entry.name, finals) if is_file_like else None
        owned = base is not None or _is_owned_staging_directory(entry)
        if not owned:
            if entry.name != "config.toml":
                keep.append(entry)
        elif not (keep_cadd and _is_cadd_name(entry, base)):
            remove.append(entry)
    return remove, keep


def _entry_size(path: Path) -> int:
    if path.is_symlink():
        return path.lstat().st_size  # only the link is removed
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
    help=(
        "Keep the CADD cache (cadd.sqlite, the largest opt-in download, ~5.8 GB). "
        "CADD full-mode files are kept either way: allelix did not create them."
    ),
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
        "Accepted for scripts written for 2.6.0 and earlier; it has no effect. "
        "db clean removes only names allelix creates, so it needs no override."
    ),
)
def db_clean(
    data_dir: Path | None,
    dry_run: bool,
    keep_cadd: bool,
    skip_confirm: bool,
    force: bool,
) -> None:
    """Remove the reference database caches allelix downloaded and built.

    Only entries allelix creates in the data directory are removed: the
    caches and their side files, GWAS's kept table, the ClinVar build log,
    staging entries a killed update left, and earlier releases' downloads.
    A directory is removed only when it is one of allelix's staging
    directories holding nothing but names allelix writes there: downloads,
    the partial files of a download or of a cache being built, an empty
    extracted/. One holding anything else (a backup, a finished cache, an
    extracted table) is kept and listed. A link is removed, never followed.
    Everything allelix did
    not create stays, and is listed: CADD full-mode files, a genome file or
    reports kept there, a directory under a cache's name. Your config.toml
    (license assertions, source toggles) is always kept, and is not listed.

    Do not run db clean while a database update is running in the same data
    directory, including an analyze without --no-update, whose refresh can
    build a cache there: it would remove that update's staging files.

    `allelix db update` downloads and builds the caches again. The ClinVar
    build log (clinvar-build.log, the sources and counts of past builds)
    cannot be fetched again; copy it elsewhere first to keep it.

    Pass `--keep-cadd` to keep the CADD cache (largest opt-in download,
    ~5.8 GB; expensive to re-download). Pass `--dry-run` to preview the
    deletion without acting. Pass `--yes` to skip the confirmation prompt
    for scripted use.
    """
    import os
    import shutil

    del force  # no guard to bypass: nothing allelix did not create is removed
    # resolve_data_dir creates the dir if missing — existence is guaranteed
    # but the dir may be empty (matched as a no-op below).
    resolved = resolve_data_dir(data_dir)
    entries, kept = _clean_plan(resolved, keep_cadd=keep_cadd)
    kept_note = ""
    # A staging-shaped folder that is kept may be an interrupted update's, or
    # the user's own: allelix cannot tell, so it says so (2.6.1).
    staging_kept = [entry for entry in kept if _has_staging_shape(entry)]
    others = [entry for entry in kept if entry not in staging_kept]
    notes = []
    if others:
        shown = ", ".join(literal(entry.name) for entry in others[:6])
        more = f" and {len(others) - 6} more" if len(others) > 6 else ""
        notes.append(f"Kept {len(others)} item(s) allelix did not create: {shown}{more}.")
    for entry in staging_kept:
        notes.append(
            f"Kept {literal(entry.name)}: may be an interrupted update's leftovers; "
            "delete it by hand if so."
        )
    kept_note = "\n".join(notes)

    if not entries:
        console.print(
            f"[dim]Nothing to clean at {literal(str(resolved))}"
            f"{' (CADD preserved via --keep-cadd)' if keep_cadd else ''}[/dim]"
        )
        if kept_note:
            console.print(f"[dim]{kept_note}[/dim]")
        return

    table = Table(
        title=f"{'Would delete' if dry_run else 'Will delete'} ({literal(str(resolved))})"
    )
    table.add_column("Item", style="cyan", no_wrap=True)
    table.add_column("Size", justify="right")
    total = 0
    for entry in entries:
        size = _entry_size(entry)
        total += size
        table.add_row(literal(entry.name), _human_bytes(size))
    table.add_row("[bold]total[/bold]", f"[bold]{_human_bytes(total)}[/bold]")
    console.print(table)

    if keep_cadd:
        console.print(
            "[dim]CADD cache preserved (--keep-cadd). Re-running `db update --cadd` "
            "would otherwise re-download ~5.8 GB.[/dim]"
        )
    if kept_note:
        console.print(f"[dim]{kept_note}[/dim]")

    if dry_run:
        console.print("[dim]Dry run — no files removed.[/dim]")
        return

    if not skip_confirm and not click.confirm(
        f"Delete {len(entries)} item(s) ({_human_bytes(total)})?",
        default=False,
    ):
        console.print("[yellow]Aborted.[/yellow]")
        return

    # Nothing is deleted from a directory that cannot be changed (2.5.2).
    if not os.access(resolved, os.W_OK | os.X_OK):
        raise click.ClickException(
            f"Cannot delete from the data directory {resolved}: it is not writable."
        )
    removed = 0
    try:
        for entry in entries:
            if entry.is_file() or entry.is_symlink():
                entry.unlink()
            else:
                shutil.rmtree(entry)
            removed += 1
    except OSError as exc:
        raise click.ClickException(
            f"Cannot delete {exc.filename or resolved}: {exc.strerror or exc}. "
            f"{removed} item(s) were removed first."
        ) from exc
    from allelix.databases.manager import CLINVAR_BUILD_LOG

    log_note = (
        " The ClinVar build log cannot be fetched again."
        if any(entry.name == CLINVAR_BUILD_LOG for entry in entries)
        else ""
    )
    console.print(
        f"[green]✓ Removed {removed} item(s) ({_human_bytes(total)}).[/green] "
        f"[bold]allelix db update[/bold] downloads and builds the caches again.{log_note}"
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

    \b
        ALLELIX_DATA=$(allelix db path)
        du -sh "$ALLELIX_DATA"

    With `--check`, additionally verify the path is writable; exit
    non-zero (and print a diagnostic to stderr) if it is not.
    Existence is guaranteed — the path is created on resolution.
    """  # noqa: D301 - plain string: click's \b marker is the backspace character
    import sys

    from allelix.databases.cache_state import directory_writable

    resolved = resolve_data_dir(data_dir)
    # Plain print (not rich console.print) so the output is shell-safe
    # for `$(allelix db path)` capture — no styling escape sequences.
    print(str(resolved))

    if check and not directory_writable(resolved, require_search=False):
        print(f"error: path is not writable: {resolved}", file=sys.stderr)
        sys.exit(1)


@db.command("status")
@_DATA_DIR_OPT
def db_status(data_dir: Path | None) -> None:
    """Show installed reference databases: readiness, version and record count."""
    from allelix.config import load_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)
    table = Table(title=f"Reference Databases ({literal(str(resolved))})")
    table.add_column("Annotator", style="cyan", no_wrap=True)
    table.add_column("Ready", justify="center")
    table.add_column("Version")
    table.add_column("Records", justify="right")
    builds = _helpers.present_clinvar_builds(resolved)
    skip_lines: list[str] = []
    for annotator in get_annotators(resolved, cadd_full=cfg.cadd_full, clinvar_builds=builds):
        with annotator:
            status = annotator.cache_status()
            if status.state in (CacheState.WAL, CacheState.READONLY_WAL):
                # ADR-0048: a WAL cache is not opened, even for its version.
                table.add_row(
                    annotator.display_name,
                    "[red]no[/red]",
                    "WAL mode from an older Allelix; run `allelix db update` to convert it",
                    "—",
                )
                continue
            ready = status.ready
            ready_marker = "[green]yes[/green]" if ready else "[red]no[/red]"
            version = annotator.version() or "—"
            sig = getattr(annotator, "cached_remote_signal", lambda: None)()
            if isinstance(annotator, PharmGKBAnnotator) and signal_lacks_cpic(sig):
                version += " (no CPIC)"
            records = "—"
            count_fn = getattr(annotator, "record_count", None)
            if callable(count_fn):
                count = count_fn()
                if count is not None:
                    records = f"{count:,}"
            table.add_row(annotator.display_name, ready_marker, literal(version), records)
            skip_line = _source_skip_line(annotator)
            if skip_line:
                skip_lines.append(skip_line)
    console.print(table)
    for line in skip_lines:
        console.print(line)

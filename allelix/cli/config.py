# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The `config` subcommand group — persistent configuration management."""

from __future__ import annotations

from typing import TYPE_CHECKING

import click
from rich.table import Table

from allelix.cli import main
from allelix.cli._helpers import console
from allelix.cli._options import _DATA_DIR_OPT
from allelix.databases import resolve_data_dir

if TYPE_CHECKING:
    from pathlib import Path


@main.group()
def config() -> None:
    """Manage persistent configuration (source toggles, license mode)."""


@config.command("show")
@_DATA_DIR_OPT
def config_show(data_dir: Path | None) -> None:
    """Display current configuration."""
    from allelix.annotators import _ANNOTATOR_CLASSES
    from allelix.annotators.base import Permission
    from allelix.config import load_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)

    table = Table(title=f"Configuration ({resolved / 'config.toml'})")
    table.add_column("Source", style="cyan", no_wrap=True)
    table.add_column("Enabled", justify="center")
    table.add_column("Note", style="dim")
    for name, enabled in sorted(cfg.sources.items()):
        perm = cfg.permission_for(name, _ANNOTATOR_CLASSES)
        note = ""
        if perm is Permission.BLOCK_PURCHASABLE:
            # cls is guaranteed non-None — permission_for only returns
            # BLOCK_PURCHASABLE when an annotator class was resolved.
            cls = _ANNOTATOR_CLASSES[name]
            marker = "[red]no[/red]"
            note = f"requires commercial license — purchase: {cls.license.purchase_url}"
        elif perm is Permission.BLOCK_FINAL:
            marker = "[red]no[/red]"
            note = "no commercial license is available"
        elif enabled:
            marker = "[green]yes[/green]"
        else:
            marker = "[red]no[/red]"
        table.add_row(name, marker, note)
    console.print(table)
    mode = "[yellow]commercial[/yellow]" if cfg.commercial else "[green]personal[/green]"
    console.print(f"License mode: {mode}")


@config.command("get")
@_DATA_DIR_OPT
@click.argument("key", required=False, default=None)
def config_get(data_dir: Path | None, key: str | None) -> None:
    r"""Get a configuration value (or dump entire config).

    \b
    Keys:
      sources.<name>       Show if a source is enabled
      license.commercial   Show commercial mode
      license.<source>     Show if a license is asserted for <source>
      options.cadd_full    Show full CADD tabix mode

    \b
    Examples:
      allelix config get                     # dump entire config
      allelix config get sources.cadd        # true
      allelix config get license.cadd        # false
      allelix config get options.cadd_full   # false
    """
    from allelix.config import _serialize, load_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)

    if key is None:
        console.print(f"[dim]Config: {resolved / 'config.toml'}[/dim]")
        click.echo(_serialize(cfg))
        return

    if key.startswith("sources."):
        source_name = key[len("sources.") :]
        val = cfg.sources.get(source_name)
        if val is None:
            raise click.ClickException(
                f"Unknown source {source_name!r}. Known sources: {', '.join(sorted(cfg.sources))}"
            )
        click.echo(str(val).lower())
    elif key == "license.commercial":
        click.echo(str(cfg.commercial).lower())
    elif key.startswith("license."):
        source_name = key[len("license.") :]
        click.echo(str(cfg.license_held(source_name)).lower())
    elif key == "options.cadd_full":
        click.echo(str(cfg.cadd_full).lower())
    else:
        raise click.ClickException(
            f"Unknown key {key!r}. Use 'sources.<name>', 'license.commercial', "
            "'license.<source>', or 'options.cadd_full'."
        )


@config.command("set")
@_DATA_DIR_OPT
@click.argument("key")
@click.argument("value")
def config_set(data_dir: Path | None, key: str, value: str) -> None:
    r"""Set a configuration value.

    \b
    Keys:
      sources.<name>       Enable/disable a source (true/false)
      license.commercial   Set commercial mode (true/false)
      license.<source>     Assert you hold a commercial license for <source>
      options.cadd_full    Use full CADD tabix file instead of cache (true/false)

    \b
    Examples:
      allelix config set sources.snpedia false
      allelix config set license.commercial true
      allelix config set license.cadd true
      allelix config set options.cadd_full true
    """
    from allelix.config import load_config, save_config

    resolved = resolve_data_dir(data_dir)
    cfg = load_config(resolved)

    val_lower = value.strip().lower()
    if val_lower not in ("true", "false"):
        raise click.ClickException(f"Value must be 'true' or 'false', got {value!r}")
    bool_val = val_lower == "true"

    if key.startswith("sources."):
        source_name = key[len("sources.") :]
        cfg.sources[source_name] = bool_val
    elif key == "license.commercial":
        cfg.commercial = bool_val
    elif key.startswith("license."):
        source_name = key[len("license.") :]
        if bool_val:
            from allelix.annotators import get_annotator_class

            cls = get_annotator_class(source_name)
            if cls is not None and not cls.license.licensable:
                raise click.ClickException(
                    f"{source_name} is not commercially licensable. "
                    f"This assertion has no effect and cannot be set."
                )
            cfg.license_overrides[source_name] = True
        else:
            from allelix.annotators import get_annotator_class

            if (
                get_annotator_class(source_name) is None
                and source_name not in cfg.license_overrides
            ):
                console.print(f"[yellow]Warning: unknown source {source_name!r}[/yellow]")
            cfg.license_overrides.pop(source_name, None)
    elif key == "options.cadd_full":
        cfg.cadd_full = bool_val
    else:
        raise click.ClickException(
            f"Unknown key {key!r}. Use 'sources.<name>', 'license.commercial', "
            "'license.<source>', or 'options.cadd_full'."
        )

    save_config(resolved, cfg)
    console.print(f"[dim]Config: {resolved / 'config.toml'}[/dim]")
    console.print(f"[green]Set {key} = {val_lower}[/green]")

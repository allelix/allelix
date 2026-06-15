# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Allelix command-line interface.

The CLI is organized as a package: `main` (the click group) is defined
here, and each subcommand module decorates against it. Importing this
package triggers all subcommand registrations as a side effect.

Public API: ``from allelix.cli import main`` — this is what the package
entry point ``allelix = "allelix.cli:main"`` calls.
"""

from __future__ import annotations

import click

from allelix import __version__


@click.group()
@click.version_option(version=__version__, prog_name="allelix")
def main() -> None:
    """Allelix: open-source genotype analysis toolkit."""


# Side-effect imports: each module decorates against `main` at import
# time. Order matters only for --help listing, not correctness.
from allelix.cli import (  # noqa: E402,F401
    analyze,
    config,
    db,
    focused,
    utility,
)

__all__ = ["main"]

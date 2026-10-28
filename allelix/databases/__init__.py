# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Local database cache management. See ADR-0004 for the offline-first model."""

from __future__ import annotations

import os
from pathlib import Path

ENV_VAR = "ALLELIX_DATA_DIR"


def default_data_dir() -> Path:
    """Resolve the default cache location.

    Precedence: ``$ALLELIX_DATA_DIR`` > ``$XDG_DATA_HOME/allelix`` >
    ``~/.local/share/allelix``. See ADR-0006.
    """
    if override := os.environ.get(ENV_VAR):
        return Path(override).expanduser()
    if xdg := os.environ.get("XDG_DATA_HOME"):
        return Path(xdg).expanduser() / "allelix"
    return Path.home() / ".local" / "share" / "allelix"


def resolve_data_dir(override: Path | None = None) -> Path:
    """Return the data directory to use, creating it if it doesn't exist.

    A directory that cannot be created raises ``UnwritableDestinationError``
    naming it (2.5.2).
    """
    path = override.expanduser() if override else default_data_dir()
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        from allelix.reports import UnwritableDestinationError

        msg = f"Cannot create the data directory {path}: {exc.strerror or exc}."
        raise UnwritableDestinationError(msg) from exc
    return path

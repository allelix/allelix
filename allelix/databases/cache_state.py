# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Read-only readiness of a cache file (ADR-0048).

``cache_state`` is the one readiness truth: every annotator's
``is_ready()`` is ``cache_status().state is CacheState.READY``, and the
CLI names the state of each source it skips. It opens the file with a
``mode=ro`` URI and runs only SELECT and PRAGMA, so checking a cache
writes nothing to the data directory.
"""

from __future__ import annotations

import contextlib
import enum
import os
import sqlite3
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from pathlib import Path

_SQLITE_MAGIC = b"SQLite format 3\x00"
# Header bytes 18 and 19 (write and read format): 1 rollback journal, 2 WAL.
_WAL_FORMAT = 2


class CacheState(enum.Enum):
    """Why a cache can or cannot be used."""

    READY = "ready"
    MISSING = "missing"
    OUTDATED = "outdated"
    UNREADABLE = "unreadable"
    CORRUPT = "corrupt"
    WAL = "wal"
    READONLY_WAL = "readonly_wal"


@dataclass(frozen=True)
class CacheStatus:
    """A cache's state, the file it concerns, and a detail for the console."""

    state: CacheState
    path: Path | None = None
    detail: str = ""

    @property
    def ready(self) -> bool:
        """True when the cache can be used."""
        return self.state is CacheState.READY


READY = CacheStatus(CacheState.READY)


def connect_readonly(path: Path) -> sqlite3.Connection:
    """Open ``path`` read-only (``mode=ro``); never creates or changes the file.

    A ``mode=ro`` open of a WAL file creates ``-shm`` and ``-wal`` beside
    it and leaves them, so callers check ``is_wal`` first and never open
    one (ADR-0048).
    """
    return sqlite3.connect(f"{path.absolute().as_uri()}?mode=ro", uri=True)


def is_wal(path: Path) -> bool:
    """True when the file's header says WAL journal mode (bytes 18 and 19; no SQLite open)."""
    try:
        with path.open("rb") as fh:
            header = fh.read(20)
    except OSError:
        return False
    return len(header) == 20 and _WAL_FORMAT in (header[18], header[19])


def cache_state(
    path: Path,
    name: str | None,
    *,
    tag: str | None = None,
    tables: Iterable[str] = (),
    columns: Mapping[str, Iterable[str]] | None = None,
) -> CacheStatus:
    """Return the state of the cache at ``path``, record ``name`` in ``database_versions``.

    With ``name`` None only the file and ``tables`` are checked (a cache
    whose readiness the caller decides from its rows, as SNPedia's).

    ``tag`` is the ``local_version_tag`` the code needs (None: any).
    ``tables`` must exist besides ``database_versions``; the first is the
    main table, which must hold a row. ``columns`` maps a table to the
    columns it must have. A missing table or an empty main table is
    ``corrupt``; a missing column, version row or tag is ``outdated``.
    """
    tables = tuple(tables)
    if not path.exists():
        return CacheStatus(CacheState.MISSING, path)
    if not os.access(path, os.R_OK):
        return CacheStatus(CacheState.UNREADABLE, path)
    try:
        with path.open("rb") as fh:
            header = fh.read(len(_SQLITE_MAGIC))
    except OSError as exc:
        return CacheStatus(CacheState.UNREADABLE, path, str(exc))
    if header != _SQLITE_MAGIC:
        return CacheStatus(CacheState.CORRUPT, path, "not an SQLite database")
    if is_wal(path):
        # Opening it read-only would leave -shm and -wal beside it, so a
        # WAL cache is skipped until the update step converts it.
        if not os.access(path.parent, os.W_OK | os.X_OK):
            return CacheStatus(CacheState.READONLY_WAL, path)
        return CacheStatus(CacheState.WAL, path)
    try:
        with contextlib.closing(connect_readonly(path)) as conn:
            present = {
                row[0]
                for row in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
            }
            needed_tables = {*tables} if name is None else {"database_versions", *tables}
            absent = sorted(needed_tables - present)
            if absent:
                return CacheStatus(CacheState.CORRUPT, path, f"missing table {', '.join(absent)}")
            if tables and conn.execute(f"SELECT 1 FROM {tables[0]} LIMIT 1").fetchone() is None:
                return CacheStatus(CacheState.CORRUPT, path, f"no rows in {tables[0]}")
            for table, needed in (columns or {}).items():
                have = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
                lacking = sorted(set(needed) - have)
                if lacking:
                    return CacheStatus(
                        CacheState.OUTDATED, path, f"{table} lacks {', '.join(lacking)}"
                    )
            if name is None:
                return READY
            if tag is None:
                row = conn.execute(
                    "SELECT 1 FROM database_versions WHERE name = ?", (name,)
                ).fetchone()
                return READY if row else CacheStatus(CacheState.OUTDATED, path, "no version row")
            have_tag_column = any(
                row[1] == "local_version_tag"
                for row in conn.execute("PRAGMA table_info(database_versions)")
            )
            if not have_tag_column:
                return CacheStatus(CacheState.OUTDATED, path, f"no version tag, needs {tag}")
            row = conn.execute(
                "SELECT local_version_tag FROM database_versions WHERE name = ?", (name,)
            ).fetchone()
    except sqlite3.DatabaseError as exc:
        text = str(exc)
        if "not a database" in text or "malformed" in text:
            return CacheStatus(CacheState.CORRUPT, path, text)
        return CacheStatus(CacheState.UNREADABLE, path, text)
    if row is None:
        return CacheStatus(CacheState.OUTDATED, path, "no version row")
    if row[0] != tag:
        built = row[0] or "no version tag"
        return CacheStatus(CacheState.OUTDATED, path, f"built for {built}, needs {tag}")
    return READY

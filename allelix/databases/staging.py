# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Staging paths that belong to one operation.

A cache writer stages its bytes beside the target and ``os.replace``s them
onto it once they are complete. A staging name derived from the target alone
(``gwas.sqlite.tmp``) is shared by every writer of that target: one writer's
start-up or failure cleanup could delete another's bytes, and its rename
could publish them. Each operation therefore gets its own name, and removes
only that name, on success, failure or interruption: SIGTERM and SIGHUP end
a command through that cleanup too (2.6.1). A file left behind by a process
killed outright (SIGKILL, a crash) is never removed by another run, since it
cannot be told from a live operation's; `db clean` removes it. On Windows a
forced end (TerminateProcess, SIGTERM there), closing the console and
Ctrl-Break do not run the SIGTERM/SIGHUP cleanup either.
"""

from __future__ import annotations

import contextlib
import logging
import os
import re
import secrets
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

_ATTEMPTS = 100

# An ``owned_staging`` name with a suffix its callers give it, or a side file
# SQLite leaves beside a staging database.
_STAGING_FILE = re.compile(
    r"\.(?P<target>.+)\.[0-9a-f]{16}(?:\.tmp|\.part|\.partial)(?:-wal|-shm|-journal)?"
)
# ``tempfile`` names a directory ``<prefix><8 characters>`` from this alphabet.
_STAGING_DIRECTORY = re.compile(r"\.(?P<prefix>[a-z_]+)-[a-z0-9_]{8}")


def staging_file_target(name: str) -> str | None:
    """The target name of an ``owned_staging`` file name (or of its side file), else None."""
    match = _STAGING_FILE.fullmatch(name)
    return match["target"] if match else None


def staging_directory_prefix(name: str) -> str | None:
    """The prefix of an ``owned_directory`` name, else None."""
    match = _STAGING_DIRECTORY.fullmatch(name)
    return match["prefix"] if match else None


@contextlib.contextmanager
def owned_staging(target: Path, suffix: str = ".tmp") -> Iterator[Path]:
    """A new, empty file beside ``target`` that only this operation uses.

    The name is hidden and unique (``.<target name>.<random><suffix>``). The
    file is created with the default mode the process's umask gives, not
    ``mkstemp``'s 0600, so a cache published from it can be read as any other
    file the user writes there (a shared read-only data directory). It is
    removed on exit unless it was already moved onto ``target``.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    for _ in range(_ATTEMPTS):
        path = target.parent / f".{target.name}.{secrets.token_hex(8)}{suffix}"
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
        except FileExistsError:
            continue
        os.close(fd)
        break
    else:
        raise FileExistsError(f"no free staging name beside {target}")
    try:
        yield path
    finally:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            logger.warning("Could not remove staging file %s", path)


@contextlib.contextmanager
def owned_directory(data_dir: Path, prefix: str) -> Iterator[Path]:
    """A new directory in ``data_dir`` for one operation's downloads, removed on exit.

    On the cache's own filesystem, so what is staged there is published by
    rename.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{prefix}-", dir=data_dir) as name:
        yield Path(name)

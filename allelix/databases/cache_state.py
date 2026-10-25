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
import ctypes
import enum
import logging
import os
import sqlite3
import sys
from ctypes import wintypes
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from pathlib import Path

_SQLITE_MAGIC = b"SQLite format 3\x00"
# Header bytes 18 and 19 (write and read format): 1 rollback journal, 2 WAL.
_WAL_FORMAT = 2
_ERROR_INSUFFICIENT_BUFFER = 122
_ERROR_NO_TOKEN = 1008
_FILE_ADD_FILE = 2
_TOKEN_QUERY_DUPLICATE = 10
_SECURITY_IMPERSONATION = 2
_OWNER_GROUP_DACL = 7
_BUFFER_RETRIES = 3
logger = logging.getLogger(__name__)


class _FileGenericMapping(ctypes.Structure):
    _fields_ = [
        (name, wintypes.DWORD)
        for name in (
            "GenericRead",
            "GenericWrite",
            "GenericExecute",
            "GenericAll",
        )
    ]


@lru_cache(maxsize=1)
def _windows_security_api() -> tuple[ctypes.CDLL, ctypes.CDLL]:
    """Bind the narrow Windows descriptor/token APIs, preserving last-error values."""
    load = getattr(ctypes, "WinDLL", None)
    if load is None:
        raise OSError("Windows security APIs are unavailable")
    kernel = load("kernel32", use_last_error=True)
    security = load("advapi32", use_last_error=True)
    kernel.GetCurrentProcess.argtypes = []
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    kernel.GetCurrentThread.argtypes = []
    kernel.GetCurrentThread.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel.CloseHandle.restype = wintypes.BOOL
    security.GetFileSecurityW.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD),
    ]
    security.GetFileSecurityW.restype = wintypes.BOOL
    security.OpenThreadToken.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        wintypes.BOOL,
        ctypes.POINTER(wintypes.HANDLE),
    ]
    security.OpenThreadToken.restype = wintypes.BOOL
    security.OpenProcessToken.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        ctypes.POINTER(wintypes.HANDLE),
    ]
    security.OpenProcessToken.restype = wintypes.BOOL
    security.DuplicateToken.argtypes = [
        wintypes.HANDLE,
        ctypes.c_int,
        ctypes.POINTER(wintypes.HANDLE),
    ]
    security.DuplicateToken.restype = wintypes.BOOL
    security.AccessCheck.argtypes = [
        wintypes.LPVOID,
        wintypes.HANDLE,
        wintypes.DWORD,
        ctypes.POINTER(_FileGenericMapping),
        wintypes.LPVOID,
        ctypes.POINTER(wintypes.DWORD),
        ctypes.POINTER(wintypes.DWORD),
        ctypes.POINTER(wintypes.BOOL),
    ]
    security.AccessCheck.restype = wintypes.BOOL
    return kernel, security


def _windows_last_error() -> int:
    get_error = getattr(ctypes, "get_last_error", None)
    if get_error is None:
        raise OSError("Windows last-error API is unavailable")
    return int(get_error())


def _windows_error(operation: str, code: int | None = None) -> OSError:
    code = _windows_last_error() if code is None else code
    return OSError(code, f"{operation} failed (Windows error {code})")


def _file_security_descriptor(security: ctypes.CDLL, path: Path) -> ctypes.Array[ctypes.c_char]:
    needed = wintypes.DWORD()
    initial = security.GetFileSecurityW(
        str(path), _OWNER_GROUP_DACL, None, 0, ctypes.byref(needed)
    )
    error = _windows_last_error()
    if initial or error != _ERROR_INSUFFICIENT_BUFFER or not needed.value:
        raise _windows_error("GetFileSecurityW size", error)
    size = needed.value
    for _attempt in range(_BUFFER_RETRIES):
        descriptor = ctypes.create_string_buffer(size)
        needed = wintypes.DWORD()
        if security.GetFileSecurityW(
            str(path), _OWNER_GROUP_DACL, descriptor, size, ctypes.byref(needed)
        ):
            if needed.value > size:
                raise OSError("GetFileSecurityW returned an incomplete descriptor")
            return descriptor
        error = _windows_last_error()
        if error != _ERROR_INSUFFICIENT_BUFFER or needed.value <= size:
            raise _windows_error("GetFileSecurityW", error)
        size = needed.value
    raise OSError("GetFileSecurityW descriptor size kept changing")


def _windows_directory_writable(path: Path) -> bool:
    """Check FILE_ADD_FILE for the effective token without writing or adjusting privileges."""
    kernel, security = _windows_security_api()
    descriptor = _file_security_descriptor(security, path.resolve())
    acquired: list[wintypes.HANDLE] = []
    try:
        token = wintypes.HANDLE()
        if not security.OpenThreadToken(
            kernel.GetCurrentThread(), _TOKEN_QUERY_DUPLICATE, True, ctypes.byref(token)
        ):
            error = _windows_last_error()
            if error != _ERROR_NO_TOKEN:
                raise _windows_error("OpenThreadToken", error)
            if not security.OpenProcessToken(
                kernel.GetCurrentProcess(), _TOKEN_QUERY_DUPLICATE, ctypes.byref(token)
            ):
                raise _windows_error("OpenProcessToken")
        acquired.append(token)
        duplicate = wintypes.HANDLE()
        if not security.DuplicateToken(token, _SECURITY_IMPERSONATION, ctypes.byref(duplicate)):
            raise _windows_error("DuplicateToken")
        acquired.append(duplicate)
        mapping = _FileGenericMapping(0x120089, 0x120116, 0x1200A0, 0x1F01FF)
        privileges = ctypes.create_string_buffer(128)
        for _attempt in range(_BUFFER_RETRIES):
            size = wintypes.DWORD(ctypes.sizeof(privileges))
            granted, allowed = wintypes.DWORD(), wintypes.BOOL()
            if security.AccessCheck(
                descriptor,
                duplicate,
                _FILE_ADD_FILE,
                ctypes.byref(mapping),
                privileges,
                ctypes.byref(size),
                ctypes.byref(granted),
                ctypes.byref(allowed),
            ):
                return bool(allowed.value)
            error = _windows_last_error()
            if error != _ERROR_INSUFFICIENT_BUFFER or size.value <= ctypes.sizeof(privileges):
                raise _windows_error("AccessCheck", error)
            privileges = ctypes.create_string_buffer(size.value)
        raise OSError("AccessCheck privilege buffer size kept changing")
    finally:
        primary_error = sys.exc_info()[0] is not None
        cleanup_error: OSError | None = None
        for handle in reversed(acquired):
            try:
                if not kernel.CloseHandle(handle):
                    close_error = _windows_error("CloseHandle")
                    if cleanup_error is None:
                        cleanup_error = close_error
            except OSError as exc:
                if cleanup_error is None:
                    cleanup_error = exc
        if cleanup_error is not None and not primary_error:
            raise cleanup_error


def directory_writable(path: Path, *, require_search: bool = True) -> bool:
    """Check directory-create permission; native uncertainty is conservatively unavailable.

    Preserve the callers' POSIX access masks. Windows metadata access checks
    ignore ACL denials, so use the effective token and the physical directory's
    descriptor there. This is a snapshot, not a guarantee of later I/O.
    """
    if os.name != "nt":
        return os.access(path, os.W_OK | (os.X_OK if require_search else 0))
    try:
        return _windows_directory_writable(path)
    except (OSError, RuntimeError) as exc:
        logger.warning("Cannot check write access to %s: %s", path, exc)
        return False


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
    # Windows SQLite can otherwise find the main file through a symlink but
    # look for its WAL/SHM beside the link, silently missing committed rows.
    return sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)


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
        try:
            parent = path.resolve().parent
        except (OSError, RuntimeError) as exc:
            return CacheStatus(CacheState.UNREADABLE, path, str(exc))
        if not directory_writable(parent):
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

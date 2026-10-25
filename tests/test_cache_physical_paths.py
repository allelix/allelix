# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""SQLite and its sidecars must share the installed physical database path."""

from __future__ import annotations

import contextlib
import hashlib
import os
import shutil
import sqlite3
from pathlib import Path

import pytest

from allelix.databases import cache_state, loader_utils


def _live_wal(tmp_path: Path) -> tuple[Path, Path]:
    build, real, aliases = (tmp_path / name for name in ("build", "real files", "aliases"))
    for directory in (build, real, aliases):
        directory.mkdir()
    filename = "cache #50%.sqlite"
    source = build / filename
    with contextlib.closing(sqlite3.connect(source)) as conn:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA wal_autocheckpoint=0")
        conn.execute("CREATE TABLE t (x INTEGER PRIMARY KEY)")
        conn.execute("INSERT INTO t VALUES (0)")
        conn.commit()
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        conn.executemany("INSERT INTO t VALUES (?)", [(i,) for i in range(1, 100)])
        conn.commit()
        for suffix in ("", "-wal"):
            shutil.copyfile(build / (filename + suffix), real / (filename + suffix))
    target = real / filename
    link = aliases / filename
    link.symlink_to(os.path.relpath(target, aliases))
    assert cache_state.is_wal(target)
    assert (real / (filename + "-wal")).stat().st_size > 0
    return target, link


def _rows(path: Path) -> list[int]:
    with contextlib.closing(cache_state.connect_readonly(path)) as conn:
        return [row[0] for row in conn.execute("SELECT x FROM t ORDER BY x")]


@pytest.mark.parametrize("relative", [False, True])
def test_readonly_link_uses_physical_uri_and_every_committed_row(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative: bool,
) -> None:
    target, link = _live_wal(tmp_path)
    connect = sqlite3.connect
    calls = []

    def capture(database, *args, **kwargs):
        calls.append((database, kwargs))
        return connect(database, *args, **kwargs)

    monkeypatch.setattr(cache_state.sqlite3, "connect", capture)
    if relative:
        monkeypatch.chdir(tmp_path)
        link = link.relative_to(tmp_path)
    assert _rows(link) == list(range(100))
    assert calls[0] == (target.resolve().as_uri() + "?mode=ro", {"uri": True})
    assert "%23" in calls[0][0] and "%25" in calls[0][0] and "%20" in calls[0][0]
    assert "immutable" not in calls[0][0]
    with (
        contextlib.closing(cache_state.connect_readonly(link)) as conn,
        pytest.raises(sqlite3.OperationalError, match="readonly"),
    ):
        conn.execute("INSERT INTO t VALUES (999)")
    assert _rows(target) == list(range(100))


@pytest.mark.parametrize("relative", [False, True])
def test_conversion_uses_physical_database_and_sidecars_without_losing_wal_rows(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relative: bool,
) -> None:
    target, link = _live_wal(tmp_path)
    connect = sqlite3.connect
    calls = []
    cleaned = []
    remove = loader_utils.remove_side_files

    def capture(database, *args, **kwargs):
        calls.append(database)
        return connect(database, *args, **kwargs)

    def cleanup(database):
        assert not cache_state.is_wal(database), "Never clean live WAL sidecars"
        cleaned.append(database)
        remove(database)

    monkeypatch.setattr(loader_utils.sqlite3, "connect", capture)
    monkeypatch.setattr(loader_utils, "remove_side_files", cleanup)
    if relative:
        monkeypatch.chdir(tmp_path)
        link = link.relative_to(tmp_path)
    assert loader_utils.convert_to_rollback_journal(link)
    assert calls == [target.resolve()]
    assert cleaned == [target.resolve()]
    assert not cache_state.is_wal(target)
    assert sorted(p.name for p in target.parent.iterdir()) == [target.name]
    assert _rows(target) == _rows(link) == list(range(100))
    before = hashlib.sha256(target.read_bytes()).hexdigest()
    assert not loader_utils.convert_to_rollback_journal(link)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == before


def test_busy_conversion_never_cleans_sidecars_or_reports_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target, link = _live_wal(tmp_path)
    with contextlib.closing(sqlite3.connect(target)) as reader:
        reader.execute("BEGIN")
        assert reader.execute("SELECT COUNT(*) FROM t").fetchone()[0] == 100
        monkeypatch.setattr(
            loader_utils, "remove_side_files", lambda _p: pytest.fail("busy cleanup")
        )
        assert not loader_utils.convert_to_rollback_journal(link)
        assert cache_state.is_wal(target)
        assert reader.execute("SELECT COUNT(*) FROM t").fetchone()[0] == 100
    assert _rows(target) == list(range(100))


def test_unreadable_header_never_opens_sqlite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target, link = _live_wal(tmp_path)
    original = Path.open

    def denied(path, *args, **kwargs):
        if path.resolve() == target.resolve():
            raise PermissionError("synthetic header denial")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", denied)
    monkeypatch.setattr(
        loader_utils.sqlite3, "connect", lambda *_a, **_kw: pytest.fail("SQLite open")
    )
    assert not loader_utils.convert_to_rollback_journal(link)


def test_missing_and_non_database_paths_are_noops(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bad = tmp_path / "not.sqlite"
    bad.write_bytes(b"not an SQLite file")
    monkeypatch.setattr(
        loader_utils.sqlite3, "connect", lambda *_a, **_kw: pytest.fail("SQLite open")
    )
    assert not loader_utils.convert_to_rollback_journal(bad)
    assert not loader_utils.convert_to_rollback_journal(tmp_path / "missing.sqlite")


@pytest.mark.parametrize("cycle", [False, True], ids=["dangling", "cycle"])
def test_unresolvable_links_preserve_conversion_noop(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cycle: bool,
) -> None:
    first, second = tmp_path / "a.sqlite", tmp_path / "b.sqlite"
    first.symlink_to(second.name)
    if cycle:
        second.symlink_to(first.name)
    monkeypatch.setattr(
        loader_utils.sqlite3, "connect", lambda *_a, **_kw: pytest.fail("SQLite open")
    )
    monkeypatch.setattr(loader_utils, "remove_side_files", lambda _p: pytest.fail("link cleanup"))
    assert not loader_utils.convert_to_rollback_journal(first)


def test_successful_conversion_does_not_swallow_cleanup_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _target, link = _live_wal(tmp_path)

    def denied(_path):
        raise PermissionError("synthetic cleanup denial")

    monkeypatch.setattr(loader_utils, "remove_side_files", denied)
    with pytest.raises(PermissionError, match="synthetic cleanup denial"):
        loader_utils.convert_to_rollback_journal(link)


@pytest.mark.parametrize("deny_target", [False, True], ids=["denied-alias", "denied-target"])
def test_wal_readiness_checks_the_physical_parent_without_opening_sqlite(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    deny_target: bool,
) -> None:
    from test_readonly_caches import _read_only

    target, link = _live_wal(tmp_path)
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in target.parent.iterdir()}
    monkeypatch.setattr(
        cache_state.sqlite3, "connect", lambda *_a, **_kw: pytest.fail("SQLite open")
    )
    directory = target.parent if deny_target else link.parent
    with _read_only(directory):
        status = cache_state.cache_state(link, None)
    expected = cache_state.CacheState.READONLY_WAL if deny_target else cache_state.CacheState.WAL
    assert status.state is expected
    assert status.path == link
    assert {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in target.parent.iterdir()
    } == before


def test_conversion_abstains_when_directory_check_denies(tmp_path, monkeypatch):
    target, link = _live_wal(tmp_path)
    monkeypatch.setattr(cache_state, "directory_writable", lambda *_a, **_kw: False)
    assert not loader_utils.convert_to_rollback_journal(link)
    assert cache_state.is_wal(target)

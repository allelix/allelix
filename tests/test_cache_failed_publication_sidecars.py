"""A failed cache replacement must retain the old database's recoverable state."""

import contextlib
import shutil
import sqlite3

import pytest

from allelix.databases import loader_utils


def test_failed_replace_keeps_the_hot_journal_needed_by_the_old_cache(tmp_path, monkeypatch):
    old = tmp_path / "cache.sqlite"
    journal = old.with_name(old.name + "-journal")
    with contextlib.closing(sqlite3.connect(old)) as conn:
        conn.execute("CREATE TABLE records (id INTEGER PRIMARY KEY, value TEXT)")
        conn.executemany(
            "INSERT INTO records VALUES (?,?)", [(i, "old" * 100) for i in range(1000)]
        )
        conn.commit()
        conn.execute("PRAGMA cache_size=5")
        conn.execute("UPDATE records SET value=?", ("uncommitted" * 100,))
        assert journal.exists()
        damaged_bytes = old.read_bytes()
        recovery_bytes = journal.read_bytes()
        conn.rollback()
    old.write_bytes(damaged_bytes)
    journal.write_bytes(recovery_bytes)
    # Prove the on-disk pair really can restore the committed old contents.
    control = tmp_path / "control.sqlite"
    shutil.copy2(old, control)
    shutil.copy2(journal, control.with_name(control.name + "-journal"))
    with contextlib.closing(sqlite3.connect(control)) as conn:
        assert conn.execute("SELECT DISTINCT value FROM records").fetchall() == [("old" * 100,)]

    staged = tmp_path / "new.sqlite"
    with contextlib.closing(sqlite3.connect(staged)) as conn:
        conn.execute("CREATE TABLE new_records (id INTEGER)")
        conn.commit()

    def refused(*args, **kwargs):
        raise PermissionError("injected replacement failure")

    monkeypatch.setattr(loader_utils.os, "replace", refused)
    with pytest.raises(PermissionError, match="injected replacement failure"):
        loader_utils.publish_cache(staged, old)

    # A failure is allowed, but not loss of the old cache's committed rows.
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT DISTINCT value FROM records").fetchall() == [("old" * 100,)]

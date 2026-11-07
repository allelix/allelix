import contextlib
import os
import sqlite3

from allelix.databases import loader_utils


def test_reader_cannot_replay_old_journal_into_new_cache(tmp_path, monkeypatch):
    old = tmp_path / "cache.sqlite"
    journal = tmp_path / "cache.sqlite-journal"
    with contextlib.closing(sqlite3.connect(old)) as conn:
        conn.execute("CREATE TABLE records (value TEXT)")
        conn.executemany("INSERT INTO records VALUES (?)", [("old" * 100,)] * 1000)
        conn.commit()
        conn.execute("PRAGMA cache_size=5")
        conn.execute("UPDATE records SET value=?", ("uncommitted" * 100,))
        old_bytes, journal_bytes = old.read_bytes(), journal.read_bytes()
        conn.rollback()
    old.write_bytes(old_bytes)
    journal.write_bytes(journal_bytes)
    new = tmp_path / "new.sqlite"
    with contextlib.closing(sqlite3.connect(new)) as conn:
        conn.execute("CREATE TABLE records (value TEXT)")
        conn.execute("INSERT INTO records VALUES (?)", ("new",))
        conn.commit()
    original_connect = sqlite3.connect
    original_replace = os.replace
    recovery_denied = []
    after_replace = []

    def connect(path, *args, **kwargs):
        if str(path) == str(old) and not after_replace:
            recovery_denied.append(True)
            raise sqlite3.OperationalError("database is locked")
        return original_connect(path, *args, **kwargs)

    def publish_then_read(source, target):
        original_replace(source, target)
        after_replace.append(True)
        # A reader can run after the atomic rename and before sidecar cleanup.
        with contextlib.closing(original_connect(old)) as conn:
            conn.execute("SELECT COUNT(*) FROM records").fetchone()

    monkeypatch.setattr(loader_utils.sqlite3, "connect", connect)
    monkeypatch.setattr(loader_utils.os, "replace", publish_then_read)
    try:
        loader_utils.publish_cache(new, old)
    except (OSError, sqlite3.DatabaseError, RuntimeError):
        assert not after_replace, (
            "failure must precede exposing a new cache beside a hot old journal"
        )
        assert journal.read_bytes() == journal_bytes
        return
    assert recovery_denied
    with contextlib.closing(original_connect(old)) as conn:
        assert conn.execute("SELECT DISTINCT value FROM records").fetchall() == [("new",)]

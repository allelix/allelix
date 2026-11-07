"""Old WAL pages must never be applied to a replacement cache."""

import contextlib
import os
import sqlite3

from allelix.databases import loader_utils


def test_old_wal_cannot_be_replayed_into_new_cache(tmp_path, monkeypatch):
    old = tmp_path / "cache.sqlite"
    with contextlib.closing(sqlite3.connect(old)) as c:
        c.execute("PRAGMA journal_mode=WAL")
        c.execute("CREATE TABLE records (value TEXT)")
        c.commit()
        c.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        c.execute("INSERT INTO records VALUES ('old')")
        c.commit()
        saved = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    for name, content in saved.items():
        (tmp_path / name).write_bytes(content)
    new = tmp_path / "new.sqlite"
    with contextlib.closing(sqlite3.connect(new)) as c:
        c.execute("CREATE TABLE records (value TEXT)")
        c.execute("INSERT INTO records VALUES ('new')")
        c.commit()
    real_replace = os.replace
    published = []

    def replace(source, target):
        real_replace(source, target)
        published.append(True)
        with contextlib.closing(sqlite3.connect(old)) as c:
            c.execute("SELECT value FROM records").fetchall()

    monkeypatch.setattr(loader_utils.os, "replace", replace)
    try:
        loader_utils.publish_cache(new, old)
    except (OSError, RuntimeError, sqlite3.DatabaseError):
        assert not published
        with contextlib.closing(sqlite3.connect(old)) as c:
            assert c.execute("SELECT value FROM records").fetchall() == [("old",)]
        return
    with contextlib.closing(sqlite3.connect(old)) as c:
        assert c.execute("SELECT value FROM records").fetchall() == [("new",)]

"""RC2 cleanup may remove abandoned bytes, never a live journal or user TSV.

Synthetic and offline. The small SQLite transaction creates a real zero-header
journal before its commit; that is not proof the journal is abandoned.
"""

import contextlib
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.databases.loader_utils import UnrecoveredJournalError, publish_cache


def _database(path, value):
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.execute("CREATE TABLE t (v TEXT)")
        conn.execute("INSERT INTO t VALUES (?)", (value,))
        conn.commit()


@pytest.mark.parametrize("transaction", ["IMMEDIATE", "EXCLUSIVE"])
def test_live_zero_header_journal_is_not_unlinked(tmp_path, transaction):
    old, staged = tmp_path / "cache.sqlite", tmp_path / "new.sqlite"
    _database(old, "old")
    _database(staged, "new")
    journal = old.with_name(old.name + "-journal")
    with contextlib.closing(sqlite3.connect(old)) as writer:
        writer.execute(f"BEGIN {transaction}")
        writer.execute("UPDATE t SET v='pending'")
        before = journal.read_bytes()
        assert before[:8] == bytes(8), "fixture must have a real non-hot journal"
        old_bytes = old.read_bytes()
        try:
            with pytest.raises(UnrecoveredJournalError):
                publish_cache(staged, old)
            assert old.read_bytes() == old_bytes
            assert journal.read_bytes() == before
            assert staged.exists()
        finally:
            writer.rollback()
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT v FROM t").fetchone() == ("old",)


@pytest.mark.parametrize("contents", [b"", bytes(512)])
def test_abandoned_non_hot_journal_still_allows_publication(tmp_path, contents):
    old, staged = tmp_path / "cache.sqlite", tmp_path / "new.sqlite"
    _database(old, "old")
    _database(staged, "new")
    journal = old.with_name(old.name + "-journal")
    journal.write_bytes(contents)
    publish_cache(staged, old)
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT v FROM t").fetchone() == ("new",)
    assert not journal.exists()


@pytest.mark.parametrize("prefix", ["gwas", "clinvar", "snpedia"])
@pytest.mark.parametrize("force", [False, True])
def test_user_tsv_in_extracted_backup_is_kept(tmp_path, prefix, force):
    backup = tmp_path / f".{prefix}-20261101" / "extracted"
    backup.mkdir(parents=True)
    user = backup / "family-genotypes.tsv"
    user.write_text("synthetic user backup\n")
    (tmp_path / "gnomad.sqlite").write_bytes(b"owned cache")
    args = ["db", "clean", "--data-dir", str(tmp_path), "--yes"]
    if force:
        args.append("--force")
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert user.exists(), "a .tsv extension is not proof of an Allelix extraction"
    assert user.read_text() == "synthetic user backup\n"
    assert not (tmp_path / "gnomad.sqlite").exists()

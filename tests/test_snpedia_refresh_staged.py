# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A SNPedia refresh parses its own staged copy and publishes only a parsed cache.

The four ``test_r02_*`` tests and their fixtures are reviewer C's (outside
review of 2.6.0 RC4, finding R4-02), unchanged but for line wrapping and
one unused name (``_downloads``): an interrupted, failed,
concurrent or killed refresh must not leave SNPedia unusable or report
success. The tests after them cover the update step's re-parse of an
outdated cache, and the stale journal beside any installed cache.
Synthetic pages only.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import shutil
import sqlite3
import urllib.request
from typing import TYPE_CHECKING

import mwparserfromhell
import pytest
from click.testing import CliRunner

from allelix.annotators import snpedia as snpedia_mod
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.cli import db as db_cli
from allelix.cli import main
from allelix.databases import snpedia_parser
from allelix.databases.cache_state import CacheState
from allelix.databases.loader_utils import install_prebuilt_gz_cache

if TYPE_CHECKING:
    from pathlib import Path

_REAL_URLOPEN = urllib.request.urlopen


@pytest.fixture(autouse=True)
def _offline(monkeypatch):
    """No test may reach the network."""

    def refuse(*args, **kwargs):
        raise AssertionError("network use in an offline regression test")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)


_SNP_PAGE = "{{Rsnum|rsid=1801133|Gene=MTHFR|StabilizedOrientation=plus}}"


def _snpedia_raw_gz(folder: Path, pages: int = 3000) -> bytes:
    """The published SNPedia asset's shape: raw wiki pages, parsed after install.

    Large enough (~26 MB parsed) that a parse's dirty pages exceed SQLite's
    default page cache, as the real ~100k-page parse does: the parse then
    writes its journal and pages to the published file before it commits.
    """
    raw = folder / "snpedia-raw.sqlite"
    raw.unlink(missing_ok=True)
    with contextlib.closing(sqlite3.connect(raw)) as conn:
        conn.execute(
            "CREATE TABLE pages (title TEXT PRIMARY KEY, category TEXT, content TEXT, "
            "scraped_at TEXT)"
        )
        conn.execute(
            "INSERT INTO pages VALUES ('Rs1801133', 'snp', ?, '2026-01-01')", (_SNP_PAGE,)
        )
        for i in range(pages):
            page = (
                f"{{{{Genotype|rsid={i}|allele1=C|allele2=T|magnitude=2.2|repute=Bad|"
                f"summary={'s' * 1500}}}}}"
            )
            conn.execute(
                "INSERT INTO pages VALUES (?, 'genotype', ?, '2026-01-01')", (f"Rs{i}(C;T)", page)
            )
        conn.commit()
    return gzip.compress(raw.read_bytes())


@pytest.fixture
def snpedia_only(tmp_path, monkeypatch):
    """`db update` managing SNPedia alone, with a fake verified download."""
    payload = _snpedia_raw_gz(tmp_path)
    downloads: list[str] = []

    def fake_download(url, dest):
        downloads.append(url)
        dest.write_bytes(payload)

    monkeypatch.setattr(snpedia_mod, "download", fake_download)
    monkeypatch.setattr(snpedia_mod, "verify_file_hash", lambda *a, **k: None)
    monkeypatch.setattr(db_cli, "get_annotators", lambda folder, **kw: [SNPediaAnnotator(folder)])
    data = tmp_path / "data"
    data.mkdir()
    first = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
    assert first.exit_code == 0, first.output
    assert SNPediaAnnotator(data).is_ready()
    return data, downloads


def test_r02_interrupted_snpedia_refresh_keeps_the_cache_you_had(snpedia_only, monkeypatch):
    """README: "a failed or interrupted update keeps the cache you had". SNPedia's
    refresh replaces the parsed cache with raw pages before its multi-minute parse,
    so a Ctrl-C during the parse leaves no usable SNPedia."""
    data, _ = snpedia_only
    rows = SNPediaAnnotator(data).record_count()
    real = mwparserfromhell.parse
    calls = [0]

    def interrupted(text, *args, **kwargs):
        calls[0] += 1
        if calls[0] == 50:
            raise KeyboardInterrupt
        return real(text, *args, **kwargs)

    monkeypatch.setattr(snpedia_parser.mwparserfromhell, "parse", interrupted)
    CliRunner().invoke(main, ["db", "update", "--data-dir", str(data), "--force"])
    assert calls[0] >= 50, "the refresh must have reached its parse"
    after = SNPediaAnnotator(data)
    assert after.is_ready(), f"interrupted refresh lost the parsed cache: {after.cache_status()}"
    assert after.record_count() == rows


def test_r02_failed_snpedia_parse_is_not_reported_as_success(snpedia_only, monkeypatch):
    """A parse that fails with an SQLite error (disk full) is swallowed by
    prepare_cache; db update prints "✓ snpedia refreshed" and exits 0 although
    the cache it replaced is gone and the new one is not usable."""
    data, _ = snpedia_only

    def disk_full(*args, **kwargs):
        raise sqlite3.OperationalError("database or disk is full")

    monkeypatch.setattr(snpedia_parser, "parse_raw_pages", disk_full)
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data), "--force"])
    ready = SNPediaAnnotator(data).is_ready()
    assert ready or result.exit_code != 0, (
        f"db update reported success with an unusable SNPedia cache:\n{result.output}"
    )


@pytest.mark.skipif(not hasattr(os, "fork"), reason="a second process is started with fork")
def test_r02_concurrent_snpedia_updates_leave_a_parsed_cache(snpedia_only, monkeypatch, tmp_path):
    """CHANGELOG 2.6.0 / README: two updates running at once cannot install each
    other's files. Two SNPedia refreshes that both report success must leave a
    usable cache. The second process runs to completion at a fixed point of the
    first one's parse (no timing): on RC4 its raw install lands under the first
    one's in-place parse, both print "✓ snpedia refreshed", and SNPedia is left
    unparsed."""
    data, _ = snpedia_only
    real = mwparserfromhell.parse
    state = {"calls": 0, "child": None}
    report = tmp_path / "second-update.txt"

    def parse(text, *args, **kwargs):
        state["calls"] += 1
        # After two 1,000-row insert batches the parse has outgrown SQLite's page
        # cache, so its journal and pages are on disk, as in any real parse.
        if state["calls"] == 2500 and state["child"] is None:
            pid = os.fork()
            if pid == 0:  # the other `allelix db update --force`, in its own process
                code = 1
                try:
                    result = CliRunner().invoke(
                        main, ["db", "update", "--data-dir", str(data), "--force"]
                    )
                    report.write_text(result.output)
                    code = result.exit_code
                finally:
                    os._exit(code)
            _, status = os.waitpid(pid, 0)
            state["child"] = os.waitstatus_to_exitcode(status)
        return real(text, *args, **kwargs)

    monkeypatch.setattr(snpedia_parser.mwparserfromhell, "parse", parse)
    first = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data), "--force"])
    assert state["child"] is not None, (
        "the second update must have run during the first one's parse"
    )
    both_succeeded = first.exit_code == 0 and state["child"] == 0
    after = SNPediaAnnotator(data)
    assert after.is_ready() or not both_succeeded, (
        f"two refreshes reported success and left SNPedia {after.cache_status()}:\n"
        f"first:\n{first.output}\nsecond:\n{report.read_text() if report.exists() else ''}"
    )


def test_r02_killed_snpedia_parse_does_not_trap_db_update(snpedia_only, tmp_path):
    """A parse killed without Python cleanup (SIGTERM, SIGKILL, closed terminal)
    leaves a hot journal beside the published snpedia.sqlite. RC4 never recovers
    it: SNPedia reads as "cannot read" and every db update re-downloads it,
    prints "✓ snpedia ready" and exits 0 while it stays unusable."""
    data, _downloads = snpedia_only
    # Reproduce the on-disk state of a killed in-place parse: the published file
    # mid-transaction with its journal, and no process holding its locks.
    crash = tmp_path / "crash"
    crash.mkdir()
    shutil.copy2(data / "snpedia.sqlite", crash / "snpedia.sqlite")
    conn = sqlite3.connect(crash / "snpedia.sqlite")
    try:
        conn.execute("PRAGMA cache_size=5")
        conn.execute("DELETE FROM snpedia_genotypes")
        conn.executemany(
            "INSERT OR IGNORE INTO snpedia_genotypes (rsid, allele1, allele2, summary) "
            "VALUES (?, 'A', 'G', ?)",
            [(f"rs{i}", "x" * 300) for i in range(2000)],
        )
        assert (crash / "snpedia.sqlite-journal").exists()
        shutil.copy2(crash / "snpedia.sqlite", data / "snpedia.sqlite")
        shutil.copy2(crash / "snpedia.sqlite-journal", data / "snpedia.sqlite-journal")
    finally:
        conn.rollback()
        conn.close()
    assert not SNPediaAnnotator(data).is_ready()

    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
    after = SNPediaAnnotator(data)
    assert after.is_ready() or result.exit_code != 0, (
        f"db update reported success but SNPedia is {after.cache_status()}:\n{result.output}"
    )


def _make_outdated(path: Path) -> None:
    """The parsed cache as an older parser left it: the update step re-parses it."""
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.execute(
            "UPDATE database_versions SET local_version_tag = 'pv:0' WHERE name='snpedia'"
        )
        conn.commit()


def test_interrupted_reparse_keeps_the_outdated_cache_you_had(snpedia_only, monkeypatch):
    """The update step's re-parse (parser version changed) parses a copy: a Ctrl-C
    during it leaves the cache as it was, rows and all."""
    data, _ = snpedia_only
    rows = SNPediaAnnotator(data).record_count()
    _make_outdated(data / "snpedia.sqlite")
    before = (data / "snpedia.sqlite").read_bytes()
    real = mwparserfromhell.parse
    calls = [0]

    def interrupted(text, *args, **kwargs):
        calls[0] += 1
        if calls[0] == 50:
            raise KeyboardInterrupt
        return real(text, *args, **kwargs)

    monkeypatch.setattr(snpedia_parser.mwparserfromhell, "parse", interrupted)
    with pytest.raises(KeyboardInterrupt):
        SNPediaAnnotator(data).prepare_cache()

    assert calls[0] >= 50
    assert (data / "snpedia.sqlite").read_bytes() == before
    assert not (data / "snpedia.sqlite-journal").exists()
    status = SNPediaAnnotator(data).cache_status()
    assert status.state is CacheState.OUTDATED, status
    assert SNPediaAnnotator(data).record_count() == rows


def test_failed_reparse_keeps_the_outdated_cache_and_says_so(snpedia_only, monkeypatch, capsys):
    data, _ = snpedia_only
    _make_outdated(data / "snpedia.sqlite")
    before = (data / "snpedia.sqlite").read_bytes()

    def disk_full(*args, **kwargs):
        raise sqlite3.OperationalError("database or disk is full")

    monkeypatch.setattr(snpedia_parser, "parse_raw_pages", disk_full)
    assert SNPediaAnnotator(data).prepare_cache() is False

    assert (data / "snpedia.sqlite").read_bytes() == before
    assert "database or disk is full" in capsys.readouterr().out


def test_reparse_publishes_a_parsed_cache(snpedia_only):
    data, _ = snpedia_only
    rows = SNPediaAnnotator(data).record_count()
    _make_outdated(data / "snpedia.sqlite")

    assert SNPediaAnnotator(data).prepare_cache() is True

    after = SNPediaAnnotator(data)
    assert after.is_ready(), after.cache_status()
    assert after.record_count() == rows
    assert sorted(p.name for p in data.iterdir() if p.name.startswith(".")) == []


def test_refresh_leaves_no_raw_pages_published_mid_parse(snpedia_only, monkeypatch):
    """While a refresh parses, the published file is still the parsed cache."""
    data, _ = snpedia_only
    seen: list[CacheState] = []
    real = snpedia_parser.parse_raw_pages

    def watching(path, **kwargs):
        seen.append(SNPediaAnnotator(data).cache_status().state)
        return real(path, **kwargs)

    monkeypatch.setattr(snpedia_parser, "parse_raw_pages", watching)
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data), "--force"])

    assert result.exit_code == 0, result.output
    assert seen == [CacheState.READY]
    assert SNPediaAnnotator(data).is_ready()


def test_installing_a_cache_removes_a_stale_journal_beside_it(tmp_path):
    """A journal left beside the old cache is settled before anything is replaced:
    SQLite removes one that is not hot when it opens the old file, so no side
    file is left and the install goes ahead."""
    source = tmp_path / "source.sqlite"
    with contextlib.closing(sqlite3.connect(source)) as conn:
        conn.execute("CREATE TABLE t (x)")
        conn.execute(
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, version TEXT, "
            "downloaded_at TEXT, source_url TEXT, record_count INTEGER, remote_signal TEXT)"
        )
        conn.commit()
    gz = tmp_path / "cache.sqlite.gz"
    gz.write_bytes(gzip.compress(source.read_bytes()))
    data = tmp_path / "data"
    data.mkdir()
    (data / "gnomad.sqlite").write_bytes(b"old")
    (data / "gnomad.sqlite-journal").write_bytes(b"stale journal")

    install_prebuilt_gz_cache(gz, data / "gnomad.sqlite", "gnomad")

    assert sorted(p.name for p in data.iterdir()) == ["gnomad.sqlite"]
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        assert conn.execute("SELECT COUNT(*) FROM t").fetchone() == (0,)


def test_publishing_over_a_hot_journal_leaves_the_new_cache_intact(tmp_path):
    """The old cache is recovered from its journal before the new one is moved in,
    so no journal is left to replay into the new file."""
    from allelix.databases.loader_utils import publish_cache

    old = tmp_path / "snpedia.sqlite"
    journal = tmp_path / "snpedia.sqlite-journal"
    with contextlib.closing(sqlite3.connect(old)) as conn:
        conn.execute("CREATE TABLE records (id INTEGER PRIMARY KEY, value TEXT)")
        conn.executemany(
            "INSERT INTO records VALUES (?, ?)", [(i, "old" * 100) for i in range(1000)]
        )
        conn.commit()
        conn.execute("PRAGMA cache_size=5")
        conn.execute("UPDATE records SET value = ?", ("uncommitted" * 100,))
        damaged, recovery = old.read_bytes(), journal.read_bytes()
        conn.rollback()
    old.write_bytes(damaged)
    journal.write_bytes(recovery)
    staged = tmp_path / "staged.sqlite"
    with contextlib.closing(sqlite3.connect(staged)) as conn:
        conn.execute("CREATE TABLE new_records (id INTEGER)")
        conn.execute("INSERT INTO new_records VALUES (7)")
        conn.commit()

    publish_cache(staged, old)

    assert not journal.exists()
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT id FROM new_records").fetchall() == [(7,)]
        assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)


def test_a_hot_journal_without_its_cache_is_removed_before_publication(tmp_path):
    """A journal whose cache file is gone holds nothing committed of anything:
    it is removed before the new cache moves in, so it cannot be replayed into it."""
    from allelix.databases.loader_utils import publish_cache

    old = tmp_path / "gnomad.sqlite"
    journal = tmp_path / "gnomad.sqlite-journal"
    journal.write_bytes(bytes.fromhex("d9d505f920a163d7") + b"\x00" * 504)
    staged = tmp_path / "staged.sqlite"
    with contextlib.closing(sqlite3.connect(staged)) as conn:
        conn.execute("CREATE TABLE t (x)")
        conn.execute("INSERT INTO t VALUES (7)")
        conn.commit()

    publish_cache(staged, old)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["gnomad.sqlite"]
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT x FROM t").fetchall() == [(7,)]


def test_wal_frames_without_their_cache_are_removed_before_publication(tmp_path):
    """WAL frames (and their -shm) whose cache file is gone are orphans: they are
    removed before the new cache moves in, so none is applied to it."""
    from allelix.databases.loader_utils import publish_cache

    old = tmp_path / "cadd.sqlite"
    (tmp_path / "cadd.sqlite-wal").write_bytes(bytes.fromhex("377f0682") + b"\x00" * 28)
    (tmp_path / "cadd.sqlite-shm").write_bytes(b"\x00" * 32768)
    staged = tmp_path / "staged.sqlite"
    with contextlib.closing(sqlite3.connect(staged)) as conn:
        conn.execute("CREATE TABLE t (x)")
        conn.execute("INSERT INTO t VALUES (7)")
        conn.commit()

    publish_cache(staged, old)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["cadd.sqlite"]
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT x FROM t").fetchall() == [(7,)]


def test_a_side_file_sqlite_leaves_in_place_is_refused_before_publication(tmp_path):
    """A stray -shm beside a rollback-mode cache is not removed by opening it, so
    it is still there after the settle step: nothing is replaced."""
    from allelix.databases.loader_utils import UnrecoveredJournalError, publish_cache

    old = tmp_path / "gwas.sqlite"
    with contextlib.closing(sqlite3.connect(old)) as conn:
        conn.execute("CREATE TABLE records (value TEXT)")
        conn.execute("INSERT INTO records VALUES ('old')")
        conn.commit()
    shm = tmp_path / "gwas.sqlite-shm"
    shm.write_bytes(b"\x00" * 32768)
    staged = tmp_path / "staged.sqlite"
    with contextlib.closing(sqlite3.connect(staged)) as conn:
        conn.execute("CREATE TABLE records (value TEXT)")
        conn.execute("INSERT INTO records VALUES ('new')")
        conn.commit()

    with pytest.raises(UnrecoveredJournalError, match=r"gwas\.sqlite-shm"):
        publish_cache(staged, old)

    assert shm.exists()
    with contextlib.closing(sqlite3.connect(old)) as conn:
        assert conn.execute("SELECT value FROM records").fetchall() == [("old",)]

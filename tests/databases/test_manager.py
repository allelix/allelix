# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for VCF parsing and SQLite cache loading."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

from allelix.databases.manager import (
    get_database_info,
    stamp_existing_clinvar_cache,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestGetDatabaseInfo:
    def test_missing_file_returns_none(self, tmp_path: Path):
        assert get_database_info(tmp_path / "nope.sqlite", "clinvar") is None

    def test_unknown_database_returns_none(self, tmp_path: Path):
        """A populated cache returns None for a name that isn't stored."""
        db = tmp_path / "clinvar.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                ("clinvar", "url", "20260101", "2026-06-01", 100, "md5:x", "iv:3"),
            )
            conn.commit()
        assert get_database_info(db, "pharmgkb") is None

    def test_garbage_file_returns_none(self, tmp_path: Path):
        f = tmp_path / "garbage.sqlite"
        f.write_text("not a database", encoding="utf-8")
        assert get_database_info(f, "clinvar") is None

    def test_legacy_v041_schema_returns_none_remote_signal(self, tmp_path: Path):
        """A v0.4.1 cache lacks the `remote_signal` column.

        get_database_info must fall back gracefully and report
        remote_signal=None so the next `db update` triggers a refresh
        (because remote != cached==None) and writes a v0.4.2 row.
        """
        db = tmp_path / "legacy.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            # Recreate the v0.4.1 schema verbatim — no remote_signal column.
            conn.executescript(
                """
                CREATE TABLE database_versions (
                    name TEXT PRIMARY KEY,
                    source_url TEXT NOT NULL,
                    version TEXT,
                    downloaded_at TEXT NOT NULL,
                    record_count INTEGER NOT NULL
                );
                """
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?)",
                ("clinvar", "old://url", "20240101", "2024-01-01T00:00:00", 100),
            )
            conn.commit()

        info = get_database_info(db, "clinvar")
        assert info is not None
        assert info["version"] == "20240101"
        assert info["record_count"] == 100
        assert info["remote_signal"] is None
        assert info["local_version_tag"] is None

    def test_pre_v150_schema_lazily_adds_local_version_tag(self, tmp_path: Path):
        """A pre-v1.5.0 cache has remote_signal but no local_version_tag.

        get_database_info lazily adds the column so that all caches
        (including gnomAD/AlphaMissense which don't use version tags)
        have a consistent schema after any db status or is_ready() call.
        """
        db = tmp_path / "pre150.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                """
                CREATE TABLE database_versions (
                    name TEXT PRIMARY KEY,
                    source_url TEXT NOT NULL,
                    version TEXT,
                    downloaded_at TEXT NOT NULL,
                    record_count INTEGER NOT NULL,
                    remote_signal TEXT
                );
                """
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?)",
                ("gnomad", "hf://url", "4.1", "2026-01-01T00:00:00", 16000000, "etag:abc"),
            )
            conn.commit()

        info = get_database_info(db, "gnomad")
        assert info is not None
        assert info["remote_signal"] == "etag:abc"
        assert info["local_version_tag"] is None

        # Column was lazily added — second read uses the 7-col path.
        with contextlib.closing(sqlite3.connect(db)) as conn:
            cols = {row[1] for row in conn.execute("PRAGMA table_info(database_versions)")}
            assert "local_version_tag" in cols

    def test_v041_schema_unknown_name_returns_none(self, tmp_path: Path):
        """4-column fallback returns None when the row doesn't match."""
        db = tmp_path / "legacy.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?)",
                ("clinvar", "url", "20240101", "2024-01-01", 100),
            )
            conn.commit()
        assert get_database_info(db, "missing_name") is None

    def test_pre_v150_unknown_name_returns_none(self, tmp_path: Path):
        """5-column fallback returns None when the row doesn't match."""
        db = tmp_path / "pre150.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?)",
                ("clinvar", "url", "20240101", "2024-01-01", 100, "etag:x"),
            )
            conn.commit()
        assert get_database_info(db, "missing_name") is None


class TestStampExistingClinvarCache:
    def test_nonexistent_db_returns_false(self, tmp_path: Path):
        assert stamp_existing_clinvar_cache(tmp_path / "nope.sqlite") is False

    def test_no_clinvar_rows_returns_false(self, tmp_path: Path):
        db = tmp_path / "empty.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
        assert stamp_existing_clinvar_cache(db) is False

    def test_no_database_versions_table_returns_false(self, tmp_path: Path):
        db = tmp_path / "bare.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE other (x TEXT)")
            conn.commit()
        assert stamp_existing_clinvar_cache(db) is False

    def test_stale_tag_returns_false(self, tmp_path: Path):
        db = tmp_path / "stale.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                ("clinvar.GRCh37", "url", "20240101", "2024-01-01", 100, "md5:x", "iv:999"),
            )
            conn.commit()
        assert stamp_existing_clinvar_cache(db) is False

    def test_null_tag_unknown_version_not_promoted(self, tmp_path: Path):
        """PR-1: a cache with NULL local_version_tag AND no baked `|iv:N`
        marker is an unknown-version legacy cache. The pre-PR-1 behavior
        auto-stamped it to the current interpreter version — silent
        promotion across the iv:2→iv:3 format boundary would serve old
        single-row VCF data labeled as fresh per-SCV TSV data. The fix:
        return False so db update reingests."""
        db = tmp_path / "pre_mechanism.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                ("clinvar.GRCh37", "url", "20240101", "2024-01-01", 100, "md5:x", None),
            )
            conn.commit()

        # PR-1 safety fix: reingest, do not silent-promote.
        assert stamp_existing_clinvar_cache(db) is False

        # Cache row left untouched — tag stays NULL so is_ready() also
        # returns False on the next call (the proper reingest signal).
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT local_version_tag, remote_signal "
                "FROM database_versions WHERE name='clinvar.GRCh37'"
            ).fetchone()
        assert row[0] is None
        assert row[1] == "md5:x"

    def test_baked_older_version_not_promoted(self, tmp_path: Path):
        """PR-1: NULL tag with a baked `|iv:N` marker for a NON-current
        N means the cache was built by a prior interpreter version.
        Promoting it to the current tag without reingesting would serve
        old-format data as fresh — return False so db update reingests."""
        db = tmp_path / "baked_old.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            # iv:1 baked — pre-v2.0.1 era. Strictly less than the
            # current CLINVAR_INTERPRETER_VERSION (3 as of PR 8b).
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                ("clinvar.GRCh37", "url", "20240101", "2024-01-01", 100, "md5:abc|iv:1", None),
            )
            conn.commit()

        assert stamp_existing_clinvar_cache(db) is False

        # Untouched: tag stays NULL, marker NOT scrubbed from signal.
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT local_version_tag, remote_signal "
                "FROM database_versions WHERE name='clinvar.GRCh37'"
            ).fetchone()
        assert row[0] is None
        assert row[1] == "md5:abc|iv:1"

    def test_legacy_iv_in_remote_signal_migrated(self, tmp_path: Path):
        """#79 missing-branch coverage: an old-format cache has the
        interpreter tag baked into remote_signal as ``...|iv:N``. The
        migration path moves the tag into local_version_tag and cleans
        the remote_signal."""
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        db = tmp_path / "legacy.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    "clinvar.GRCh38",
                    "url",
                    "20240101",
                    "2024-01-01",
                    100,
                    f"md5:abcdef|iv:{CLINVAR_INTERPRETER_VERSION}",
                    None,
                ),
            )
            conn.commit()

        assert stamp_existing_clinvar_cache(db) is True

        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT local_version_tag, remote_signal "
                "FROM database_versions WHERE name='clinvar.GRCh38'"
            ).fetchone()
        assert row[0] == f"iv:{CLINVAR_INTERPRETER_VERSION}"
        assert row[1] == "md5:abcdef"  # the |iv:N suffix was scrubbed

    def test_already_current_tag_is_idempotent(self, tmp_path: Path):
        """#79 missing-branch coverage: a cache already stamped with the
        current interpreter tag is a no-op success — return True without
        writing."""
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        current_tag = f"iv:{CLINVAR_INTERPRETER_VERSION}"
        db = tmp_path / "current.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                ("clinvar.GRCh37", "url", "20240101", "2024-01-01", 100, "md5:x", current_tag),
            )
            conn.commit()

        assert stamp_existing_clinvar_cache(db) is True

        # remote_signal must be untouched (no |iv:N scrub on a current tag)
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT local_version_tag, remote_signal "
                "FROM database_versions WHERE name='clinvar.GRCh37'"
            ).fetchone()
        assert row[0] == current_tag
        assert row[1] == "md5:x"

    def test_dual_build_caches_both_self_heal_when_baked_current(self, tmp_path: Path):
        """Both GRCh37 and GRCh38 rows self-heal in one call when both
        carry a baked `|iv:CURRENT` marker. Post-PR-1 the self-heal path
        requires positive version evidence on every row."""
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        marker = f"|iv:{CLINVAR_INTERPRETER_VERSION}"
        db = tmp_path / "dual.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.executemany(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                [
                    ("clinvar.GRCh37", "url37", "v37", "2024-01-01", 50, f"md5:a{marker}", None),
                    ("clinvar.GRCh38", "url38", "v38", "2024-01-01", 60, f"md5:b{marker}", None),
                ],
            )
            conn.commit()

        assert stamp_existing_clinvar_cache(db) is True

        expected_tag = f"iv:{CLINVAR_INTERPRETER_VERSION}"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            rows = list(
                conn.execute(
                    "SELECT name, local_version_tag, remote_signal "
                    "FROM database_versions ORDER BY name"
                )
            )
        tags = {name: tag for name, tag, _ in rows}
        signals = {name: sig for name, _, sig in rows}
        assert tags["clinvar.GRCh37"] == expected_tag
        assert tags["clinvar.GRCh38"] == expected_tag
        # The |iv: marker was scrubbed from both signals.
        assert signals["clinvar.GRCh37"] == "md5:a"
        assert signals["clinvar.GRCh38"] == "md5:b"

    def test_dual_builds_one_unknown_blocks_both(self, tmp_path: Path):
        """If ANY row fails the self-heal check the whole cache is treated
        as untrustworthy and reingest. Mixed state isn't safe — partial
        promotion would land per-build dispatch on inconsistent data.
        """
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        marker_current = f"|iv:{CLINVAR_INTERPRETER_VERSION}"
        db = tmp_path / "mixed.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                "CREATE TABLE database_versions ("
                "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
                "version TEXT, downloaded_at TEXT NOT NULL, "
                "record_count INTEGER NOT NULL, remote_signal TEXT, "
                "local_version_tag TEXT);"
            )
            conn.executemany(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?, ?)",
                [
                    # First row has the marker, would self-heal alone…
                    (
                        "clinvar.GRCh37",
                        "url37",
                        "v37",
                        "2024-01-01",
                        50,
                        f"md5:a{marker_current}",
                        None,
                    ),
                    # …but second row has no marker, so the whole call returns False.
                    ("clinvar.GRCh38", "url38", "v38", "2024-01-01", 60, "md5:b", None),
                ],
            )
            conn.commit()

        assert stamp_existing_clinvar_cache(db) is False

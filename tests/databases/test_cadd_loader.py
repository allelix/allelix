# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for `allelix.databases.cadd_loader.install_prebuilt_cache`.

The other three pre-built cache loaders (gnomad / alphamissense / snpedia)
are exercised end-to-end by their respective annotator tests, so the
shared loader_utils body sees real input under tests. cadd_loader is the
exception — every cadd test fakes `install_prebuilt_cache`, so the
facade body itself never runs (cadd_loader.py line 42 missed in
coverage, audit #79 batch 2).

This file exercises the cadd_loader facade directly, end-to-end:
construct a synthetic gzipped SQLite cache, point install_prebuilt_cache
at it, verify the decompressed cache lands with the right schema_version
tag and (when supplied) remote_signal.

Naming note: `test_cadd_loader.py` (this file) tests the cadd_loader
module. `test_build_cadd_cache.py` (separate, pre-existing) tests
`scripts/build_cadd_cache.py` — the standalone cache-builder script.
The two are different layers and the audit (#79) flagged the naming
collision as misleading; this file establishes the explicit
`test_<loader>.py` convention for future loader tests.
"""

from __future__ import annotations

import contextlib
import gzip
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.databases._versions import CADD_SCHEMA_VERSION
from allelix.databases.cadd_loader import (
    CADD_DB_FILENAME,
    install_prebuilt_cache,
)
from allelix.databases.schema import CADD_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path


def _make_synthetic_gz_cache(tmp_path: Path) -> Path:
    """Build a tiny gzipped SQLite file shaped like the real CADD cache.

    Returns the gz path. The CADD_SCHEMA matches what the production
    builder emits; one synthetic row gives install_prebuilt_cache
    something non-empty to land.
    """
    src_db = tmp_path / "src.sqlite"
    with contextlib.closing(sqlite3.connect(src_db)) as conn:
        for stmt in CADD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.execute(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
            ("1", 12345, "A", "G", 25.5),
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count)"
            " VALUES (?, ?, ?, ?, ?)",
            ("cadd", "test://prebuilt", "v1.7", "2026-01-01", 1),
        )
        conn.commit()

    gz_path = tmp_path / "cadd.sqlite.gz"
    with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        f_out.write(f_in.read())
    return gz_path


class TestInstallPrebuiltCache:
    """End-to-end tests for the cadd_loader facade.

    `install_prebuilt_cache` is a 4-line delegation to
    `loader_utils.install_prebuilt_gz_cache` with the cadd-specific
    record_name and schema_version_tag. Indirect coverage via the
    annotator tests bypasses this delegation (the cadd annotator tests
    monkeypatch `allelix.annotators.cadd.install_prebuilt_cache`), so
    these tests are the only place the cadd_loader body actually runs.
    """

    def test_decompresses_and_stamps_schema_version(self, tmp_path: Path) -> None:
        """Happy path: schema_version_tag must be stamped to `sv:N`."""
        gz_path = _make_synthetic_gz_cache(tmp_path)
        dest_db = tmp_path / "dest" / CADD_DB_FILENAME
        dest_db.parent.mkdir()

        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:cadd-test")

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal, local_version_tag, source_url "
                "FROM database_versions WHERE name = 'cadd'"
            ).fetchone()
        assert row[0] == "etag:cadd-test"
        assert row[1] == f"sv:{CADD_SCHEMA_VERSION}"

    def test_decompresses_without_remote_signal(self, tmp_path: Path) -> None:
        """remote_signal=None (the default) must not write a signal but must
        still stamp the schema_version_tag."""
        gz_path = _make_synthetic_gz_cache(tmp_path)
        dest_db = tmp_path / "dest" / CADD_DB_FILENAME
        dest_db.parent.mkdir()

        install_prebuilt_cache(gz_path, dest_db)

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal, local_version_tag "
                "FROM database_versions WHERE name = 'cadd'"
            ).fetchone()
        # Pre-existing signal from the synthetic source survives the
        # decompression since loader_utils.stamp_remote_signal is only
        # called when remote_signal is provided.
        assert row[0] is None or row[0] == ""
        assert row[1] == f"sv:{CADD_SCHEMA_VERSION}"

    def test_atomic_replace_uses_tmp_path(self, tmp_path: Path) -> None:
        """If a .tmp file is left over from a prior aborted run, it's cleared."""
        gz_path = _make_synthetic_gz_cache(tmp_path)
        dest_db = tmp_path / "dest" / CADD_DB_FILENAME
        dest_db.parent.mkdir()

        # Plant a stale tmp file matching the loader_utils naming convention.
        stale_tmp = dest_db.parent / f"{dest_db.name}.tmp"
        stale_tmp.write_bytes(b"garbage from previous aborted run")
        assert stale_tmp.exists()

        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:atomic")

        # The stale tmp is gone (replaced and then renamed to dest)
        assert not stale_tmp.exists()
        # The real dest is the freshly decompressed cache
        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_overwrites_existing_dest(self, tmp_path: Path) -> None:
        """An existing cache at the destination must be replaced atomically."""
        gz_path = _make_synthetic_gz_cache(tmp_path)
        dest_db = tmp_path / "dest" / CADD_DB_FILENAME
        dest_db.parent.mkdir()

        # Plant an old cache at the destination.
        dest_db.write_bytes(b"old-cache-bytes")
        old_size = dest_db.stat().st_size

        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:replace")

        assert dest_db.exists()
        new_size = dest_db.stat().st_size
        assert new_size != old_size  # was replaced
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            tables = {
                row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
            }
        # New schema present → it's actually the decompressed cache, not
        # the planted bytes.
        assert "cadd_scores" in tables

    def test_insufficient_disk_space_raises(self, tmp_path: Path, monkeypatch) -> None:
        """The pre-decompression disk-space check refuses with a clear OSError.

        Exercises loader_utils' DISK_SPACE_MULTIPLIER guard from this
        facade. Without the patch the guard never fires on a normal dev
        machine — testing it via the cadd facade keeps the assertion at
        the loader boundary where the user sees it.
        """
        from collections import namedtuple

        from allelix.databases import loader_utils

        gz_path = _make_synthetic_gz_cache(tmp_path)
        dest_db = tmp_path / "dest" / CADD_DB_FILENAME
        dest_db.parent.mkdir()

        FakeUsage = namedtuple("FakeUsage", ["total", "used", "free"])

        def fake_disk_usage(_path: Path) -> FakeUsage:
            return FakeUsage(total=1, used=1, free=1)  # 1 byte free

        monkeypatch.setattr(loader_utils.shutil, "disk_usage", fake_disk_usage)

        with pytest.raises(OSError, match="Not enough disk space"):
            install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:nospace")

        # Destination must not be created on failure.
        assert not dest_db.exists()

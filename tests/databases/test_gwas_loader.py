# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the GWAS Catalog SQLite loader."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.databases import gwas_loader
from allelix.databases.gwas_loader import (
    _CATEGORIZER_VERSION,
    classify_gwas_trait,
    load_gwas_tsv,
    schema_is_current,
)
from allelix.databases.schema import GWAS_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path


class TestBatchedInsert:
    """Exercise the batch-flush path in load_gwas_tsv."""

    def test_batched_insert_flushes(
        self, tmp_path: Path, mock_gwas_tsv: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(gwas_loader, "INSERT_BATCH_SIZE", 3)
        executemany_payloads: list[int] = []
        real_connect = sqlite3.connect

        class _SpyConn:
            def __init__(self, real: sqlite3.Connection) -> None:
                self._real = real

            def executemany(self, sql: str, seq: list) -> sqlite3.Cursor:
                seq_list = list(seq)
                executemany_payloads.append(len(seq_list))
                return self._real.executemany(sql, seq_list)

            def __getattr__(self, name: str) -> object:
                return getattr(self._real, name)

        def spying_connect(*args: object, **kwargs: object) -> _SpyConn:
            return _SpyConn(real_connect(*args, **kwargs))

        monkeypatch.setattr(sqlite3, "connect", spying_connect)

        db = tmp_path / "gwas.sqlite"
        count = load_gwas_tsv(mock_gwas_tsv, db, source_url="test://batch")
        assert count == 8
        # 8 records / batch_size 3 = 2 full batches + 1 remainder.
        assert executemany_payloads == [3, 3, 2]


class TestClassifierDiseaseTrait:
    """classify_gwas_trait uses DISEASE/TRAIT when MAPPED_TRAIT is empty."""

    def test_uses_disease_trait_when_mapped_trait_empty(self) -> None:
        assert (
            classify_gwas_trait(
                mapped_trait="",
                mapped_trait_uri="http://purl.obolibrary.org/obo/HP_0000924",
                disease_trait="Impedance of whole body (UKB data field 23106)",
            )
            == "body_measurement"
        )

    def test_uses_disease_trait_for_arm_impedance(self) -> None:
        assert (
            classify_gwas_trait(
                "",
                "http://purl.obolibrary.org/obo/HP_0000924",
                disease_trait="Impedance of arm left (UKB data field 23110)",
            )
            == "body_measurement"
        )
        assert (
            classify_gwas_trait(
                "",
                "http://purl.obolibrary.org/obo/HP_0000924",
                disease_trait="Impedance of arm right (UKB data field 23109)",
            )
            == "body_measurement"
        )

    def test_uses_mapped_trait_when_disease_trait_empty(self) -> None:
        assert (
            classify_gwas_trait(
                mapped_trait="Whole body water mass",
                mapped_trait_uri="",
            )
            == "body_measurement"
        )


class TestSchemaIsCurrentEdge:
    """schema_is_current defensive branches."""

    def test_missing_file_returns_false(self, tmp_path: Path) -> None:
        assert not schema_is_current(tmp_path / "nonexistent.sqlite")

    def test_empty_db_returns_false(self, tmp_path: Path) -> None:
        db = tmp_path / "empty.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(GWAS_SCHEMA)
        assert not schema_is_current(db)

    def test_missing_gwas_table_returns_false(self, tmp_path: Path) -> None:
        db = tmp_path / "no_table.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript("""
                CREATE TABLE database_versions (
                    name TEXT PRIMARY KEY,
                    source_url TEXT NOT NULL,
                    version TEXT,
                    downloaded_at TEXT NOT NULL,
                    record_count INTEGER NOT NULL,
                    remote_signal TEXT,
                    local_version_tag TEXT
                );
            """)
            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES ('gwas', 'http://x', '2026', '2026-01-01T00:00:00Z', "
                "0, 'etag:abc', ?)",
                (f"cv:{_CATEGORIZER_VERSION}",),
            )
            conn.commit()
        assert not schema_is_current(db)


class TestCategorizerVersion:
    """_CATEGORIZER_VERSION marker in local_version_tag for cache invalidation."""

    def test_schema_is_current_rejects_cache_without_cv_tag(self, tmp_path: Path) -> None:
        db = tmp_path / "stale.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(GWAS_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES ('gwas', 'http://x', '2026-05-19', "
                "'2026-05-19T00:00:00Z', 0, 'etag:abc', NULL)",
            )
            conn.commit()
        assert not schema_is_current(db)

    def test_schema_is_current_accepts_matching_cv_tag(self, tmp_path: Path) -> None:
        db = tmp_path / "fresh.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(GWAS_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES ('gwas', 'http://x', '2026-05-19', "
                "'2026-05-19T00:00:00Z', 0, 'etag:abc', ?)",
                (f"cv:{_CATEGORIZER_VERSION}",),
            )
            conn.commit()
        assert schema_is_current(db)

    def test_schema_is_current_rejects_old_cv_tag(self, tmp_path: Path) -> None:
        db = tmp_path / "old.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(GWAS_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES ('gwas', 'http://x', '2026-05-19', "
                "'2026-05-19T00:00:00Z', 0, 'etag:abc', 'cv:1')",
            )
            conn.commit()
        assert not schema_is_current(db)

    def test_load_stamps_categorizer_version(self, tmp_path: Path, mock_gwas_tsv: Path) -> None:
        db = tmp_path / "gwas.sqlite"
        load_gwas_tsv(mock_gwas_tsv, db, source_url="test://cv", remote_signal="etag:xyz")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT remote_signal, local_version_tag FROM database_versions WHERE name='gwas'"
            ).fetchone()
        assert row is not None
        assert row[0] == "etag:xyz"
        assert row[1] == f"cv:{_CATEGORIZER_VERSION}"


class TestMinRowsFloor:
    """GH #19: post-load row-count sanity gate catches mid-stream truncation.

    Production callers pass ``min_rows=GWAS_MIN_ROWS``. A truncated download
    (chunked transfer without Content-Length, connection drop mid-stream)
    produces a parseable-but-short TSV that previously committed silently
    to the cache. The floor check raises before ``os.replace`` so the cache
    isn't overwritten with the short file.
    """

    def test_below_floor_raises_and_keeps_cache_clean(self, mock_gwas_tsv, tmp_path) -> None:
        from allelix.databases.gwas_loader import load_gwas_tsv

        db = tmp_path / "gwas.sqlite"
        # Mock fixture has a handful of rows; set min_rows above that.
        with pytest.raises(OSError, match="rows ingested"):
            load_gwas_tsv(mock_gwas_tsv, db, source_url="test://truncated", min_rows=999_999)
        # tmp_path cleaned, real db not created
        assert not db.exists()
        assert not (tmp_path / "gwas.sqlite.tmp").exists()

    def test_default_min_rows_zero_loads_small_fixtures(self, mock_gwas_tsv, tmp_path) -> None:
        """Default ``min_rows=0`` lets test fixtures load regardless of size."""
        from allelix.databases.gwas_loader import load_gwas_tsv

        db = tmp_path / "gwas.sqlite"
        count = load_gwas_tsv(mock_gwas_tsv, db, source_url="test://small")
        assert count > 0
        assert db.exists()

    def test_at_or_above_floor_succeeds(self, mock_gwas_tsv, tmp_path) -> None:
        """``min_rows`` ≤ actual count → load completes."""
        from allelix.databases.gwas_loader import load_gwas_tsv

        db = tmp_path / "gwas.sqlite"
        # Pick a floor at 1 — the mock fixture has more than that.
        count = load_gwas_tsv(mock_gwas_tsv, db, source_url="test://at-floor", min_rows=1)
        assert count >= 1
        assert db.exists()

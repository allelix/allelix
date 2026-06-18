# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the AlphaMissense enrichment annotator."""

from __future__ import annotations

import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.databases.schema import ALPHAMISSENSE_SCHEMA
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


def _build_db(tmp_path: Path) -> Path:
    """Build a minimal AlphaMissense SQLite cache for testing."""
    db_path = tmp_path / "alphamissense.sqlite"
    conn = sqlite3.connect(db_path)
    conn.executescript(ALPHAMISSENSE_SCHEMA)
    conn.executemany(
        "INSERT INTO alphamissense_scores "
        "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
        "protein_variant, am_pathogenicity, am_class) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            (
                "1",
                100,
                "A",
                "G",
                "rs1001",
                "P12345",
                "ENST001",
                "A100G",
                0.95,
                "likely_pathogenic",
            ),
            ("1", 200, "C", "T", "rs1002", "P12345", "ENST001", "C200T", 0.20, "likely_benign"),
            ("2", 300, "G", "A", "rs2001", "P67890", "ENST002", "G300A", 0.45, "ambiguous"),
            ("3", 400, "T", "C", None, "P99999", "ENST003", "T400C", 0.80, "likely_pathogenic"),
        ],
    )
    conn.execute(
        "INSERT INTO database_versions (name, source_url, version, downloaded_at, record_count) "
        "VALUES ('alphamissense', 'https://test', '2023.1', '2026-06-08', 4)",
    )
    conn.commit()
    conn.close()
    return db_path


class TestAlphaMissenseAnnotator:
    def test_is_ready(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.is_ready()

    def test_not_ready_without_db(self, tmp_path: Path) -> None:
        am = AlphaMissenseAnnotator(tmp_path)
        assert not am.is_ready()

    def test_version(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.version() == "2023.1"

    def test_record_count(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.record_count() == 4

    def test_annotate_returns_empty(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        v = Variant(rsid="rs1001", chromosome="1", position=100, allele1="A", allele2="G")
        assert am.annotate(v) == []

    def test_lookup_found(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.lookup("rs1001")
        assert result is not None
        score, cls = result
        assert score == 0.95
        assert cls == "likely_pathogenic"
        am.close()

    def test_lookup_not_found(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.lookup("rs9999") is None
        am.close()

    def test_bulk_lookup(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs1001", "rs1002", "rs9999"})
        assert len(result) == 2
        assert result["rs1001"] == (0.95, "likely_pathogenic")
        assert result["rs1002"] == (0.20, "likely_benign")
        assert "rs9999" not in result
        am.close()

    def test_bulk_lookup_empty(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.bulk_lookup(set()) == {}
        am.close()

    def test_close(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        am.lookup("rs1001")
        am.close()
        assert am._conn is None

    def test_am_class_thresholds(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs1001", "rs1002", "rs2001"})
        assert result["rs1001"][1] == "likely_pathogenic"
        assert result["rs1002"][1] == "likely_benign"
        assert result["rs2001"][1] == "ambiguous"
        am.close()


class TestCheckGnomadVersion:
    """Warning paths in _check_gnomad_version.

    The AlphaMissense cache stamps which gnomAD version it was built
    against (`alphamissense_gnomad_source` row). On first connection, the
    annotator compares that stamp against the installed gnomAD cache and
    warns if they diverge. Three branches:

    1. `stamped == "no_gnomad"` — cache was built without gnomAD rsID
       resolution; rsID lookups will silently miss. Warn loudly.
    2. `stamped != installed` — cache was built against a different
       gnomAD version; mappings may be stale. Warn with rebuild hint.
    3. `stamped is None` or `installed is None` — no comparison possible;
       silently return.
    """

    def _build_with_gnomad_source(self, tmp_path: Path, gnomad_source: str) -> Path:
        """Build an AM cache stamped with a specific alphamissense_gnomad_source value."""
        db_path = _build_db(tmp_path)
        conn = sqlite3.connect(db_path)
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES (?, ?, ?, ?, ?)",
            ("alphamissense_gnomad_source", "internal://stamp", gnomad_source, "2026-06-08", 0),
        )
        conn.commit()
        conn.close()
        return db_path

    def _install_gnomad_cache(self, tmp_path: Path, version: str) -> None:
        """Stage a fake gnomad.sqlite at the expected path with a version row."""
        from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME

        gnomad_db = tmp_path / GNOMAD_DB_FILENAME
        conn = sqlite3.connect(gnomad_db)
        conn.execute(
            "CREATE TABLE database_versions ("
            "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
            "version TEXT, downloaded_at TEXT NOT NULL, "
            "record_count INTEGER NOT NULL, remote_signal TEXT, "
            "local_version_tag TEXT)"
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES (?, ?, ?, ?, ?)",
            ("gnomad", "test://gnomad", version, "2026-06-08", 0),
        )
        conn.commit()
        conn.close()

    def test_no_gnomad_stamp_warns_loudly(self, tmp_path: Path, caplog) -> None:
        """Branch 1: stamp == 'no_gnomad'. Warn that rsID lookups will miss."""
        self._build_with_gnomad_source(tmp_path, "no_gnomad")
        am = AlphaMissenseAnnotator(tmp_path)
        with caplog.at_level("WARNING"):
            am._connection()  # triggers _check_gnomad_version
        am.close()
        joined = " ".join(rec.message for rec in caplog.records)
        assert "no_gnomad" in joined or "without gnomAD" in joined
        assert "rsID lookups" in joined

    def test_version_mismatch_warns_with_rebuild_hint(self, tmp_path: Path, caplog) -> None:
        """Branch 2: stamp != installed. Warn with the rebuild command."""
        self._build_with_gnomad_source(tmp_path, "4.0")
        self._install_gnomad_cache(tmp_path, "4.1")  # different version
        am = AlphaMissenseAnnotator(tmp_path)
        with caplog.at_level("WARNING"):
            am._connection()
        am.close()
        joined = " ".join(rec.message for rec in caplog.records)
        assert "4.0" in joined
        assert "4.1" in joined
        assert "stale" in joined.lower() or "rebuild" in joined.lower()
        assert "build_alphamissense_cache" in joined

    def test_version_match_does_not_warn(self, tmp_path: Path, caplog) -> None:
        """Stamp == installed. Silent — no warning emitted."""
        self._build_with_gnomad_source(tmp_path, "4.1")
        self._install_gnomad_cache(tmp_path, "4.1")
        am = AlphaMissenseAnnotator(tmp_path)
        with caplog.at_level("WARNING"):
            am._connection()
        am.close()
        warnings = [
            rec for rec in caplog.records if "gnomAD" in rec.message or "gnomad" in rec.message
        ]
        assert warnings == []

    def test_no_gnomad_cache_installed_returns_silently(self, tmp_path: Path, caplog) -> None:
        """Branch 3: stamp present but no gnomad.sqlite at the expected path.
        get_database_info returns None and the check silently exits."""
        self._build_with_gnomad_source(tmp_path, "4.1")
        # Deliberately do NOT install a gnomad cache.
        am = AlphaMissenseAnnotator(tmp_path)
        with caplog.at_level("WARNING"):
            am._connection()
        am.close()
        warnings = [rec for rec in caplog.records if "stale" in rec.message.lower()]
        assert warnings == []

    def test_no_stamp_row_returns_silently(self, tmp_path: Path, caplog) -> None:
        """No `alphamissense_gnomad_source` row at all — older cache; silent."""
        _build_db(tmp_path)  # base build, no gnomad_source stamp
        am = AlphaMissenseAnnotator(tmp_path)
        with caplog.at_level("WARNING"):
            am._connection()
        am.close()
        warnings = [
            rec for rec in caplog.records if "gnomAD" in rec.message or "gnomad" in rec.message
        ]
        assert warnings == []


class TestBulkLookupByAlt:
    """Exact (rsid, alt) lookup for multi-allelic enrichment."""

    def test_exact_match(self, tmp_path: Path) -> None:
        db_path = tmp_path / "alphamissense.sqlite"
        conn = sqlite3.connect(db_path)
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores "
            "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
            "protein_variant, am_pathogenicity, am_class) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("1", 100, "A", "G", "rs5000", "P1", "E1", "A1G", 0.20, "likely_benign"),
                ("1", 100, "A", "T", "rs5000", "P1", "E1", "A1T", 0.92, "likely_pathogenic"),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES ('alphamissense', 'test', '1.0', '2026-01-01', 2)",
        )
        conn.commit()
        conn.close()

        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup_by_alt({("rs5000", "G"), ("rs5000", "T")})
        assert result[("rs5000", "G")] == (0.20, "likely_benign")
        assert result[("rs5000", "T")] == (0.92, "likely_pathogenic")
        am.close()

    def test_miss_returns_empty(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup_by_alt({("rs1001", "C")})
        assert result == {}
        am.close()

    def test_empty_input(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        assert am.bulk_lookup_by_alt(set()) == {}
        am.close()

    def test_mixed_hit_and_miss(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup_by_alt({("rs1001", "G"), ("rs1001", "X")})
        assert ("rs1001", "G") in result
        assert ("rs1001", "X") not in result
        am.close()


class TestInstallPrebuiltCache:
    """install_prebuilt_cache decompresses and stamps signal."""

    def test_decompress_and_stamp(self, tmp_path: Path) -> None:
        import contextlib
        import gzip

        from allelix.databases.alphamissense_loader import (
            ALPHAMISSENSE_DB_FILENAME,
            install_prebuilt_cache,
        )

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            conn.executescript(ALPHAMISSENSE_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("alphamissense", "test://prebuilt", "2023.1", "2026-01-01", 100),
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / "dest" / ALPHAMISSENSE_DB_FILENAME
        dest_db.parent.mkdir()
        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:am123")

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal, local_version_tag "
                "FROM database_versions WHERE name = 'alphamissense'"
            ).fetchone()
        assert row[0] == "etag:am123"
        assert row[1] == "sv:1"

    def test_decompress_without_signal(self, tmp_path: Path) -> None:
        import contextlib
        import gzip

        from allelix.databases.alphamissense_loader import (
            ALPHAMISSENSE_DB_FILENAME,
            install_prebuilt_cache,
        )

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            conn.executescript(ALPHAMISSENSE_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("alphamissense", "test://prebuilt", "2023.1", "2026-01-01", 100),
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / ALPHAMISSENSE_DB_FILENAME
        install_prebuilt_cache(gz_path, dest_db)

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal, local_version_tag "
                "FROM database_versions WHERE name = 'alphamissense'"
            ).fetchone()
        assert row[0] is None
        assert row[1] == "sv:1"

    def test_disk_space_check(self, tmp_path: Path) -> None:
        import gzip
        from unittest.mock import patch

        from allelix.databases.alphamissense_loader import install_prebuilt_cache

        gz_path = tmp_path / "tiny.gz"
        with gzip.open(gz_path, "wb") as f:
            f.write(b"x" * 100)

        dest_db = tmp_path / "out.sqlite"
        fake_usage = type("Usage", (), {"free": 1})()
        target = "allelix.databases.loader_utils.shutil.disk_usage"
        with (
            patch(target, return_value=fake_usage),
            pytest.raises(OSError, match="Not enough disk space"),
        ):
            install_prebuilt_cache(gz_path, dest_db)

    def test_replaces_existing_tmp(self, tmp_path: Path) -> None:
        import contextlib
        import gzip

        from allelix.databases.alphamissense_loader import (
            ALPHAMISSENSE_DB_FILENAME,
            install_prebuilt_cache,
        )

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            conn.executescript(ALPHAMISSENSE_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("alphamissense", "test://prebuilt", "2023.1", "2026-01-01", 100),
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / ALPHAMISSENSE_DB_FILENAME
        stale_tmp = tmp_path / f"{ALPHAMISSENSE_DB_FILENAME}.tmp"
        stale_tmp.write_text("stale")

        install_prebuilt_cache(gz_path, dest_db)
        assert dest_db.exists()
        assert not stale_tmp.exists()

    def test_stamps_signal_when_versions_table_missing(self, tmp_path: Path) -> None:
        """Pre-built cache without database_versions must not crash."""
        import contextlib
        import gzip

        from allelix.databases.alphamissense_loader import (
            ALPHAMISSENSE_DB_FILENAME,
            install_prebuilt_cache,
        )

        src_db = tmp_path / "source.sqlite"
        with contextlib.closing(sqlite3.connect(src_db)) as conn:
            conn.execute(
                "CREATE TABLE alphamissense_scores ("
                "chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, rsid TEXT,"
                " uniprot_id TEXT, transcript_id TEXT, protein_variant TEXT,"
                " am_pathogenicity REAL NOT NULL, am_class TEXT NOT NULL,"
                " PRIMARY KEY (chrom, pos, ref, alt))"
            )
            conn.commit()

        gz_path = tmp_path / "test.sqlite.gz"
        with src_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            f_out.write(f_in.read())

        dest_db = tmp_path / ALPHAMISSENSE_DB_FILENAME
        install_prebuilt_cache(gz_path, dest_db, remote_signal="etag:no-table")

        assert dest_db.exists()
        with contextlib.closing(sqlite3.connect(dest_db)) as conn:
            row = conn.execute(
                "SELECT remote_signal FROM database_versions WHERE name = 'alphamissense'"
            ).fetchone()
        assert row[0] == "etag:no-table"


class TestMultiAllelicMax:
    """MAX(am_pathogenicity) aggregation for multi-allelic sites sharing an rsID."""

    def test_bulk_lookup_returns_max_score(self, tmp_path: Path) -> None:
        db_path = tmp_path / "alphamissense.sqlite"
        conn = sqlite3.connect(db_path)
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores "
            "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
            "protein_variant, am_pathogenicity, am_class) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("1", 100, "A", "G", "rs5000", "P1", "E1", "A1G", 0.20, "likely_benign"),
                ("1", 100, "A", "T", "rs5000", "P1", "E1", "A1T", 0.92, "likely_pathogenic"),
                ("1", 100, "A", "C", "rs5000", "P1", "E1", "A1C", 0.45, "ambiguous"),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES ('alphamissense', 'test', '1.0', '2026-01-01', 3)",
        )
        conn.commit()
        conn.close()

        am = AlphaMissenseAnnotator(tmp_path)
        result = am.bulk_lookup({"rs5000"})
        assert result["rs5000"][0] == 0.92
        am.close()

    def test_lookup_returns_max_score(self, tmp_path: Path) -> None:
        db_path = tmp_path / "alphamissense.sqlite"
        conn = sqlite3.connect(db_path)
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores "
            "(chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
            "protein_variant, am_pathogenicity, am_class) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("2", 200, "C", "T", "rs6000", "P2", "E2", "C2T", 0.10, "likely_benign"),
                ("2", 200, "C", "A", "rs6000", "P2", "E2", "C2A", 0.88, "likely_pathogenic"),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name, source_url, version, downloaded_at, record_count) "
            "VALUES ('alphamissense', 'test', '1.0', '2026-01-01', 2)",
        )
        conn.commit()
        conn.close()

        am = AlphaMissenseAnnotator(tmp_path)
        result = am.lookup("rs6000")
        assert result is not None
        assert result[0] == 0.88
        am.close()


class TestBulkLookupByPosition:
    """Position-keyed fallback for rsIDs whose AlphaMissense rsid index is sparse.

    Same shape as gnomAD: resolution via ClinVar (GH #8) populates a variant's
    rsid, but AlphaMissense's rsid index may not list it. PK lookup recovers
    the score in a second enrichment pass.
    """

    def test_known_position(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        try:
            result = am.bulk_lookup_by_position({("1", 100, "A", "G")})
            assert result == {("1", 100, "A", "G"): (0.95, "likely_pathogenic")}
        finally:
            am.close()

    def test_unknown_position(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        try:
            assert am.bulk_lookup_by_position({("1", 99_999, "A", "G")}) == {}
        finally:
            am.close()

    def test_empty_input(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        try:
            assert am.bulk_lookup_by_position(set()) == {}
        finally:
            am.close()

    def test_finds_row_with_null_rsid(self, tmp_path: Path) -> None:
        """A row that's null in the rsid column is still findable by PK.

        This is the case the position fallback exists to handle.
        """
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        try:
            result = am.bulk_lookup_by_position({("3", 400, "T", "C")})
            assert result == {("3", 400, "T", "C"): (0.80, "likely_pathogenic")}
        finally:
            am.close()

    def test_mixed_hits_and_misses(self, tmp_path: Path) -> None:
        _build_db(tmp_path)
        am = AlphaMissenseAnnotator(tmp_path)
        try:
            keys = {
                ("1", 100, "A", "G"),  # hit
                ("2", 300, "G", "A"),  # hit
                ("9", 999, "T", "C"),  # miss
            }
            result = am.bulk_lookup_by_position(keys)
            assert len(result) == 2
            assert ("1", 100, "A", "G") in result
            assert ("2", 300, "G", "A") in result
            assert ("9", 999, "T", "C") not in result
        finally:
            am.close()

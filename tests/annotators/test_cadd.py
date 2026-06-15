# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the CADD variant deleteriousness annotator."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.cadd import (
    _BULK_BATCH_SIZE,
    CADD_FULL_FILENAME,
    CADD_INDEL_FILENAME,
    CaddAnnotator,
)
from allelix.databases._versions import CADD_SCHEMA_VERSION
from allelix.databases.cadd_loader import CADD_DB_FILENAME
from allelix.databases.schema import CADD_SCHEMA
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def cadd_db(tmp_path: Path) -> Path:
    """Create a minimal CADD SQLite cache with test data."""
    db_path = tmp_path / CADD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in CADD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.executemany(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
            [
                ("1", 11796321, "G", "A", 24.3),
                ("22", 19963748, "G", "A", 15.7),
                ("19", 44908684, "T", "C", 33.0),
                ("17", 43093517, "A", "G", 28.5),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count,"
            "  local_version_tag)"
            " VALUES (?, ?, ?, ?, ?, ?)",
            (
                "cadd",
                "test://mock",
                "v1.7",
                "2026-01-01T00:00:00Z",
                4,
                f"sv:{CADD_SCHEMA_VERSION}",
            ),
        )
        conn.commit()
    return tmp_path


class TestSetupAndStatus:
    """Annotator lifecycle: ready, version, record_count, close."""

    def test_unconfigured_is_not_ready(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path)
        assert not annotator.is_ready()

    def test_configured_is_ready(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        assert annotator.is_ready()

    def test_version_returns_string(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        assert annotator.version() == "v1.7"

    def test_version_returns_none_when_missing(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path)
        assert annotator.version() is None

    def test_record_count(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        assert annotator.record_count() == 4

    def test_record_count_returns_none_when_missing(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path)
        assert annotator.record_count() is None

    def test_close_is_idempotent(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        annotator.close()
        annotator.close()

    def test_close_after_lookup(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        annotator.lookup("1", 11796321, "G", "A")
        annotator.close()
        assert annotator._conn is None

    def test_fetch_remote_signal_returns_none(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        assert annotator.fetch_remote_signal() is None

    def test_cached_remote_signal_returns_none(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        assert annotator.cached_remote_signal() is None


class TestAnnotateReturnsEmpty:
    """CADD does not participate in the per-variant annotation loop."""

    def test_annotate_returns_empty_list(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        v = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="G",
            allele2="A",
            build="GRCh38",
        )
        assert annotator.annotate(v) == []


class TestLookup:
    """Single-coordinate PHRED lookup."""

    def test_known_variant(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            assert annotator.lookup("1", 11796321, "G", "A") == pytest.approx(24.3)
        finally:
            annotator.close()

    def test_unknown_variant(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            assert annotator.lookup("1", 99999999, "A", "T") is None
        finally:
            annotator.close()

    def test_missing_db_raises(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path)
        with pytest.raises(FileNotFoundError, match="CADD cache not found"):
            annotator.lookup("1", 11796321, "G", "A")


class TestBulkLookup:
    """Batched coordinate lookup."""

    def test_all_found(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            keys = {("1", 11796321, "G", "A"), ("22", 19963748, "G", "A")}
            result = annotator.bulk_lookup(keys)
            assert result[("1", 11796321, "G", "A")] == pytest.approx(24.3)
            assert result[("22", 19963748, "G", "A")] == pytest.approx(15.7)
        finally:
            annotator.close()

    def test_partial_match(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            keys = {("1", 11796321, "G", "A"), ("1", 99999999, "A", "T")}
            result = annotator.bulk_lookup(keys)
            assert ("1", 11796321, "G", "A") in result
            assert ("1", 99999999, "A", "T") not in result
        finally:
            annotator.close()

    def test_empty_input(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            assert annotator.bulk_lookup(set()) == {}
        finally:
            annotator.close()

    def test_no_matches(self, cadd_db: Path) -> None:
        annotator = CaddAnnotator(cadd_db)
        try:
            result = annotator.bulk_lookup({("99", 1, "X", "Y")})
            assert result == {}
        finally:
            annotator.close()

    def test_batch_size_constant(self) -> None:
        """Pin: batch size stays within SQLite's variable limit."""
        assert _BULK_BATCH_SIZE == 900


class TestFullMode:
    """Full mode (tabix) status checks without pysam."""

    def test_full_mode_not_ready_without_files(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert not annotator.is_ready()

    def test_full_mode_not_ready_with_gz_only(self, tmp_path: Path) -> None:
        (tmp_path / CADD_FULL_FILENAME).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert not annotator.is_ready()

    def test_full_mode_ready_with_both_files(self, tmp_path: Path) -> None:
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert annotator.is_ready()

    def test_full_mode_version(self, tmp_path: Path) -> None:
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert annotator.version() == "v1.7 (full)"

    def test_full_mode_version_not_ready(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert annotator.version() is None

    def test_full_mode_record_count_returns_none(self, tmp_path: Path) -> None:
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        assert annotator.record_count() is None

    def test_full_mode_missing_pysam_raises(self, tmp_path: Path, monkeypatch) -> None:
        """Without pysam installed, _open_tabix raises ImportError."""
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)

        bi = __builtins__
        original_import = bi.__import__ if hasattr(bi, "__import__") else __import__

        def mock_import(name, *args, **kwargs):
            if name == "pysam":
                raise ImportError("No module named 'pysam'")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", mock_import)
        with pytest.raises(ImportError, match="pysam"):
            annotator._open_tabix()


class TestLicenseDescriptor:
    """License metadata on the annotator class."""

    def test_commercial_ok_is_false(self) -> None:
        assert CaddAnnotator.license.commercial_ok is False

    def test_spdx_is_cadd(self) -> None:
        assert CaddAnnotator.license.spdx == "LicenseRef-CADD"

    def test_license_url(self) -> None:
        assert "cadd.gs.washington.edu" in CaddAnnotator.license.license_url


class TestSetup:
    """Download + install path (mocked network)."""

    def test_setup_downloads_and_installs(self, tmp_path: Path, monkeypatch) -> None:
        """setup() downloads, optionally verifies hash, installs, and cleans up."""
        calls: dict[str, list] = {"download": [], "verify": [], "install": []}

        def fake_download(url, dest):
            calls["download"].append((url, dest))
            dest.touch()

        def fake_verify(path, algo, expected):
            calls["verify"].append((path, algo, expected))

        def fake_install(gz, db, *, source_url="", remote_signal=None):
            calls["install"].append((gz, db))
            db.touch()

        monkeypatch.setattr("allelix.annotators.cadd.download", fake_download)
        monkeypatch.setattr("allelix.annotators.cadd.verify_file_hash", fake_verify)
        monkeypatch.setattr("allelix.annotators.cadd.install_prebuilt_cache", fake_install)
        monkeypatch.setattr("allelix.annotators.cadd.CADD_EXPECTED_SHA256", "abc123")

        annotator = CaddAnnotator(tmp_path)
        annotator.setup()

        assert len(calls["download"]) == 1
        assert len(calls["verify"]) == 1
        assert len(calls["install"]) == 1

    def test_setup_always_verifies_hash(self, tmp_path: Path, monkeypatch) -> None:
        """setup() always calls verify_file_hash unconditionally."""
        calls: dict[str, list] = {"verify": []}

        def fake_download(url, dest):
            dest.touch()

        def fake_verify(path, algo, expected):
            calls["verify"].append(True)

        def fake_install(gz, db, *, source_url="", remote_signal=None):
            db.touch()

        monkeypatch.setattr("allelix.annotators.cadd.download", fake_download)
        monkeypatch.setattr("allelix.annotators.cadd.verify_file_hash", fake_verify)
        monkeypatch.setattr("allelix.annotators.cadd.install_prebuilt_cache", fake_install)

        annotator = CaddAnnotator(tmp_path)
        annotator.setup()
        assert len(calls["verify"]) == 1

    def test_setup_tolerates_unlink_failure(self, tmp_path: Path, monkeypatch) -> None:
        """setup() logs warning if staged gz can't be removed."""

        def fake_download(url, dest):
            dest.touch()

        def fake_install(gz, db, *, source_url="", remote_signal=None):
            db.touch()

        def fake_unlink(self):
            raise OSError("permission denied")

        monkeypatch.setattr("allelix.annotators.cadd.download", fake_download)
        monkeypatch.setattr("allelix.annotators.cadd.verify_file_hash", lambda *a: None)
        monkeypatch.setattr("allelix.annotators.cadd.install_prebuilt_cache", fake_install)
        monkeypatch.setattr("pathlib.Path.unlink", fake_unlink)

        annotator = CaddAnnotator(tmp_path)
        annotator.setup()


class TestFullModeTabix:
    """Full mode tabix query paths using a mock tabix object."""

    @pytest.fixture()
    def full_annotator(self, tmp_path: Path) -> CaddAnnotator:
        """CaddAnnotator in full mode with mocked tabix."""
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)

        class MockTabixFile:
            def fetch(self, chrom, start, end):
                if chrom == "1" and start == 11796320 and end == 11796321:
                    return ["1\t11796321\tG\tA\t1.23\t24.3"]
                if chrom == "99":
                    raise ValueError("invalid contig")
                return []

            def close(self):
                pass

        annotator._tabix = MockTabixFile()
        return annotator

    def test_tabix_lookup_hit(self, full_annotator: CaddAnnotator) -> None:
        score = full_annotator.lookup("1", 11796321, "G", "A")
        assert score == pytest.approx(24.3)

    def test_tabix_lookup_miss(self, full_annotator: CaddAnnotator) -> None:
        assert full_annotator.lookup("2", 999, "A", "T") is None

    def test_tabix_lookup_wrong_alleles(self, full_annotator: CaddAnnotator) -> None:
        assert full_annotator.lookup("1", 11796321, "G", "T") is None

    def test_tabix_lookup_strips_chr_prefix(self, full_annotator: CaddAnnotator) -> None:
        score = full_annotator.lookup("chr1", 11796321, "G", "A")
        assert score == pytest.approx(24.3)

    def test_tabix_lookup_value_error(self, full_annotator: CaddAnnotator) -> None:
        assert full_annotator.lookup("99", 100, "A", "G") is None

    def test_tabix_bulk_lookup(self, full_annotator: CaddAnnotator) -> None:
        keys = {("1", 11796321, "G", "A"), ("2", 999, "A", "T")}
        result = full_annotator.bulk_lookup(keys)
        assert ("1", 11796321, "G", "A") in result
        assert result[("1", 11796321, "G", "A")] == pytest.approx(24.3)
        assert ("2", 999, "A", "T") not in result

    def test_tabix_bulk_lookup_empty(self, full_annotator: CaddAnnotator) -> None:
        assert full_annotator.bulk_lookup(set()) == {}

    def test_tabix_routes_indel_to_indel_file(self, tmp_path: Path) -> None:
        """Indel lookups query the indel tabix, not the SNV tabix."""
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        (tmp_path / CADD_INDEL_FILENAME).touch()
        (tmp_path / (CADD_INDEL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)

        class MockSNVTabix:
            def fetch(self, chrom, start, end):
                return []

            def close(self):
                pass

        class MockIndelTabix:
            def fetch(self, chrom, start, end):
                if chrom == "1" and start == 99 and end == 100:
                    return ["1\t100\tAC\tA\t0.8\t12.5"]
                return []

            def close(self):
                pass

        annotator._tabix = MockSNVTabix()
        annotator._indel_tabix = MockIndelTabix()

        assert annotator.lookup("1", 100, "AC", "A") == pytest.approx(12.5)
        assert annotator.lookup("1", 100, "A", "G") is None

    def test_indel_without_indel_file_returns_none(self, full_annotator: CaddAnnotator) -> None:
        """Indel lookup with no indel tabix file returns None."""
        assert full_annotator.lookup("1", 100, "AC", "A") is None

    def test_close_with_tabix(self, full_annotator: CaddAnnotator) -> None:
        assert full_annotator._tabix is not None
        full_annotator.close()
        assert full_annotator._tabix is None

    def test_close_clears_indel_tabix(self, tmp_path: Path) -> None:
        """close() clears both SNV and indel tabix handles."""
        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)

        class MockTabix:
            def close(self):
                pass

        annotator._tabix = MockTabix()
        annotator._indel_tabix = MockTabix()
        annotator.close()
        assert annotator._tabix is None
        assert annotator._indel_tabix is None

    def test_open_tabix_missing_file_raises(self, tmp_path: Path, monkeypatch) -> None:
        """Tabix file not found raises FileNotFoundError."""
        import types

        fake_pysam = types.ModuleType("pysam")
        monkeypatch.setitem(__import__("sys").modules, "pysam", fake_pysam)

        (tmp_path / CADD_FULL_FILENAME).touch()
        (tmp_path / (CADD_FULL_FILENAME + ".tbi")).touch()
        annotator = CaddAnnotator(tmp_path, full_mode=True)
        (tmp_path / CADD_FULL_FILENAME).unlink()
        with pytest.raises(FileNotFoundError, match="CADD tabix file not found"):
            annotator._open_tabix()


class TestContextManager:
    """CaddAnnotator supports the context manager protocol."""

    def test_context_manager(self, cadd_db: Path) -> None:
        with CaddAnnotator(cadd_db) as annotator:
            assert annotator.lookup("1", 11796321, "G", "A") == pytest.approx(24.3)
        assert annotator._conn is None

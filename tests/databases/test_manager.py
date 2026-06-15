# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for VCF parsing and SQLite cache loading."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.databases import manager
from allelix.databases.manager import (
    _parse_info,
    _pick,
    get_database_info,
    iter_clinvar_records,
    load_clinvar_vcf,
    parse_clinvar_version,
    stamp_existing_clinvar_cache,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestPick:
    """Pin the defensive empty-list guard. Current callers never hit it (they
    pass `[""]` instead of `[]`), but the guard exists so a future caller that
    drops the convention degrades gracefully rather than IndexError-ing."""

    def test_empty_list_returns_empty_string(self):
        assert _pick([], 0) == ""
        assert _pick([], 5) == ""

    def test_in_range_returns_value(self):
        assert _pick(["a", "b", "c"], 1) == "b"

    def test_out_of_range_pads_with_last(self):
        assert _pick(["a", "b"], 5) == "b"


class TestParseInfo:
    def test_basic_kv_pairs(self):
        out = _parse_info("RS=12345;CLNSIG=Pathogenic;CLNDN=Some_disease")
        assert out["RS"] == "12345"
        assert out["CLNSIG"] == "Pathogenic"
        assert out["CLNDN"] == "Some_disease"

    def test_flag_without_value(self):
        out = _parse_info("RS=12345;SOMATIC;ALLELEID=999")
        assert out["RS"] == "12345"
        assert "SOMATIC" in out
        assert out["ALLELEID"] == "999"

    def test_empty_value(self):
        out = _parse_info("RS=12345;CLNSIG=")
        assert out["CLNSIG"] == ""


class TestIterClinvarRecords:
    def test_yields_expected_records(self, mock_clinvar_vcf: Path):
        records = list(iter_clinvar_records(mock_clinvar_vcf))
        # 11 single-allele records + 2 from one multi-allelic row = 13
        # (Round 23 / ADR-0021: rs104894490 NIPA1 added for strand-inversion pin.)
        assert len(records) == 13
        rsids = {r["rsid"] for r in records}
        assert "rs1801133" in rsids
        assert "rs80357906" in rsids
        assert "rs113993960" in rsids  # CFTR indel
        assert "rs104894490" in rsids  # NIPA1 strand-inversion regression

    def test_extracts_clinical_significance(self, mock_clinvar_vcf: Path):
        records = list(iter_clinvar_records(mock_clinvar_vcf))
        mthfr = next(r for r in records if r["rsid"] == "rs1801133")
        assert mthfr["clinical_significance"] == "Pathogenic"
        assert mthfr["gene"] == "MTHFR"

    def test_unescapes_underscores_in_condition(self, mock_clinvar_vcf: Path):
        records = list(iter_clinvar_records(mock_clinvar_vcf))
        brca = next(r for r in records if r["rsid"] == "rs80357906")
        assert "Hereditary breast" in brca["condition"]
        assert "_" not in brca["condition"]

    def test_skips_comment_lines(self, tmp_path: Path):
        f = tmp_path / "tiny.vcf"
        f.write_text(
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tA\tT\t.\t.\tRS=1;CLNSIG=Benign\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert len(records) == 1
        assert records[0]["rsid"] == "rs1"

    def test_skips_records_without_rs(self, tmp_path: Path):
        f = tmp_path / "no_rs.vcf"
        f.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tA\tT\t.\t.\tCLNSIG=Benign\n"
            "1\t200\t2\tC\tG\t.\t.\tRS=42;CLNSIG=Pathogenic\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert len(records) == 1
        assert records[0]["rsid"] == "rs42"

    def test_skips_short_row(self, tmp_path: Path):
        """A VCF row with fewer than 8 columns is malformed; skip with a warning."""
        f = tmp_path / "short.vcf"
        f.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tA\n"  # 4 columns
            "1\t200\t2\tC\tG\t.\t.\tRS=42;CLNSIG=Pathogenic\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert [r["rsid"] for r in records] == ["rs42"]

    def test_skips_invalid_position(self, tmp_path: Path):
        f = tmp_path / "bad_pos.vcf"
        f.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\tBAD\t1\tA\tT\t.\t.\tRS=1;CLNSIG=Benign\n"
            "1\t200\t2\tC\tG\t.\t.\tRS=42;CLNSIG=Pathogenic\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert [r["rsid"] for r in records] == ["rs42"]

    def test_multi_allelic_split(self, tmp_path: Path):
        """C-2: One multi-allelic ALT row yields one record per ALT, paired by index."""
        f = tmp_path / "multi.vcf"
        f.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tG\tA,T\t.\t.\tRS=42;CLNSIG=Pathogenic|Benign;"
            "CLNDN=Cancer|Healthy;ALLELEID=900|901;GENEINFO=GENE1:1\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert len(records) == 2
        assert records[0]["alt"] == "A"
        assert records[0]["clinical_significance"] == "Pathogenic"
        assert records[0]["condition"] == "Cancer"
        assert records[0]["allele_id"] == 900
        assert records[1]["alt"] == "T"
        assert records[1]["clinical_significance"] == "Benign"
        assert records[1]["condition"] == "Healthy"
        assert records[1]["allele_id"] == 901

    def test_multi_allelic_pads_when_clnsig_short(self, tmp_path: Path):
        """If CLNSIG has fewer tokens than ALTs, extras inherit the last value."""
        f = tmp_path / "padded.vcf"
        f.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tG\tA,T,C\t.\t.\tRS=42;CLNSIG=Pathogenic\n",
            encoding="utf-8",
        )
        records = list(iter_clinvar_records(f))
        assert [r["alt"] for r in records] == ["A", "T", "C"]
        assert all(r["clinical_significance"] == "Pathogenic" for r in records)


class TestParseClinvarVersion:
    def test_extracts_filedate(self, mock_clinvar_vcf: Path):
        version = parse_clinvar_version(mock_clinvar_vcf)
        assert version == "20260101"

    def test_returns_none_when_absent(self, tmp_path: Path):
        f = tmp_path / "no_date.vcf"
        f.write_text(
            "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            encoding="utf-8",
        )
        assert parse_clinvar_version(f) is None

    def test_returns_none_when_only_meta_lines(self, tmp_path: Path):
        """m-version-fallback: cover the fall-through return at end-of-file."""
        f = tmp_path / "only_meta.vcf"
        f.write_text(
            "##fileformat=VCFv4.1\n##source=test\n##reference=GRCh37\n",
            encoding="utf-8",
        )
        assert parse_clinvar_version(f) is None


class TestLoadClinvarVcf:
    def test_populates_sqlite(self, tmp_path: Path, mock_clinvar_vcf: Path):
        db = tmp_path / "clinvar.sqlite"
        count = load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://mock")
        assert count == 13
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT clinical_significance, gene FROM clinvar_variants WHERE rsid = ?",
                ("rs1801133",),
            ).fetchone()
        assert row == ("Pathogenic", "MTHFR")

    def test_records_database_version(self, tmp_path: Path, mock_clinvar_vcf: Path):
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://mock")
        info = get_database_info(db, "clinvar")
        assert info is not None
        assert info["source_url"] == "test://mock"
        assert info["record_count"] == 13
        # m-6: version comes from VCF ##fileDate, not today's date
        assert info["version"] == "20260101"
        assert info["downloaded_at"]

    def test_reload_replaces_existing_cache(self, tmp_path: Path, mock_clinvar_vcf: Path):
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://mock")
        load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://mock-v2")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0]
            assert count == 13  # not 26 — reload replaces, not appends.
            versions = conn.execute("SELECT COUNT(*) FROM database_versions").fetchone()[0]
            assert versions == 1

    def test_batched_insert_flushes(
        self, tmp_path: Path, mock_clinvar_vcf: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """W-1: pin that mid-iteration flushes actually happen, not only the trailing flush.

        With INSERT_BATCH_SIZE=3 and 13 records, executemany is called 5 times:
        4 mid-iteration flushes of size 3 and a trailing flush of size 1.
        A mutation that disables mid-iteration flushing collapses this to a
        single trailing call.

        We can't monkey-patch `sqlite3.Connection.executemany` directly (it's a
        C slot, read-only). Wrap the connection in a delegating proxy instead.
        """
        monkeypatch.setattr(manager, "INSERT_BATCH_SIZE", 3)
        executemany_payloads: list[int] = []
        real_connect = sqlite3.connect

        class _SpyConn:
            def __init__(self, real):
                self._real = real

            def executemany(self, sql, seq):
                seq_list = list(seq)
                executemany_payloads.append(len(seq_list))
                return self._real.executemany(sql, seq_list)

            def __getattr__(self, name):
                return getattr(self._real, name)

        def spying_connect(*args, **kwargs):
            return _SpyConn(real_connect(*args, **kwargs))

        monkeypatch.setattr(sqlite3, "connect", spying_connect)

        db = tmp_path / "clinvar.sqlite"
        count = load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://batch")
        assert count == 13
        # 13 records / batch_size 3 = 4 mid-iteration flushes of size 3 + a
        # trailing flush of size 1.
        assert executemany_payloads == [3, 3, 3, 3, 1]

    def test_atomic_load_preserves_old_cache_on_failure(
        self, tmp_path: Path, mock_clinvar_vcf: Path, monkeypatch: pytest.MonkeyPatch
    ):
        """M-1: If a reload fails mid-parse, the previous cache must survive."""
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://v1")
        info_before = get_database_info(db, "clinvar")
        assert info_before is not None
        original = manager.iter_clinvar_records

        def failing_iter(path):
            for i, record in enumerate(original(path)):
                if i >= 2:
                    raise RuntimeError("simulated parse failure")
                yield record

        monkeypatch.setattr(manager, "iter_clinvar_records", failing_iter)

        with pytest.raises(RuntimeError, match="simulated parse failure"):
            load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://v2")

        info_after = get_database_info(db, "clinvar")
        assert info_after == info_before
        # No leftover .tmp file
        assert not (tmp_path / "clinvar.sqlite.tmp").exists()


class TestGetDatabaseInfo:
    def test_missing_file_returns_none(self, tmp_path: Path):
        assert get_database_info(tmp_path / "nope.sqlite", "clinvar") is None

    def test_unknown_database_returns_none(self, tmp_path: Path, mock_clinvar_vcf: Path):
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_vcf(mock_clinvar_vcf, db, source_url="test://mock")
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

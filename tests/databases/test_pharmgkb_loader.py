# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for ClinPGx TSV parsing and SQLite cache loading."""

from __future__ import annotations

import contextlib
import sqlite3
import zipfile
from typing import TYPE_CHECKING

import pytest

from allelix.databases import pharmgkb_loader
from allelix.databases.manager import get_database_info
from allelix.databases.pharmgkb_loader import (
    FUNCTION_CLASS_DECREASED,
    FUNCTION_CLASS_NO_FUNCTION,
    FUNCTION_CLASS_NORMAL,
    FUNCTION_CLASS_UNKNOWN,
    _is_single_rsid,
    _normalize_genotype,
    _safe_float,
    classify_function,
    is_nonfinding,
    is_nonfinding_by_allele_lookup,
    is_nonfinding_for_row,
    iter_pharmgkb_records,
    load_pharmgkb_tsv,
    schema_is_current,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestNormalizeGenotype:
    def test_two_letter_sorted(self):
        assert _normalize_genotype("AG") == "AG"
        assert _normalize_genotype("GA") == "AG"
        assert _normalize_genotype("TC") == "CT"

    def test_strips_separators(self):
        assert _normalize_genotype("A:G") == "AG"
        assert _normalize_genotype("A;G") == "AG"
        assert _normalize_genotype("A/G") == "AG"

    def test_lowercase_normalized_to_upper(self):
        assert _normalize_genotype("ag") == "AG"

    def test_indel_returns_none(self):
        assert _normalize_genotype("CTT/C") is None
        assert _normalize_genotype("ins") is None

    def test_star_allele_returns_none(self):
        assert _normalize_genotype("*1/*2") is None

    def test_invalid_letter_returns_none(self):
        assert _normalize_genotype("AX") is None
        assert _normalize_genotype("NN") is None


class TestStructuredFunctionClassifier:
    """ADR-0016: classification reads ClinPGx's structured `Allele Function`
    column, never regex on annotation text. The earlier regex-based filter
    was architecturally wrong and is now in a test safety net only.
    """

    def test_normal_function(self):
        assert classify_function("Normal function") == FUNCTION_CLASS_NORMAL
        assert classify_function("normal function") == FUNCTION_CLASS_NORMAL
        assert classify_function("  Normal function  ") == FUNCTION_CLASS_NORMAL

    def test_decreased_function(self):
        assert classify_function("Decreased function") == FUNCTION_CLASS_DECREASED

    def test_no_function(self):
        assert classify_function("No function") == "no_function"

    def test_increased_function(self):
        assert classify_function("Increased function") == "increased"

    def test_empty_field_is_unknown(self):
        assert classify_function("") == FUNCTION_CLASS_UNKNOWN
        assert classify_function(None) == FUNCTION_CLASS_UNKNOWN

    def test_unrecognized_value_is_unknown(self):
        assert classify_function("Function not assigned") == FUNCTION_CLASS_UNKNOWN
        assert classify_function("garbage") == FUNCTION_CLASS_UNKNOWN


class TestIsNonfindingStructured:
    """`is_nonfinding(function_class)` is the pure structured check.

    Used in tests and back-compat shims. Production code uses
    `is_nonfinding_for_row` which also handles the empty-field fallback.
    """

    def test_normal_is_nonfinding(self):
        assert is_nonfinding(FUNCTION_CLASS_NORMAL) is True

    def test_decreased_is_finding(self):
        assert is_nonfinding(FUNCTION_CLASS_DECREASED) is False

    def test_no_function_is_finding(self):
        assert is_nonfinding("no_function") is False

    def test_increased_is_finding(self):
        assert is_nonfinding("increased") is False

    def test_unknown_is_finding(self):
        assert is_nonfinding(FUNCTION_CLASS_UNKNOWN) is False


class TestIsNonfindingByAlleleLookup:
    """ADR-0020 per-allele structured classifier (the filter, as a join)."""

    def test_homozygous_reference_is_nonfinding(self):
        # The exact production leakers from v0.7.0 + v0.8.0.
        lookup = {
            ("rs1800559", "C"): FUNCTION_CLASS_NORMAL,
            ("rs1800559", "T"): FUNCTION_CLASS_DECREASED,
            ("rs116855232", "C"): FUNCTION_CLASS_NORMAL,
            ("rs116855232", "T"): FUNCTION_CLASS_NO_FUNCTION,
        }
        assert is_nonfinding_by_allele_lookup("rs1800559", "CC", lookup) is True
        assert is_nonfinding_by_allele_lookup("rs116855232", "CC", lookup) is True

    def test_heterozygous_carrier_is_finding(self):
        lookup = {
            ("rs1801133", "G"): FUNCTION_CLASS_NORMAL,
            ("rs1801133", "A"): FUNCTION_CLASS_DECREASED,
        }
        assert is_nonfinding_by_allele_lookup("rs1801133", "AG", lookup) is False

    def test_homozygous_variant_is_finding(self):
        lookup = {
            ("rs1801133", "G"): FUNCTION_CLASS_NORMAL,
            ("rs1801133", "A"): FUNCTION_CLASS_DECREASED,
        }
        assert is_nonfinding_by_allele_lookup("rs1801133", "AA", lookup) is False

    def test_no_function_allele_in_genotype_is_finding(self):
        lookup = {
            ("rs116855232", "C"): FUNCTION_CLASS_NORMAL,
            ("rs116855232", "T"): FUNCTION_CLASS_NO_FUNCTION,
        }
        assert is_nonfinding_by_allele_lookup("rs116855232", "CT", lookup) is False
        assert is_nonfinding_by_allele_lookup("rs116855232", "TT", lookup) is False

    def test_dpyd_cluster_reference_homozygotes(self):
        # The v0.6.1 reviewer's cited leakers — must be non-findings.
        lookup = {
            ("rs115232898", "T"): FUNCTION_CLASS_NORMAL,
            ("rs115232898", "C"): FUNCTION_CLASS_DECREASED,
            ("rs1801266", "G"): FUNCTION_CLASS_NORMAL,
            ("rs1801266", "A"): FUNCTION_CLASS_NO_FUNCTION,
            ("rs3918290", "C"): FUNCTION_CLASS_NORMAL,
            ("rs3918290", "T"): FUNCTION_CLASS_NO_FUNCTION,
        }
        assert is_nonfinding_by_allele_lookup("rs115232898", "TT", lookup) is True
        assert is_nonfinding_by_allele_lookup("rs1801266", "GG", lookup) is True
        assert is_nonfinding_by_allele_lookup("rs3918290", "CC", lookup) is True

    def test_rsid_absent_returns_none(self):
        """Tier 2 abstains for rsids not in CPIC; caller treats as finding."""
        lookup = {("rs1800559", "C"): FUNCTION_CLASS_NORMAL}
        assert is_nonfinding_by_allele_lookup("rs_unknown", "AA", lookup) is None

    def test_unclassified_base_in_known_rsid_is_finding(self):
        """ADR-0020: an allele missing from the lookup at an rsid that HAS
        entries is uncharacterized — never silently suppressed.
        """
        lookup = {("rs1", "C"): FUNCTION_CLASS_NORMAL}
        # User has CA: C is Normal, A is unclassified → emit as finding.
        assert is_nonfinding_by_allele_lookup("rs1", "AC", lookup) is False


class TestIsNonfindingForRow:
    """ADR-0020 (v0.9.0) filter: structured Allele Function → CPIC lookup → emit."""

    def test_structured_normal_wins_regardless_of_text(self):
        assert is_nonfinding_for_row("Normal function") is True

    def test_structured_decreased_wins(self):
        assert is_nonfinding_for_row("Decreased function") is False

    def test_lookup_decides_when_structured_field_empty(self):
        lookup = {
            ("rs1800559", "C"): FUNCTION_CLASS_NORMAL,
            ("rs1800559", "T"): FUNCTION_CLASS_DECREASED,
        }
        assert (
            is_nonfinding_for_row(
                "", rsid="rs1800559", genotype="CC", allele_function_lookup=lookup
            )
            is True
        )
        assert (
            is_nonfinding_for_row(
                "", rsid="rs1800559", genotype="CT", allele_function_lookup=lookup
            )
            is False
        )

    def test_no_data_emits(self):
        """When both tiers have no data, the row emits (never silently suppress)."""
        assert is_nonfinding_for_row("") is False
        assert is_nonfinding_for_row(None) is False
        assert (
            is_nonfinding_for_row("", rsid="rs999", genotype="AA", allele_function_lookup={})
            is False
        )


class TestSchemaIsCurrent:
    """ADR-0016 schema migration: pre-v0.6.0 caches lack `function_class`."""

    def test_returns_false_for_missing_file(self, tmp_path: Path):
        assert not schema_is_current(tmp_path / "missing.sqlite")

    def test_returns_true_for_v060_cache(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        db = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test")
        assert schema_is_current(db)

    def test_returns_false_for_v05x_cache(self, tmp_path: Path):
        """A v0.5.x cache has `is_nonfinding` + `is_somatic` but no `function_class`."""
        db = tmp_path / "legacy.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                """
                CREATE TABLE pharmgkb_annotations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    rsid TEXT NOT NULL, genotype TEXT NOT NULL, gene TEXT,
                    drugs TEXT, phenotype TEXT, phenotype_category TEXT,
                    annotation_text TEXT, level_of_evidence TEXT, score REAL,
                    pgkb_annotation_id TEXT, allele_function TEXT,
                    is_nonfinding INTEGER NOT NULL, is_somatic INTEGER NOT NULL
                );
                """
            )
            conn.commit()
        assert not schema_is_current(db)

    def test_returns_false_for_v04x_cache(self, tmp_path: Path):
        """A v0.4.x cache lacks both `function_class` AND `is_nonfinding`."""
        db = tmp_path / "legacy.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                """
                CREATE TABLE pharmgkb_annotations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    rsid TEXT NOT NULL, genotype TEXT NOT NULL, gene TEXT,
                    drugs TEXT, phenotype TEXT, phenotype_category TEXT,
                    annotation_text TEXT, level_of_evidence TEXT, score REAL,
                    pgkb_annotation_id TEXT
                );
                """
            )
            conn.commit()
        assert not schema_is_current(db)


class TestSafeFloat:
    def test_valid_float(self):
        assert _safe_float("0.85") == 0.85

    def test_empty_returns_none(self):
        assert _safe_float("") is None

    def test_garbage_returns_none(self):
        assert _safe_float("not-a-number") is None


class TestIsSingleRsid:
    def test_valid_rsid(self):
        assert _is_single_rsid("rs1801133")
        assert _is_single_rsid("rs999999999")

    def test_rejects_star_allele(self):
        assert not _is_single_rsid("CYP2D6*1/*2")

    def test_rejects_multi_rsid(self):
        assert not _is_single_rsid("rs1801133, rs1801131")

    def test_rejects_blank(self):
        assert not _is_single_rsid("")
        assert not _is_single_rsid("   ")


class TestIterPharmgkbRecords:
    def test_yields_expected_records(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        assert len(records) == 16

    def test_function_class_populated_from_allele_function(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        for r in records:
            assert r["function_class"] in {
                FUNCTION_CLASS_NORMAL,
                FUNCTION_CLASS_DECREASED,
                "no_function",
                "increased",
                FUNCTION_CLASS_UNKNOWN,
            }

    def test_is_nonfinding_uses_structured_lookup(
        self, mock_pharmgkb_dir: Path, mock_cpic_lookup: dict[tuple[str, str], str]
    ):
        """ADR-0020: is_nonfinding decided by structured allele_function_lookup."""
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir, mock_cpic_lookup))
        by_key = {(r["rsid"], r["genotype"]): bool(r["is_nonfinding"]) for r in records}
        # Reference homozygotes → non-finding.
        assert by_key[("rs1801133", "GG")] is True
        assert by_key[("rs4680", "GG")] is True
        assert by_key[("rs900000010", "GG")] is True
        # Carriers → finding.
        assert by_key[("rs1801133", "AG")] is False
        assert by_key[("rs1801133", "AA")] is False
        assert by_key[("rs4680", "AA")] is False

    def test_skips_star_allele_annotation(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        ids = {r["pgkb_annotation_id"] for r in records}
        assert "PA-005" not in ids

    def test_skips_multi_rsid_annotation(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        ids = {r["pgkb_annotation_id"] for r in records}
        assert "PA-006" not in ids

    def test_inner_normalize_skip_for_non_snv_alleles(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        ids = {r["pgkb_annotation_id"] for r in records}
        assert "PA-007" not in ids

    def test_genotypes_are_normalized(self, mock_pharmgkb_dir: Path):
        records = list(iter_pharmgkb_records(mock_pharmgkb_dir))
        rs1801133 = [r for r in records if r["rsid"] == "rs1801133"]
        genotypes = {r["genotype"] for r in rs1801133}
        assert genotypes == {"AG", "AA", "GG"}

    def test_handles_zip_input(self, mock_pharmgkb_dir: Path, tmp_path: Path):
        zip_path = tmp_path / "clinicalAnnotations.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            for f in mock_pharmgkb_dir.iterdir():
                zf.write(f, arcname=f.name)
        records = list(iter_pharmgkb_records(zip_path))
        assert len(records) == 16

    def test_missing_files_raise(self, tmp_path: Path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="missing required TSVs"):
            list(iter_pharmgkb_records(empty))


class TestLoadPharmgkbTsv:
    def test_populates_sqlite(
        self,
        tmp_path: Path,
        mock_pharmgkb_dir: Path,
        mock_cpic_lookup: dict[tuple[str, str], str],
    ):
        db = tmp_path / "pharmgkb.sqlite"
        count = load_pharmgkb_tsv(
            mock_pharmgkb_dir,
            db,
            source_url="test://mock",
            allele_function_lookup=mock_cpic_lookup,
        )
        assert count == 16
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT gene, drugs, level_of_evidence, function_class "
                "FROM pharmgkb_annotations WHERE rsid = ? AND genotype = ?",
                ("rs1801133", "AG"),
            ).fetchone()
        # SNV rows have Allele Function = "" in real ClinPGx; the
        # structured function_class column reflects that (`unknown`).
        assert row == ("MTHFR", "methotrexate", "2A", FUNCTION_CLASS_UNKNOWN)

    def test_records_database_version(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        db = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://mock", version="2026-test")
        info = get_database_info(db, "pharmgkb")
        assert info is not None
        assert info["version"] == "2026-test"
        assert info["record_count"] == 16
        assert info["source_url"] == "test://mock"

    def test_batched_insert_flushes(
        self, tmp_path: Path, mock_pharmgkb_dir: Path, monkeypatch: pytest.MonkeyPatch
    ):
        monkeypatch.setattr(pharmgkb_loader, "INSERT_BATCH_SIZE", 3)
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

        db = tmp_path / "pharmgkb.sqlite"
        count = load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://batch")
        assert count == 16
        # 16 records / batch_size 3 = 5 full batches + 1 remainder.
        assert executemany_payloads == [3, 3, 3, 3, 3, 1]

    def test_leftover_tmp_file_is_cleared(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        db = tmp_path / "pharmgkb.sqlite"
        stale_tmp = tmp_path / "pharmgkb.sqlite.tmp"
        stale_tmp.write_bytes(b"\x00\x01garbage\x02\x03")
        assert stale_tmp.exists()

        count = load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://restart")
        assert count == 16
        assert db.exists()
        assert not stale_tmp.exists()

    def test_atomic_load_preserves_old_cache_on_failure(
        self,
        tmp_path: Path,
        mock_pharmgkb_dir: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        db = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://v1")
        info_before = get_database_info(db, "pharmgkb")

        original = pharmgkb_loader.iter_pharmgkb_records

        def failing_iter(path, lookup=None):
            for i, record in enumerate(original(path, lookup)):
                if i >= 2:
                    raise RuntimeError("simulated parse failure")
                yield record

        monkeypatch.setattr(pharmgkb_loader, "iter_pharmgkb_records", failing_iter)

        with pytest.raises(RuntimeError, match="simulated parse failure"):
            load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://v2")

        info_after = get_database_info(db, "pharmgkb")
        assert info_after == info_before
        assert not (tmp_path / "pharmgkb.sqlite.tmp").exists()


class TestMinRowsFloor:
    """GH #19: post-load row-count sanity gate catches mid-stream truncation.

    Production callers pass ``min_rows=PHARMGKB_MIN_ROWS``. A truncated
    download produces a parseable-but-short cache that previously committed
    silently. The floor check raises before ``os.replace``.
    """

    def test_below_floor_raises_and_keeps_cache_clean(self, mock_pharmgkb_dir, tmp_path) -> None:
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db = tmp_path / "pharmgkb.sqlite"
        with pytest.raises(OSError, match="rows ingested"):
            load_pharmgkb_tsv(
                mock_pharmgkb_dir,
                db,
                source_url="test://truncated",
                min_rows=999_999,
            )
        assert not db.exists()
        assert not (tmp_path / "pharmgkb.sqlite.tmp").exists()

    def test_default_min_rows_zero_loads_small_fixtures(self, mock_pharmgkb_dir, tmp_path) -> None:
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db = tmp_path / "pharmgkb.sqlite"
        count = load_pharmgkb_tsv(mock_pharmgkb_dir, db, source_url="test://small")
        assert count > 0
        assert db.exists()

    def test_at_or_above_floor_succeeds(self, mock_pharmgkb_dir, tmp_path) -> None:
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db = tmp_path / "pharmgkb.sqlite"
        count = load_pharmgkb_tsv(
            mock_pharmgkb_dir,
            db,
            source_url="test://at-floor",
            min_rows=1,
        )
        assert count >= 1
        assert db.exists()

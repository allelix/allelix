# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for scripts/build_cadd_cache.py helper functions."""

from __future__ import annotations

import contextlib
import gzip
import sqlite3
import sys
from pathlib import Path

import pytest

from allelix.databases.schema import ALPHAMISSENSE_SCHEMA, CLINVAR_SCHEMA, GNOMAD_SCHEMA

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from build_cadd_cache import (  # noqa: I001
    _BATCH_SIZE,
    _load_position_set,
    _normalize_chrom,
    _pack,
    build_cache,
)


# ── Helpers ──────────────────────────────────────────────────────────


def _make_gnomad_db(path: Path, rows: list[tuple[str, int, str, str]]) -> None:
    """Create a minimal gnomAD SQLite with the given coordinate rows."""
    db_path = path / "gnomad.sqlite"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        for chrom, pos, ref, alt in rows:
            conn.execute(
                "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (chrom, pos, ref, alt, f"rs{pos}", 0.01),
            )
        conn.commit()


def _make_alphamissense_db(path: Path, rows: list[tuple[str, int, str, str]]) -> None:
    """Create a minimal AlphaMissense SQLite."""
    db_path = path / "alphamissense.sqlite"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in ALPHAMISSENSE_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        for chrom, pos, ref, alt in rows:
            conn.execute(
                "INSERT INTO alphamissense_scores"
                " (chrom, pos, ref, alt, am_pathogenicity, am_class)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (chrom, pos, ref, alt, 0.5, "ambiguous"),
            )
        conn.commit()


def _make_clinvar_db(path: Path, rows: list[tuple[str, int, str, str]]) -> None:
    """Create a minimal ClinVar GRCh38 SQLite."""
    db_path = path / "clinvar.GRCh38.sqlite"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in CLINVAR_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        for chrom, pos, ref, alt in rows:
            conn.execute(
                "INSERT INTO clinvar_variants"
                " (rsid, chromosome, position, ref, alt, clinical_significance)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (f"rs{pos}", chrom, pos, ref, alt, "Pathogenic"),
            )
        conn.commit()


_CADD_HEADER = "## CADD GRCh38-v1.7\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n"


def _make_cadd_line(
    chrom: str = "1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "G",
    raw: str = "1.234",
    phred: str = "15.5",
) -> str:
    return f"{chrom}\t{pos}\t{ref}\t{alt}\t{raw}\t{phred}\n"


def _make_cadd_gz(
    tmp_path: Path, lines: list[str], *, filename: str = "whole_genome_SNVs.tsv.gz"
) -> Path:
    """Write a gzipped CADD file from text lines."""
    gz_path = tmp_path / filename
    with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
        fh.write("".join(lines))
    return gz_path


# ── _pack ────────────────────────────────────────────────────────────


class TestPack:
    """Tests for int64 packing of SNV keys."""

    def test_known_value(self) -> None:
        result = _pack("1", 100, "A", "G")
        assert result is not None
        assert result == (1 << 34) | (100 << 4) | (0 << 2) | 2

    def test_unmapped_contig_returns_none(self) -> None:
        assert _pack("GL000220.1", 100, "A", "G") is None

    def test_non_acgt_ref_returns_none(self) -> None:
        assert _pack("1", 100, "N", "G") is None

    def test_non_acgt_alt_returns_none(self) -> None:
        assert _pack("1", 100, "A", "N") is None

    def test_mt_and_m_equivalent(self) -> None:
        assert _pack("MT", 100, "A", "G") == _pack("M", 100, "A", "G")

    def test_position_overflow(self) -> None:
        with pytest.raises(ValueError, match="exceeds 30-bit budget"):
            _pack("1", 2**30, "A", "G")

    def test_x_chrom(self) -> None:
        result = _pack("X", 500, "C", "T")
        assert result is not None
        assert result == (23 << 34) | (500 << 4) | (1 << 2) | 3

    def test_all_nucleotide_combos(self) -> None:
        for ref_nuc, ref_idx in [("A", 0), ("C", 1), ("G", 2), ("T", 3)]:
            for alt_nuc, alt_idx in [("A", 0), ("C", 1), ("G", 2), ("T", 3)]:
                result = _pack("1", 1, ref_nuc, alt_nuc)
                assert result == (1 << 34) | (1 << 4) | (ref_idx << 2) | alt_idx


# ── _normalize_chrom ─────────────────────────────────────────────────


class TestNormalizeChrom:
    """Tests for contig normalization."""

    def test_strips_chr_prefix(self) -> None:
        assert _normalize_chrom("chr1") == "1"

    def test_m_to_mt(self) -> None:
        assert _normalize_chrom("M") == "MT"

    def test_chr_m_to_mt(self) -> None:
        assert _normalize_chrom("chrM") == "MT"

    def test_passthrough(self) -> None:
        assert _normalize_chrom("22") == "22"


# ── _load_position_set ───────────────────────────────────────────────


class TestLoadPositionSet:
    """Tests for multi-source position loading."""

    def test_loads_gnomad_snvs(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G"), ("22", 200, "C", "T")])
        snv_keys, indel_keys = _load_position_set(tmp_path)
        assert len(snv_keys) == 2
        assert _pack("1", 100, "A", "G") in snv_keys
        assert len(indel_keys) == 0

    def test_loads_indels_from_gnomad(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "AC", "A")])
        snv_keys, indel_keys = _load_position_set(tmp_path)
        assert len(snv_keys) == 0
        assert ("1", 100, "AC", "A") in indel_keys

    def test_merges_alphamissense(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        _make_alphamissense_db(tmp_path, [("2", 300, "T", "C")])
        snv_keys, _ = _load_position_set(tmp_path)
        assert _pack("1", 100, "A", "G") in snv_keys
        assert _pack("2", 300, "T", "C") in snv_keys

    def test_merges_clinvar(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        _make_clinvar_db(tmp_path, [("3", 400, "G", "A")])
        snv_keys, _ = _load_position_set(tmp_path)
        assert _pack("3", 400, "G", "A") in snv_keys

    def test_missing_gnomad_exits(self, tmp_path: Path) -> None:
        with pytest.raises(SystemExit):
            _load_position_set(tmp_path)

    def test_missing_alphamissense_warns(self, tmp_path: Path, capsys) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        _load_position_set(tmp_path)
        captured = capsys.readouterr()
        assert "alphamissense.sqlite not found" in captured.out

    def test_empty_gnomad(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [])
        snv_keys, indel_keys = _load_position_set(tmp_path)
        assert len(snv_keys) == 0
        assert len(indel_keys) == 0

    def test_m_chrom_normalized(self, tmp_path: Path) -> None:
        """M in source is normalized to MT for packing."""
        _make_gnomad_db(tmp_path, [("M", 100, "A", "G")])
        snv_keys, _ = _load_position_set(tmp_path)
        assert _pack("MT", 100, "A", "G") in snv_keys


# ── build_cache integration ──────────────────────────────────────────


class TestBuildCache:
    """Integration tests: synthetic CADD gz + database SQLites → CADD SQLite."""

    def test_filters_to_database_positions(self, tmp_path: Path) -> None:
        """Only rows matching database positions are kept."""
        _make_gnomad_db(
            tmp_path,
            [("1", 100, "A", "G"), ("22", 200, "C", "T")],
        )
        cadd_lines = [
            _CADD_HEADER,
            _make_cadd_line("1", 100, "A", "G", phred="24.3"),
            _make_cadd_line("1", 999, "T", "C", phred="5.0"),
            _make_cadd_line("22", 200, "C", "T", phred="15.7"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT chrom, pos, ref, alt, phred FROM cadd_scores ORDER BY chrom, pos"
            ).fetchall()
        assert len(rows) == 2
        assert rows[0] == ("1", 100, "A", "G", pytest.approx(24.3))
        assert rows[1] == ("22", 200, "C", "T", pytest.approx(15.7))

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            "## CADD GRCh38-v1.7\n",
            "## Some other header\n",
            "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n",
            _make_cadd_line("1", 100, "A", "G", phred="20.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_skips_short_lines(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            "1\t100\tA\n",
            _make_cadd_line("1", 100, "A", "G", phred="18.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_skips_bad_position(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            "1\tXYZ\tA\tG\t0.5\t10.0\n",
            _make_cadd_line("1", 100, "A", "G", phred="22.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT phred FROM cadd_scores").fetchall()
        assert len(rows) == 1
        assert rows[0][0] == pytest.approx(22.0)

    def test_skips_bad_phred(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            "1\t100\tA\tG\t0.5\tNA\n",
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 0

    def test_stamps_database_version(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            _make_cadd_line("1", 100, "A", "G", phred="12.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT name, version, record_count, local_version_tag"
                " FROM database_versions WHERE name='cadd'"
            ).fetchone()
        assert row is not None
        assert row[0] == "cadd"
        assert row[1] == "v1.7"
        assert row[2] == 1
        assert row[3].startswith("sv:")

    def test_empty_positions_exits(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [])
        cadd_lines = [_CADD_HEADER, _make_cadd_line()]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        with pytest.raises(SystemExit):
            build_cache(snv_gz, None, output, tmp_path)

    def test_batch_size_constant(self) -> None:
        assert _BATCH_SIZE == 100_000

    def test_atomic_write(self, tmp_path: Path) -> None:
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [_CADD_HEADER, _make_cadd_line("1", 100, "A", "G")]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        assert output.exists()
        tmp_file = tmp_path / "cadd.sqlite.tmp"
        assert not tmp_file.exists()

    def test_indel_pass(self, tmp_path: Path) -> None:
        """Indel rows from the indel file are cached when they match."""
        _make_gnomad_db(tmp_path, [("1", 100, "AC", "A")])
        indel_lines = [
            _CADD_HEADER,
            "1\t100\tAC\tA\t0.8\t12.5\n",
            "1\t200\tGG\tG\t0.3\t5.0\n",
        ]
        indel_gz = _make_cadd_gz(tmp_path, indel_lines, filename="indel.tsv.gz")
        snv_lines = [_CADD_HEADER]
        snv_gz = _make_cadd_gz(tmp_path, snv_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, indel_gz, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT chrom, pos, ref, alt, phred FROM cadd_scores").fetchall()
        assert len(rows) == 1
        assert rows[0] == ("1", 100, "AC", "A", pytest.approx(12.5))

    def test_snv_and_indel_combined_count(self, tmp_path: Path) -> None:
        """Record count includes both SNV and indel matches."""
        _make_gnomad_db(
            tmp_path,
            [("1", 100, "A", "G"), ("2", 200, "AC", "A")],
        )
        snv_lines = [
            _CADD_HEADER,
            _make_cadd_line("1", 100, "A", "G", phred="24.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, snv_lines)
        indel_lines = [
            _CADD_HEADER,
            "2\t200\tAC\tA\t0.9\t18.0\n",
        ]
        indel_gz = _make_cadd_gz(tmp_path, indel_lines, filename="indel.tsv.gz")
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, indel_gz, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT record_count FROM database_versions WHERE name='cadd'"
            ).fetchone()
            total = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert row[0] == 2
        assert total == 2

    def test_chr_prefix_normalized(self, tmp_path: Path) -> None:
        """CADD rows with chr prefix match positions stored without prefix."""
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            _make_cadd_line("chr1", 100, "A", "G", phred="20.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_unmapped_contig_skipped(self, tmp_path: Path) -> None:
        """Alt contigs like GL000220.1 are silently skipped."""
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            "GL000220.1\t100\tA\tG\t0.5\t10.0\n",
            _make_cadd_line("1", 100, "A", "G", phred="22.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_non_acgt_allele_skipped(self, tmp_path: Path) -> None:
        """Non-ACGT alleles in CADD rows are silently skipped."""
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        cadd_lines = [
            _CADD_HEADER,
            "1\t100\tN\tG\t0.5\t10.0\n",
            _make_cadd_line("1", 100, "A", "G", phred="22.0"),
        ]
        snv_gz = _make_cadd_gz(tmp_path, cadd_lines)
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM cadd_scores").fetchone()[0]
        assert count == 1

    def test_no_indel_file(self, tmp_path: Path) -> None:
        """build_cache works with indel_path=None."""
        _make_gnomad_db(tmp_path, [("1", 100, "A", "G")])
        snv_gz = _make_cadd_gz(
            tmp_path,
            [_CADD_HEADER, _make_cadd_line("1", 100, "A", "G")],
        )
        output = tmp_path / "cadd.sqlite"

        build_cache(snv_gz, None, output, tmp_path)

        assert output.exists()

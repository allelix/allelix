# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for scripts/build_alphamissense_cache.py helper functions."""

from __future__ import annotations

import contextlib
import gzip
import io
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from build_alphamissense_cache import (  # noqa: I001
    _load_gnomad_rsid_map,
    _open_tsv,
    build_cache,
)

from allelix.databases.schema import GNOMAD_SCHEMA


# ── Test TSV data ─────────────────────────────────────────────────────

_TSV_HEADER = (
    "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\t"
    "protein_variant\tam_pathogenicity\tam_class\n"
)


def _make_tsv_line(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "G",
    genome: str = "hg38",
    uniprot_id: str = "P12345",
    transcript_id: str = "ENST001",
    protein_variant: str = "A100G",
    am_pathogenicity: str = "0.95",
    am_class: str = "likely_pathogenic",
) -> str:
    """Build one AlphaMissense TSV line."""
    return (
        f"{chrom}\t{pos}\t{ref}\t{alt}\t{genome}\t{uniprot_id}\t"
        f"{transcript_id}\t{protein_variant}\t{am_pathogenicity}\t{am_class}\n"
    )


def _make_gzipped_tsv(lines: list[str]) -> bytes:
    """Produce a gzip-compressed TSV from text lines."""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write("".join(lines).encode("utf-8"))
    return buf.getvalue()


def _build_gnomad_db(path: Path, rows: list[tuple]) -> Path:
    """Build a minimal gnomAD SQLite cache for rsID join testing."""
    conn = sqlite3.connect(path)
    conn.executescript(GNOMAD_SCHEMA)
    conn.executemany(
        "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.execute(
        "INSERT INTO database_versions (name, source_url, version, downloaded_at, record_count) "
        "VALUES ('gnomad', 'test', '4.1', '2026-01-01', ?)",
        (len(rows),),
    )
    conn.commit()
    conn.close()
    return path


# ── _load_gnomad_rsid_map ────────────────────────────────────────────


class TestLoadGnomadRsidMap:
    """Tests for the coordinate → rsID map loader."""

    def test_loads_mappings(self, tmp_path: Path) -> None:
        db = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [("1", 100, "A", "G", "rs100", 0.05), ("2", 200, "C", "T", "rs200", 0.10)],
        )
        result = _load_gnomad_rsid_map(db)
        assert result[("1", 100, "A", "G")] == "rs100"
        assert result[("2", 200, "C", "T")] == "rs200"

    def test_skips_null_rsids(self, tmp_path: Path) -> None:
        db = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [("1", 100, "A", "G", "rs100", 0.05), ("1", 200, "C", "T", None, 0.10)],
        )
        result = _load_gnomad_rsid_map(db)
        assert len(result) == 1
        assert ("1", 200, "C", "T") not in result

    def test_missing_db_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="gnomAD cache not found"):
            _load_gnomad_rsid_map(tmp_path / "nonexistent.sqlite")

    def test_multi_allelic_entries(self, tmp_path: Path) -> None:
        db = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [
                ("1", 100, "A", "G", "rs100", 0.05),
                ("1", 100, "A", "T", "rs100", 0.02),
            ],
        )
        result = _load_gnomad_rsid_map(db)
        assert result[("1", 100, "A", "G")] == "rs100"
        assert result[("1", 100, "A", "T")] == "rs100"


# ── _open_tsv ─────────────────────────────────────────────────────────


class TestOpenTsv:
    """Tests for the local file opener."""

    def test_reads_plain_tsv(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line())
        with _open_tsv(tsv) as fh:
            lines = fh.readlines()
        assert len(lines) == 2

    def test_reads_gzipped_tsv(self, tmp_path: Path) -> None:
        gz = tmp_path / "data.tsv.gz"
        content = _TSV_HEADER + _make_tsv_line()
        with gzip.open(gz, "wt", encoding="utf-8") as f:
            f.write(content)
        with _open_tsv(gz) as fh:
            lines = fh.readlines()
        assert len(lines) == 2

    def test_none_would_stream(self) -> None:
        """Passing None returns the Zenodo streaming context manager (not called here)."""
        ctx = _open_tsv(None)
        assert hasattr(ctx, "__enter__")


# ── TSV parsing ───────────────────────────────────────────────────────


class TestTsvParsing:
    """Tests for TSV column parsing via build_cache."""

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            "# license header line 1\n# license header line 2\n" + _TSV_HEADER + _make_tsv_line()
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
        assert count == 1

    def test_skips_short_lines(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + "chr1\t100\tA\n"  # too few columns
            + _make_tsv_line(chrom="chr2", pos=200)
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT chrom, pos FROM alphamissense_scores").fetchall()
        assert len(rows) == 1
        assert rows[0] == ("2", 200)

    def test_skips_invalid_position(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + "chr1\tBAD\tA\tG\thg38\tP1\tE1\tA1G\t0.50\tambiguous\n"
            + _make_tsv_line(chrom="chr1", pos=500)
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
        assert count == 1

    def test_skips_invalid_pathogenicity(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + "chr1\t100\tA\tG\thg38\tP1\tE1\tA1G\tN/A\tambiguous\n"
            + _make_tsv_line(chrom="chr1", pos=500)
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
        assert count == 1

    def test_chr_prefix_stripped(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line(chrom="chr17", pos=300))
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute("SELECT chrom FROM alphamissense_scores").fetchone()
        assert row[0] == "17"

    def test_empty_fields_become_null(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + "chr1\t100\tA\tG\thg38\t\t\t\t0.50\tambiguous\n")
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT uniprot_id, transcript_id, protein_variant FROM alphamissense_scores"
            ).fetchone()
        assert row == (None, None, None)


# ── gnomAD rsID join ──────────────────────────────────────────────────


class TestGnomadJoin:
    """Tests for the coordinate → rsID join during build."""

    def test_rsid_populated_from_gnomad(self, tmp_path: Path) -> None:
        gnomad = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [("1", 100, "A", "G", "rs1001", 0.05)],
        )
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line(chrom="chr1", pos=100, ref="A", alt="G"))

        output = tmp_path / "am.sqlite"
        build_cache(output, gnomad, tsv_path=tsv)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute("SELECT rsid FROM alphamissense_scores").fetchone()
        assert row[0] == "rs1001"

    def test_rsid_null_when_not_in_gnomad(self, tmp_path: Path) -> None:
        gnomad = _build_gnomad_db(tmp_path / "gnomad.sqlite", [])
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line(chrom="chr1", pos=100))

        output = tmp_path / "am.sqlite"
        build_cache(output, gnomad, tsv_path=tsv)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute("SELECT rsid FROM alphamissense_scores").fetchone()
        assert row[0] is None

    def test_skip_gnomad_all_null_rsids(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + _make_tsv_line(chrom="chr1", pos=100)
            + _make_tsv_line(chrom="chr1", pos=200, ref="C", alt="T")
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT rsid FROM alphamissense_scores").fetchall()
        assert all(row[0] is None for row in rows)

    def test_chr_prefix_normalized_for_join(self, tmp_path: Path) -> None:
        """gnomAD stores '1', AM source has 'chr1' — join must match."""
        gnomad = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [("1", 100, "A", "G", "rs42", 0.01)],
        )
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line(chrom="chr1", pos=100, ref="A", alt="G"))

        output = tmp_path / "am.sqlite"
        build_cache(output, gnomad, tsv_path=tsv)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute("SELECT rsid FROM alphamissense_scores").fetchone()
        assert row[0] == "rs42"


# ── Multi-allelic / PK behavior ───────────────────────────────────────


class TestMultiAllelic:
    """Composite PK (chrom, pos, ref, alt) preserves multi-allelic sites."""

    def test_both_alleles_stored(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + _make_tsv_line(chrom="chr1", pos=100, ref="A", alt="G", am_pathogenicity="0.20")
            + _make_tsv_line(chrom="chr1", pos=100, ref="A", alt="T", am_pathogenicity="0.92")
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT alt, am_pathogenicity FROM alphamissense_scores ORDER BY alt"
            ).fetchall()
        assert len(rows) == 2
        assert rows[0] == ("G", pytest.approx(0.20))
        assert rows[1] == ("T", pytest.approx(0.92))

    def test_insert_or_replace_on_duplicate_pk(self, tmp_path: Path) -> None:
        """Exact same (chrom,pos,ref,alt) replaces — not an error."""
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + _make_tsv_line(
                chrom="chr1",
                pos=100,
                ref="A",
                alt="G",
                am_pathogenicity="0.50",
                am_class="ambiguous",
            )
            + _make_tsv_line(
                chrom="chr1",
                pos=100,
                ref="A",
                alt="G",
                am_pathogenicity="0.95",
                am_class="likely_pathogenic",
            )
        )
        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT am_pathogenicity, am_class FROM alphamissense_scores"
            ).fetchall()
        assert len(rows) == 1
        assert rows[0] == (pytest.approx(0.95), "likely_pathogenic")


# ── Integration: build_cache end-to-end ───────────────────────────────


class TestBuildCacheIntegration:
    """End-to-end integration test: local TSV → build_cache → SQLite."""

    def test_full_build_from_local_tsv(self, tmp_path: Path) -> None:
        gnomad = _build_gnomad_db(
            tmp_path / "gnomad.sqlite",
            [
                ("1", 100, "A", "G", "rs1001", 0.05),
                ("2", 300, "G", "A", "rs2001", 0.20),
            ],
        )
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            _TSV_HEADER
            + _make_tsv_line(
                chrom="chr1",
                pos=100,
                ref="A",
                alt="G",
                uniprot_id="P12345",
                transcript_id="ENST001",
                protein_variant="A100G",
                am_pathogenicity="0.95",
                am_class="likely_pathogenic",
            )
            + _make_tsv_line(
                chrom="chr1",
                pos=200,
                ref="C",
                alt="T",
                uniprot_id="P12345",
                transcript_id="ENST001",
                protein_variant="C200T",
                am_pathogenicity="0.20",
                am_class="likely_benign",
            )
            + _make_tsv_line(
                chrom="chr2",
                pos=300,
                ref="G",
                alt="A",
                uniprot_id="P67890",
                transcript_id="ENST002",
                protein_variant="G300A",
                am_pathogenicity="0.45",
                am_class="ambiguous",
            )
        )

        output = tmp_path / "am.sqlite"
        build_cache(output, gnomad, tsv_path=tsv)

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT chrom, pos, ref, alt, rsid, uniprot_id, transcript_id, "
                "protein_variant, am_pathogenicity, am_class "
                "FROM alphamissense_scores ORDER BY chrom, pos"
            ).fetchall()

        assert len(rows) == 3
        r0 = rows[0]
        assert r0[:5] == ("1", 100, "A", "G", "rs1001")
        assert r0[5:8] == ("P12345", "ENST001", "A100G")
        assert r0[8] == pytest.approx(0.95)
        assert r0[9] == "likely_pathogenic"
        r1 = rows[1]
        assert r1[:5] == ("1", 200, "C", "T", None)
        assert r1[8] == pytest.approx(0.20)
        assert r1[9] == "likely_benign"
        r2 = rows[2]
        assert r2[:5] == ("2", 300, "G", "A", "rs2001")
        assert r2[8] == pytest.approx(0.45)
        assert r2[9] == "ambiguous"

    def test_database_versions_stamped(self, tmp_path: Path) -> None:
        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line())

        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT name, source_url, version, record_count "
                "FROM database_versions WHERE name='alphamissense'"
            ).fetchone()
        assert row is not None
        assert row[0] == "alphamissense"
        assert "zenodo" in row[1]
        assert row[2] == "2023.2"
        assert row[3] == 1

    def test_atomic_replace(self, tmp_path: Path) -> None:
        """Output file is atomically replaced via tmp → rename."""
        output = tmp_path / "am.sqlite"
        output.write_text("old data")

        tsv = tmp_path / "data.tsv"
        tsv.write_text(_TSV_HEADER + _make_tsv_line())
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
        assert count == 1

    def test_gzipped_local_tsv(self, tmp_path: Path) -> None:
        line2 = _make_tsv_line(chrom="chr2", pos=200, ref="C", alt="T")
        content = _TSV_HEADER + _make_tsv_line() + line2
        gz_path = tmp_path / "data.tsv.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as f:
            f.write(content)

        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=gz_path, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
        assert count == 2

    def test_batched_insert(self, tmp_path: Path) -> None:
        """Records exceeding _BATCH_SIZE are committed in multiple batches."""
        import build_alphamissense_cache

        original_batch = build_alphamissense_cache._BATCH_SIZE
        build_alphamissense_cache._BATCH_SIZE = 3

        try:
            tsv = tmp_path / "data.tsv"
            lines = [_TSV_HEADER]
            for i in range(7):
                lines.append(_make_tsv_line(chrom=f"chr{i + 1}", pos=100 + i, ref="A", alt="G"))
            tsv.write_text("".join(lines))

            output = tmp_path / "am.sqlite"
            build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

            with contextlib.closing(sqlite3.connect(output)) as conn:
                count = conn.execute("SELECT COUNT(*) FROM alphamissense_scores").fetchone()[0]
            assert count == 7
        finally:
            build_alphamissense_cache._BATCH_SIZE = original_batch

    def test_record_count_in_versions(self, tmp_path: Path) -> None:
        """database_versions.record_count matches actual inserted rows."""
        tsv = tmp_path / "data.tsv"
        lines = [_TSV_HEADER]
        for i in range(5):
            lines.append(_make_tsv_line(chrom=f"chr{i + 1}", pos=100 + i))
        tsv.write_text("".join(lines))

        output = tmp_path / "am.sqlite"
        build_cache(output, tmp_path / "gnomad.sqlite", tsv_path=tsv, skip_gnomad=True)

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT record_count FROM database_versions WHERE name='alphamissense'"
            ).fetchone()
        assert row[0] == 5

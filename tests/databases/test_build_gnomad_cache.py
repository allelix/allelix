# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for scripts/build_gnomad_cache.py helper functions."""

from __future__ import annotations

import contextlib
import gzip
import io
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from build_gnomad_cache import (  # noqa: I001
    GNOMAD_FILENAME,
    _iter_vcf_records,
    _load_manifest,
    _parse_info_field,
    _read_local_vcf_gz,
    _safe_float,
    build_cache,
)


# ── _parse_info_field ──────────────────────────────────────────────────


class TestParseInfoField:
    """Tests for VCF INFO field parsing."""

    def test_key_value_pairs(self) -> None:
        info = "AF=0.5;AF_grpmax=0.6;grpmax=nfe"
        result = _parse_info_field(info)
        assert result == {"AF": "0.5", "AF_grpmax": "0.6", "grpmax": "nfe"}

    def test_flag_without_value(self) -> None:
        info = "AF=0.5;PASS;grpmax=nfe"
        result = _parse_info_field(info)
        assert result["AF"] == "0.5"
        assert result["PASS"] == ""
        assert result["grpmax"] == "nfe"

    def test_single_entry(self) -> None:
        assert _parse_info_field("AF=0.1") == {"AF": "0.1"}

    def test_value_with_equals(self) -> None:
        result = _parse_info_field("VQSR=a=b;AF=0.1")
        assert result["VQSR"] == "a=b"
        assert result["AF"] == "0.1"


# ── _safe_float ────────────────────────────────────────────────────────


class TestSafeFloat:
    """Tests for safe float conversion."""

    def test_valid_float(self) -> None:
        assert _safe_float("0.123") == pytest.approx(0.123)

    def test_scientific_notation(self) -> None:
        assert _safe_float("1.5e-4") == pytest.approx(1.5e-4)

    def test_none_input(self) -> None:
        assert _safe_float(None) is None

    def test_empty_string(self) -> None:
        assert _safe_float("") is None

    def test_dot(self) -> None:
        assert _safe_float(".") is None

    def test_non_numeric(self) -> None:
        assert _safe_float("NA") is None

    def test_zero(self) -> None:
        assert _safe_float("0") == 0.0


# ── _load_manifest ─────────────────────────────────────────────────────


class TestLoadManifest:
    """Tests for rsID manifest loading."""

    def test_loads_rsids(self, tmp_path: Path) -> None:
        manifest = tmp_path / "manifest.txt"
        manifest.write_text("rs123\nrs456\nrs789\n")
        result = _load_manifest(manifest)
        assert result == {"rs123", "rs456", "rs789"}

    def test_skips_blanks_and_non_rs(self, tmp_path: Path) -> None:
        manifest = tmp_path / "manifest.txt"
        manifest.write_text("rs123\n\n# comment\nchr1:100\nrs456\n")
        result = _load_manifest(manifest)
        assert result == {"rs123", "rs456"}

    def test_strips_whitespace(self, tmp_path: Path) -> None:
        manifest = tmp_path / "manifest.txt"
        manifest.write_text("  rs123  \nrs456\t\n")
        result = _load_manifest(manifest)
        assert result == {"rs123", "rs456"}

    def test_empty_file(self, tmp_path: Path) -> None:
        manifest = tmp_path / "manifest.txt"
        manifest.write_text("")
        result = _load_manifest(manifest)
        assert result == set()


# ── _iter_vcf_records ──────────────────────────────────────────────────

_VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _make_vcf_line(
    rsid: str = "rs123",
    af: str = "0.05",
    af_grpmax: str = "0.08",
    grpmax: str = "nfe",
    af_afr: str = "0.01",
    af_amr: str = "0.02",
    af_asj: str = "0.03",
    af_eas: str = "0.04",
    af_fin: str = "0.06",
    af_nfe: str = "0.08",
    af_sas: str = "0.07",
) -> str:
    info = (
        f"AF={af};AF_grpmax={af_grpmax};grpmax={grpmax}"
        f";AF_afr={af_afr};AF_amr={af_amr};AF_asj={af_asj}"
        f";AF_eas={af_eas};AF_fin={af_fin};AF_nfe={af_nfe}"
        f";AF_sas={af_sas}"
    )
    return f"chr1\t100\t{rsid}\tA\tG\t.\tPASS\t{info}\n"


class TestIterVcfRecords:
    """Tests for VCF record iterator."""

    def test_yields_records(self) -> None:
        vcf = _VCF_HEADER + _make_vcf_line(rsid="rs100", af="0.25")
        records = list(_iter_vcf_records(io.StringIO(vcf)))
        assert len(records) == 1
        assert records[0][0] == "1"
        assert records[0][1] == 100
        assert records[0][2] == "A"
        assert records[0][3] == "G"
        assert records[0][4] == "rs100"
        assert records[0][5] == pytest.approx(0.25)

    def test_skips_header_lines(self) -> None:
        vcf = _VCF_HEADER + _make_vcf_line()
        records = list(_iter_vcf_records(io.StringIO(vcf)))
        assert len(records) == 1

    def test_skips_dot_rsid(self) -> None:
        line = "chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.1\n"
        records = list(_iter_vcf_records(io.StringIO(line)))
        assert records == []

    def test_skips_non_rs_id(self) -> None:
        line = "chr1\t100\tgnomAD_123\tA\tG\t.\tPASS\tAF=0.1\n"
        records = list(_iter_vcf_records(io.StringIO(line)))
        assert records == []

    def test_rsid_filter(self) -> None:
        vcf = (
            _VCF_HEADER
            + _make_vcf_line(rsid="rs100")
            + _make_vcf_line(rsid="rs200")
            + _make_vcf_line(rsid="rs300")
        )
        records = list(_iter_vcf_records(io.StringIO(vcf), rsid_filter={"rs100", "rs300"}))
        assert len(records) == 2
        assert records[0][4] == "rs100"
        assert records[1][4] == "rs300"

    def test_missing_info_fields(self) -> None:
        line = "chr1\t100\trs999\tA\tG\t.\tPASS\tAF=0.1\n"
        records = list(_iter_vcf_records(io.StringIO(line)))
        assert len(records) == 1
        assert records[0][0] == "1"
        assert records[0][1] == 100
        assert records[0][2] == "A"
        assert records[0][3] == "G"
        assert records[0][4] == "rs999"
        assert records[0][5] == pytest.approx(0.1)
        assert records[0][6] is None  # af_popmax
        assert records[0][7] is None  # popmax
        assert records[0][8] is None  # af_afr

    def test_all_population_frequencies(self) -> None:
        vcf = _VCF_HEADER + _make_vcf_line(
            rsid="rs500",
            af="0.05",
            af_afr="0.01",
            af_amr="0.02",
            af_asj="0.03",
            af_eas="0.04",
            af_fin="0.06",
            af_nfe="0.08",
            af_sas="0.07",
        )
        records = list(_iter_vcf_records(io.StringIO(vcf)))
        row = records[0]
        assert row[0] == "1"
        assert row[1] == 100
        assert row[2] == "A"
        assert row[3] == "G"
        assert row[4] == "rs500"
        assert row[5] == pytest.approx(0.05)  # af
        assert row[8] == pytest.approx(0.01)  # afr
        assert row[9] == pytest.approx(0.02)  # amr
        assert row[10] == pytest.approx(0.03)  # asj
        assert row[11] == pytest.approx(0.04)  # eas
        assert row[12] == pytest.approx(0.06)  # fin
        assert row[13] == pytest.approx(0.08)  # nfe
        assert row[14] == pytest.approx(0.07)  # sas

    def test_short_line_skipped(self) -> None:
        line = "chr1\t100\trs1\n"
        records = list(_iter_vcf_records(io.StringIO(line)))
        assert records == []

    def test_multi_allelic_same_rsid(self) -> None:
        """Two VCF lines with the same rsID but different ALT alleles."""
        vcf = (
            _VCF_HEADER
            + "chr1\t100\trs100\tA\tG\t.\tPASS\tAF=0.05\n"
            + "chr1\t100\trs100\tA\tT\t.\tPASS\tAF=0.02\n"
        )
        records = list(_iter_vcf_records(io.StringIO(vcf)))
        assert len(records) == 2
        assert records[0][4] == "rs100"
        assert records[0][3] == "G"
        assert records[0][5] == pytest.approx(0.05)
        assert records[1][4] == "rs100"
        assert records[1][3] == "T"
        assert records[1][5] == pytest.approx(0.02)


class TestBuildCacheMultiAllelic:
    """Composite PK preserves multi-allelic sites with the same rsID."""

    def test_both_alleles_stored(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Two lines sharing an rsID but different ALTs both persist."""
        vcf_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t100\trs100\tA\tG\t.\tPASS\tAF=0.05\n",
            "chr1\t100\trs100\tA\tT\t.\tPASS\tAF=0.02\n",
        ]
        gz_data = _make_gzipped_vcf(vcf_lines)

        import urllib.request

        class FakeResponse(io.BytesIO):
            """Simulates an HTTP response returning gzipped VCF data."""

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                pass

        def fake_urlopen(_req, **_kw):
            return FakeResponse(gz_data)

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        output = tmp_path / "gnomad.sqlite"
        build_cache(output, full=True, chromosomes=["1"])

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT chrom, pos, ref, alt, rsid, af FROM gnomad_frequencies ORDER BY alt"
            ).fetchall()
        assert len(rows) == 2
        assert rows[0] == ("1", 100, "A", "G", "rs100", pytest.approx(0.05))
        assert rows[1] == ("1", 100, "A", "T", "rs100", pytest.approx(0.02))


# ── Integration: gzipped VCF → SQLite ─────────────────────────────────


def _make_gzipped_vcf(lines: list[str]) -> bytes:
    """Produce a gzip-compressed VCF from text lines."""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write("".join(lines).encode("utf-8"))
    return buf.getvalue()


class TestBuildCacheIntegration:
    """Integration test: synthetic gzipped VCF → build_cache → SQLite."""

    def test_build_cache_from_manifest(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Build a cache from a manifest-filtered synthetic VCF."""
        vcf_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t100\trs100\tA\tG\t.\tPASS\tAF=0.05;AF_grpmax=0.08;grpmax=nfe;AF_afr=0.01;AF_amr=0.02;AF_asj=0.03;AF_eas=0.04;AF_fin=0.06;AF_nfe=0.08;AF_sas=0.07\n",
            "chr1\t200\trs200\tC\tT\t.\tPASS\tAF=0.10;AF_nfe=0.12\n",
            "chr1\t300\trs300\tG\tA\t.\tPASS\tAF=0.50\n",
        ]
        gz_data = _make_gzipped_vcf(vcf_lines)

        manifest = tmp_path / "manifest.txt"
        manifest.write_text("rs100\nrs300\n")

        import urllib.request

        class FakeResponse(io.BytesIO):
            """Simulates an HTTP response returning gzipped VCF data."""

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                pass

        def fake_urlopen(_req, **_kw):
            return FakeResponse(gz_data)

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=False,
            manifest_path=manifest,
            chromosomes=["1"],
        )

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute(
                "SELECT rsid, chrom, pos, ref, alt, af, af_nfe"
                " FROM gnomad_frequencies ORDER BY rsid"
            ).fetchall()
        assert len(rows) == 2
        assert rows[0][0] == "rs100"
        assert rows[0][1] == "1"
        assert rows[0][2] == 100
        assert rows[0][3] == "A"
        assert rows[0][4] == "G"
        assert rows[0][5] == pytest.approx(0.05)
        assert rows[0][6] == pytest.approx(0.08)
        assert rows[1][0] == "rs300"
        assert rows[1][1] == "1"
        assert rows[1][2] == 300
        assert rows[1][3] == "G"
        assert rows[1][4] == "A"
        assert rows[1][5] == pytest.approx(0.50)
        assert rows[1][6] is None

    def test_build_cache_full_mode(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Full mode inserts all records (no rsid filtering)."""
        vcf_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr22\t100\trs1001\tA\tG\t.\tPASS\tAF=0.01\n",
            "chr22\t200\trs1002\tC\tT\t.\tPASS\tAF=0.02\n",
        ]
        gz_data = _make_gzipped_vcf(vcf_lines)

        import urllib.request

        class FakeResponse(io.BytesIO):
            """Simulates an HTTP response returning gzipped VCF data."""

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                pass

        def fake_urlopen(_req, **_kw):
            return FakeResponse(gz_data)

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=True,
            chromosomes=["22"],
        )

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT rsid FROM gnomad_frequencies ORDER BY rsid").fetchall()
        assert len(rows) == 2
        assert rows[0][0] == "rs1001"
        assert rows[1][0] == "rs1002"

    def test_streaming_decompression(self, tmp_path: Path) -> None:
        """Verify gzip.GzipFile + TextIOWrapper yields records from compressed data."""
        vcf_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t100\trs42\tA\tG\t.\tPASS\tAF=0.33;AF_nfe=0.40\n",
        ]
        gz_data = _make_gzipped_vcf(vcf_lines)
        decompressor = gzip.GzipFile(fileobj=io.BytesIO(gz_data))
        text_stream = io.TextIOWrapper(decompressor, encoding="utf-8")
        records = list(_iter_vcf_records(text_stream))
        assert len(records) == 1
        assert records[0][0] == "1"
        assert records[0][1] == 100
        assert records[0][2] == "A"
        assert records[0][3] == "G"
        assert records[0][4] == "rs42"
        assert records[0][5] == pytest.approx(0.33)
        assert records[0][13] == pytest.approx(0.40)  # af_nfe


# ── Local file support ────────────────────────────────────────────────


_LOCAL_VCF_LINES = [
    "##fileformat=VCFv4.2\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    "chr1\t100\trs10\tA\tG\t.\tPASS\tAF=0.10;AF_nfe=0.12\n",
    "chr1\t200\trs20\tC\tT\t.\tPASS\tAF=0.20\n",
    "chr1\t300\trs30\tG\tA\t.\tPASS\tAF=0.30;AF_afr=0.35\n",
]


def _write_local_bgz(directory: Path, chrom: str, lines: list[str]) -> Path:
    """Write a gzipped VCF to the expected gnomAD filename in directory."""
    path = directory / GNOMAD_FILENAME.format(chrom=chrom)
    with gzip.open(path, "wb") as gz:
        gz.write("".join(lines).encode("utf-8"))
    return path


class TestReadLocalVcfGz:
    """Tests for _read_local_vcf_gz."""

    def test_reads_all_records(self, tmp_path: Path) -> None:
        path = _write_local_bgz(tmp_path, "1", _LOCAL_VCF_LINES)
        records = list(_read_local_vcf_gz(path))
        assert len(records) == 3
        assert records[0][0] == "1"
        assert records[0][1] == 100
        assert records[0][2] == "A"
        assert records[0][3] == "G"
        assert records[0][4] == "rs10"
        assert records[1][0] == "1"
        assert records[1][1] == 200
        assert records[1][2] == "C"
        assert records[1][3] == "T"
        assert records[1][4] == "rs20"
        assert records[2][0] == "1"
        assert records[2][1] == 300
        assert records[2][2] == "G"
        assert records[2][3] == "A"
        assert records[2][4] == "rs30"

    def test_reads_with_rsid_filter(self, tmp_path: Path) -> None:
        path = _write_local_bgz(tmp_path, "1", _LOCAL_VCF_LINES)
        records = list(_read_local_vcf_gz(path, rsid_filter={"rs10", "rs30"}))
        assert len(records) == 2
        assert records[0][4] == "rs10"
        assert records[1][4] == "rs30"

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "nonexistent.vcf.bgz"
        with pytest.raises(FileNotFoundError):
            list(_read_local_vcf_gz(path))


class TestBuildCacheLocal:
    """Integration tests for build_cache with --local-dir."""

    def test_local_dir_full(self, tmp_path: Path) -> None:
        """Full mode from local files inserts all records."""
        vcf_dir = tmp_path / "vcfs"
        vcf_dir.mkdir()
        _write_local_bgz(vcf_dir, "22", _LOCAL_VCF_LINES)

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=True,
            chromosomes=["22"],
            local_dir=vcf_dir,
        )

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT rsid, af FROM gnomad_frequencies ORDER BY rsid").fetchall()
        assert len(rows) == 3
        assert rows[0] == ("rs10", pytest.approx(0.10))
        assert rows[1] == ("rs20", pytest.approx(0.20))
        assert rows[2] == ("rs30", pytest.approx(0.30))

    def test_local_dir_with_manifest(self, tmp_path: Path) -> None:
        """Manifest filtering works with local files."""
        vcf_dir = tmp_path / "vcfs"
        vcf_dir.mkdir()
        _write_local_bgz(vcf_dir, "1", _LOCAL_VCF_LINES)

        manifest = tmp_path / "manifest.txt"
        manifest.write_text("rs10\nrs30\n")

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=False,
            manifest_path=manifest,
            chromosomes=["1"],
            local_dir=vcf_dir,
        )

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT rsid FROM gnomad_frequencies ORDER BY rsid").fetchall()
        assert len(rows) == 2
        assert rows[0][0] == "rs10"
        assert rows[1][0] == "rs30"

    def test_local_dir_missing_chrom_skips(self, tmp_path: Path) -> None:
        """Missing chromosome file is skipped, not fatal."""
        vcf_dir = tmp_path / "vcfs"
        vcf_dir.mkdir()

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=True,
            chromosomes=["22"],
            local_dir=vcf_dir,
        )

        assert output.exists()
        with contextlib.closing(sqlite3.connect(output)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM gnomad_frequencies").fetchone()[0]
        assert count == 0

    def test_local_dir_stamps_source(self, tmp_path: Path) -> None:
        """database_versions.source_url records the local directory path."""
        vcf_dir = tmp_path / "vcfs"
        vcf_dir.mkdir()
        _write_local_bgz(vcf_dir, "22", _LOCAL_VCF_LINES)

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=True,
            chromosomes=["22"],
            local_dir=vcf_dir,
        )

        with contextlib.closing(sqlite3.connect(output)) as conn:
            row = conn.execute(
                "SELECT source_url FROM database_versions WHERE name='gnomad'"
            ).fetchone()
        assert row is not None
        assert str(vcf_dir) in row[0]

    def test_local_dir_multi_chrom(self, tmp_path: Path) -> None:
        """Multiple local chromosome files are processed together."""
        vcf_dir = tmp_path / "vcfs"
        vcf_dir.mkdir()

        chr1_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t100\trs111\tA\tG\t.\tPASS\tAF=0.11\n",
        ]
        chr22_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr22\t100\trs222\tC\tT\t.\tPASS\tAF=0.22\n",
        ]
        _write_local_bgz(vcf_dir, "1", chr1_lines)
        _write_local_bgz(vcf_dir, "22", chr22_lines)

        output = tmp_path / "gnomad.sqlite"
        build_cache(
            output,
            full=True,
            chromosomes=["1", "22"],
            local_dir=vcf_dir,
        )

        with contextlib.closing(sqlite3.connect(output)) as conn:
            rows = conn.execute("SELECT rsid, af FROM gnomad_frequencies ORDER BY rsid").fetchall()
        assert len(rows) == 2
        assert rows[0] == ("rs111", pytest.approx(0.11))
        assert rows[1] == ("rs222", pytest.approx(0.22))

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Cut or damaged gzip input: stop with a named error, or warn when the cut cannot be seen.

Every file is built in ``tmp_path`` from a committed synthetic fixture.
BGZF blocks are written here with the standard layout (gzip members with
the ``BC`` extra subfield), split so the file has two data blocks.
"""

from __future__ import annotations

import gzip
import json
import logging
import struct
import zlib
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import PARSERS, DamagedFileError, _helpers, detect_parser
from allelix.parsers._helpers import BGZF_EOF_BLOCK, open_genotype_text
from allelix.parsers.vcf import VcfParser

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"
VCF = FIXTURES_DIR / "mock_vcf.vcf"
FTDNA = FIXTURES_DIR / "mock_ftdna.csv"
TRAILING_HTML = b"\r\n\r\n<!DOCTYPE html><html></html>\n"


def _bgzf_block(data: bytes) -> bytes:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = compressor.compress(data) + compressor.flush()
    bsize = 18 + len(cdata) + 8 - 1
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff" + struct.pack("<H", 6)
    extra = b"BC" + struct.pack("<HH", 2, bsize)
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + extra + cdata + trailer


def _two_blocks(source: Path) -> tuple[bytes, bytes]:
    """The file as two BGZF blocks; the second holds the last 8 lines."""
    lines = source.read_bytes().splitlines(keepends=True)
    split = len(lines) - 8
    return _bgzf_block(b"".join(lines[:split])), _bgzf_block(b"".join(lines[split:]))


@pytest.fixture(autouse=True)
def _fresh_warnings(monkeypatch: pytest.MonkeyPatch) -> None:
    # Warnings are once per file per process; each test starts clean.
    monkeypatch.setattr(_helpers, "_warned", set())


@pytest.fixture
def intact(tmp_path: Path) -> Path:
    first, second = _two_blocks(VCF)
    path = tmp_path / "intact.vcf.gz"
    path.write_bytes(first + second + BGZF_EOF_BLOCK)
    return path


@pytest.fixture
def cut_at_boundary(tmp_path: Path) -> Path:
    first, _ = _two_blocks(VCF)
    path = tmp_path / "boundary.vcf.gz"
    path.write_bytes(first)
    return path


@pytest.fixture
def cut_mid_block(tmp_path: Path) -> Path:
    first, second = _two_blocks(VCF)
    path = tmp_path / "midblock.vcf.gz"
    path.write_bytes(first + second[: len(second) // 2])
    return path


def test_block_writer_matches_the_standard_eof_block() -> None:
    assert _bgzf_block(b"") == BGZF_EOF_BLOCK


class TestBgzfEndBlock:
    def test_intact_file_reads_in_full_without_warning(
        self, intact: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            got = list(VcfParser().parse(intact))
        assert got == list(VcfParser().parse(VCF))
        assert "end-of-file block" not in caplog.text

    def test_cut_at_block_boundary_warns_and_keeps_the_records_before_the_cut(
        self, cut_at_boundary: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            got = list(VcfParser().parse(cut_at_boundary))
        full = list(VcfParser().parse(VCF))
        assert 0 < len(got) < len(full)
        assert got == full[: len(got)]
        assert "boundary.vcf.gz does not end with the BGZF end-of-file block" in caplog.text

    def test_warning_appears_once_across_repeated_reads(
        self, cut_at_boundary: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            for _ in range(3):
                list(VcfParser().parse(cut_at_boundary))
        assert caplog.text.count("end-of-file block") == 1

    def test_detection_alone_does_not_warn(
        self, cut_at_boundary: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            assert detect_parser(cut_at_boundary).name == "vcf"
        assert "end-of-file block" not in caplog.text

    def test_plain_gzip_is_not_held_to_the_bgzf_rule(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        path = tmp_path / "plain.vcf.gz"
        path.write_bytes(gzip.compress(VCF.read_bytes()))
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            assert list(VcfParser().parse(path)) == list(VcfParser().parse(VCF))
        assert "end-of-file block" not in caplog.text


class TestCutInsideABlock:
    def test_parse_raises_a_named_error(self, cut_mid_block: Path) -> None:
        with pytest.raises(
            DamagedFileError, match=r"midblock\.vcf\.gz looks truncated or damaged"
        ):
            list(VcfParser().parse(cut_mid_block))

    def test_detection_still_works(self, cut_mid_block: Path) -> None:
        # The header block is intact, so the file is recognised; the damage
        # is reported when the records are read.
        assert detect_parser(cut_mid_block).name == "vcf"

    def test_corrupt_block_raises_a_named_error(self, tmp_path: Path) -> None:
        first, second = _two_blocks(VCF)
        damaged = bytearray(second)
        damaged[len(damaged) // 2] ^= 0xFF
        path = tmp_path / "corrupt.vcf.gz"
        path.write_bytes(first + bytes(damaged) + BGZF_EOF_BLOCK)
        with pytest.raises(DamagedFileError, match=r"corrupt\.vcf\.gz looks truncated or damaged"):
            list(VcfParser().parse(path))

    def test_cut_array_file_raises_a_named_error(self, tmp_path: Path) -> None:
        data = gzip.compress(FTDNA.read_bytes())
        path = tmp_path / "cut.csv.gz"
        # Keep the header and first rows, lose the end of the stream.
        path.write_bytes(data[: len(data) - 12])
        parser = detect_parser(FTDNA)
        with pytest.raises(DamagedFileError, match=r"cut\.csv\.gz looks truncated or damaged"):
            list(parser.parse(path))

    @pytest.mark.parametrize("command", ["stats", "extract"])
    def test_cli_prints_the_error_and_exits_1(self, cut_mid_block: Path, command: str) -> None:
        args = [command, str(cut_mid_block)]
        if command == "extract":
            args += ["--snps", "rs1"]
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 1
        assert "midblock.vcf.gz looks truncated or damaged" in result.output
        assert "Aborted" not in result.output
        assert not isinstance(result.exception, DamagedFileError)


class TestBytesAfterTheLastMember:
    def test_trailing_non_gzip_bytes_are_ignored_with_a_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        path = tmp_path / "trailing.csv.gz"
        path.write_bytes(gzip.compress(FTDNA.read_bytes()) + TRAILING_HTML)
        parser = detect_parser(path)
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            got = list(parser.parse(path))
        assert got == list(parser.parse(FTDNA))
        assert f"ignored {len(TRAILING_HTML)} bytes after the end of the compressed data" in (
            caplog.text
        )

    def test_several_members_are_all_read(self, tmp_path: Path) -> None:
        lines = FTDNA.read_bytes().splitlines(keepends=True)
        path = tmp_path / "members.csv.gz"
        path.write_bytes(gzip.compress(b"".join(lines[:6])) + gzip.compress(b"".join(lines[6:])))
        parser = detect_parser(FTDNA)
        assert list(parser.parse(path)) == list(parser.parse(FTDNA))

    def test_non_gzip_after_the_magic_bytes_is_damage(self, tmp_path: Path) -> None:
        path = tmp_path / "fake.gz"
        path.write_bytes(b"\x1f\x8b" + b"\x00" * 40)
        with pytest.raises(DamagedFileError), open_genotype_text(path) as fh:
            fh.read()

    @pytest.mark.parametrize("parser", PARSERS, ids=lambda p: p.name)
    def test_damaged_file_is_not_claimed(self, tmp_path: Path, parser: object) -> None:
        path = tmp_path / "fake.gz"
        path.write_bytes(b"\x1f\x8b" + b"\x00" * 40)
        assert parser.can_parse(path) is False  # type: ignore[attr-defined]


class TestNoticesOnTheConsole:
    def test_stats_shows_the_warning_and_does_not_count_it_as_a_skipped_row(
        self, cut_at_boundary: Path
    ) -> None:
        result = CliRunner().invoke(main, ["stats", str(cut_at_boundary)])
        assert result.exit_code == 0, result.output
        assert (
            result.output.count("boundary.vcf.gz does not end with the BGZF end-of-file block")
            == 1
        )
        assert "Skipped (malformed)" not in result.output

    def test_trailing_bytes_notice_is_not_a_skipped_row(self, tmp_path: Path) -> None:
        path = tmp_path / "trailing.csv.gz"
        path.write_bytes(gzip.compress(FTDNA.read_bytes()) + TRAILING_HTML)
        result = CliRunner().invoke(main, ["stats", str(path)])
        assert result.exit_code == 0, result.output
        assert "ignored" in result.output
        assert "Skipped (malformed)" not in result.output

    def test_tiny_file_detection_stays_quiet(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        # Detection can read a small file to its end; the notice waits for the parse.
        path = tmp_path / "tiny.csv.gz"
        path.write_bytes(gzip.compress(FTDNA.read_bytes()) + TRAILING_HTML)
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            parser = detect_parser(path)
        assert caplog.text == ""
        with caplog.at_level(logging.WARNING, logger="allelix.parsers"):
            list(parser.parse(path))
        assert "ignored" in caplog.text


def test_header_build_scan_on_a_damaged_file_is_none(tmp_path: Path) -> None:
    """The header build scan answers "no build stated" for a damaged file; the parse reports it."""
    from allelix.parsers.base import read_header_build

    path = tmp_path / "fake.gz"
    path.write_bytes(b"\x1f\x8b" + b"\x00" * 40)
    assert read_header_build(path) is None


class TestAnalyze:
    """The reviewer's case: analyze on a gzipped FTDNA file with trailing bytes, and a cut one."""

    def test_trailing_bytes_are_ignored_with_a_warning(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        path = tmp_path / "trailing.csv.gz"
        path.write_bytes(gzip.compress(FTDNA.read_bytes()) + TRAILING_HTML)
        out = tmp_path / "report.json"
        result = CliRunner().invoke(
            main,
            ["analyze", str(path), "--data-dir", str(clinvar_data_dir), "--no-update",
             "--build", "grch37", "--output", str(out)],
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert result.output.count(f"ignored {len(TRAILING_HTML)} bytes") == 1
        assert out.exists()

    def test_cut_file_stops_with_a_named_error_and_writes_nothing(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        data = gzip.compress(FTDNA.read_bytes())
        path = tmp_path / "cut.csv.gz"
        path.write_bytes(data[: len(data) - 12])
        out = tmp_path / "report.json"
        result = CliRunner().invoke(
            main,
            ["analyze", str(path), "--data-dir", str(clinvar_data_dir), "--no-update",
             "--build", "grch37", "--output", str(out)],
        )  # fmt: skip
        assert result.exit_code == 1
        assert "cut.csv.gz looks truncated or damaged" in " ".join(result.output.split())
        assert "Traceback" not in result.output
        assert not isinstance(result.exception, (DamagedFileError, EOFError))
        assert not out.exists()


class TestNoticesInTheReport:
    """Input notices are report warnings too (ADR-0049 ``report_warnings``), once each."""

    def _analyze(self, path: Path, data_dir: Path, out: Path, *extra: str) -> str:
        result = CliRunner().invoke(
            main,
            ["analyze", str(path), "--data-dir", str(data_dir), "--no-update",
             "--output", str(out), *extra],
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        return result.output

    def test_trailing_bytes_in_json_and_html(self, tmp_path: Path, clinvar_data_dir: Path) -> None:
        path = tmp_path / "trailing.csv.gz"
        path.write_bytes(gzip.compress(FTDNA.read_bytes()) + TRAILING_HTML)
        out = tmp_path / "report.json"
        self._analyze(path, clinvar_data_dir, out, "--build", "grch37")
        warnings = json.loads(out.read_text())["report_warnings"]
        notices = [w for w in warnings if w["code"] == "input_trailing_bytes"]
        assert len(notices) == 1
        assert f"ignored {len(TRAILING_HTML)} bytes" in notices[0]["message"]
        html = tmp_path / "report.html"
        _helpers._warned.clear()
        self._analyze(path, clinvar_data_dir, html, "--build", "grch37")
        assert f"ignored {len(TRAILING_HTML)} bytes" in html.read_text()

    def test_cut_at_a_block_boundary_in_json(
        self, cut_at_boundary: Path, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        out = tmp_path / "report.json"
        self._analyze(cut_at_boundary, clinvar_data_dir, out)
        codes = [w["code"] for w in json.loads(out.read_text())["report_warnings"]]
        assert codes.count("input_possibly_truncated") == 1

    def test_an_intact_file_has_no_input_notice(
        self, intact: Path, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        out = tmp_path / "report.json"
        self._analyze(intact, clinvar_data_dir, out)
        codes = [w["code"] for w in json.loads(out.read_text())["report_warnings"]]
        assert not [c for c in codes if c.startswith("input_")]

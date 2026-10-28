# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Zip downloads: one genotype file in a zip is read in place; anything else says why not.

Every zip is built in ``tmp_path`` from a committed synthetic fixture. The
line-ending matrix (``test_line_endings.py``) covers each format zipped; this
file covers what a zip can hold besides one plain file.
"""

from __future__ import annotations

import io
import struct
import zipfile
import zlib
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import (
    DamagedFileError,
    InputContainerError,
    ParserNotFoundError,
    _helpers,
    detect_parser,
)
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser

FIXTURE = Path(__file__).parent.parent / "fixtures" / "mock_23andme.txt"


def _zip(path: Path, members: dict[str, bytes], method: int = zipfile.ZIP_DEFLATED) -> Path:
    with zipfile.ZipFile(path, "w", method) as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    return path


def _understated(path: Path, data: bytes, declared: int, method: int) -> Path:
    """A zip whose entry declares ``declared`` bytes, with a CRC written for that many."""
    _zip(path, {"genome.txt": data}, method)
    raw = bytearray(path.read_bytes())
    central = raw.index(b"PK\x01\x02")
    for crc_at, size_at in ((14, 22), (central + 16, central + 24)):
        struct.pack_into("<I", raw, crc_at, zlib.crc32(data[:declared]))
        struct.pack_into("<I", raw, size_at, declared)
    path.write_bytes(bytes(raw))
    return path


def _variants(path: Path) -> list:
    return list(TwentyThreeAndMeParser().parse(path))


class TestOneFile:
    def test_stored_and_deflated_read_like_the_plain_file(self, tmp_path: Path) -> None:
        data = FIXTURE.read_bytes()
        for method in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED, zipfile.ZIP_BZIP2):
            path = _zip(tmp_path / f"m{method}.zip", {"genome.txt": data}, method)
            assert detect_parser(path).name == "23andme"
            assert _variants(path) == _variants(FIXTURE)

    def test_extension_does_not_matter(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "genome.txt", {"genome.txt": FIXTURE.read_bytes()})
        assert detect_parser(path).name == "23andme"

    def test_folder_and_macos_debris_are_not_files(self, tmp_path: Path) -> None:
        path = tmp_path / "finder.zip"
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("export/", b"")
            archive.writestr("export/genome.txt", FIXTURE.read_bytes())
            archive.writestr("__MACOSX/export/._genome.txt", b"\x00\x05\x16\x07")
            archive.writestr("export/.DS_Store", b"\x00\x00\x00\x01Bud1")
        assert detect_parser(path).name == "23andme"
        assert _variants(path) == _variants(FIXTURE)

    def test_byte_order_mark_in_the_member_is_dropped(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "bom.zip", {"genome.txt": b"\xef\xbb\xbf" + FIXTURE.read_bytes()})
        assert detect_parser(path).name == "23andme"
        assert _variants(path) == _variants(FIXTURE)


class TestRefused:
    def test_several_files_are_refused_with_the_count(self, tmp_path: Path) -> None:
        path = _zip(
            tmp_path / "bundle.zip",
            {"a.vcf": b"##fileformat=VCFv4.2\n", "a.bed": b"chrY\t1\t2\n", "readme.txt": b"x"},
        )
        with pytest.raises(InputContainerError) as exc:
            detect_parser(path)
        message = str(exc.value)
        assert "bundle.zip is a zip file holding 3 files (a.vcf, a.bed, readme.txt)" in message
        assert "unzip this one and pass the genotype file itself" in message
        # A forced format reads through the same opener and gets the same error.
        with pytest.raises(InputContainerError, match="holding 3 files"):
            _variants(path)

    def test_long_member_lists_are_cut(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "many.zip", {f"f{n}.txt": b"x" for n in range(8)})
        with pytest.raises(InputContainerError, match=r"\(f0.txt, .*, f4.txt and 3 more\)"):
            detect_parser(path)

    def test_an_empty_zip_is_refused(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "empty.zip", {"folder/": b""})
        with pytest.raises(InputContainerError, match="no file in it"):
            detect_parser(path)

    def test_a_zip_with_no_entries_at_all_is_refused(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "bare.zip", {})
        assert path.read_bytes()[:4] == b"PK\x05\x06"  # only the end-of-directory record
        with pytest.raises(InputContainerError, match="no file in it"):
            detect_parser(path)

    def test_a_zip_inside_a_zip_is_refused(self, tmp_path: Path) -> None:
        inner = _zip(tmp_path / "inner.zip", {"genome.txt": FIXTURE.read_bytes()}).read_bytes()
        path = _zip(tmp_path / "outer.zip", {"inner.zip": inner})
        with pytest.raises(InputContainerError, match="holding another zip file"):
            detect_parser(path)
        with pytest.raises(InputContainerError, match="holding another zip file"):
            _variants(path)

    def test_password_protected(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "locked.zip", {"genome.txt": FIXTURE.read_bytes()})
        data = bytearray(path.read_bytes())
        # Set the "encrypted" flag in the local and central headers; zipfile
        # will not write it, and the data need not really be encrypted.
        data[6] |= 0x1
        data[data.index(b"PK\x01\x02") + 8] |= 0x1
        path.write_bytes(bytes(data))
        with pytest.raises(InputContainerError, match="password-protected"):
            detect_parser(path)

    def test_unreadable_method_is_named(self, tmp_path: Path) -> None:
        path = _zip(tmp_path / "d64.zip", {"genome.txt": FIXTURE.read_bytes()}, zipfile.ZIP_STORED)
        data = bytearray(path.read_bytes())
        central = data.index(b"PK\x01\x02")
        data[central + 10 : central + 12] = (9).to_bytes(2, "little")  # Deflate64
        path.write_bytes(bytes(data))
        with pytest.raises(InputContainerError, match="compressed with Deflate64"):
            detect_parser(path)

    def test_size_cap_on_the_declared_size(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A member that declares more than the cap stops before it is read."""
        path = _zip(tmp_path / "big.zip", {"genome.txt": FIXTURE.read_bytes()})
        monkeypatch.setattr(_helpers, "ZIP_MAX_MEMBER_BYTES", 100)
        message = r"big\.zip holds a file over .* GiB unzipped.* Unzip it, or gzip it instead"
        with pytest.raises(InputContainerError, match=message):
            detect_parser(path)
        with pytest.raises(InputContainerError, match=message):
            _variants(path)

    def test_size_cap_on_the_bytes_read(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A member that declares less than it holds still stops at the cap."""
        data = FIXTURE.read_bytes() * 200
        path = _understated(tmp_path / "big.zip", data, 50, zipfile.ZIP_DEFLATED)
        monkeypatch.setattr(_helpers, "ZIP_MAX_MEMBER_BYTES", 1000)
        with pytest.raises(InputContainerError, match=r"over .* GiB unzipped"):
            _variants(path)


class TestDamaged:
    def test_missing_directory(self, tmp_path: Path) -> None:
        whole = _zip(tmp_path / "whole.zip", {"genome.txt": FIXTURE.read_bytes()}).read_bytes()
        cut = tmp_path / "cut.zip"
        cut.write_bytes(whole[: len(whole) // 2])
        with pytest.raises(DamagedFileError, match=r"cut\.zip looks truncated or damaged"):
            detect_parser(cut)
        with pytest.raises(DamagedFileError, match="zip directory is missing or corrupt"):
            _variants(cut)

    def test_corrupt_member_data(self, tmp_path: Path) -> None:
        data = FIXTURE.read_bytes()
        path = _zip(tmp_path / "crc.zip", {"genome.txt": data}, zipfile.ZIP_STORED)
        raw = bytearray(path.read_bytes())
        last_row = raw.rindex(data.splitlines()[-1])
        raw[last_row] ^= 0x01  # one flipped bit in the last row: the CRC no longer matches
        path.write_bytes(bytes(raw))
        with pytest.raises(DamagedFileError, match="corrupt zip member"):
            _variants(path)

    @pytest.mark.parametrize("method", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
    @pytest.mark.parametrize("newline", [None, b"\n", b"\r\n"], ids=["fixture", "LF", "CRLF"])
    @pytest.mark.parametrize(
        ("copies", "message"),
        [
            # One read holds it all: the CRC (written for the declared part) fails.
            (1, "Bad CRC-32"),
            # Over several reads, the count passes the declared size first.
            (200, "holds more than the {declared:,} bytes its entry declares"),
        ],
    )
    def test_more_bytes_than_declared(
        self, tmp_path: Path, method: int, copies: int, message: str, newline: bytes | None
    ) -> None:
        """zipfile alone would stop at the declared size and drop the rest."""
        fixture = FIXTURE.read_bytes()
        if newline is not None:
            fixture = fixture.replace(b"\r\n", b"\n").replace(b"\r", b"\n").replace(b"\n", newline)
        data = fixture * copies
        declared = len(data) // 2
        path = _understated(tmp_path / "long.zip", data, declared, method)
        message = message.format(declared=declared)
        with pytest.raises(DamagedFileError, match=message):
            _variants(path)

    @pytest.mark.parametrize("method", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
    def test_fewer_bytes_than_declared(self, tmp_path: Path, method: int) -> None:
        data = FIXTURE.read_bytes()
        path = _zip(tmp_path / "short.zip", {"genome.txt": data}, method)
        raw = bytearray(path.read_bytes())
        for at in (22, raw.index(b"PK\x01\x02") + 24):  # local and central sizes
            struct.pack_into("<I", raw, at, len(data) + 1000)
        path.write_bytes(bytes(raw))
        with pytest.raises(DamagedFileError, match=r"ends after .* bytes its entry declares"):
            _variants(path)

    def test_not_a_zip_is_not_a_container(self, tmp_path: Path) -> None:
        path = tmp_path / "notes.txt"
        path.write_bytes(b"PK is how this note starts, not a zip\n")
        with pytest.raises(ParserNotFoundError):
            detect_parser(path)


def test_cli_prints_a_plain_error(tmp_path: Path) -> None:
    path = _zip(tmp_path / "bundle.zip", {"a.txt": b"x", "b.txt": b"y"})
    result = CliRunner().invoke(main, ["stats", str(path)])
    assert result.exit_code == 1
    assert "Error: bundle.zip is a zip file holding 2 files (a.txt, b.txt)" in result.output
    assert "Traceback" not in result.output


def test_member_is_streamed_not_extracted(tmp_path: Path) -> None:
    path = _zip(tmp_path / "dl.zip", {"genome.txt": FIXTURE.read_bytes()})
    before = sorted(p.name for p in tmp_path.iterdir())
    with _helpers.open_genotype_text(path) as fh:
        # A text stream over the zip member, behind the stray-prefix filter.
        assert isinstance(fh._inner, io.TextIOWrapper)  # type: ignore[attr-defined]
        fh.readline()
    assert sorted(p.name for p in tmp_path.iterdir()) == before

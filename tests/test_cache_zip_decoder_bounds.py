"""ZIP metadata must not permit unbounded internal decompression allocations."""

import shutil
import struct
import zipfile
import zlib

import pytest

from allelix.annotators import gwas
from allelix.databases import pharmgkb_loader

ANN = pharmgkb_loader.CLINICAL_ANN_FILENAME
ALLELES = pharmgkb_loader.CLINICAL_ANN_ALLELES_FILENAME
PATHS = ("clinpgx_hash", "clinpgx_extract", "gwas_extract")


def _source(path, compression, forged):
    with zipfile.ZipFile(path, "w", compression=compression) as archive:
        archive.writestr(ANN, b"A" * (16 * 1024 * 1024) if forged else b"A")
        archive.writestr(ALLELES, b"B")
    if forged:
        with zipfile.ZipFile(path) as archive:
            central = archive.start_dir
            local = archive.getinfo(ANN).header_offset
        raw = bytearray(path.read_bytes())
        # Both headers and the CRC describe the byte zipfile actually returns.
        # No fake ZipInfo or decoder is used: the archive bytes are the input.
        for offset in (central + 16, local + 14):
            struct.pack_into("<I", raw, offset, zlib.crc32(b"A"))
        for offset in (central + 24, local + 22):
            struct.pack_into("<I", raw, offset, 1)
        path.write_bytes(raw)
    return path


def _invoke(kind, source, root, monkeypatch):
    if kind == "clinpgx_hash":
        assert pharmgkb_loader.pharmgkb_content_hash(source, {}) is not None
    elif kind == "clinpgx_extract":
        directory, cleanup = pharmgkb_loader._open_directory(source)
        try:
            assert (directory / ANN).read_bytes() == b"A"
            assert (directory / ALLELES).read_bytes() == b"B"
        finally:
            if cleanup is not None:
                cleanup.cleanup()
    else:
        monkeypatch.setattr(gwas, "download", lambda url, dest: shutil.copyfile(source, dest))
        loaded = []
        monkeypatch.setattr(
            gwas, "load_gwas_tsv", lambda path, *a, **kw: loaded.append(path.read_bytes())
        )
        annotator = gwas.GWASCatalogAnnotator(root)
        monkeypatch.setattr(annotator, "fetch_remote_signal", lambda: "etag:synthetic")
        try:
            annotator.setup()
        finally:
            annotator.close()
        assert loaded == [b"A"]


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("compression", (zipfile.ZIP_BZIP2, zipfile.ZIP_LZMA))
def test_false_small_zip_size_cannot_allocate_a_whole_member(
    tmp_path, monkeypatch, kind, compression
):
    source = _source(tmp_path / "source.zip", compression, forged=True)
    with zipfile.ZipFile(source) as archive:
        assert sum(info.file_size for info in archive.infolist()) == 2
    assert source.stat().st_size < 4096
    real_decompressor = zipfile._get_decompressor
    allocations = []

    class ObservedDecoder:
        def __init__(self, decoder):
            self.decoder = decoder

        def __getattr__(self, name):
            return getattr(self.decoder, name)

        def decompress(self, *args, **kwargs):
            data = self.decoder.decompress(*args, **kwargs)
            allocations.append(len(data))
            return data

    monkeypatch.setattr(
        zipfile, "_get_decompressor", lambda method: ObservedDecoder(real_decompressor(method))
    )
    error = None
    try:
        _invoke(kind, source, tmp_path / "data", monkeypatch)
    except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as exc:
        error = exc
    assert error is not None or allocations, "The source must be refused or actually read"
    assert max(allocations, default=0) <= 8 * 1024 * 1024, (
        "A tiny advertised size must not allocate the whole hidden member inside zipfile; "
        "reject unsupported compression before reading, or bound the decoder itself",
        allocations,
    )


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("compression", (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED))
def test_small_stored_and_deflated_sources_remain_usable(tmp_path, monkeypatch, kind, compression):
    source = _source(tmp_path / "source.zip", compression, forged=False)
    _invoke(kind, source, tmp_path / "data", monkeypatch)

"""Downloaded cache ZIPs are bounded before member bytes are decompressed."""

import contextlib
import hashlib
import json
import shutil
import zipfile

import pytest

from allelix.annotators import gwas
from allelix.databases import pharmgkb_loader

ANN = pharmgkb_loader.CLINICAL_ANN_FILENAME
ALLELES = pharmgkb_loader.CLINICAL_ANN_ALLELES_FILENAME
HUGE = 1 << 60  # Metadata only: no large file or allocation is made.


def _archive(path, members):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, body in members.items():
            archive.writestr(name, body)
    return path


def _instrument(monkeypatch, oversized):
    opened = []
    original_init, original_open = zipfile.ZipFile.__init__, zipfile.ZipFile.open

    def init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        if self.mode == "r" and oversized is not None and oversized in self.NameToInfo:
            self.NameToInfo[oversized].file_size = HUGE

    def open_member(self, name, *args, **kwargs):
        opened.append(name.filename if isinstance(name, zipfile.ZipInfo) else name)
        return original_open(self, name, *args, **kwargs)

    monkeypatch.setattr(zipfile.ZipFile, "__init__", init)
    monkeypatch.setattr(zipfile.ZipFile, "open", open_member)
    return opened


@pytest.mark.parametrize("target", (ANN, ALLELES, "extra.txt"))
def test_clinpgx_unpack_refuses_large_metadata_before_any_member_read(
    tmp_path, monkeypatch, target
):
    source = _archive(
        tmp_path / "clinical.zip", {ANN: b"tiny", ALLELES: b"tiny", "extra.txt": b"tiny"}
    )
    opened = _instrument(monkeypatch, target)
    cleanup = None
    error = None
    try:
        _, cleanup = pharmgkb_loader._open_directory(source)
    except (OSError, ValueError, RuntimeError) as exc:
        error = exc
    finally:
        if cleanup is not None:
            cleanup.cleanup()
    assert not opened, "Inspect every extracted member's size before decompressing any member"
    assert error is not None and any(
        word in str(error).lower() for word in ("size", "large", "limit")
    )


@pytest.mark.parametrize("target", (ANN, ALLELES))
def test_clinpgx_hash_cannot_decompress_an_oversized_tsv_before_unpack_guard(
    tmp_path, monkeypatch, target
):
    source = _archive(tmp_path / "clinical.zip", {ANN: b"one", ALLELES: b"two"})
    opened = _instrument(monkeypatch, target)
    result = None
    with contextlib.suppress(OSError, ValueError, RuntimeError):
        result = pharmgkb_loader.pharmgkb_content_hash(source, {})
    assert target not in opened, "The pre-load content hash must not bypass the archive bound"
    assert result is None, "An oversized source must not produce an accepted content hash"


@pytest.mark.parametrize("member", ("associations.tsv", "nested/associations.tsv"))
def test_gwas_setup_refuses_large_member_before_unpack_or_cache_publication(
    tmp_path, monkeypatch, member
):
    source = _archive(tmp_path / "source.zip", {member: b"tiny table"})
    data = tmp_path / "data"
    data.mkdir()
    db = data / gwas.GWAS_DB_FILENAME
    table = data / gwas.GWAS_TSV_FILENAME
    db.write_bytes(b"previous cache")
    table.write_bytes(b"previous table")
    before = {p.name: (p.read_bytes(), p.stat().st_ino) for p in data.iterdir()}
    monkeypatch.setattr(gwas, "download", lambda url, dest: shutil.copyfile(source, dest))
    loads = []
    monkeypatch.setattr(gwas, "load_gwas_tsv", lambda *args, **kwargs: loads.append(args))
    annotator = gwas.GWASCatalogAnnotator(data)
    monkeypatch.setattr(annotator, "fetch_remote_signal", lambda: "etag:synthetic")
    opened = _instrument(monkeypatch, member)
    error = None
    try:
        annotator.setup()
    except (OSError, ValueError, RuntimeError) as exc:
        error = exc
    finally:
        annotator.close()
    assert not opened, "The selected TSV's declared size must be checked before extraction"
    assert not loads and error is not None
    assert any(word in str(error).lower() for word in ("size", "large", "limit"))
    assert {p.name: (p.read_bytes(), p.stat().st_ino) for p in data.iterdir()} == before


@pytest.mark.parametrize("compression", (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED))
def test_clinpgx_hash_streams_bounded_reads_and_preserves_existing_digest(
    tmp_path, monkeypatch, compression
):
    members = {ANN: b"abc\tdef\n" * 16384, ALLELES: b"genotype\nAA\n" * 8192}
    lookup = {("rs2", "A"): "normal", ("rs1", "G"): "decreased"}
    source = tmp_path / "clinical.zip"
    with zipfile.ZipFile(source, "w", compression=compression) as archive:
        for name, body in members.items():
            archive.writestr(name, body)
    digest = hashlib.sha256()
    for name in (ANN, ALLELES):
        digest.update(f"{name}\0{len(members[name])}\0".encode())
        digest.update(members[name])
    digest.update(json.dumps(sorted(lookup.items())).encode())
    original = zipfile.ZipExtFile.read
    reads = []

    def read(self, n=-1):
        reads.append(n)
        return original(self, n)

    monkeypatch.setattr(zipfile.ZipExtFile, "read", read)
    actual = pharmgkb_loader.pharmgkb_content_hash(source, lookup)
    assert actual == f"sha256:{digest.hexdigest()}"
    assert reads and all(0 < n <= 8 * 1024 * 1024 for n in reads), (
        "Hashing must read bounded chunks, not a whole expanded TSV",
        reads,
    )


def test_small_clinpgx_archive_and_directory_keep_the_same_hash_and_contents(tmp_path):
    members = {ANN: b"one\trow\n", ALLELES: b"second\trow\n"}
    source = _archive(tmp_path / "clinical.zip", members)
    directory, cleanup = pharmgkb_loader._open_directory(source)
    try:
        assert {name: (directory / name).read_bytes() for name in members} == members
        assert pharmgkb_loader.pharmgkb_content_hash(
            source, {}
        ) == pharmgkb_loader.pharmgkb_content_hash(directory, {})
    finally:
        cleanup.cleanup()


def test_missing_clinpgx_member_still_has_no_content_hash(tmp_path):
    source = _archive(tmp_path / "missing.zip", {ANN: b"one"})
    assert pharmgkb_loader.pharmgkb_content_hash(source, {}) is None


@pytest.mark.parametrize("member", ("associations.tsv", "nested/associations.tsv"))
def test_small_gwas_member_still_reaches_loader_and_is_kept(tmp_path, monkeypatch, member):
    source = _archive(tmp_path / "source.zip", {member: b"small table\n"})
    data = tmp_path / "data"
    data.mkdir()
    monkeypatch.setattr(gwas, "download", lambda url, dest: shutil.copyfile(source, dest))
    loaded = []
    monkeypatch.setattr(
        gwas, "load_gwas_tsv", lambda path, *args, **kwargs: loaded.append(path.read_bytes())
    )
    annotator = gwas.GWASCatalogAnnotator(data)
    monkeypatch.setattr(annotator, "fetch_remote_signal", lambda: "etag:synthetic")
    try:
        annotator.setup()
    finally:
        annotator.close()
    assert loaded == [b"small table\n"]
    assert (data / gwas.GWAS_TSV_FILENAME).read_bytes() == b"small table\n"
    assert not list(data.glob(".gwas-*"))

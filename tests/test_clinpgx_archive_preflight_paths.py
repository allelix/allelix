"""Validate the whole downloaded ClinPGx ZIP before reading even the hash tables."""

import shutil
import zipfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.annotators import pharmgkb as pa
from allelix.cli import db as db_cli
from allelix.cli import main
from allelix.config import load_config
from allelix.databases import loader_utils
from allelix.databases import pharmgkb_loader as pl

FIXTURES_DIR = Path(__file__).parent / "fixtures"
EXTRA = "nested/README.txt"


def _archive(path, method=zipfile.ZIP_DEFLATED, extra_first=False):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        if extra_first:
            zf.writestr(EXTRA, b"synthetic extra", compress_type=method)
        for source in sorted((FIXTURES_DIR / "mock_pharmgkb").iterdir()):
            zf.write(source, source.name)
        if not extra_first:
            zf.writestr(EXTRA, b"synthetic extra", compress_type=method)
    return path


def _spy(monkeypatch, bad):
    opened = []
    init, open_member = zipfile.ZipFile.__init__, zipfile.ZipFile.open

    def initialize(zf, *args, **kwargs):
        init(zf, *args, **kwargs)
        if zf.mode != "r":
            return
        if bad == "oversized":
            zf.getinfo(EXTRA).file_size = loader_utils.ZIP_EXTRACT_LIMIT_BYTES + 1
        elif bad == "total":
            # Individually in bounds; the whole archive exceeds the limit.
            for name in (EXTRA, pl.CLINICAL_ANN_FILENAME):
                zf.getinfo(name).file_size = loader_utils.ZIP_EXTRACT_LIMIT_BYTES // 2 + 1

    def read(zf, member, *args, **kwargs):
        opened.append(member.filename if isinstance(member, zipfile.ZipInfo) else member)
        return open_member(zf, member, *args, **kwargs)

    monkeypatch.setattr(zipfile.ZipFile, "__init__", initialize)
    monkeypatch.setattr(zipfile.ZipFile, "open", read)
    return opened


def _bad_archive(tmp_path, bad, order=False):
    method = {"lzma": zipfile.ZIP_LZMA, "bzip2": zipfile.ZIP_BZIP2}.get(bad, zipfile.ZIP_DEFLATED)
    return _archive(tmp_path / "clinical.zip", method, order)


def _snapshot(data):
    return {
        p.name: (p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns)
        for p in data.iterdir()
        if p.is_file()
    }


def _prepare(tmp_path, monkeypatch, archive, prior, lookup):
    data = tmp_path / "data"
    data.mkdir()
    if prior != "new":
        digest = pl.pharmgkb_content_hash(FIXTURES_DIR / "mock_pharmgkb", lookup)
        pl.load_pharmgkb_tsv(
            FIXTURES_DIR / "mock_pharmgkb",
            data / "pharmgkb.sqlite",
            source_url="test://synthetic",
            allele_function_lookup=lookup,
            remote_signal='pgkb:etag:"old"|cpic:lastchange:synthetic',
            content_hash=digest if prior == "unchanged" else "sha256:old-content",
        )
    monkeypatch.setattr(pa, "PHARMGKB_MIN_ROWS", 0)
    monkeypatch.setattr(pa, "download", lambda url, dest: shutil.copyfile(archive, dest))
    monkeypatch.setattr(pa, "fetch_cpic_allele_functions", lambda *a, **kw: dict(lookup))
    monkeypatch.setattr(
        pa.PharmGKBAnnotator,
        "fetch_remote_signal",
        lambda self: 'pgkb:etag:"new"|cpic:lastchange:synthetic',
    )
    return data


@pytest.mark.parametrize("bad", ["lzma", "bzip2", "oversized", "total"])
@pytest.mark.parametrize("extra_first", [False, True])
def test_hash_refuses_bad_extra_before_reading_either_table(
    tmp_path, monkeypatch, bad, extra_first, mock_cpic_lookup
):
    archive = _bad_archive(tmp_path, bad, extra_first)
    opened = _spy(monkeypatch, bad)
    error = None
    try:
        pl.pharmgkb_content_hash(archive, mock_cpic_lookup)
    except ValueError as exc:
        error = exc
    assert not opened, "Even hashing must first validate every member, regardless of archive order"
    assert error is not None and any(s in str(error).lower() for s in ("compression", "limit"))


@pytest.mark.parametrize("bad", ["lzma", "bzip2", "oversized", "total"])
@pytest.mark.parametrize("prior", ["new", "changed", "unchanged"])
@pytest.mark.parametrize("cli", [False, True])
def test_setup_and_update_refuse_before_content_shortcut_and_preserve_cache(
    tmp_path, monkeypatch, bad, prior, cli, mock_cpic_lookup
):
    archive = _bad_archive(tmp_path, bad)
    data = _prepare(tmp_path, monkeypatch, archive, prior, mock_cpic_lookup)
    downloads = []
    if cli:
        # db update does not filter most downloads by the source config toggles.
        # Select the real annotator at its factory boundary; do not stub setup.
        monkeypatch.setattr(
            db_cli, "get_annotators", lambda data_dir, **kwargs: [pa.PharmGKBAnnotator(data_dir)]
        )
        load_config(data)
        local_download = pa.download

        def download(url, dest):
            downloads.append(url)
            return local_download(url, dest)

        monkeypatch.setattr(pa, "download", download)
    before = _snapshot(data)
    opened = _spy(monkeypatch, bad)
    if cli:
        run = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data), "--force"])
        assert downloads == [pa.PHARMGKB_CLINICAL_URL], run.output
        assert run.exit_code != 2, run.output
        assert run.exit_code == 1, run.output
        assert "content unchanged; cache kept" not in run.output
        assert "pharmgkb" in run.output and ("compression" in run.output or "limit" in run.output)
    else:
        with pa.PharmGKBAnnotator(data) as annotator:
            error = None
            try:
                annotator.setup(force=True)
            except ValueError as exc:
                error = exc
            assert error is not None, (
                "An unchanged-table shortcut must not accept an invalid archive"
            )
    assert not opened, "No table may be read before the invalid extra member is refused"
    assert _snapshot(data) == before, "Refusal must not publish data or update its freshness stamp"
    assert not [p for p in data.iterdir() if p.is_dir()], "Owned staging must be cleaned"


@pytest.mark.parametrize("method", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
@pytest.mark.parametrize("prior", ["new", "changed", "unchanged"])
def test_valid_extra_and_same_tables_keep_hash_and_refresh_semantics(
    tmp_path, monkeypatch, method, prior, mock_cpic_lookup
):
    archive = _archive(tmp_path / "clinical.zip", method)
    expected = pl.pharmgkb_content_hash(FIXTURES_DIR / "mock_pharmgkb", mock_cpic_lookup)
    assert pl.pharmgkb_content_hash(archive, mock_cpic_lookup) == expected
    data = _prepare(tmp_path, monkeypatch, archive, prior, mock_cpic_lookup)
    with pa.PharmGKBAnnotator(data) as annotator:
        annotator.setup(force=True)
        assert annotator.is_ready()
        assert annotator.content_unchanged == (prior == "unchanged")
    assert pl.stored_content_hash(data / "pharmgkb.sqlite") == expected

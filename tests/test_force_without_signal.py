# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""`db update --force` for sources whose freshness signal no checksum depends on.

GWAS Catalog and ClinPGx publish no checksum: their downloads are checked by
TLS, Content-Length and the row floor, with or without a freshness signal.
So ``--force`` refreshes them when the signal probe fails, stores no signal,
and the next ordinary update refreshes once more. Without ``--force`` the
cache is kept. Failures come from the probe itself (``head_request_headers``).
"""

from __future__ import annotations

import importlib
import zipfile
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators import gwas as gwas_module
from allelix.annotators import pharmgkb as pharmgkb_module
from allelix.annotators.gwas import GWAS_DB_FILENAME, GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.cli import main
from tests.annotators.test_pharmgkb import _fake_download

if TYPE_CHECKING:
    from pathlib import Path

_HEADERS = {"ETag": '"v1"'}


class _Source:
    """One source's module, class, cache file and fake download."""

    def __init__(self, name: str, mock_gwas_tsv: Path, mock_pharmgkb_dir: Path) -> None:
        self.name = name
        if name == "gwas":
            self.module, self.cls, self.cache = gwas_module, GWASCatalogAnnotator, GWAS_DB_FILENAME

            def fake(url: str, dest: Path) -> None:
                with zipfile.ZipFile(dest, "w") as archive:
                    archive.write(mock_gwas_tsv, arcname="gwas_catalog_associations.tsv")

        else:
            self.module, self.cls, self.cache = (
                pharmgkb_module,
                PharmGKBAnnotator,
                "pharmgkb.sqlite",
            )

            def fake(url: str, dest: Path) -> None:
                _fake_download(mock_pharmgkb_dir, url, dest)

        self.fake = fake


@pytest.fixture(params=["gwas", "pharmgkb"])
def source(request, monkeypatch, mock_gwas_tsv: Path, mock_pharmgkb_dir: Path) -> _Source:
    src = _Source(request.param, mock_gwas_tsv, mock_pharmgkb_dir)
    db_module = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db_module, "get_annotators", lambda data, **kw: [src.cls(data)])
    monkeypatch.setattr(pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:x")
    monkeypatch.setattr(
        pharmgkb_module, "fetch_cpic_allele_functions", lambda: {("rs1801133", "G"): "normal"}
    )
    return src


def _update(data: Path, monkeypatch, src: _Source, *, headers, download, force=False):
    calls: list[str] = []

    def counted(url: str, dest: Path) -> None:
        calls.append(url)
        download(url, dest)

    monkeypatch.setattr(src.module, "head_request_headers", lambda url, *a, **k: headers)
    monkeypatch.setattr(src.module, "download", counted)
    args = ["db", "update", "--data-dir", str(data), *(["--force"] if force else [])]
    result = CliRunner().invoke(main, args)
    return result, " ".join(result.output.split()), calls


def _built(tmp_path: Path, monkeypatch, src: _Source) -> tuple[Path, bytes]:
    data = tmp_path / "data"
    data.mkdir()
    result, output, calls = _update(data, monkeypatch, src, headers=_HEADERS, download=src.fake)
    assert result.exit_code == 0, output
    assert len(calls) == 1
    return data, (data / src.cache).read_bytes()


def test_without_force_the_cache_is_kept(tmp_path, monkeypatch, source) -> None:
    data, before = _built(tmp_path, monkeypatch, source)
    result, output, calls = _update(data, monkeypatch, source, headers=None, download=source.fake)
    assert result.exit_code == 0, output
    assert "freshness can't be verified" in output and "Pass --force to refresh anyway" in output
    assert calls == []
    assert (data / source.cache).read_bytes() == before


def test_force_failure_keeps_the_cache(tmp_path, monkeypatch, source) -> None:
    data, before = _built(tmp_path, monkeypatch, source)

    def broken(url: str, dest: Path) -> None:
        raise OSError("connection reset")

    result, output, _ = _update(
        data, monkeypatch, source, headers=None, download=broken, force=True
    )
    assert result.exit_code == 1
    assert "connection reset" in output
    assert (data / source.cache).read_bytes() == before


def test_force_refreshes_without_a_signal_then_once_more(tmp_path, monkeypatch, source) -> None:
    data, _ = _built(tmp_path, monkeypatch, source)
    result, output, calls = _update(
        data, monkeypatch, source, headers=None, download=source.fake, force=True
    )
    assert result.exit_code == 0, output
    assert len(calls) == 1
    assert "refreshed without a freshness signal; none is stored" in output
    with source.cls(data) as ann:
        assert ann.is_ready() and ann.cached_remote_signal() is None
    # The signal is back: one more refresh stores it, then the cache is current.
    result, output, calls = _update(
        data, monkeypatch, source, headers=_HEADERS, download=source.fake
    )
    assert result.exit_code == 0, output
    assert "predates the freshness signal" in output and len(calls) == 1
    result, output, calls = _update(
        data, monkeypatch, source, headers=_HEADERS, download=source.fake
    )
    assert "already current" in output and calls == []


def test_clinpgx_force_without_signal_and_without_cpic(tmp_path, monkeypatch, mock_pharmgkb_dir):
    """N1 still holds in the double failure: no ClinPGx signal and no CPIC data.

    The stamp keeps only ``cpic:unavailable``, so status shows "(no CPIC)" and
    the next ordinary update, with both back, refreshes.
    """
    import urllib.error

    src = _Source("pharmgkb", tmp_path / "unused.tsv", mock_pharmgkb_dir)
    db_module = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db_module, "get_annotators", lambda data, **kw: [PharmGKBAnnotator(data)])

    def cpic_down():
        raise urllib.error.URLError("CPIC is down")

    monkeypatch.setattr(pharmgkb_module, "fetch_cpic_allele_functions", cpic_down)
    monkeypatch.setattr(pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:x")
    data = tmp_path / "data"
    data.mkdir()
    result, output, calls = _update(
        data, monkeypatch, src, headers=None, download=src.fake, force=True
    )
    assert result.exit_code == 0, output
    assert len(calls) == 1
    assert "refreshed without a freshness signal; none is stored" in output
    assert "CPIC unavailable" in output
    with PharmGKBAnnotator(data) as ann:
        assert ann.is_ready() and ann.cached_remote_signal() == "cpic:unavailable"
    status = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])
    assert "(no CPIC)" in status.output

    monkeypatch.setattr(
        pharmgkb_module, "fetch_cpic_allele_functions", lambda: {("rs1801133", "G"): "normal"}
    )
    result, output, calls = _update(data, monkeypatch, src, headers=_HEADERS, download=src.fake)
    assert result.exit_code == 0, output
    assert "remote signal changed; refreshing" in output and len(calls) == 1
    with PharmGKBAnnotator(data) as ann:
        assert ann.cached_remote_signal() == 'pgkb:etag:"v1"|cpic:lastchange:x'

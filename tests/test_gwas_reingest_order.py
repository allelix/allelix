"""O5: a GWAS re-ingest copies the kept table before it reads the cache's stamp.

A refresh that finishes between the copy and the read is seen in the read: its
hash does not match the copied table, so nothing is rebuilt and its own cache
and signal stay. Production setup, parser, loader and hashing; only transport
and scheduling are substituted.
"""

import contextlib
import os
import shutil
import sqlite3
import zipfile
from pathlib import Path

import pytest

from allelix.annotators import gwas
from tests.test_cache_publication_consistency import _gwas_snapshot


@pytest.mark.parametrize("hash_state", ["missing_column", "null_hash", "recorded_hash"])
@pytest.mark.parametrize("copy_only", [False, True])
def test_a_refresh_between_the_copy_and_the_stamp_read_is_not_rebuilt_over(
    tmp_path, monkeypatch, hash_state, copy_only
):
    fixture = Path(__file__).resolve().parent / "fixtures/mock_gwas_catalog.tsv"
    lines = fixture.read_text().splitlines()
    header, fields = lines[0], lines[1].split("\t")
    trait_column = header.split("\t").index("DISEASE/TRAIT")
    label = "A"
    db = tmp_path / "gwas.sqlite"

    def download(url, path):
        row = fields.copy()
        row[trait_column] = "Synthetic " + label
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("rows.tsv", header + "\n" + "\t".join(row) + "\n")

    monkeypatch.setattr(gwas, "download", download)
    monkeypatch.setattr(gwas, "GWAS_MIN_ROWS", 1)
    monkeypatch.setattr(
        gwas.GWASCatalogAnnotator, "fetch_remote_signal", lambda self: "etag:" + label
    )
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        annotator.setup()
    with contextlib.closing(sqlite3.connect(db)) as conn, conn:
        if hash_state == "missing_column":
            conn.execute("ALTER TABLE database_versions DROP COLUMN content_hash")
        elif hash_state == "null_hash":
            conn.execute("UPDATE database_versions SET content_hash=NULL")
        conn.execute("UPDATE database_versions SET local_version_tag='cv:0'")

    fired = False

    def concurrent_refresh():
        nonlocal label, fired
        fired = True
        label = "B"
        with gwas.GWASCatalogAnnotator(tmp_path) as concurrent:
            concurrent.setup()

    real_link, real_copy = os.link, shutil.copyfile

    def link(src, dst, *args, **kwargs):
        if copy_only:
            raise OSError("synthetic filesystem without hard links")
        real_link(src, dst, *args, **kwargs)
        if Path(src).name == gwas.GWAS_TSV_FILENAME and not fired:
            concurrent_refresh()

    def copyfile(src, dst, *args, **kwargs):
        result = real_copy(src, dst, *args, **kwargs)
        if Path(src).name == gwas.GWAS_TSV_FILENAME and not fired:
            concurrent_refresh()
        return result

    monkeypatch.setattr(os, "link", link)
    monkeypatch.setattr(shutil, "copyfile", copyfile)
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        assert annotator.prepare_cache() is False
    assert fired
    assert _gwas_snapshot(db) == ([("Synthetic B",)], "etag:B")
    assert not list(tmp_path.glob(".gwas-*"))

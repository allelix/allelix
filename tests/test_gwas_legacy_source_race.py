"""O5 r2: legacy source metadata must not permit cross-refresh re-ingestion."""

import contextlib
import os
import sqlite3
import zipfile
from pathlib import Path

import pytest

from allelix.annotators import gwas
from tests.test_cache_publication_consistency import _gwas_snapshot


@pytest.mark.parametrize("hash_state", ["missing_column", "null_hash", "recorded_hash"])
@pytest.mark.parametrize("copy_only", [False, True])
@pytest.mark.parametrize("first,second", [("A", "B"), ("B", "A")])
def test_legacy_reingest_never_labels_a_new_table_with_the_old_signal(
    tmp_path, monkeypatch, hash_state, copy_only, first, second
):
    fixture = Path(gwas.__file__).resolve().parents[2] / "tests/fixtures/mock_gwas_catalog.tsv"
    lines = fixture.read_text().splitlines()
    header, fields = lines[0], lines[1].split("\t")
    trait_column = header.split("\t").index("DISEASE/TRAIT")
    label = first
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
    before = db.read_bytes()
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        assert annotator.is_ready(), "missing optional metadata must not invalidate a cache"
    assert db.read_bytes() == before
    with contextlib.closing(sqlite3.connect(db)) as conn, conn:
        conn.execute("UPDATE database_versions SET local_version_tag='cv:0'")

    if copy_only:

        def no_hard_links(*args, **kwargs):
            raise OSError("synthetic filesystem without hard links")

        monkeypatch.setattr(os, "link", no_hard_links)
    read_stamp = gwas.retained_source_stamp
    fired = False

    def interleaved_stamp(path):
        nonlocal label, fired
        value = read_stamp(path)
        if not fired:
            fired = True
            label = second
            with gwas.GWASCatalogAnnotator(tmp_path) as concurrent:
                concurrent.setup()
        return value

    monkeypatch.setattr(gwas, "retained_source_stamp", interleaved_stamp)
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        annotator.prepare_cache()
    assert fired
    # A safe refusal may retain the concurrent cache; a completed re-ingest
    # may publish its own pinned older source. Either must keep rows and
    # signal from ONE refresh, including when the old cache had no hash.
    traits, signal = _gwas_snapshot(db)
    assert signal in ("etag:" + first, "etag:" + second)
    assert traits == [("Synthetic " + signal.removeprefix("etag:"),)], (
        "legacy re-ingest paired another refresh's table with the old signal"
    )
    assert not list(tmp_path.glob(".gwas-*"))

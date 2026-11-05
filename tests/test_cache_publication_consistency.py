"""O5 r1: source bytes, installed rows and freshness metadata stay together.

Only transport and scheduling are substituted. Production parsers, loaders,
content hashing and SQLite publication run on tiny synthetic repository fixtures.
"""

import contextlib
import os
import sqlite3
import zipfile
from pathlib import Path

import pytest

from allelix.annotators import gwas, pharmgkb
from tests.test_clinpgx_content_refresh import Remote, _rows, _stamp


def _gwas_snapshot(db):
    with contextlib.closing(sqlite3.connect(db)) as conn:
        traits = conn.execute("SELECT trait FROM gwas_associations ORDER BY trait").fetchall()
        signal = conn.execute(
            "SELECT remote_signal FROM database_versions WHERE name='gwas'"
        ).fetchone()[0]
    return traits, signal


@pytest.mark.parametrize("outer,inner", [("A", "B"), ("B", "A")])
@pytest.mark.parametrize("nested_after_publication", [False, True])
def test_gwas_offline_reingest_never_pairs_another_refreshs_rows_and_stamp(
    tmp_path, monkeypatch, outer, inner, nested_after_publication
):
    fixture = Path(gwas.__file__).resolve().parents[2] / "tests/fixtures/mock_gwas_catalog.tsv"
    lines = fixture.read_text().splitlines()
    header, fields = lines[0], lines[1].split("\t")
    trait_column = header.split("\t").index("DISEASE/TRAIT")
    label = outer
    fired = False
    db = tmp_path / "gwas.sqlite"
    real_replace = os.replace

    def download(url, path):
        row = fields.copy()
        row[trait_column] = "Synthetic " + label
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("rows.tsv", header + "\n" + "\t".join(row) + "\n")

    def nested():
        nonlocal label
        label = inner
        try:
            with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
                annotator.setup()
        finally:
            label = outer

    def replace(src, dst):
        nonlocal fired
        if Path(dst) == db and not fired:
            fired = True
            if nested_after_publication:
                real_replace(src, dst)
                nested()
                return
            nested()
        real_replace(src, dst)

    monkeypatch.setattr(gwas, "download", download)
    monkeypatch.setattr(gwas, "GWAS_MIN_ROWS", 1)
    monkeypatch.setattr(
        gwas.GWASCatalogAnnotator, "fetch_remote_signal", lambda self: "etag:" + label
    )
    monkeypatch.setattr(os, "replace", replace)
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        annotator.setup()
    assert fired
    winner = inner if nested_after_publication else outer
    expected = ([("Synthetic " + winner,)], "etag:" + winner)
    assert _gwas_snapshot(db) == expected

    # A later interpreter change must read the source belonging to this cache.
    # Refusing an unavailable/mismatched retained source is also safe: it must
    # never install another source's rows under this cache's freshness signal.
    with contextlib.closing(sqlite3.connect(db)) as conn, conn:
        conn.execute("UPDATE database_versions SET local_version_tag='cv:0'")
    with gwas.GWASCatalogAnnotator(tmp_path) as annotator:
        annotator.prepare_cache()
    assert _gwas_snapshot(db) == expected, (
        "re-ingest combined one refresh's source with another's stamp"
    )


@pytest.mark.parametrize("changed_input", ["cpic", "annotations"])
@pytest.mark.parametrize("inner_fails", [False, True])
def test_clinpgx_unchanged_content_cannot_restamp_a_concurrent_refresh(
    tmp_path, monkeypatch, changed_input, inner_fails
):
    remote = Remote()
    monkeypatch.setattr(pharmgkb, "head_request_headers", remote.headers)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_remote_signal", lambda: remote.cpic_date)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", remote.cpic)
    monkeypatch.setattr(pharmgkb, "download", remote.download)
    monkeypatch.setattr(pharmgkb, "PHARMGKB_MIN_ROWS", 1)
    with pharmgkb.PharmGKBAnnotator(tmp_path) as annotator:
        annotator.setup()
    original_rows = _rows(tmp_path)
    original_hash = _stamp(tmp_path)["content_hash"]
    remote.etag = '"repacked-A"'  # New signal, but exactly the installed input content.
    real_hash = pharmgkb.stored_content_hash
    fired = False
    winner_rows = None
    winner_stamp = None

    def stored_hash(path):
        nonlocal fired, winner_rows, winner_stamp
        value = real_hash(path)
        if not fired:
            fired = True
            remote.etag = '"changed-B"'
            if changed_input == "cpic":
                remote.lookup[("rs1801133", "A")] = "normal"
                remote.cpic_date = "lastchange:2020-02-01"
            else:
                name = "clinical_ann_alleles.tsv"
                remote.files[name] = remote.files[name].replace(
                    b"Heterozygous carriers", b"Synthetic changed annotation for carriers"
                )
            with pharmgkb.PharmGKBAnnotator(tmp_path) as inner:
                if inner_fails:
                    with monkeypatch.context() as scope:
                        scope.setattr(pharmgkb, "PHARMGKB_MIN_ROWS", 999_999)
                        with pytest.raises(OSError, match="floor"):
                            inner.setup()
                else:
                    inner.setup()
            winner_rows = _rows(tmp_path)
            winner_stamp = _stamp(tmp_path)
        return value

    monkeypatch.setattr(pharmgkb, "stored_content_hash", stored_hash)
    with pharmgkb.PharmGKBAnnotator(tmp_path) as annotator:
        annotator.setup()
    assert fired
    actual = _stamp(tmp_path)
    if inner_fails:
        assert _rows(tmp_path) == original_rows
        assert actual["content_hash"] == original_hash
        assert 'pgkb:etag:"repacked-A"' in actual["remote_signal"]
    else:
        assert winner_rows != original_rows, "the concurrent build must change actual rows"
        assert winner_stamp["content_hash"] != original_hash
        assert _rows(tmp_path) == winner_rows
        assert actual["content_hash"] == winner_stamp["content_hash"]
        assert actual["remote_signal"] == winner_stamp["remote_signal"], (
            "unchanged-content shortcut stamped a different writer's installed content"
        )

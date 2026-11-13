"""Every cache publication settles the previous database before exposing a new one."""

import contextlib
import gzip
import importlib.util
import os
import sqlite3
from pathlib import Path

import pytest
from test_cache_staging_ownership import harness

from allelix.databases.rsid_identity_loader import write_identity_db

KINDS = (
    "gwas",
    "pharmgkb",
    "clinvar",
    "identity",
    "prebuilt",
    "offline-gnomad",
    "offline-alphamissense",
    "offline-cadd",
)
SIDES = ("-journal", "-wal", "-shm")


def _read(db):
    with contextlib.closing(sqlite3.connect(db)) as conn:
        assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        return {
            row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }


def _builder(root, monkeypatch, kind):
    if kind in ("gwas", "pharmgkb", "clinvar", "prebuilt"):
        run, _stamp, db = harness(root, monkeypatch, kind)
        return lambda: run("A"), db
    root.mkdir()
    db = root / "cache.sqlite"
    if kind == "identity":
        return lambda: write_identity_db(
            [("rs900000001", "1", 100, "G", "A", 200, "G", "A", 0)],
            db,
            version="synthetic",
            source_url="test://synthetic",
            downloaded_at="test",
        ), db
    name = kind.removeprefix("offline-")
    script = Path(__file__).resolve().parents[1] / "scripts" / f"build_{name}_cache.py"
    spec = importlib.util.spec_from_file_location(f"publication_{name}", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if name == "gnomad":
        source = root / module.GNOMAD_FILENAME.format(chrom="1")
        with gzip.open(source, "wt") as f:
            f.write(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t100\trs100\tA\tG\t.\tPASS\tAF=0.05\n"
            )
        return lambda: module.build_cache(db, full=True, chromosomes=["1"], local_dir=root), db
    if name == "alphamissense":
        from tests.databases.test_build_alphamissense_cache import _TSV_HEADER, _make_tsv_line

        source = root / "source.tsv"
        source.write_text(_TSV_HEADER + _make_tsv_line())
        return lambda: module.build_cache(
            db, root / "unused", tsv_path=source, skip_gnomad=True
        ), db
    from tests.databases.test_build_cadd_cache import (
        _CADD_HEADER,
        _make_cadd_gz,
        _make_cadd_line,
        _make_gnomad_db,
    )

    _make_gnomad_db(root, [("1", 100, "A", "G")])
    source = _make_cadd_gz(root, [_CADD_HEADER, _make_cadd_line("1", 100, "A", "G")])
    return lambda: module.build_cache(source, None, db, root), db


def _old_state(db, kind):
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("CREATE TABLE publication_old_only (value TEXT)")
        conn.executemany(
            "INSERT INTO publication_old_only VALUES (?)", [("committed" * 60,)] * 500
        )
        conn.commit()
        if kind == "wal":
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            conn.execute("UPDATE publication_old_only SET value='committed'")
            conn.commit()
        else:
            conn.execute("PRAGMA cache_size=5")
            conn.execute("UPDATE publication_old_only SET value=?", ("uncommitted" * 60,))
        saved = {db: db.read_bytes()}
        for suffix in SIDES:
            side = Path(str(db) + suffix)
            if side.exists():
                saved[side] = side.read_bytes()
        assert Path(str(db) + ("-wal" if kind == "wal" else "-journal")) in saved
        if kind != "wal":
            conn.rollback()
    for path, content in saved.items():
        path.write_bytes(content)
    # Independently prove recovery of exactly this old pair/triple.
    control = db.with_name("recovery-control.sqlite")
    for path, content in saved.items():
        suffix = str(path)[len(str(db)) :]
        Path(str(control) + suffix).write_bytes(content)
    assert "publication_old_only" in _read(control)
    with contextlib.closing(sqlite3.connect(control)) as conn:
        assert all(
            not value.startswith("uncommitted")
            for (value,) in conn.execute("SELECT DISTINCT value FROM publication_old_only")
        )
    control.unlink()


def _sides(db):
    return [p for suffix in SIDES if (p := Path(str(db) + suffix)).exists()]


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("state", ("hot", "wal", "recovery-refused", "orphan"))
def test_install_never_exposes_new_cache_beside_old_recovery_bytes(
    tmp_path, monkeypatch, kind, state
):
    run, db = _builder(tmp_path / "cache", monkeypatch, kind)
    run()
    _old_state(db, "wal" if state == "wal" else "hot")
    if state == "orphan":
        db.unlink()
    original_connect, original_replace = sqlite3.connect, os.replace
    publications = []
    exposed = []

    def connect(path, *args, **kwargs):
        if state == "recovery-refused" and str(path) == str(db):
            raise sqlite3.OperationalError("synthetic recovery unavailable")
        return original_connect(path, *args, **kwargs)

    def replace(src, dst, *args, **kwargs):
        original_replace(src, dst, *args, **kwargs)
        if Path(dst) == db:
            publications.append(True)
            exposed.extend(_sides(db))  # A reader can run before post-rename cleanup.

    error = None
    with monkeypatch.context() as patch:
        patch.setattr(sqlite3, "connect", connect)
        patch.setattr(os, "replace", replace)
        try:
            run()
        except (OSError, RuntimeError, sqlite3.Error) as exc:
            error = exc
    assert not exposed, "Old journal/WAL bytes were exposed beside the new cache"
    if error is not None:
        assert not publications, "Refusal must happen before publication"
        assert "publication_old_only" in _read(db), (
            "Failed publication must keep committed old data"
        )
    else:
        assert state != "recovery-refused", "Failed recovery must not report installation success"
        assert publications and not _sides(db)
        assert "publication_old_only" not in _read(db)


@pytest.mark.parametrize("kind", KINDS)
def test_live_writer_is_refused_without_replacing_its_database(tmp_path, monkeypatch, kind):
    run, db = _builder(tmp_path / "cache", monkeypatch, kind)
    run()
    inode = db.stat().st_ino
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("CREATE TABLE publication_old_only (value TEXT)")
        conn.commit()
        conn.execute("INSERT INTO publication_old_only VALUES ('uncommitted')")
        assert _sides(db)
        error = None
        try:
            run()
        except (OSError, RuntimeError, sqlite3.Error) as exc:
            error = exc
        try:
            assert error is not None, "An active old writer must not be replaced as a success"
            assert db.stat().st_ino == inode
        finally:
            conn.rollback()
    assert "publication_old_only" in _read(db)


@pytest.mark.parametrize("kind", KINDS)
def test_failed_final_move_keeps_recoverable_old_data(tmp_path, monkeypatch, kind):
    run, db = _builder(tmp_path / "cache", monkeypatch, kind)
    run()
    _old_state(db, "hot")
    original_replace = os.replace

    def replace(src, dst, *args, **kwargs):
        if Path(dst) == db:
            raise PermissionError("synthetic final move refused")
        return original_replace(src, dst, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(os, "replace", replace)
        with pytest.raises((OSError, RuntimeError, sqlite3.Error)):
            run()
    assert "publication_old_only" in _read(db)
    assert not list(db.parent.glob(".cache.sqlite.*.tmp"))


@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize("state", ("hot", "wal", "recovery-refused"))
def test_clinvar_setup_final_publication_has_the_same_recovery_rule(
    tmp_path, monkeypatch, build, state
):
    import hashlib

    from click.testing import CliRunner
    from test_clinvar_damaged_cache_replaced import _clinvar_sources

    from allelix.annotators import clinvar as module
    from allelix.cli import main
    from allelix.databases import manager

    payloads = _clinvar_sources(tmp_path, build)
    monkeypatch.setattr(
        module,
        "_fetch_checksum",
        lambda url: hashlib.md5(payloads[url.rsplit("/", 1)[-1]]).hexdigest(),
    )
    monkeypatch.setattr(
        module, "download", lambda url, dest: dest.write_bytes(payloads[url.rsplit("/", 1)[-1]])
    )
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", {"GRCh37": 1, "GRCh38": 1})
    data = tmp_path / "data"
    data.mkdir()
    with module.ClinVarAnnotator(data, builds=(build,)) as annotator:
        annotator.setup()
    db = data / module.clinvar_db_filename(build)
    original_build = module.ClinVarAnnotator._setup_one_from_tsv
    original_connect, original_replace = sqlite3.connect, os.replace
    injected, publications, exposed = [], [], []

    def build_candidate(self, *args, **kwargs):
        original_build(self, *args, **kwargs)
        # A previous user can leave recovery bytes after the baseline count was
        # read, while a new candidate is built. Exercise the final move too.
        _old_state(db, "wal" if state == "wal" else "hot")
        injected.append(True)

    def connect(path, *args, **kwargs):
        if injected and state == "recovery-refused" and str(path) == str(db):
            raise sqlite3.OperationalError("synthetic recovery unavailable")
        return original_connect(path, *args, **kwargs)

    def replace(src, dst, *args, **kwargs):
        original_replace(src, dst, *args, **kwargs)
        if Path(dst) == db:
            publications.append(True)
            exposed.extend(_sides(db))

    with monkeypatch.context() as patch:
        patch.setattr(module.ClinVarAnnotator, "_setup_one_from_tsv", build_candidate)
        patch.setattr(sqlite3, "connect", connect)
        patch.setattr(os, "replace", replace)
        result = CliRunner().invoke(
            main,
            [
                "db",
                "update",
                "--data-dir",
                str(data),
                "--only",
                "clinvar",
                "--build",
                build.lower(),
                "--force",
            ],
        )
    assert injected, (result.output, result.exception)
    assert not exposed, "ClinVar's final candidate move bypassed recovery"
    if result.exit_code == 0:
        assert state != "recovery-refused"
        assert publications and "publication_old_only" not in _read(db)
    else:
        assert not publications
        assert "publication_old_only" in _read(db)
        assert "error" in result.output.lower() or "fail" in result.output.lower()

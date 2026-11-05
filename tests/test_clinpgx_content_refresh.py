# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Refresh decisions depend on every build input, not the ZIP's packaging."""

from __future__ import annotations

import contextlib
import hashlib
import importlib
import shutil
import sqlite3
import zipfile
from pathlib import Path

import pytest
from click.testing import CliRunner

import allelix
from allelix.annotators import pharmgkb
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.cli import main

_FIXTURES = Path(allelix.__file__).resolve().parent.parent / "tests/fixtures/mock_pharmgkb"
_LEGACY_FIELDS = {
    "name",
    "source_url",
    "version",
    "downloaded_at",
    "record_count",
    "remote_signal",
    "local_version_tag",
}
_LOOKUP = {("rs1801133", "A"): "decreased", ("rs1801133", "G"): "normal"}


class Remote:
    def __init__(self):
        self.etag = '"first"'
        self.cpic_date = "lastchange:2020-01-01"
        self.lookup = dict(_LOOKUP)
        self.files = {p.name: p.read_bytes() for p in _FIXTURES.glob("*.tsv")}
        self.reverse = False
        self.compression = zipfile.ZIP_STORED
        self.year = 2020
        self.downloads = 0
        self.builds = 0
        self.heads = 0
        self.cpic_fetches = 0

    def headers(self, *args, **kwargs):
        self.heads += 1
        return {"ETag": self.etag, "Last-Modified": "Wed, 01 Jan 2020 00:00:00 GMT"}

    def download(self, url, dest):
        self.downloads += 1
        with zipfile.ZipFile(dest, "w") as archive:
            for name in sorted(self.files, reverse=self.reverse):
                entry = zipfile.ZipInfo(name, date_time=(self.year, 1, 1, 0, 0, 0))
                entry.compress_type = self.compression
                archive.writestr(entry, self.files[name])

    def cpic(self):
        self.cpic_fetches += 1
        return dict(self.lookup)

    def reset_counts(self):
        self.downloads = self.builds = self.heads = self.cpic_fetches = 0

    @property
    def signal(self):
        return f"pgkb:etag:{self.etag}|cpic:{self.cpic_date}"


@pytest.fixture
def remote(monkeypatch):
    remote = Remote()
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db, "get_annotators", lambda data, **kwargs: [PharmGKBAnnotator(data)])
    monkeypatch.setattr(pharmgkb, "head_request_headers", remote.headers)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_remote_signal", lambda: remote.cpic_date)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", remote.cpic)
    monkeypatch.setattr(pharmgkb, "download", remote.download)
    monkeypatch.setattr(pharmgkb, "PHARMGKB_MIN_ROWS", 1)
    real_load = pharmgkb.load_pharmgkb_tsv

    def build(*args, **kwargs):
        remote.builds += 1
        return real_load(*args, **kwargs)

    monkeypatch.setattr(pharmgkb, "load_pharmgkb_tsv", build)
    return remote


def _assert_signal_parts(actual, expected):
    """Keep the source and CPIC evidence exact, allowing another HTTP field."""
    assert actual is not None
    parts = actual.split("|")
    source, cpic = expected.split("|")
    prefix = ":".join(source.split(":")[:2]) + ":"
    assert [part for part in parts if part.startswith(prefix)] == [source]
    assert [part for part in parts if part.startswith("cpic:")] == [cpic]


def _update(data):
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
    assert result.exit_code == 0, (result.output, result.exception)
    with PharmGKBAnnotator(data) as ann:
        assert ann.is_ready()
    return " ".join(result.output.split()).lower()


def _install(tmp_path, remote):
    data = tmp_path / "data"
    _update(data)
    assert remote.downloads == remote.builds == 1
    remote.reset_counts()
    return data


def _snapshot(data):
    return {
        p.name: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
        for p in data.iterdir()
        if p.is_file()
    }


def _stamp(data):
    with contextlib.closing(
        sqlite3.connect(f"file:{data / 'pharmgkb.sqlite'}?mode=ro", uri=True)
    ) as conn:
        conn.row_factory = sqlite3.Row
        return dict(
            conn.execute("SELECT * FROM database_versions WHERE name='pharmgkb'").fetchone()
        )


def _digest_stamp(data):
    # HTTP freshness metadata is optional too; it is not the content digest.
    value = _stamp(data).get("content_hash")
    assert value, "No content hash persisted in the stamp"
    return value


def _rows(data):
    with contextlib.closing(
        sqlite3.connect(f"file:{data / 'pharmgkb.sqlite'}?mode=ro", uri=True)
    ) as conn:
        return (
            sorted(conn.execute("SELECT * FROM pharmgkb_annotations").fetchall()),
            sorted(conn.execute("SELECT * FROM pharmgkb_allele_function").fetchall()),
        )


def _legacy_stamp(data):
    """A published cache from before the optional hash existed, with its tag kept."""
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.execute("ALTER TABLE database_versions RENAME TO saved_versions")
        conn.execute(
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT, "
            "version TEXT, downloaded_at TEXT, record_count INTEGER, "
            "remote_signal TEXT, local_version_tag TEXT)"
        )
        columns = ", ".join(sorted(_LEGACY_FIELDS))
        conn.execute(
            f"INSERT INTO database_versions ({columns}) SELECT {columns} FROM saved_versions"
        )
        conn.execute("DROP TABLE saved_versions")
        conn.commit()


def _forbid_copy(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("An unchanged cache was copied")

    for name in ("copy", "copy2", "copyfile", "copytree"):
        monkeypatch.setattr(shutil, name, forbidden)


@pytest.mark.parametrize("legacy", [False, True])
def test_unchanged_signal_only_checks_remote_metadata(tmp_path, monkeypatch, remote, legacy):
    data = _install(tmp_path, remote)
    if legacy:
        _legacy_stamp(data)
    before = _snapshot(data)
    _forbid_copy(monkeypatch)
    output = _update(data)
    assert "already current" in output
    assert remote.heads == 1
    assert (remote.downloads, remote.cpic_fetches, remote.builds) == (0, 0, 0)
    assert _snapshot(data) == before


@pytest.mark.parametrize("packaging", ["identical", "zip-metadata", "unused-file", "cpic-order"])
def test_changed_etag_same_inputs_downloads_once_without_rebuild_or_copy(
    tmp_path, monkeypatch, remote, packaging
):
    data = _install(tmp_path, remote)
    rows = _rows(data)
    prior_digest = _stamp(data).get("content_hash")
    inode = (data / "pharmgkb.sqlite").stat().st_ino
    remote.etag = '"second"'
    if packaging == "zip-metadata":
        remote.reverse = True
        remote.compression = zipfile.ZIP_DEFLATED
        remote.year = 2021
    elif packaging == "unused-file":
        remote.files["unused.tsv"] = b"This file is not a build input.\n"
    elif packaging == "cpic-order":
        remote.lookup = dict(reversed(list(remote.lookup.items())))
    _forbid_copy(monkeypatch)
    output = _update(data)
    assert remote.downloads == 1
    assert remote.builds == 0, "Same build inputs caused a rebuild"
    assert "content unchanged" in output
    assert _rows(data) == rows
    assert (data / "pharmgkb.sqlite").stat().st_ino == inode
    _assert_signal_parts(_stamp(data)["remote_signal"], remote.signal)
    # A fresh annotator must see the new signal; repeated runs cannot redownload it.
    digest = _digest_stamp(data)
    assert digest == prior_digest
    remote.reset_counts()
    assert "already current" in _update(data)
    assert remote.downloads == remote.builds == remote.cpic_fetches == 0
    assert _digest_stamp(data) == digest


@pytest.mark.parametrize("part", ["annotations", "alleles", "cpic-only", "etag-and-cpic"])
def test_every_changed_build_input_rebuilds_and_persists_its_hash(tmp_path, remote, part):
    data = _install(tmp_path, remote)
    old_digest = _digest_stamp(data)
    before = _rows(data)
    if part in ("annotations", "alleles"):
        remote.etag = '"second"'
        filename = (
            "clinical_annotations.tsv" if part == "annotations" else "clinical_ann_alleles.tsv"
        )
        old = (
            b"Increased toxicity risk"
            if part == "annotations"
            else b"Heterozygous carriers may have elevated toxicity risk"
        )
        assert old in remote.files[filename]
        remote.files[filename] = remote.files[filename].replace(old, b"Changed synthetic response")
    else:
        remote.cpic_date = "lastchange:2020-02-01"
        remote.lookup[("rs1801133", "A")] = "normal"
        if part == "etag-and-cpic":
            remote.etag = '"second"'
    _update(data)
    assert (remote.downloads, remote.builds) == (1, 1)
    assert _rows(data) != before
    assert _digest_stamp(data) != old_digest
    _assert_signal_parts(_stamp(data)["remote_signal"], remote.signal)
    # Prove the replacement's digest is actually used by the next refresh.
    remote.etag = '"third"'
    remote.reset_counts()
    assert "content unchanged" in _update(data)
    assert (remote.downloads, remote.builds) == (1, 0)


def test_old_stamp_is_kept_then_learns_hash_on_a_real_content_refresh(tmp_path, remote):
    data = _install(tmp_path, remote)
    _legacy_stamp(data)
    before = _snapshot(data)
    with PharmGKBAnnotator(data) as ann:
        assert ann.is_ready(), "Missing optional hash invalidated a published cache"
    assert "already current" in _update(data)
    assert _snapshot(data) == before
    assert remote.downloads == remote.builds == 0
    remote.etag = '"second"'
    old = b"Heterozygous carriers may have elevated toxicity risk"
    name = "clinical_ann_alleles.tsv"
    remote.files[name] = remote.files[name].replace(old, b"Changed synthetic response")
    _update(data)
    assert remote.downloads == remote.builds == 1
    _digest_stamp(data)
    remote.reset_counts()
    remote.etag = '"third"'
    assert "content unchanged" in _update(data)
    assert (remote.downloads, remote.builds) == (1, 0)


def test_same_bytes_do_not_keep_an_outdated_interpreter(tmp_path, remote):
    data = _install(tmp_path, remote)
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.execute("UPDATE database_versions SET local_version_tag='iv:0'")
        conn.commit()
    with PharmGKBAnnotator(data) as ann:
        assert not ann.is_ready()
    # Content equality is not permission to retain an obsolete interpretation.
    _update(data)
    assert (remote.downloads, remote.builds) == (1, 1)

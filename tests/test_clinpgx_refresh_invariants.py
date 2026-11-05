# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Freshness invariants independent of how an ETag-only change is decided."""

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
from allelix.config import load_config
from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

_FIXTURES = Path(allelix.__file__).resolve().parent.parent / "tests" / "fixtures" / "mock_pharmgkb"
_OLD_DATE = "lastchange:2020-01-01"
_NEW_DATE = "lastchange:2020-02-01"
_OLD_LM = "Wed, 01 Jan 2020 00:00:00 GMT"
_NEW_LM = "Sat, 01 Feb 2020 00:00:00 GMT"
_LOOKUP = {("rs1801133", "A"): "decreased", ("rs1801133", "G"): "normal"}


def _signal(kind, *, changed=False, cpic_changed=False):
    pg = (
        ('etag:"new"' if changed else 'etag:"old"')
        if kind == "etag"
        else ("lm:" + (_NEW_LM if changed else _OLD_LM))
    )
    return f"pgkb:{pg}|cpic:{_NEW_DATE if cpic_changed else _OLD_DATE}"


def _built(tmp_path, kind):
    data = tmp_path / "data"
    data.mkdir()
    load_pharmgkb_tsv(
        _FIXTURES,
        data / "pharmgkb.sqlite",
        version="synthetic-old",
        remote_signal=_signal(kind),
        allele_function_lookup=_LOOKUP,
        min_rows=1,
    )
    load_config(data, create=True)
    return data


def _snapshot(data):
    return {
        p.name: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
        for p in data.iterdir()
        if p.is_file()
    }


def _wire(monkeypatch, kind, *, changed=False, cpic_changed=False):
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db, "get_annotators", lambda data, **kw: [PharmGKBAnnotator(data)])
    headers = {"Last-Modified": _NEW_LM if changed else _OLD_LM}
    if kind == "etag":
        headers["ETag"] = '"new"' if changed else '"old"'
    monkeypatch.setattr(pharmgkb, "head_request_headers", lambda *a, **kw: headers)
    monkeypatch.setattr(
        pharmgkb, "fetch_cpic_remote_signal", lambda: _NEW_DATE if cpic_changed else _OLD_DATE
    )
    monkeypatch.setattr(pharmgkb, "PHARMGKB_MIN_ROWS", 1)


def _assert_signal_parts(actual, expected):
    """Keep the source and CPIC evidence exact, allowing another HTTP field."""
    assert actual is not None
    parts = actual.split("|")
    source, cpic = expected.split("|")
    prefix = ":".join(source.split(":")[:2]) + ":"
    assert [part for part in parts if part.startswith(prefix)] == [source]
    assert [part for part in parts if part.startswith("cpic:")] == [cpic]


def _update(data):
    result = CliRunner(env={"COLUMNS": "500"}).invoke(
        main, ["db", "update", "--data-dir", str(data)]
    )
    assert result.exit_code == 0, (result.output, result.exception)
    with PharmGKBAnnotator(data) as ann:
        assert ann.is_ready()
        return result.output, ann.cached_remote_signal()


@pytest.mark.parametrize("kind", ["etag", "last-modified"])
def test_unchanged_signal_does_not_download_copy_or_rewrite(tmp_path, monkeypatch, kind):
    data = _built(tmp_path, kind)
    _wire(monkeypatch, kind)
    before = _snapshot(data)

    def forbidden(*a, **kw):
        pytest.fail("Unchanged source attempted a download, data fetch or cache copy")

    monkeypatch.setattr(pharmgkb, "download", forbidden)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", forbidden)
    for name in ("copy", "copy2", "copyfile", "copytree"):
        monkeypatch.setattr(shutil, name, forbidden)
    output, stored = _update(data)
    assert "already current" in output, output
    assert stored == _signal(kind)
    assert _snapshot(data) == before


@pytest.mark.parametrize("kind", ["etag", "last-modified"])
@pytest.mark.parametrize("changed_part", ["clinical-annotations", "cpic"])
def test_changed_content_is_ingested_and_stamped(tmp_path, monkeypatch, kind, changed_part):
    data = _built(tmp_path, kind)
    pg_changed = changed_part == "clinical-annotations"
    _wire(monkeypatch, kind, changed=pg_changed, cpic_changed=not pg_changed)
    downloads = []

    def download(url, dest):
        downloads.append(url)
        with zipfile.ZipFile(dest, "w") as archive:
            for path in sorted(_FIXTURES.glob("*.tsv")):
                text = path.read_text()
                if pg_changed and path.name == "clinical_ann_alleles.tsv":
                    text = text.replace(
                        "Heterozygous carriers may have elevated toxicity risk",
                        "Changed synthetic clinical response",
                    )
                archive.writestr(path.name, text)

    new_lookup = dict(_LOOKUP)
    if not pg_changed:
        new_lookup[("rs1801133", "A")] = "normal"
    monkeypatch.setattr(pharmgkb, "download", download)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", lambda: new_lookup)
    _, stored = _update(data)
    assert len(downloads) == 1
    _assert_signal_parts(stored, _signal(kind, changed=pg_changed, cpic_changed=not pg_changed))
    with contextlib.closing(
        sqlite3.connect(f"file:{data / 'pharmgkb.sqlite'}?mode=ro", uri=True)
    ) as conn:
        if pg_changed:
            assert (
                conn.execute(
                    "SELECT COUNT(*) FROM pharmgkb_annotations WHERE "
                    "annotation_text='Changed synthetic clinical response'"
                ).fetchone()[0]
                > 0
            )
        else:
            assert (
                conn.execute(
                    "SELECT function_class FROM pharmgkb_allele_function "
                    "WHERE rsid='rs1801133' AND allele='A'"
                ).fetchone()[0]
                == "normal"
            )

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""HTTP freshness is case-insensitive and compares shared evidence fields."""

from __future__ import annotations

import contextlib
import importlib
import sqlite3
import zipfile
from email.message import Message
from pathlib import Path

import pytest
from click.testing import CliRunner

import allelix
from allelix.annotators import gwas, pharmgkb
from allelix.cli import main
from allelix.config import load_config
from allelix.databases import manager
from allelix.databases.gwas_loader import load_gwas_tsv
from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

FIXTURES = Path(allelix.__file__).resolve().parent.parent / "tests/fixtures"
LM = "Wed, 01 Jan 2020 00:00:00 GMT"
NEW_LM = "Sat, 01 Feb 2020 00:00:00 GMT"
ETAG = '"epoch-A"'
NEW_ETAG = '"epoch-a"'
CPIC = "lastchange:2020-01-01"
LOOKUP = {("rs1801133", "A"): "decreased", ("rs1801133", "G"): "normal"}
LEGACY_COLUMNS = (
    "name, source_url, version, downloaded_at, record_count, remote_signal, local_version_tag"
)


class Response:
    def __init__(self, headers):
        self.headers = Message()
        for key, value in headers.items():
            self.headers[key] = value

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class Remote:
    def __init__(self, source):
        self.source = source
        self.module = gwas if source == "gwas" else pharmgkb
        self.cls = gwas.GWASCatalogAnnotator if source == "gwas" else pharmgkb.PharmGKBAnnotator
        self.headers = {"ETag": ETAG, "Last-Modified": LM}
        self.downloads = self.heads = self.cpic_fetches = 0
        self.revision = 0
        self.cpic = CPIC

    def open(self, request, *args, **kwargs):
        assert request.get_method() == "HEAD", "Freshness must not fetch a data payload"
        self.heads += 1
        return Response(self.headers)

    def download(self, url, dest):
        self.downloads += 1
        with zipfile.ZipFile(dest, "w") as archive:
            if self.source == "gwas":
                text = (FIXTURES / "mock_gwas_catalog.tsv").read_text()
                if self.revision:
                    # Change a loaded trait, not just ZIP metadata.
                    lines = text.splitlines()
                    header = lines[0].split("\t")
                    index = header.index("DISEASE/TRAIT")
                    row = lines[1].split("\t")
                    row[index] += " changed"
                    lines[1] = "\t".join(row)
                    text = "\n".join(lines) + "\n"
                archive.writestr("associations.tsv", text)
            else:
                for file in sorted((FIXTURES / "mock_pharmgkb").glob("*.tsv")):
                    text = file.read_text()
                    if self.revision and file.name == "clinical_ann_alleles.tsv":
                        text = text.replace(
                            "Heterozygous carriers may have elevated toxicity risk",
                            "Changed synthetic response",
                        )
                    archive.writestr(file.name, text)

    def cpic_data(self):
        self.cpic_fetches += 1
        lookup = dict(LOOKUP)
        if self.cpic != CPIC:
            lookup[("rs1801133", "A")] = "normal"
        return lookup

    def reset(self):
        self.downloads = self.heads = self.cpic_fetches = 0


@pytest.fixture(params=["gwas", "pharmgkb"])
def remote(request, monkeypatch):
    remote = Remote(request.param)
    monkeypatch.setattr(manager.urllib.request, "urlopen", remote.open)
    monkeypatch.setattr(remote.module, "download", remote.download)
    monkeypatch.setattr(gwas, "GWAS_MIN_ROWS", 1)
    monkeypatch.setattr(pharmgkb, "PHARMGKB_MIN_ROWS", 1)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_remote_signal", lambda: remote.cpic)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", remote.cpic_data)
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db, "get_annotators", lambda data, **kwargs: [remote.cls(data)])
    return remote


def _update(data, remote):
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
    assert result.exit_code == 0, (result.output, result.exception)
    with remote.cls(data) as annotator:
        assert annotator.is_ready(), result.output
    return " ".join(result.output.lower().split())


def _stamp(data, source):
    with contextlib.closing(
        sqlite3.connect(f"file:{data / (source + '.sqlite')}?mode=ro", uri=True)
    ) as c:
        c.row_factory = sqlite3.Row
        return dict(
            c.execute("SELECT * FROM database_versions WHERE name=?", (source,)).fetchone()
        )


def _payload_state(data, source):
    path = data / (source + ".sqlite")
    with contextlib.closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as conn:
        names = [
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
            if row[0] != "database_versions" and not row[0].startswith("sqlite_")
        ]
        rows = {
            name: conn.execute(f'SELECT * FROM "{name}" ORDER BY rowid').fetchall()
            for name in names
        }
    stamp = _stamp(data, source)
    # Optional HTTP metadata may be learned in place. The cache itself and
    # its interpretation/version must remain; no rebuild or replacement.
    return (
        path.stat().st_ino,
        rows,
        tuple(stamp[k] for k in ("version", "local_version_tag", "record_count")),
    )


def _install(tmp_path, remote):
    data = tmp_path / "data"
    _update(data, remote)
    assert remote.downloads == 1
    remote.reset()
    return data


def _legacy(tmp_path, remote, kind):
    data = tmp_path / "data"
    data.mkdir()
    signal = f"etag:{ETAG}" if kind == "etag" else f"lm:{LM}"
    if remote.source == "gwas":
        load_gwas_tsv(
            FIXTURES / "mock_gwas_catalog.tsv",
            data / "gwas.sqlite",
            remote_signal=signal,
            min_rows=1,
        )
    else:
        load_pharmgkb_tsv(
            FIXTURES / "mock_pharmgkb",
            data / "pharmgkb.sqlite",
            version="synthetic",
            remote_signal=f"pgkb:{signal}|cpic:{CPIC}",
            allele_function_lookup=LOOKUP,
            min_rows=1,
        )
    # Recreate the published seven-column stamp, even after new fields are added.
    with contextlib.closing(sqlite3.connect(data / (remote.source + ".sqlite"))) as conn:
        conn.execute("ALTER TABLE database_versions RENAME TO saved_versions")
        conn.execute(
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT, "
            "version TEXT, downloaded_at TEXT, record_count INTEGER, "
            "remote_signal TEXT, local_version_tag TEXT)"
        )
        conn.execute(
            f"INSERT INTO database_versions ({LEGACY_COLUMNS}) "
            f"SELECT {LEGACY_COLUMNS} FROM saved_versions"
        )
        conn.execute("DROP TABLE saved_versions")
        conn.commit()
    load_config(data, create=True)
    remote.reset()
    return data


@pytest.mark.parametrize("etag_name", ["ETag", "Etag", "etag", "eTaG"])
@pytest.mark.parametrize("lm_name", ["Last-Modified", "last-modified"])
def test_shared_head_helper_makes_source_signal_independent_of_header_case(
    tmp_path, remote, etag_name, lm_name
):
    with remote.cls(tmp_path) as annotator:
        canonical = annotator.fetch_remote_signal()
    remote.headers = {etag_name: ETAG, lm_name: LM}
    with remote.cls(tmp_path) as annotator:
        actual = annotator.fetch_remote_signal()
    assert actual == canonical and actual is not None
    assert "epoch-A" in actual, "ETag takes precedence when both fields exist"
    assert remote.heads == 2


@pytest.mark.parametrize(
    "lm_name", ["Last-Modified", "Last-modified", "last-modified", "LAST-MODIFIED"]
)
def test_last_modified_alone_is_case_insensitive(tmp_path, remote, lm_name):
    remote.headers = {"Last-Modified": LM}
    with remote.cls(tmp_path) as annotator:
        canonical = annotator.fetch_remote_signal()
    remote.headers = {lm_name: LM}
    with remote.cls(tmp_path) as annotator:
        actual = annotator.fetch_remote_signal()
    assert actual == canonical and actual is not None


@pytest.mark.parametrize("kind", ["etag", "lm"])
@pytest.mark.parametrize("lowercase", [False, True])
def test_old_stamp_with_shared_unchanged_field_needs_no_download(
    tmp_path, remote, kind, lowercase
):
    data = _legacy(tmp_path, remote, kind)
    if lowercase:
        remote.headers = {"etag": ETAG, "last-modified": LM}
    before = _payload_state(data, remote.source)
    output = _update(data, remote)
    assert remote.downloads == remote.cpic_fetches == 0, output
    assert "already current" in output, output
    assert _payload_state(data, remote.source) == before, (
        "An old valid cache must not be rebuilt or replaced"
    )


@pytest.mark.parametrize("lowercase", [False, True])
def test_build_persists_both_header_values_for_the_next_process(tmp_path, remote, lowercase):
    if lowercase:
        remote.headers = {"etag": ETAG, "last-modified": LM}
    data = _install(tmp_path, remote)
    stored = str(_stamp(data, remote.source))
    assert "epoch-A" in stored and LM in stored, stored
    # The absent ETag is not a change. _update creates a new annotator instance.
    remote.headers = {"last-modified": LM}
    output = _update(data, remote)
    assert remote.downloads == 0, output
    assert "already current" in output, output


@pytest.mark.parametrize(
    "stored_kind,remote_kind,etag,lm,changed",
    [
        ("both", "lm", ETAG, LM, False),
        ("both", "lm", ETAG, NEW_LM, True),
        ("lm", "both", ETAG, LM, False),
        ("lm", "both", ETAG, NEW_LM, True),
        ("both", "both", ETAG, NEW_LM, False),
        ("both", "both", NEW_ETAG, LM, True),
        ("etag", "both", ETAG, NEW_LM, False),
        ("both", "etag", ETAG, LM, False),
    ],
)
def test_compare_shared_fields_etag_first(
    tmp_path, remote, stored_kind, remote_kind, etag, lm, changed
):
    if stored_kind == "lm":
        remote.headers = {"Last-Modified": LM}
    elif stored_kind == "etag":
        remote.headers = {"ETag": ETAG}
    data = _install(tmp_path, remote)
    remote.headers = {}
    if remote_kind in ("etag", "both"):
        remote.headers["ETag"] = etag
    if remote_kind in ("lm", "both"):
        remote.headers["Last-Modified"] = lm
    remote.revision = int(changed)
    before = _payload_state(data, remote.source)
    output = _update(data, remote)
    assert remote.downloads == int(changed), output
    if not changed:
        assert "already current" in output
        assert _payload_state(data, remote.source) == before
    else:
        assert "refreshed" in output
        assert _payload_state(data, remote.source) != before
        remote.reset()
        output = _update(data, remote)
        assert remote.downloads == 0 and "already current" in output


@pytest.mark.parametrize("old_kind", ["etag", "lm"])
def test_disjoint_header_kinds_alone_do_not_claim_a_change(tmp_path, remote, old_kind):
    data = _legacy(tmp_path, remote, old_kind)
    remote.headers = {"Last-Modified": LM} if old_kind == "etag" else {"ETag": ETAG}
    before = _payload_state(data, remote.source)
    output = _update(data, remote)
    assert remote.downloads == remote.cpic_fetches == 0, output
    assert "remote signal changed" not in output
    assert _payload_state(data, remote.source) == before


def test_clinpgx_cpic_change_is_not_hidden_by_equal_http_fields(tmp_path, remote):
    if remote.source != "pharmgkb":
        # A GWAS control: no CPIC dependency or fetch in its refresh path.
        data = _install(tmp_path, remote)
        remote.cpic = "lastchange:2020-02-01"
        output = _update(data, remote)
        assert remote.downloads == remote.cpic_fetches == 0, output
        return
    data = _install(tmp_path, remote)
    remote.headers = {"Last-Modified": LM}
    remote.cpic = "lastchange:2020-02-01"
    output = _update(data, remote)
    assert remote.downloads == remote.cpic_fetches == 1, output
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        value = conn.execute(
            "SELECT function_class FROM pharmgkb_allele_function "
            "WHERE rsid='rs1801133' AND allele='A'"
        ).fetchone()[0]
    assert value == "normal"

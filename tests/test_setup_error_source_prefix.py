# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A source failure has one CLI source prefix and keeps the useful reason."""

from __future__ import annotations

import contextlib
import importlib
import io
import re
import sqlite3
import zipfile

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.annotators import get_annotators
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.cli import _helpers, main


@pytest.fixture
def cli_source(tmp_path, monkeypatch):
    db = importlib.import_module("allelix.cli.db")
    data = tmp_path / "data"
    data.mkdir()

    def select(name, *, refresh=False):
        registry = get_annotators(data, clinvar_builds=("GRCh37",))
        ann = next(a for a in registry if a.name == name)
        for other in registry:
            if other is not ann:
                other.close()
        monkeypatch.setattr(db, "get_annotators", lambda *a, **kw: [ann])
        if refresh:
            monkeypatch.setattr(ann, "is_ready", lambda: True)
            monkeypatch.setattr(ann, "fetch_remote_signal", lambda: "changed")
            monkeypatch.setattr(ann, "cached_remote_signal", lambda: "old")
            # Commit-pinned sources use the forced update path.
        monkeypatch.setattr(db, "_confirm_cadd_license", lambda **kw: True)
        return ann

    return data, select


def _invoke(data, name, *, force=False):
    args = ["db", "update", "--data-dir", str(data)]
    if force:
        args += ["--force"]
    if name == "cadd":
        args += ["--cadd"]
    result = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
    assert result.exit_code == 1, (result.output, result.exception)
    assert "Traceback" not in result.output
    return " ".join(result.output.split())


def _once(text, name, reason):
    # The progress and final summary may name the source again. Count only
    # the error sentence, and only source prefixes, not names inside paths.
    match = re.search(rf"(?:{name}:\s*)+[^\n]*{re.escape(reason)}", text, re.I)
    assert match, text
    segment = match.group(0)
    assert not re.search(rf"{name}:\s*{name}:", segment, re.I), text


@pytest.mark.parametrize("refresh", [False, True])
@pytest.mark.parametrize(
    "case", ["submission-checksum", "variant-checksum", "old-count", "publish"]
)
def test_clinvar_actual_setup_errors_have_one_prefix(cli_source, monkeypatch, refresh, case):
    data, select = cli_source
    ann = select("clinvar", refresh=refresh)
    mod = importlib.import_module("allelix.annotators.clinvar")
    monkeypatch.setattr(ann, "_fetch_remote_signal_for_tsv", lambda: "md5:" + "0" * 32)
    monkeypatch.setattr(mod, "_fetch_checksum", lambda url: "1" * 32)
    if case == "submission-checksum":
        monkeypatch.setattr(ann, "_fetch_remote_signal_for_tsv", lambda: None)
        reason = "checksum for submission_summary"
    elif case == "variant-checksum":
        monkeypatch.setattr(mod, "_fetch_checksum", lambda url: None)
        reason = "variant_summary checksum"
    elif case == "old-count":
        with contextlib.closing(sqlite3.connect(data / "clinvar.GRCh37.sqlite")) as conn, conn:
            conn.execute("CREATE TABLE mock (n INTEGER)")
        reason = "existing cache record count"
    else:
        monkeypatch.setattr(mod, "download", lambda url, dest: dest.write_bytes(b"synthetic"))
        monkeypatch.setattr(mod, "verify_file_hash", lambda *a, **kw: None)
        monkeypatch.setattr(ann, "_setup_one_from_tsv", lambda vs, ss, dest, *a: dest.touch())
        original = mod.os.replace

        def replace(src, dest):
            if dest == data / "clinvar.GRCh37.sqlite":
                raise PermissionError("synthetic publication denied")
            return original(src, dest)

        monkeypatch.setattr(mod.os, "replace", replace)
        reason = "publication failed"
    text = _invoke(data, "clinvar", force=refresh)
    _once(text, "clinvar", reason)
    assert "kept" in text or "no cache installed" in text, text


@pytest.mark.parametrize("name", ["gwas", "pharmgkb", "rsid_identity"])
@pytest.mark.parametrize("refresh", [False, True])
def test_actual_cannot_verify_errors_have_one_prefix(cli_source, monkeypatch, name, refresh):
    data, select = cli_source
    ann = select(name, refresh=refresh)
    mod = importlib.import_module(f"allelix.annotators.{name}")
    if name == "rsid_identity":
        monkeypatch.setattr(mod, "RSID_IDENTITY_EXPECTED_SHA256", "")
        reason = "no published cache is pinned"
    else:
        # A refresh first sees a changed remote signal; the subsequent setup
        # probe fails. A fresh install goes straight to the failing probe.
        replies = iter(["new", None] if refresh else [None])
        monkeypatch.setattr(ann, "fetch_remote_signal", lambda: next(replies))
        reason = "cannot verify remote freshness signal"
    text = _invoke(data, name, force=refresh and name == "rsid_identity")
    _once(text, name, reason)


@pytest.mark.parametrize(
    "name",
    ["clinvar", "gwas", "pharmgkb", "rsid_identity", "snpedia", "gnomad", "alphamissense", "cadd"],
)
@pytest.mark.parametrize("refresh", [False, True])
def test_propagated_download_failure_has_one_source_prefix(cli_source, monkeypatch, name, refresh):
    data, select = cli_source
    ann = select(name, refresh=refresh)
    mod = importlib.import_module(f"allelix.annotators.{name}")
    if name == "clinvar":
        monkeypatch.setattr(ann, "_fetch_remote_signal_for_tsv", lambda: "md5:" + "0" * 32)
        monkeypatch.setattr(mod, "_fetch_checksum", lambda url: "1" * 32)
    elif name in ("gwas", "pharmgkb"):
        monkeypatch.setattr(ann, "fetch_remote_signal", lambda: "new")

    def fail(*a, **kw):
        raise OSError("synthetic connection interrupted")

    monkeypatch.setattr(mod, "download", fail)
    text = _invoke(data, name, force=refresh)
    _once(text, name, "synthetic connection interrupted")


def test_gwas_malformed_archive_keeps_its_reason(cli_source, monkeypatch):
    data, select = cli_source
    ann = select("gwas")
    mod = importlib.import_module("allelix.annotators.gwas")
    monkeypatch.setattr(ann, "fetch_remote_signal", lambda: "new")

    def empty_archive(url, dest):
        with zipfile.ZipFile(dest, "w") as zf:
            zf.writestr("readme.txt", "Synthetic archive without a table")

    monkeypatch.setattr(mod, "download", empty_archive)
    _once(_invoke(data, "gwas"), "gwas", "No .tsv file found")


def test_automatic_setup_uses_the_same_single_prefix(tmp_path, monkeypatch):
    """The shared helper also serves analysis' automatic setup, not just db update."""
    ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
    monkeypatch.setattr(ann, "_fetch_remote_signal_for_tsv", lambda: None)
    stream = io.StringIO()
    monkeypatch.setattr(_helpers, "console", Console(file=stream, width=300, color_system=None))
    assert not _helpers._run_setup(ann)
    _once(stream.getvalue(), "clinvar", "checksum for submission_summary")

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An offline report discloses an unavailable CPIC normal-function filter."""

from __future__ import annotations

import contextlib
import hashlib
import html
import json
import re
import socket
import sqlite3
import urllib.request

import pytest
from click.testing import CliRunner

from allelix.annotators import pharmgkb
from allelix.cli import main
from tests.test_annotation_parity_and_named_disclosures import _clinvar, _pharm, _vcf


@pytest.fixture
def cache(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    _pharm(data, "rs900260201", "MOCK", ("AG",))
    for build in ("GRCh37", "GRCh38"):
        _clinvar(
            data, build, [("rs900260201", "1", 12345, "G", "A", "MOCK", "Benign", "Benign", 1, 1)]
        )
    return data


def _state(data, state):
    signals = {
        "unavailable": "pgkb:etag:mock|cpic:unavailable",
        "marker-only": "cpic:unavailable",
        "legacy-empty": "pgkb:etag:mock|cpic:lastchange:2020-01-01",
        "unstamped-empty": None,
        "healthy": "pgkb:etag:mock|cpic:lastchange:2020-01-01",
    }
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.execute("UPDATE database_versions SET remote_signal=?", (signals[state],))
        conn.execute("DELETE FROM pharmgkb_allele_function")
        if state == "healthy":
            conn.execute(
                "INSERT INTO pharmgkb_allele_function (rsid,allele,function_class) "
                "VALUES ('rs900260201','A','decreased')"
            )
        conn.commit()


def _snapshot(data):
    return {
        p.name: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
        for p in data.iterdir()
        if p.is_file()
    }


def _offline(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Report generation attempted network access or a CPIC refresh")

    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(socket, "getaddrinfo", forbidden)
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    monkeypatch.setattr(pharmgkb, "head_request_headers", forbidden)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_remote_signal", forbidden)
    monkeypatch.setattr(pharmgkb, "fetch_cpic_allele_functions", forbidden)
    monkeypatch.setattr(pharmgkb, "download", forbidden)


def _report(tmp_path, cache, monkeypatch, *, command, surface, build, extra=()):
    _offline(monkeypatch)
    vcf = _vcf(tmp_path / "mock.vcf", [("rs900260201", "1", 12345, "G", "A", "0/1", ".")])
    output = tmp_path / f"report.{surface}"
    args = [
        command,
        str(vcf),
        "--build",
        build.lower(),
        "--data-dir",
        str(cache),
        "--no-update",
        "--min-magnitude",
        "0",
        *extra,
    ]
    if surface != "terminal":
        args += ["--output", str(output)]
    before = _snapshot(cache)
    files = list(cache.iterdir())
    for path in files:
        path.chmod(0o444)
    cache.chmod(0o555)
    try:
        result = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
        assert result.exit_code == 0, (result.output, result.exception)
        assert _snapshot(cache) == before, "Analysis changed its source cache"
    finally:
        cache.chmod(0o755)
        for path in files:
            path.chmod(0o644)
    if surface == "json":
        report = json.loads(output.read_text())
        notices = [w for w in report["report_warnings"] if w["code"] == "cpic_unavailable"]
        return [w["message"] for w in notices], result.output
    if surface == "html":
        raw = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", output.read_text(), flags=re.S | re.I)
        return [html.unescape(re.sub(r"<[^>]+>", " ", raw))], result.output
    return [result.output], result.output


@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics", "methylation"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "state", ["unavailable", "marker-only", "legacy-empty", "unstamped-empty"]
)
def test_missing_cpic_is_disclosed_without_repair(
    tmp_path, cache, monkeypatch, command, surface, build, state
):
    _state(cache, state)
    notices, console = _report(
        tmp_path, cache, monkeypatch, command=command, surface=surface, build=build
    )
    text = " ".join(" ".join(notices).split())
    assert "CPIC" in text, text
    assert "normal" in text.lower() and "function" in text.lower(), text
    assert "allelix db update" in text, text
    assert any(word in text.lower() for word in ("shown", "show", "appear")), text
    assert "CPIC" in console and "allelix db update" in " ".join(console.split()), console
    if surface == "json":
        assert len(notices) == 1, notices


@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics", "methylation"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
@pytest.mark.parametrize("mode", ["healthy", "disabled-config", "absent"])
def test_unused_or_healthy_cpic_has_no_degraded_filter_notice(
    tmp_path, cache, monkeypatch, command, surface, mode
):
    _state(cache, "healthy" if mode == "healthy" else "unavailable")
    extra = ()
    if mode == "disabled-config":
        (cache / "config.toml").write_text("[sources]\npharmgkb = false\n")
    elif mode == "absent":
        (cache / "pharmgkb.sqlite").unlink()
    notices, console = _report(
        tmp_path, cache, monkeypatch, command=command, surface=surface, build="GRCh38", extra=extra
    )
    if surface == "json":
        assert notices == [], notices
    for text in (*notices, console):
        assert not re.search(r"CPIC.{0,80}(?:unavailable|missing|not available)", text, re.I), text


@pytest.mark.parametrize(
    "state", ["unavailable", "marker-only", "legacy-empty", "unstamped-empty"]
)
def test_notice_does_not_change_cached_findings(tmp_path, cache, monkeypatch, state):
    """Adding the notice does not reclassify rows already stored in this cache."""
    _state(cache, "healthy")
    _report(tmp_path, cache, monkeypatch, command="analyze", surface="json", build="GRCh38")
    before = json.loads((tmp_path / "report.json").read_text())["annotations"]
    assert any(a["source"] == "pharmgkb" for a in before)
    _state(cache, state)
    _report(tmp_path, cache, monkeypatch, command="analyze", surface="json", build="GRCh38")
    after = json.loads((tmp_path / "report.json").read_text())["annotations"]
    assert after == before

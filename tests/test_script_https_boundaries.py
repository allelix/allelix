"""Offline source builders use the same HTTPS redirect boundary as updates."""

import ast
import gzip
import importlib.util
import io
import shlex
import sys
import urllib.error
import urllib.request
from types import SimpleNamespace

import pytest
from test_publish_release_gates import ROOT

from allelix.databases import manager

SCRIPTS = ("build_gnomad_cache", "build_alphamissense_cache", "scrape_snpedia")


def _load(name, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "probe_" + name, ROOT / "scripts" / (name + ".py")
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("name", SCRIPTS)
@pytest.mark.parametrize("refused", [False, True])
def test_python_source_downloads_use_the_shared_opener(monkeypatch, name, refused):
    seen = []
    bare = []

    def open_shared(request, timeout):
        seen.append((request.full_url, timeout))
        if refused:
            raise urllib.error.URLError("Refusing non-HTTPS redirect")
        if name == "scrape_snpedia":
            return io.BytesIO(b'{"query": {"synthetic": true}}')
        payload = b"# synthetic header\n"
        if name == "build_gnomad_cache":
            payload += b"1\t1000\trs123\tG\tA\t.\tPASS\tAF=0.1\n"
        return io.BytesIO(gzip.compress(payload))

    def open_bare(*args, **kwargs):
        bare.append(True)
        # No connection occurs, even on the unfixed tree.
        raise urllib.error.URLError("unprotected opener intercepted by test")

    monkeypatch.setattr(manager, "open_url", open_shared)
    monkeypatch.setattr(urllib.request, "urlopen", open_bare)
    module = _load(name, monkeypatch)
    if name == "scrape_snpedia":
        monkeypatch.setattr(module, "MAX_RETRIES", 0)

    def run():
        if name == "build_gnomad_cache":
            return list(module._stream_vcf_gz_filtered("https://source.invalid/test.gz"))
        if name == "build_alphamissense_cache":
            with module._stream_zenodo() as stream:
                return stream.read()
        return module.api_get({"action": "query", "format": "json"})

    error = None
    result = None
    try:
        result = run()
    except urllib.error.URLError as exc:
        error = exc
    assert not bare, f"{name} bypassed manager.open_url"
    assert len(seen) == 1 and seen[0][0].startswith("https://") and seen[0][1] > 0
    if refused:
        assert error is not None and "non-HTTPS redirect" in str(error)
    else:
        assert error is None and result


def _curl_bound(command):
    tokens = shlex.split(command)
    assert "curl" in tokens, command
    i = tokens.index("curl")
    tokens = tokens[i:]
    assert "--proto-redir" in tokens, "curl redirect targets must be restricted before contact"
    pos = tokens.index("--proto-redir")
    assert tokens[pos + 1] == "=https", "Do not leave default redirect protocols enabled"


@pytest.mark.parametrize("line_index", [0, 1])
def test_fixture_fetcher_binds_each_curl_redirect(line_index):
    lines = [
        line.strip()
        for line in (ROOT / "scripts/fetch_testdata.sh").read_text().splitlines()
        if line.lstrip().startswith("curl ")
    ]
    assert len(lines) == 2
    _curl_bound(lines[line_index])


@pytest.mark.parametrize("remote", [False, True])
def test_identity_builder_bounds_remote_redirects_but_keeps_local_input(
    tmp_path, monkeypatch, remote
):
    module = _load("build_rsid_identity_cache", monkeypatch)
    commands = []

    class Process:
        returncode = 0

        def __init__(self, args, **kwargs):
            commands.append(args[-1])
            self.stdout = io.StringIO("")

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.stdout.close()

    monkeypatch.setattr(module.subprocess, "Popen", Process)
    monkeypatch.setattr(module.shutil, "which", lambda name: "/bin/bash")
    # The empty input never queries the synthetic Bloom object.
    module.filter_source(
        "https://source.invalid/db.gz" if remote else "/synthetic/db.gz",
        SimpleNamespace(),
        tmp_path / "out.tsv",
    )
    assert len(commands) == 1
    if remote:
        _curl_bound(commands[0])
    else:
        assert shlex.split(commands[0]) == ["gzip", "-dc", "/synthetic/db.gz"]


def test_no_remaining_bare_urlopen_in_python_scripts():
    found = []
    for path in sorted((ROOT / "scripts").glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and (
                (isinstance(node.func, ast.Attribute) and node.func.attr == "urlopen")
                or (isinstance(node.func, ast.Name) and node.func.id == "urlopen")
            ):
                found.append(f"{path.name}:{node.lineno}")
    assert not found, f"Unprotected source download sites: {found}"

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The builder launches the PATH-selected shell and retains bounded stream/error behavior."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


@pytest.fixture
def builder(monkeypatch: pytest.MonkeyPatch):
    # The script prepends its repository to sys.path; keep that test-local.
    monkeypatch.setattr(sys, "path", list(sys.path))
    path = Path(__file__).resolve().parents[1] / "scripts" / "build_rsid_identity_cache.py"
    spec = importlib.util.spec_from_file_location("identity_builder_stream_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _LineStream:
    def __iter__(self):
        return iter(["synthetic row one\n", "synthetic row two\n"])

    def read(self, *_args):
        pytest.fail("Builder must iterate the input, not read it all")


class _Process:
    def __init__(self, exit_code=0):
        self.stdout = _LineStream()
        self.returncode = None
        self.exit_code = exit_code
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.returncode = self.exit_code
        self.closed = True


@pytest.mark.parametrize(
    "source,command",
    [
        (r"C:\synthetic files\dbsnp.vcf.gz", r"gzip -dc 'C:\synthetic files\dbsnp.vcf.gz'"),
        (
            "https://example.invalid/dbsnp.vcf.gz?x=1&y=2",
            "curl -sSfL --retry 5 'https://example.invalid/dbsnp.vcf.gz?x=1&y=2' | gzip -dc",
        ),
    ],
)
def test_local_and_remote_streams_use_the_resolved_shell(
    builder,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source,
    command,
) -> None:
    selected = str(tmp_path / "Git Bash" / "bash.exe")
    process = _Process()
    calls = []
    bloom = object()
    monkeypatch.setattr(
        builder.shutil, "which", lambda executable: selected if executable == "bash" else None
    )

    def launch(argv, **kwargs):
        calls.append((argv, kwargs))
        return process

    def filter_rows(stream, actual_bloom, out):
        assert actual_bloom is bloom
        count = 0
        for row in stream:
            out.write(row)
            count += 1
        return count, count

    monkeypatch.setattr(builder.subprocess, "Popen", launch)
    monkeypatch.setattr(builder, "filter_dbsnp_stream", filter_rows)
    out = tmp_path / "placements.tsv"
    builder.filter_source(source, bloom, out)
    assert calls == [
        (
            [selected, "-o", "pipefail", "-c", command],
            {"stdout": builder.subprocess.PIPE, "text": True, "bufsize": 1 << 20},
        )
    ]
    assert out.read_text(encoding="utf-8") == "synthetic row one\nsynthetic row two\n"
    assert process.closed


def test_missing_shell_is_actionable_and_does_not_spawn(
    builder,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(builder.shutil, "which", lambda _exe: None)
    monkeypatch.setattr(
        builder.subprocess, "Popen", lambda *_a, **_kw: pytest.fail("missing shell spawned")
    )
    out = tmp_path / "placements.tsv"
    with pytest.raises(SystemExit, match=r"[Bb]ash.*PATH"):
        builder.filter_source("input.gz", object(), out)
    assert not out.exists()


@pytest.mark.parametrize("source", ["input.gz", "https://example.invalid/input.gz"])
def test_child_stream_failure_is_not_reported_as_success(
    builder,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source,
) -> None:
    process = _Process(exit_code=7)
    monkeypatch.setattr(builder.shutil, "which", lambda _exe: "/synthetic/bash")
    monkeypatch.setattr(builder.subprocess, "Popen", lambda *_a, **_kw: process)
    monkeypatch.setattr(
        builder, "filter_dbsnp_stream", lambda stream, _b, out: (len(list(stream)), 0)
    )
    with pytest.raises(SystemExit, match="failed with exit code 7"):
        builder.filter_source(source, object(), tmp_path / "placements.tsv")
    assert process.closed


def test_shell_launch_error_still_propagates(
    builder,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(builder.shutil, "which", lambda _exe: "/synthetic/bash")

    def refused(*_args, **_kwargs):
        raise PermissionError("synthetic shell launch denied")

    monkeypatch.setattr(builder.subprocess, "Popen", refused)
    with pytest.raises(PermissionError, match="synthetic shell launch denied"):
        builder.filter_source("input.gz", object(), tmp_path / "placements.tsv")

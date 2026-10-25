# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""PLINK temporary ownership and per-file publication guarantees."""

from __future__ import annotations

import os

import pytest

from allelix.exporters import plink
from allelix.models import Variant


def _variants(fail=False):
    yield Variant("rs1", "1", 100, "A", "G", build="GRCh37")
    if fail:
        raise RuntimeError("synthetic iterator failure")


@pytest.mark.parametrize("failure", [False, True])
def test_plink_keeps_unrelated_legacy_temps(tmp_path, failure):
    prefix = tmp_path / "release.v1"
    paths = [prefix.with_suffix(ext) for ext in (".fam", ".bim", ".bed")]
    for path in paths:
        path.write_bytes(b"prior")
        path.with_name(path.name + ".tmp").write_bytes(b"unrelated legacy")
    if failure:
        with pytest.raises(RuntimeError, match="iterator failure"):
            plink.export_plink(_variants(True), prefix, "GRCh37")
    else:
        plink.export_plink(_variants(), prefix, "GRCh37")
    for path in paths:
        assert path.with_name(path.name + ".tmp").read_bytes() == b"unrelated legacy"
        if failure:
            assert path.read_bytes() == b"prior"
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("failed_rename", [1, 2, 3])
def test_plink_publication_is_atomic_per_file(tmp_path, monkeypatch, failed_rename):
    prefix = tmp_path / "out"
    paths = [prefix.with_suffix(ext) for ext in (".fam", ".bim", ".bed")]
    for path in paths:
        path.write_bytes(b"prior")
    calls = []
    replace = os.replace

    def publish(source, destination):
        calls.append(destination)
        if len(calls) == failed_rename:
            raise OSError("synthetic rename failure")
        replace(source, destination)

    monkeypatch.setattr(os, "replace", publish)
    with pytest.raises(OSError, match="rename failure"):
        plink.export_plink(_variants(), prefix, "GRCh37")
    for index, path in enumerate(paths):
        assert (path.read_bytes() != b"prior") == (index < failed_rename - 1)
    assert not list(tmp_path.glob(".*.tmp"))


def test_plink_all_content_is_closed_before_first_publication(tmp_path, monkeypatch):
    import contextlib

    from allelix import reports

    stages = []
    original = reports.staged_output
    replace = os.replace
    published = []

    @contextlib.contextmanager
    def staged(path):
        with original(path) as output:
            stages.append(output)
            yield output

    def publish(source, destination):
        assert len(stages) == 3
        assert all(stage.stream is not None and stage.stream.closed for stage in stages)
        if not published:
            assert all(stage.temporary.exists() for stage in stages)
        published.append(destination.suffix)
        replace(source, destination)

    monkeypatch.setattr(reports, "staged_output", staged)
    monkeypatch.setattr(os, "replace", publish)
    prefix = tmp_path / "out"
    counts = plink.export_plink(_variants(), prefix, "GRCh37")
    assert counts.written == 1
    assert published == [".fam", ".bim", ".bed"]
    assert (
        prefix.with_suffix(".fam").read_bytes() == ("0\tSAMPLE\t0\t0\t0\t-9" + os.linesep).encode()
    )
    assert (
        prefix.with_suffix(".bim").read_bytes() == ("1\trs1\t0\t100\tA\tG" + os.linesep).encode()
    )
    assert prefix.with_suffix(".bed").read_bytes() == b"\x6c\x1b\x01\x02"

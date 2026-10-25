# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Bounded synthetic failures exercise descriptor ownership and publication boundaries."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from allelix import reports
from allelix.exporters.plink import export_plink
from allelix.models import Variant
from allelix.reports.vcf import render_vcf

_VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
)


def _writer(tmp_path, kind, prior):
    source = tmp_path / "input.vcf"
    source.write_text(_VCF, encoding="utf-8")
    if kind == "plink":
        prefix = tmp_path / "out"
        outputs = [prefix.with_suffix(ext) for ext in (".fam", ".bim", ".bed")]

        def run():
            return export_plink(
                iter([Variant("rs1", "1", 100, "A", "G", build="GRCh37")]), prefix, "GRCh37"
            )
    else:
        outputs = [tmp_path / ("out.vcf" if kind == "vcf" else "out.txt")]

        def run():
            if kind == "vcf":
                return render_vcf(
                    input_path=source,
                    output_path=outputs[0],
                    build="GRCh37",
                    annotators_used=[],
                    annotations_by_key={},
                )
            return reports.atomic_write_text(outputs[0], "complete content\n")

    if prior:
        for path in outputs:
            path.write_bytes(b"prior")
    return run, outputs


def _track_allocations(monkeypatch):
    created = []
    allocate = reports.tempfile.mkstemp

    def tracked(**kwargs):
        descriptor, name = allocate(**kwargs)
        created.append((descriptor, Path(name)))
        assert Path(name).parent == Path(kwargs["dir"])
        assert Path(name).name.startswith(kwargs["prefix"])
        assert Path(name).name.endswith(".tmp")
        if os.name != "nt":
            assert stat.S_IMODE(os.fstat(descriptor).st_mode) == 0o600
        return descriptor, name

    monkeypatch.setattr(reports.tempfile, "mkstemp", tracked)
    return created


@pytest.mark.parametrize("kind", ["text", "vcf", "plink"])
@pytest.mark.parametrize("phase", ["fdopen", "write", "close"])
@pytest.mark.parametrize("failure", [OSError, KeyboardInterrupt])
@pytest.mark.parametrize("prior", [False, True])
def test_prepublication_failures_close_descriptors_and_preserve_outputs(
    tmp_path, monkeypatch, kind, phase, failure, prior
):
    run, outputs = _writer(tmp_path, kind, prior)
    created = _track_allocations(monkeypatch)
    fdopen = reports.os.fdopen
    primary = failure("synthetic primary failure")
    selected = []

    class Faulty:
        def __init__(self, raw):
            self.raw = raw
            self.failed = False

        @property
        def closed(self):
            return self.raw.closed

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.close()

        def write(self, value):
            if phase == "write":
                raise primary
            return self.raw.write(value)

        def close(self):
            self.raw.close()
            if phase == "close" and not self.failed:
                self.failed = True
                raise primary

    def opened(descriptor, *args, **kwargs):
        if not selected:
            selected.append(descriptor)
            if phase == "fdopen":
                raise primary
            return Faulty(fdopen(descriptor, *args, **kwargs))
        return fdopen(descriptor, *args, **kwargs)

    monkeypatch.setattr(reports.os, "fdopen", opened)
    with pytest.raises(failure, match="synthetic primary failure") as caught:
        run()
    assert caught.value is primary
    for path in outputs:
        if prior:
            assert path.read_bytes() == b"prior"
        else:
            assert not path.exists()
    assert created
    for descriptor, path in created:
        assert not path.exists()
        with pytest.raises(OSError):
            os.fstat(descriptor)


@pytest.mark.parametrize("failure_at", [1, 2, 3])
def test_partial_plink_allocation_removes_only_owned_temps(tmp_path, monkeypatch, failure_at):
    run, outputs = _writer(tmp_path, "plink", True)
    original = reports.tempfile.mkstemp
    calls, created = [], []

    def allocate(**kwargs):
        calls.append(kwargs)
        if len(calls) == failure_at:
            raise OSError("synthetic allocation failure")
        result = original(**kwargs)
        created.append(result)
        return result

    monkeypatch.setattr(reports.tempfile, "mkstemp", allocate)
    with pytest.raises(OSError, match="allocation failure"):
        run()
    assert [path.read_bytes() for path in outputs] == [b"prior"] * 3
    for descriptor, name in created:
        assert not Path(name).exists()
        with pytest.raises(OSError):
            os.fstat(descriptor)


def test_overlapping_stages_are_independent_and_closed_before_publication(tmp_path):
    destination = tmp_path / "out.txt"
    destination.write_bytes(b"prior")
    with reports.staged_output(destination) as first, reports.staged_output(destination) as second:
        assert first.temporary != second.temporary
        with first.open_text() as stream:
            stream.write("first")
        with second.open_text() as stream:
            stream.write("second")
        assert destination.read_bytes() == b"prior"
        first.publish()
        assert destination.read_bytes() == b"first"
        assert second.temporary.read_bytes() == b"second"
        second.publish()
    assert destination.read_bytes() == b"second"
    assert not list(tmp_path.glob(".out.txt.*.tmp"))


@pytest.mark.parametrize("mode", [0o400, 0o600, 0o640, 0o755])
def test_existing_ordinary_modes_preserved_or_old_file_kept_on_denial(tmp_path, mode):
    target = tmp_path / "out.txt"
    target.write_bytes(b"prior")
    target.chmod(mode)
    expected = target.stat().st_mode & 0o777
    try:
        try:
            reports.atomic_write_text(target, "complete")
        except PermissionError:
            assert os.name == "nt"
            assert target.read_bytes() == b"prior"
        else:
            assert target.read_bytes() == b"complete"
        assert target.stat().st_mode & 0o777 == expected
        assert not list(tmp_path.glob(".out.txt.*.tmp"))
    finally:
        target.chmod(stat.S_IREAD | stat.S_IWRITE)


def test_cleanup_failure_preserves_original_error_and_unrelated_files(tmp_path, monkeypatch):
    target, unrelated = tmp_path / "out.txt", tmp_path / "unrelated.txt"
    unrelated.write_bytes(b"unrelated")
    primary = OSError("primary write failure")
    unlink = Path.unlink
    owned = []

    def deny(path, *args, **kwargs):
        if path in owned:
            raise PermissionError("secondary cleanup failure")
        return unlink(path, *args, **kwargs)

    try:
        with (
            pytest.raises(OSError, match="primary write failure") as caught,
            reports.staged_output(target) as staged,
        ):
            owned.append(staged.temporary)
            monkeypatch.setattr(Path, "unlink", deny)
            raise primary
        assert caught.value is primary
        assert unrelated.read_bytes() == b"unrelated"
        assert not target.exists()
    finally:
        for path in owned:
            path.chmod(stat.S_IREAD | stat.S_IWRITE)
            unlink(path)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Direct report APIs protect inputs and own their temporary files."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from allelix import reports
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from allelix.reports.vcf import render_vcf

_VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
)


def _result(source):
    return AnalysisResult(
        file_path=source,
        parser_name="23andme",
        parser_display_name="23andMe",
        sample_id="SAMPLE",
        build="GRCh37",
        total_variants=1,
        skipped_count=0,
        annotators_used=[],
        annotations=[],
    )


@pytest.mark.parametrize("writer", [render_json, render_html])
@pytest.mark.parametrize("alias", ["direct", "symlink", "hardlink", "legacy", "input-symlink"])
def test_direct_renderer_rejects_input_alias(tmp_path, writer, alias):
    source, output = tmp_path / "input.txt", tmp_path / "report.json"
    source.write_text("synthetic input", encoding="utf-8")
    if alias == "direct":
        output = source
    elif alias == "symlink":
        output.symlink_to(source.name)
    elif alias == "hardlink":
        os.link(source, output)
    elif alias == "input-symlink":
        source.rename(output)
        source.symlink_to(output.name)
    else:
        source.rename(output.with_name(output.name + ".tmp"))
        source = output.with_name(output.name + ".tmp")
    before = source.read_bytes()
    with pytest.raises(ValueError, match=r"alias.*input|input.*alias"):
        writer(_result(source), output_path=output)
    assert source.read_bytes() == before


@pytest.mark.parametrize("alias", ["direct", "symlink", "hardlink", "input-symlink"])
def test_direct_vcf_rejects_input_alias(tmp_path, alias):
    source, output = tmp_path / "input.vcf", tmp_path / "annotated.vcf"
    source.write_text(_VCF, encoding="utf-8")
    if alias == "direct":
        output = source
    elif alias == "symlink":
        output.symlink_to(source.name)
    elif alias == "hardlink":
        os.link(source, output)
    else:
        source.rename(output)
        source.symlink_to(output.name)
    before = source.read_bytes()
    with pytest.raises(ValueError, match=r"alias.*input|input.*alias"):
        render_vcf(
            input_path=source,
            output_path=output,
            build="GRCh37",
            annotators_used=[],
            annotations_by_key={},
        )
    assert source.read_bytes() == before


@pytest.mark.parametrize("failure", [False, True])
@pytest.mark.parametrize("legacy_kind", ["file", "symlink", "hardlink"])
def test_text_writer_leaves_unrelated_fixed_temp_untouched(
    tmp_path, monkeypatch, failure, legacy_kind
):
    target, victim = tmp_path / "out.txt", tmp_path / "unrelated.txt"
    target.write_text("prior", encoding="utf-8")
    victim.write_text("unrelated", encoding="utf-8")
    legacy = target.with_name(target.name + ".tmp")
    if legacy_kind == "file":
        legacy.write_text("legacy", encoding="utf-8")
    elif legacy_kind == "symlink":
        legacy.symlink_to(victim.name)
    else:
        os.link(victim, legacy)
    before = legacy.read_bytes()
    if failure:

        def denied(*args):
            raise OSError("synthetic publish failure")

        monkeypatch.setattr(reports.os, "replace", denied)
        with pytest.raises(OSError, match="publish failure"):
            reports.atomic_write_text(target, "replacement")
    else:
        reports.atomic_write_text(target, "replacement")
    assert legacy.read_bytes() == before
    assert victim.read_text(encoding="utf-8") == "unrelated"
    assert target.read_text(encoding="utf-8") == ("prior" if failure else "replacement")


@pytest.mark.parametrize("prior", [False, True])
def test_vcf_midstream_failure_preserves_destination(tmp_path, monkeypatch, prior):
    from allelix.reports import vcf

    source, output = tmp_path / "input.vcf", tmp_path / "out.vcf"
    source.write_text(_VCF, encoding="utf-8")
    if prior:
        output.write_bytes(b"prior VCF")

    def fail(*args, **kwargs):
        raise OSError("synthetic row failure")

    monkeypatch.setattr(vcf, "_stamp_data_line", fail)
    with pytest.raises(OSError, match="row failure"):
        render_vcf(
            input_path=source,
            output_path=output,
            build="GRCh37",
            annotators_used=[],
            annotations_by_key={},
        )
    assert output.read_bytes() == b"prior VCF" if prior else not output.exists()
    assert source.read_text(encoding="utf-8") == _VCF
    assert not list(tmp_path.glob(".out.vcf.*.tmp"))


def test_preserved_mode_excludes_privileged_bits(tmp_path, monkeypatch):
    import stat
    from types import SimpleNamespace

    target = tmp_path / "out.txt"
    original_stat = Path.stat

    def mode(path, *args, **kwargs):
        if path == target:
            return SimpleNamespace(st_mode=stat.S_IFREG | 0o7640)
        return original_stat(path, *args, **kwargs)

    monkeypatch.setattr(Path, "stat", mode)
    with reports.staged_output(target) as staged:
        assert staged.mode == 0o640


def test_published_temp_name_is_no_longer_owned(tmp_path):
    target = tmp_path / "out.txt"
    with reports.staged_output(target) as staged:
        former = staged.temporary
        with staged.open_text() as stream:
            stream.write("published")
        staged.publish()
        former.write_text("new unrelated file", encoding="utf-8")
    assert former.read_text(encoding="utf-8") == "new unrelated file"
    assert target.read_text(encoding="utf-8") == "published"


def test_windows_readonly_unpublished_temp_is_removed_without_changing_final(
    tmp_path, monkeypatch
):
    import stat
    from types import SimpleNamespace

    target, unrelated = tmp_path / "out.txt", tmp_path / "unrelated.txt"
    target.write_bytes(b"prior")
    target.chmod(stat.S_IREAD)
    unrelated.write_bytes(b"unrelated")
    old_mode = stat.S_IMODE(target.stat().st_mode)
    unlink, chmod = Path.unlink, Path.chmod
    changes = []
    owned = []
    allocate = reports.tempfile.mkstemp

    def allocated(**kwargs):
        descriptor, name = allocate(**kwargs)
        owned.append(Path(name))
        return descriptor, name

    def publish(*args):
        raise PermissionError("primary publication failure")

    def remove(path, *args, **kwargs):
        if path in owned and path.exists() and not path.stat().st_mode & stat.S_IWUSR:
            raise PermissionError("Windows readonly temp")
        return unlink(path, *args, **kwargs)

    def permissions(path, mode, *args, **kwargs):
        changes.append(path)
        return chmod(path, mode, *args, **kwargs)

    monkeypatch.setattr(
        reports,
        "os",
        SimpleNamespace(name="nt", fdopen=os.fdopen, close=os.close, replace=publish),
    )
    monkeypatch.setattr(reports.tempfile, "mkstemp", allocated)
    monkeypatch.setattr(Path, "unlink", remove)
    monkeypatch.setattr(Path, "chmod", permissions)
    try:
        with pytest.raises(PermissionError, match="primary publication failure"):
            reports.atomic_write_text(target, "replacement")
        assert target.read_bytes() == b"prior"
        assert stat.S_IMODE(target.stat().st_mode) == old_mode
        assert unrelated.read_bytes() == b"unrelated"
        assert all(path in owned for path in changes)
        assert owned and all(not path.exists() for path in owned)
    finally:
        chmod(target, stat.S_IREAD | stat.S_IWRITE)
        for path in owned:
            if path.exists():
                chmod(path, stat.S_IREAD | stat.S_IWRITE)
                unlink(path)


def test_input_stat_error_cannot_certify_a_missing_destination_safe(tmp_path, monkeypatch):
    source, destination = tmp_path / "input.txt", tmp_path / "new.json"
    source.write_bytes(b"synthetic")
    original = Path.stat

    def denied(path, *args, **kwargs):
        if path == source:
            raise PermissionError("synthetic input stat denial")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "stat", denied)
    with pytest.raises(PermissionError, match="input stat denial"):
        reports.validate_output_path(source, destination)
    assert not destination.exists()


@pytest.mark.parametrize("operation", ["resolve", "stat", "samefile"])
def test_real_destination_identity_errors_are_not_ignored(tmp_path, monkeypatch, operation):
    source, destination = tmp_path / "input.txt", tmp_path / "out.json"
    source.write_bytes(b"input")
    destination.write_bytes(b"prior")
    original = getattr(Path, operation)

    def denied(path, *args, **kwargs):
        if path == destination:
            raise PermissionError("synthetic identity denial")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, operation, denied)
    with pytest.raises(PermissionError, match="identity denial"):
        reports.validate_output_path(source, destination)
    assert source.read_bytes() == b"input"
    assert destination.read_bytes() == b"prior"

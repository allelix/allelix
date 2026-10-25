# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ZIP resources are disposed even when construction or explicit close fails."""

from __future__ import annotations

import gc
import io
import zipfile

import pytest

import allelix.parsers._helpers as helpers


def _archive(tmp_path, monkeypatch, kind="good", compression=zipfile.ZIP_STORED):
    path = tmp_path / "synthetic.zip"
    if kind == "directory":
        path.write_bytes(b"PK\x03\x04" + bytes(40))
    else:
        with zipfile.ZipFile(path, "w", compression=compression) as archive:
            archive.writestr("synthetic.txt", "synthetic\n")
            if kind == "multiple":
                archive.writestr("second.txt", "second\n")
        if kind == "password":
            raw = bytearray(path.read_bytes())
            raw[6] |= 1
            raw[raw.index(b"PK\x01\x02") + 8] |= 1
            path.write_bytes(raw)
        if kind == "cap":
            monkeypatch.setattr(helpers, "ZIP_MAX_MEMBER_BYTES", 1)
    return path


def _assert_disposed(archives, members):
    assert all(archive.fp is None for archive in archives)
    assert all(member.closed for member in members)


@pytest.mark.parametrize(
    "kind,error,match",
    [
        ("directory", helpers.DamagedFileError, "zip directory is missing or corrupt"),
        ("multiple", helpers.InputContainerError, "holding 2 files"),
        ("cap", helpers.InputContainerError, "holds a file over"),
        ("password", helpers.InputContainerError, "password-protected"),
    ],
)
def test_failed_initialization_explicit_close_twice(tmp_path, monkeypatch, kind, error, match):
    path = _archive(tmp_path, monkeypatch, kind)
    wrapper = helpers._ZipMember.__new__(helpers._ZipMember)
    try:
        with pytest.raises(error, match=match):
            wrapper.__init__(path)
        archive = getattr(wrapper, "_archive", None)
        assert archive is None or archive.fp is None
        wrapper.close()
        wrapper.close()
        assert wrapper.closed
    finally:
        archive = getattr(wrapper, "_archive", None)
        if archive is not None:
            archive.close()
        io.RawIOBase.close(wrapper)


@pytest.mark.parametrize("kind", ["directory", "multiple", "cap", "password"])
def test_failed_construction_natural_collection_has_no_unraisable(tmp_path, monkeypatch, kind):
    path = _archive(tmp_path, monkeypatch, kind)

    def construct():
        with pytest.raises((helpers.DamagedFileError, helpers.InputContainerError)):
            helpers._ZipMember(path)

    construct()
    gc.collect()  # Existing warning-as-error filters remain active, including on 3.13.


def test_member_acquisition_failure_preserves_primary_and_closes_archive(tmp_path, monkeypatch):
    path = _archive(tmp_path, monkeypatch)
    acquired = []
    real = zipfile.ZipFile

    def acquire(path):
        archive = real(path)
        acquired.append(archive)

        def refuse(*args, **kwargs):
            raise OSError("member acquisition refused")

        monkeypatch.setattr(archive, "open", refuse)
        return archive

    monkeypatch.setattr(helpers.zipfile, "ZipFile", acquire)
    partial = helpers._ZipMember.__new__(helpers._ZipMember)
    try:
        with pytest.raises(OSError, match="member acquisition refused"):
            partial.__init__(path)
        assert acquired[0].fp is None
        partial.close()
        assert partial.closed
    finally:
        acquired[0].close()
        io.RawIOBase.close(partial)


def test_constructor_cleanup_error_does_not_replace_refusal(tmp_path, monkeypatch):
    path = _archive(tmp_path, monkeypatch, "multiple")
    real = zipfile.ZipFile
    archives = []

    def acquire(path):
        archive = real(path)
        archives.append(archive)
        original = archive.close

        def close():
            original()
            raise RuntimeError("secondary archive cleanup failure")

        monkeypatch.setattr(archive, "close", close)
        return archive

    monkeypatch.setattr(helpers.zipfile, "ZipFile", acquire)
    partial = helpers._ZipMember.__new__(helpers._ZipMember)
    try:
        with pytest.raises(helpers.InputContainerError, match="holding 2 files"):
            partial.__init__(path)
        assert archives[0].fp is None
        partial.close()
        assert partial.closed
    finally:
        io.RawIOBase.close(partial)
        monkeypatch.undo()
        for archive in archives:
            archive.close()


@pytest.mark.parametrize("close_first", [False, True])
def test_member_close_failure_attempts_archive_and_base_once(tmp_path, monkeypatch, close_first):
    path = _archive(tmp_path, monkeypatch)
    wrapper = helpers._ZipMember(path)
    archive, member = wrapper._archive, wrapper._member
    original = member.close
    calls = []

    def close():
        calls.append("member")
        if close_first:
            original()
        raise RuntimeError("member cleanup failed")

    monkeypatch.setattr(member, "close", close)
    try:
        with pytest.raises(RuntimeError, match="member cleanup failed"):
            wrapper.close()
        assert archive.fp is None
        assert wrapper.closed
        assert member.closed is close_first  # A failed-before-close double did not release itself.
        wrapper.close()
        assert calls == ["member"]
    finally:
        monkeypatch.setattr(member, "close", original)
        original()
        archive.close()
        io.RawIOBase.close(wrapper)


@pytest.mark.parametrize("compression", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
@pytest.mark.parametrize("early_exit", [False, True])
def test_normal_member_read_and_repeated_close(tmp_path, monkeypatch, compression, early_exit):
    path = _archive(tmp_path, monkeypatch, compression=compression)
    wrapper = helpers._ZipMember(path)
    archive, member = wrapper._archive, wrapper._member
    try:
        assert wrapper.read(1 if early_exit else -1) == (b"s" if early_exit else b"synthetic\n")
        wrapper.close()
        wrapper.close()
        _assert_disposed([archive], [member])
        assert wrapper.closed
    finally:
        member.close()
        archive.close()
        io.RawIOBase.close(wrapper)


def test_corrupt_member_read_exit_closes_resources(tmp_path, monkeypatch):
    path = _archive(tmp_path, monkeypatch)
    raw = bytearray(path.read_bytes())
    raw[30 + len("synthetic.txt")] ^= 1
    path.write_bytes(raw)
    wrapper = helpers._ZipMember(path)
    archive, member = wrapper._archive, wrapper._member
    with pytest.raises(helpers.DamagedFileError, match="corrupt zip member"), wrapper:
        wrapper.read()
    _assert_disposed([archive], [member])
    assert wrapper.closed


@pytest.mark.parametrize("member_fails", [False, True])
def test_archive_cleanup_failure_still_closes_base_and_preserves_first_error(
    tmp_path, monkeypatch, member_fails
):
    path = _archive(tmp_path, monkeypatch)
    wrapper = helpers._ZipMember(path)
    archive, member = wrapper._archive, wrapper._member
    archive_close, member_close = archive.close, member.close
    calls = []

    def close_archive():
        calls.append("archive")
        archive_close()
        raise RuntimeError("archive cleanup failed")

    def close_member():
        calls.append("member")
        member_close()
        if member_fails:
            raise RuntimeError("member cleanup failed")

    monkeypatch.setattr(archive, "close", close_archive)
    monkeypatch.setattr(member, "close", close_member)
    try:
        with pytest.raises(
            RuntimeError,
            match="member cleanup failed" if member_fails else "archive cleanup failed",
        ):
            wrapper.close()
        assert calls == ["member", "archive"]
        _assert_disposed([archive], [member])
        assert wrapper.closed
        wrapper.close()
        assert calls == ["member", "archive"]
    finally:
        monkeypatch.setattr(archive, "close", archive_close)
        monkeypatch.setattr(member, "close", member_close)
        member_close()
        archive_close()
        io.RawIOBase.close(wrapper)

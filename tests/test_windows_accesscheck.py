# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Native permission call contracts, independent of the host OS and real ACL fixtures."""

from __future__ import annotations

import ctypes
import json
import subprocess
import sys
from ctypes import wintypes
from pathlib import Path
from types import SimpleNamespace

import pytest

from allelix.databases import cache_state


def _set(pointer, kind, value):
    ctypes.cast(pointer, ctypes.POINTER(kind)).contents.value = value


class FakeSecurity:
    def __init__(self):
        self.error = 0
        self.thread = False
        self.allowed = True
        self.fail = {}
        self.descriptor_grows = 0
        self.privileges_grow = 0
        self.closed = []
        self.calls = []

    def _fails(self, name):
        if name in self.fail:
            self.error = self.fail[name]
            return True
        return False

    def GetFileSecurityW(self, path, flags, buffer, size, needed):  # noqa: N802 - native ABI name
        self.calls.append(("security", path, flags, size))
        if self._fails("security"):
            return False
        if buffer is None:
            _set(needed, wintypes.DWORD, 32)
            self.error = 122
            return False
        if self.descriptor_grows:
            self.descriptor_grows -= 1
            _set(needed, wintypes.DWORD, size * 2)
            self.error = 122
            return False
        _set(needed, wintypes.DWORD, 32)
        return True

    def GetCurrentThread(self):  # noqa: N802 - native ABI name
        return -2

    def GetCurrentProcess(self):  # noqa: N802 - native ABI name
        return -1

    def OpenThreadToken(self, handle, rights, as_self, output):  # noqa: N802 - native ABI name
        self.calls.append(("thread", handle, rights, as_self))
        if self._fails("thread"):
            return False
        if not self.thread:
            self.error = 1008
            return False
        _set(output, wintypes.HANDLE, 101)
        return True

    def OpenProcessToken(self, handle, rights, output):  # noqa: N802 - native ABI name
        self.calls.append(("process", handle, rights))
        if self._fails("process"):
            return False
        _set(output, wintypes.HANDLE, 102)
        return True

    def DuplicateToken(self, handle, level, output):  # noqa: N802 - native ABI name
        self.calls.append(("duplicate", handle.value, level))
        if self._fails("duplicate"):
            return False
        _set(output, wintypes.HANDLE, 202)
        return True

    def AccessCheck(self, sd, token, rights, mapping, privileges, size, granted, allowed):  # noqa: N802 - native ABI name
        length = ctypes.cast(size, ctypes.POINTER(wintypes.DWORD)).contents.value
        self.calls.append(("access", token.value, rights, length))
        if self._fails("access"):
            return False
        if self.privileges_grow:
            self.privileges_grow -= 1
            _set(size, wintypes.DWORD, length * 2)
            self.error = 122
            return False
        _set(granted, wintypes.DWORD, 2 if self.allowed else 0)
        _set(allowed, wintypes.BOOL, self.allowed)
        # LastError need not be clear after a successful API call.
        self.error = 5
        return True

    def CloseHandle(self, handle):  # noqa: N802 - native ABI name
        self.closed.append(handle.value)
        return not self._fails(f"close{handle.value}")


@pytest.fixture
def security(monkeypatch: pytest.MonkeyPatch):
    native = FakeSecurity()
    monkeypatch.setattr(cache_state, "_windows_security_api", lambda: (native, native))
    monkeypatch.setattr(ctypes, "get_last_error", lambda: native.error, raising=False)
    return native


@pytest.mark.parametrize("thread", [False, True])
@pytest.mark.parametrize("allowed", [False, True])
def test_effective_token_and_denied_access_status_are_distinct_from_api_failure(
    tmp_path: Path,
    security,
    thread,
    allowed,
) -> None:
    security.thread, security.allowed = thread, allowed
    assert cache_state._windows_directory_writable(tmp_path) is allowed
    assert security.closed == [202, 101 if thread else 102]
    assert ("thread", -2, 10, True) in security.calls
    assert bool([c for c in security.calls if c[0] == "process"]) is not thread
    assert ("duplicate", 101 if thread else 102, 2) in security.calls
    assert all(c[2] == 2 for c in security.calls if c[0] == "access")


def test_descriptor_and_privilege_buffers_retry_when_sizes_grow(tmp_path: Path, security) -> None:
    security.descriptor_grows = security.privileges_grow = 1
    assert cache_state._windows_directory_writable(tmp_path)
    assert len([c for c in security.calls if c[0] == "security"]) == 3
    lengths = [c[3] for c in security.calls if c[0] == "access"]
    assert len(lengths) == 2 and lengths[1] > lengths[0]
    assert security.closed == [202, 102]


@pytest.mark.parametrize(
    "operation,closed",
    [
        ("security", []),
        ("thread", []),
        ("process", []),
        ("duplicate", [102]),
        ("access", [202, 102]),
    ],
)
def test_native_api_errors_close_every_acquired_token(
    tmp_path: Path,
    security,
    operation,
    closed,
) -> None:
    security.fail[operation] = 5
    with pytest.raises(OSError):
        cache_state._windows_directory_writable(tmp_path)
    assert security.closed == closed
    if operation == "thread":
        assert not [c for c in security.calls if c[0] == "process"]


@pytest.mark.parametrize("buffer", ["descriptor_grows", "privileges_grow"])
def test_unstable_buffer_sizes_stop_with_an_error(tmp_path: Path, security, buffer) -> None:
    setattr(security, buffer, 100)
    with pytest.raises(OSError):
        cache_state._windows_directory_writable(tmp_path)
    assert len(security.calls) < 20
    if buffer == "privileges_grow":
        assert security.closed == [202, 102]


def test_zero_sized_descriptor_is_not_accepted(tmp_path: Path, security, monkeypatch) -> None:
    def zero_size(_path, _flags, _buffer, _size, needed):
        _set(needed, wintypes.DWORD, 0)
        security.error = 122
        return False

    monkeypatch.setattr(security, "GetFileSecurityW", zero_size)
    with pytest.raises(OSError):
        cache_state._windows_directory_writable(tmp_path)
    assert security.closed == []


@pytest.mark.parametrize("primary", [False, True])
def test_cleanup_attempts_both_closes_and_preserves_primary_error(
    tmp_path: Path, security, primary
) -> None:
    security.fail["close202"] = 6
    if primary:
        security.fail["access"] = 5
    with pytest.raises(OSError, match="AccessCheck" if primary else "CloseHandle"):
        cache_state._windows_directory_writable(tmp_path)
    assert security.closed == [202, 102]


def test_windows_api_error_is_conservative_without_os_access_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog,
) -> None:
    monkeypatch.setattr(cache_state, "os", SimpleNamespace(name="nt"))

    def unavailable(_path):
        raise OSError("synthetic permission service failure")

    monkeypatch.setattr(cache_state, "_windows_directory_writable", unavailable)
    assert not cache_state.directory_writable(tmp_path)
    assert "Cannot check write access" in caplog.text


@pytest.mark.parametrize("search,mask", [(False, 2), (True, 3)])
def test_posix_access_masks_are_preserved(tmp_path: Path, monkeypatch, search, mask) -> None:
    calls = []
    monkeypatch.setattr(
        cache_state,
        "os",
        SimpleNamespace(
            name="posix",
            W_OK=2,
            X_OK=1,
            access=lambda path, mode: calls.append((path, mode)) or True,
        ),
    )
    assert cache_state.directory_writable(tmp_path, require_search=search)
    assert calls == [(tmp_path, mask)]


def test_token_closes_still_execute_under_optimized_python(tmp_path: Path) -> None:
    code = """
import ctypes,json,runpy,sys
from pathlib import Path
from allelix.databases import cache_state
scope=runpy.run_path(sys.argv[1]); native=scope['FakeSecurity']()
ctypes.get_last_error=lambda:native.error
cache_state._windows_security_api=lambda:(native,native)
result=cache_state._windows_directory_writable(Path.cwd())
print(json.dumps({'allowed':result,'closed':native.closed}))
"""
    result = subprocess.run(
        [sys.executable, "-O", "-c", code, str(Path(__file__).resolve())],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {"allowed": True, "closed": [202, 102]}


def test_db_check_honors_denial_when_metadata_access_allows(tmp_path, monkeypatch):
    from click.testing import CliRunner

    from allelix.cli import main

    monkeypatch.setattr(cache_state.os, "access", lambda *_a: True)
    monkeypatch.setattr(cache_state, "directory_writable", lambda *_a, **_kw: False)
    runner = CliRunner()
    args = ["db", "path", "--data-dir", str(tmp_path)]
    assert runner.invoke(main, args).exit_code == 0
    result = runner.invoke(main, [*args, "--check"])
    assert result.exit_code == 1
    assert "not writable" in result.output

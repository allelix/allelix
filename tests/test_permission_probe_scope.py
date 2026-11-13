"""Permission guards must probe the operation and preserve ACL-backed coverage."""

import ast
import builtins
import errno
import io
import os
import stat
from pathlib import Path

import pytest
from test_publish_release_gates import ROOT

READ_FIXTURE = "require_unreadable_files"
ACL_CASES = (
    ("test_readonly_caches.py", "test_read_only_data_directory"),
    ("test_readonly_caches.py", "test_unreadable"),
    ("test_readonly_caches.py", "test_readonly_wal"),
    ("test_readonly_caches.py", "test_update_step_leaves_a_read_only_directory_alone"),
    ("test_readonly_caches.py", "test_wal_cache_that_cannot_be_converted_is_not_downloaded"),
    ("test_cli.py", "test_db_path_check_fails_on_readonly_dir"),
)


def _node(filename, function):
    tree = ast.parse((ROOT / "tests" / filename).read_text())
    return next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == function)


@pytest.mark.parametrize(
    "filename,function",
    (
        ("test_clinvar_damaged_cache_replaced.py", "test_an_unreadable_cache_is_still_refused"),
        (
            "test_early_output_checks.py",
            "test_an_unreadable_diff_baseline_is_refused_by_the_option_check",
        ),
    ),
)
def test_chmod_unreadable_cases_request_the_read_probe(filename, function):
    node = _node(filename, function)
    assert READ_FIXTURE in ast.unparse(node), "Write denial does not establish read denial"
    assert "require_readonly_filesystem" not in {a.arg for a in node.args.args}
    assert not any("require_readonly_filesystem" in ast.unparse(d) for d in node.decorator_list)


@pytest.mark.parametrize("uid", (0, 1000))
@pytest.mark.parametrize("outcome", ("allowed", "denied", "io-error"))
def test_read_probe_uses_an_actual_read_and_restores_access(
    tmp_path, monkeypatch, request, uid, outcome
):
    tree = ast.parse((ROOT / "tests/conftest.py").read_text())
    assert any(isinstance(n, ast.FunctionDef) and n.name == READ_FIXTURE for n in tree.body), (
        "Supply the SUP-approved shared read-probe sibling"
    )
    monkeypatch.setattr(os, "geteuid", lambda: uid, raising=False)
    attempts = []
    for owner in (builtins, io, os):
        original = owner.open

        def open_file(path, mode="r", *args, _original=original, **kwargs):
            if isinstance(path, (str, bytes, os.PathLike)):
                target = Path(os.fsdecode(path))
                reading = (
                    not mode & os.O_WRONLY if isinstance(mode, int) else "r" in mode or "+" in mode
                )
                if reading and target.is_relative_to(tmp_path) and target.is_file():
                    permissions = stat.S_IMODE(target.stat().st_mode)
                    if not permissions & 0o444:
                        attempts.append(target)
                        if outcome == "denied":
                            raise PermissionError(
                                errno.EACCES, "synthetic read denied", str(target)
                            )
                        if outcome == "io-error":
                            raise OSError(errno.EIO, "synthetic I/O error", str(target))
                        target.chmod(permissions | 0o600)
                        try:
                            return _original(path, mode, *args, **kwargs)
                        finally:
                            target.chmod(permissions)
            return _original(path, mode, *args, **kwargs)

        monkeypatch.setattr(owner, "open", open_file)
    if outcome == "allowed":
        with pytest.raises(pytest.skip.Exception) as caught:
            request.getfixturevalue(READ_FIXTURE)
        assert "read" in str(caught.value).lower()
    elif outcome == "io-error":
        with pytest.raises(OSError) as caught:
            request.getfixturevalue(READ_FIXTURE)
        assert caught.value.errno == errno.EIO
    else:
        request.getfixturevalue(READ_FIXTURE)
    assert attempts, "Must attempt to read a file with its read bits removed"
    assert all(not p.exists() for p in attempts), "Remove probe files even on skip or I/O error"


@pytest.mark.parametrize("filename,function", ACL_CASES)
def test_acl_cases_are_not_unconditionally_gated_by_a_chmod_only_probe(filename, function):
    node = _node(filename, function)
    guards = {a.arg for a in node.args.args}
    guards.update(ast.unparse(d) for d in node.decorator_list)
    assert not any("require_readonly_filesystem" in g or READ_FIXTURE in g for g in guards), (
        "These tests already use _denied_permissions' Windows ACL path; a chmod-only "
        "fixture must not skip them before it runs. Probe through that helper, or use "
        "a separate ACL-aware guard."
    )


@pytest.mark.parametrize("read", (False, True))
def test_acl_capable_helper_skips_a_bypassed_permission_after_restoring_mode(
    tmp_path, monkeypatch, read
):
    from test_readonly_caches import _denied_permissions

    path = tmp_path / "file" if read else tmp_path / "directory"
    if read:
        path.write_bytes(b"synthetic")
    else:
        path.mkdir()
    old_mode = stat.S_IMODE(path.stat().st_mode)
    real_open = Path.open

    def bypass(target, *args, **kwargs):
        guarded = target if read else target.parent
        mode = stat.S_IMODE(guarded.stat().st_mode)
        guarded.chmod(mode | (0o600 if read else 0o700))
        try:
            return real_open(target, *args, **kwargs)
        finally:
            guarded.chmod(mode)

    monkeypatch.setattr(Path, "open", bypass)
    with pytest.raises(pytest.skip.Exception), _denied_permissions(path, read=read):
        pytest.fail("Permission bypass must be reported as a skip")
    assert stat.S_IMODE(path.stat().st_mode) == old_mode
    assert (
        not (path / ".permission-denial-probe").exists()
        if not read
        else path.read_bytes() == b"synthetic"
    )


def test_contributing_names_the_read_capability_probe():
    text = (ROOT / "CONTRIBUTING.md").read_text()
    assert READ_FIXTURE in text

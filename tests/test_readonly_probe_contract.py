"""Permission-dependent safety tests are guarded by capability, never UID."""

import ast
import builtins
import errno
import io
import os
import stat
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = "require_readonly_filesystem"
CASES = (
    ("test_rc_unwritable_destinations.py", "test_report_permission_failure_precedes_analysis"),
    ("test_rc_unwritable_destinations.py", "test_vcf_permission_failure_precedes_report_write"),
    ("test_rc_unwritable_destinations.py", "test_plink_permission_failure_precedes_reading_calls"),
    ("test_rc_unwritable_destinations.py", "test_config_set_permission_error_preserves_config"),
    ("test_rc_unwritable_destinations.py", "test_db_clean_unwritable_directory_is_a_clean_error"),
    (
        "test_rc_unwritable_destinations.py",
        "test_db_update_unwritable_destination_precedes_download",
    ),
    (
        "test_rc_destination_edges.py",
        "test_new_data_directory_below_readonly_parent_is_a_clean_error",
    ),
    (
        "test_rc_destination_edges.py",
        "test_db_path_check_keeps_its_existing_readonly_directory_diagnostic",
    ),
    ("test_rc_destination_edges.py", "test_nonclinvar_readonly_cache_stops_before_download"),
    ("test_rc_extract_bad_index.py", "test_unreadable_vcf_is_a_clean_path_error"),
    (
        "test_rc_extract_bad_index.py",
        "test_bad_index_warns_and_returns_the_sequential_selected_sample",
    ),
)


def _definition():
    tree = ast.parse((ROOT / "tests/conftest.py").read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == FIXTURE]
    assert len(nodes) == 1, "Provide one shared capability-probe fixture in tests/conftest.py"
    return nodes[0]


def _call_probe(request):
    _definition()  # Assert the missing contract, rather than a fixture-lookup setup error.
    return request.getfixturevalue(FIXTURE)


def _simulate_write(monkeypatch, tmp_path, outcome):
    attempts = []
    for owner, name in ((io, "open"), (builtins, "open"), (os, "open")):
        original = getattr(owner, name)

        def open_file(path, mode="r", *args, _original=original, **kwargs):
            if isinstance(path, (str, bytes, os.PathLike)):
                target = Path(os.fsdecode(path))
                writing = (
                    bool(mode & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
                    if isinstance(mode, int)
                    else any(c in mode for c in "wax+")
                )
                parent = target.parent
                if writing and target.is_relative_to(tmp_path) and parent.exists():
                    permissions = stat.S_IMODE(parent.stat().st_mode)
                    if not permissions & 0o222:
                        attempts.append(target)
                        if outcome == "denied":
                            raise PermissionError(
                                errno.EACCES, "synthetic write denied", str(target)
                            )
                        if outcome == "io-error":
                            raise OSError(errno.EIO, "synthetic filesystem error", str(target))
                        parent.chmod(permissions | 0o700)
                        try:
                            return _original(path, mode, *args, **kwargs)
                        finally:
                            parent.chmod(permissions)
            return _original(path, mode, *args, **kwargs)

        monkeypatch.setattr(owner, name, open_file)
    return attempts


@pytest.mark.parametrize("uid", (0, 1000))
@pytest.mark.parametrize("outcome", ("allowed", "denied", "io-error"))
def test_probe_decides_from_a_real_write_attempt_not_user_id(
    tmp_path, monkeypatch, request, uid, outcome
):
    monkeypatch.setattr(os, "geteuid", lambda: uid, raising=False)
    attempts = _simulate_write(monkeypatch, tmp_path, outcome)
    if outcome == "allowed":
        with pytest.raises(pytest.skip.Exception) as caught:
            _call_probe(request)
        text = str(caught.value).lower()
        assert "writ" in text and ("read-only" in text or "permission" in text)
    elif outcome == "io-error":
        with pytest.raises(OSError) as caught:
            _call_probe(request)
        assert caught.value.errno == errno.EIO, "Unexpected filesystem errors must not be hidden"
    else:
        _call_probe(request)  # Permission denial means these tests can run, even with uid=0.
    assert attempts, "The shared fixture must actually attempt a write into a read-only folder"
    assert all(not path.exists() for path in attempts), "Probe files must be removed"
    assert all(
        stat.S_IMODE(path.stat().st_mode) & 0o700 == 0o700
        for path in tmp_path.rglob("*")
        if path.is_dir()
    ), "Restore probe-directory access for cleanup"


@pytest.mark.parametrize("filename,function", CASES)
def test_each_permission_dependent_test_uses_the_shared_probe(filename, function):
    tree = ast.parse((ROOT / "tests/rc_safety" / filename).read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function)
    uses = FIXTURE in {a.arg for a in node.args.args} or any(
        isinstance(n, ast.Constant) and n.value == FIXTURE for n in ast.walk(node)
    )
    module_marks = [
        n
        for n in tree.body
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "pytestmark" for t in n.targets)
    ]
    uses |= any(FIXTURE in ast.unparse(n) for n in module_marks)
    assert uses, f"{filename}::{function} still assumes chmod is enforced"


def test_probe_is_opt_in_and_corrupt_index_controls_do_not_require_permissions():
    fixture = _definition()
    assert not any(
        isinstance(n, ast.keyword)
        and n.arg == "autouse"
        and isinstance(n.value, ast.Constant)
        and n.value.value is True
        for d in fixture.decorator_list
        for n in ast.walk(d)
    ), "Do not skip unrelated tests"
    source = (ROOT / "tests/rc_safety/conftest.py").read_text()
    tree = ast.parse(source)
    case = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "rc_case")
    assert FIXTURE not in ast.unparse(case), "The common rc_case must not skip every safety test"
    tree = ast.parse((ROOT / "tests/rc_safety/test_rc_extract_bad_index.py").read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == CASES[-1][1])
    assert FIXTURE not in {a.arg for a in node.args.args}
    guarded = [
        n
        for n in ast.walk(node)
        if isinstance(n, ast.If)
        and "index_state" in ast.unparse(n.test)
        and "unreadable" in ast.unparse(n.test)
    ]
    assert any(FIXTURE in ast.unparse(n) for n in guarded), (
        "Only unreadable-index cases need the probe"
    )


def test_contributing_explains_the_capability_based_skip():
    text = (ROOT / "CONTRIBUTING.md").read_text().lower()
    lines = [line for line in text.splitlines() if "skip" in line and "writ" in line]
    assert any("read-only" in line or "permission" in line for line in lines)

"""Remaining permission tests must opt into capability checks, not UID skips."""

import ast

import pytest
from test_publish_release_gates import ROOT
from test_readonly_probe_contract import (
    test_probe_decides_from_a_real_write_attempt_not_user_id as _probe_contract,
)

FIXTURE = "require_readonly_filesystem"
FILES = (
    "test_readonly_caches.py",
    "test_cli.py",
    "annotators/test_gnomad_position_query.py",
    "test_clinvar_damaged_cache_replaced.py",
    "test_early_output_checks.py",
)
CASES = (
    (FILES[0], "test_read_only_data_directory"),
    (FILES[0], "test_read_only_file"),
    (FILES[0], "test_unreadable"),
    (FILES[0], "test_readonly_wal"),
    (FILES[0], "test_update_step_leaves_a_read_only_directory_alone"),
    (FILES[0], "test_wal_cache_that_cannot_be_converted_is_not_downloaded"),
    (FILES[1], "test_db_path_check_fails_on_readonly_dir"),
    (FILES[2], "test_resolver_runs_on_a_read_only_cache"),
    (FILES[3], "test_an_unreadable_cache_is_still_refused"),
    (FILES[4], "test_an_unreadable_diff_baseline_is_refused_by_the_option_check"),
)


def _tree(filename):
    return ast.parse((ROOT / "tests" / filename).read_text())


@pytest.mark.parametrize("filename", FILES)
def test_remaining_permission_modules_do_not_decide_by_uid(filename):
    offenders = [
        n.lineno
        for n in ast.walk(_tree(filename))
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr in ("geteuid", "getuid")
    ]
    assert not offenders, f"{filename} still decides permission coverage from UID at {offenders}"


@pytest.mark.parametrize("filename,function", CASES)
def test_each_remaining_permission_case_uses_the_shared_probe(filename, function):
    tree = _tree(filename)
    node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == function)
    read_case = filename in (FILES[3], FILES[4])
    fixture = "require_unreadable_files" if read_case else FIXTURE
    uses = fixture in {arg.arg for arg in node.args.args} or fixture in ast.unparse(node)
    # These helpers probe the operation with the platform's actual mechanism,
    # including Windows ACLs; do not demand a preceding chmod-only fixture.
    if filename in (FILES[0], FILES[1]):
        uses |= any(
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Name)
            and n.func.id in ("_read_only", "_denied_permissions")
            for n in ast.walk(node)
        )
    # A class mark is fine when all its cases need that capability; a module
    # blanket would skip unrelated read-only connection/schema controls.
    for parent in ast.walk(tree):
        if isinstance(parent, ast.ClassDef) and node in parent.body:
            uses |= any(fixture in ast.unparse(d) for d in parent.decorator_list)
    assert uses, f"{filename}::{function} has no shared permission-capability guard"


def test_readonly_cache_module_does_not_blanket_skip_unrelated_tests():
    tree = _tree(FILES[0])
    marks = [
        n
        for n in tree.body
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "pytestmark" for t in n.targets)
    ]
    assert not any("skip" in ast.unparse(n) or FIXTURE in ast.unparse(n) for n in marks), (
        "Read-only connection, schema and cache-state controls must still run on a writable probe"
    )


@pytest.mark.parametrize("uid", [0, 1000])
@pytest.mark.parametrize("outcome", ["allowed", "denied", "io-error"])
def test_shared_probe_still_uses_write_evidence_and_restores_access(
    tmp_path, monkeypatch, request, uid, outcome
):
    _probe_contract(tmp_path, monkeypatch, request, uid, outcome)

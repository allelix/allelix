# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Every sqlite3 connection that code opens is closed when the call returns.

``with sqlite3.connect(path) as conn:`` commits or rolls back on exit but
does not close the connection. The connection stays open until the garbage
collector finalizes it. Python 3.13 and later then emit ``ResourceWarning:
unclosed database``, which this suite turns into an error in whatever test
happens to be running when the collector runs. These tests do not depend on
the collector or on the Python version:

- a source scan finds no ``with sqlite3.connect(...)`` in ``allelix/`` or
  ``tests/`` (the context manager does not close);
- the test helper that builds the unknown-locus fixture, and the CLI runs on
  that fixture, are run while every ``sqlite3.connect`` is recorded, and each
  recorded connection must be closed when the call returns. A closed
  connection raises ``sqlite3.ProgrammingError`` on use.

All rows are synthetic.
"""

from __future__ import annotations

import ast
import sqlite3
from pathlib import Path

import pytest
import test_iyg_unknown_locus_contract as iyg
from click.testing import CliRunner

from allelix.cli import main

_ROOT = Path(__file__).resolve().parent.parent


def _is_sqlite_connect(node: ast.expr) -> bool:
    return (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "connect"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "sqlite3"
    )


def _with_connect_sites() -> list[str]:
    sites = []
    for top in ("allelix", "tests"):
        for path in sorted((_ROOT / top).rglob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, (ast.With, ast.AsyncWith)):
                    for item in node.items:
                        if _is_sqlite_connect(item.context_expr):
                            sites.append(f"{path.relative_to(_ROOT)}:{node.lineno}")
    return sites


def test_no_code_uses_a_sqlite_connection_as_its_own_context_manager():
    """``with sqlite3.connect(...)`` leaves the connection open: close it explicitly."""
    assert _with_connect_sites() == []


@pytest.fixture
def opened(monkeypatch):
    """Record every connection opened through ``sqlite3.connect``."""
    connections: list[sqlite3.Connection] = []
    real = sqlite3.connect

    def connect(*args, **kwargs):
        conn = real(*args, **kwargs)
        connections.append(conn)
        return conn

    monkeypatch.setattr(sqlite3, "connect", connect)
    yield connections
    # The asserts have run. Close what is still open so a leak found here is
    # not finalized later, inside another test, on Python 3.13 and later.
    for conn in connections:
        conn.close()


def _still_open(connections: list[sqlite3.Connection]) -> int:
    count = 0
    for conn in connections:
        try:
            conn.execute("SELECT 1")
        except sqlite3.ProgrammingError:
            continue
        count += 1
    return count


@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("reason", ("missing", "placed_same", "di_call"))
def test_the_unknown_locus_fixture_closes_every_connection(tmp_path, opened, reason, reverse):
    iyg._fixture(tmp_path, reason=reason, reverse=reverse)
    assert opened, "the fixture opens the synthetic caches"
    assert _still_open(opened) == 0


@pytest.mark.parametrize("suffix", ("json", "html"))
def test_analyze_closes_every_connection(tmp_path, opened, suffix):
    """Control on the product side: passes on RC1 and must keep passing."""
    data, path = iyg._fixture(tmp_path)
    # The fixture's own connections are checked in the test above.
    for conn in opened:
        conn.close()
    opened.clear()
    output = tmp_path / f"report.{suffix}"
    result = CliRunner().invoke(
        main,
        ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(output)],
    )
    assert result.exit_code == 0, result.output
    assert opened, "analyze opens the caches"
    assert _still_open(opened) == 0

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""No development-tracker reference (a ``#`` and three or four digits) in tracked text.

The repository's code, tests and documentation are published; an issue or pull
request number of the development tracker means nothing there. A CSS colour
value of a colour property or custom property, closed by ``;`` or ``}``
(``--border: #333;``), and an HTML character reference are not references.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = re.compile(r"(?<![&\w])#\d{3,4}\b")
CSS_PROPERTY = re.compile(
    r"(?:--[\w-]+|[\w-]*(?:color|background|border|fill|stroke|outline))\s*:\s*$"
)
CSS_END = re.compile(r"\s*[;}]")
SKIP_DIRS = {".git", ".venv", "venv", "build", "dist", "node_modules", "__pycache__"}


def references(line: str) -> list[str]:
    """The tracker references on ``line``, leaving out CSS colour values."""
    return [
        m.group()
        for m in REFERENCE.finditer(line)
        if not (CSS_PROPERTY.search(line[: m.start()]) and CSS_END.match(line, m.end()))
    ]


def _tracked_files() -> list[Path]:
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "-z"], capture_output=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return [
            p
            for p in ROOT.rglob("*")
            if p.is_file() and not SKIP_DIRS.intersection(p.relative_to(ROOT).parts)
        ]
    return [ROOT / name for name in out.decode().split("\0") if name]


@pytest.mark.parametrize(
    ("line", "found"),
    [
        (f"named in the report (ADR-0051, #{1234}).", [f"#{1234}"]),
        (f"# #{567}: a comment", [f"#{567}"]),
        (f"Issue: #{476}", [f"#{476}"]),
        (f"Naming: see #{477}; then", [f"#{477}"]),
        (f"--border: #{476} in prose", [f"#{476}"]),
        ("  --border: #333;", []),
        ("a { color: #333 }", []),
        ("color: #0066cc;", []),
        ("&#8217; is an apostrophe", []),
        ("rs1234 and chr1:#", []),
    ],
)
def test_reference_matcher(line: str, found: list[str]) -> None:
    assert references(line) == found


def test_no_tracker_reference_in_tracked_text() -> None:
    hits = []
    for path in _tracked_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for number, line in enumerate(text.splitlines(), 1):
            hits += [f"{path.relative_to(ROOT)}:{number}: {ref}" for ref in references(line)]
    assert not hits, "\n".join(hits)

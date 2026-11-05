#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Released ADRs are append-only (ADR-0001): the lock and its check.

The lock (``tests/data/adr_released.json``) holds, per released ADR file,
the number of lines released and a SHA-256 of those lines, plus a short
digest per line so a failure can name the first changed line. A locked
file may grow; its first ``lines`` lines may not change, with one
exception: the Status line may change to a "Superseded by ADR-NNNN" or
"Superseded in part by ADR-NNNN" form, when every ADR it names exists.
The index, ``docs/adr/README.md``, is not locked.

    scripts/adr_lock.py              refresh the lock at release: add new
                                     ADRs, grow each entry to the full
                                     file; refuses when a locked prefix
                                     changed
    scripts/adr_lock.py --covers     exit 1 unless the lock covers every
                                     ADR in full (release preflight)
    --adr-dir DIR / --lock FILE      other paths (seeding from a tag's
                                     docs/adr exported to DIR)

No git is needed: the check runs on a shallow checkout without tags.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ADR_DIR = ROOT / "docs" / "adr"
LOCK = ROOT / "tests" / "data" / "adr_released.json"
INDEX = "README.md"

# The three Status forms in the tree: "- **Status:** ...", "**Status:** ...",
# "Status: ...".
_STATUS = re.compile(r"^(?P<prefix>(?:- )?(?:\*\*)?Status:(?:\*\*)?)\s")
_SUPERSEDED = re.compile(r"^Superseded (?:in part )?by ADR-\d{4}\b")
_ADR_REF = re.compile(r"ADR-(\d{4})")


def read_lines(path: Path) -> list[str]:
    """The file's lines, without line endings; CRLF reads as LF."""
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def _sha256(lines: list[str]) -> str:
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def _line_digest(line: str) -> str:
    return hashlib.sha256(line.encode("utf-8")).hexdigest()[:12]


def adr_files(adr_dir: Path) -> list[Path]:
    """Every ADR file in ``adr_dir``; the index is not one."""
    return sorted(p for p in adr_dir.glob("*.md") if p.name != INDEX)


def _status_index(lines: list[str]) -> int | None:
    return next((i for i, line in enumerate(lines) if _STATUS.match(line)), None)


def entry(lines: list[str]) -> dict:
    """The lock entry for a file released with ``lines``."""
    status = _status_index(lines)
    return {
        "lines": len(lines),
        "sha256": _sha256(lines),
        "line_digests": [_line_digest(line) for line in lines],
        "status_line": None if status is None else status + 1,
        "status": None if status is None else lines[status],
    }


def _superseded_ok(released: str, new: str, adr_dir: Path) -> bool:
    """``new`` is a Superseded form of the ``released`` Status line, naming existing ADRs."""
    released_match, new_match = _STATUS.match(released), _STATUS.match(new)
    if released_match is None or new_match is None:
        return False
    if new_match.group("prefix") != released_match.group("prefix"):
        return False
    value = new[new_match.end() :].strip()
    if not _SUPERSEDED.match(value):
        return False
    return all(
        any(p.is_file() for p in adr_dir.glob(f"{number}-*.md"))
        for number in _ADR_REF.findall(value)
    )


def check_file(path: Path, locked: dict, adr_dir: Path) -> str | None:
    """Why ``path`` breaks its lock entry, or None when it keeps it."""
    if not path.exists():
        return f"{path.name}: released ADR deleted"
    lines = read_lines(path)
    count = locked["lines"]
    status = locked["status_line"]
    kept: list[str] = []
    # In file order, so the first difference is the one named: a deleted
    # line shows as the next line moved up, or as the end of the file.
    for i, digest in enumerate(locked["line_digests"]):
        if i >= len(lines):
            return f"{path.name}: line {i + 1} deleted (released with {count} lines)"
        line = lines[i]
        if status is not None and i == status - 1 and line != locked["status"]:
            if not _superseded_ok(locked["status"], line, adr_dir):
                return (
                    f"{path.name}: line {status} changed: the Status line may change only "
                    f"to a Superseded-by form naming an existing ADR: {line!r}"
                )
            line = locked["status"]
        if _line_digest(line) != digest:
            return f"{path.name}: line {i + 1} changed: {line!r}"
        kept.append(line)
    if _sha256(kept) != locked["sha256"]:
        return f"{path.name}: released lines changed (SHA-256 differs)"
    return None


def check(adr_dir: Path, lock: dict) -> list[str]:
    """Every way the ADRs in ``adr_dir`` break ``lock``."""
    problems = []
    for name, locked in sorted(lock["files"].items()):
        problem = check_file(adr_dir / name, locked, adr_dir)
        if problem:
            problems.append(problem)
    return problems


def refresh(adr_dir: Path, lock: dict) -> dict:
    """The lock grown to every ADR in ``adr_dir`` in full.

    Raises ValueError when a locked prefix changed.
    """
    problems = check(adr_dir, lock)
    if problems:
        raise ValueError("locked ADR lines changed:\n" + "\n".join(problems))
    return {
        "comment": lock.get("comment", _COMMENT),
        "files": {p.name: entry(read_lines(p)) for p in adr_files(adr_dir)},
    }


def uncovered(adr_dir: Path, lock: dict) -> list[str]:
    """ADRs the lock does not cover in full."""
    missing = []
    for path in adr_files(adr_dir):
        locked = lock["files"].get(path.name)
        if locked is None:
            missing.append(f"{path.name}: not locked")
        elif locked["lines"] != len(read_lines(path)):
            missing.append(
                f"{path.name}: locked {locked['lines']} of {len(read_lines(path))} lines"
            )
    return missing


_COMMENT = (
    "Released ADR lines (ADR-0001): per file, the line count released and their "
    "SHA-256. Refreshed at release by scripts/adr_lock.py; never edited by hand."
)


def load(path: Path) -> dict:
    """The lock at ``path``, or an empty one when there is none yet."""
    if not path.exists():
        return {"comment": _COMMENT, "files": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    """Refresh the lock, or with ``--covers`` check it; 0 on success."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--adr-dir", type=Path, default=ADR_DIR)
    parser.add_argument("--lock", type=Path, default=LOCK)
    parser.add_argument("--covers", action="store_true", help="check the lock covers every ADR")
    args = parser.parse_args(argv)
    lock = load(args.lock)
    if args.covers:
        problems = check(args.adr_dir, lock) + uncovered(args.adr_dir, lock)
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1 if problems else 0
    try:
        new = refresh(args.adr_dir, lock)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 1
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    args.lock.write_text(json.dumps(new, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"locked {len(new['files'])} ADRs in {args.lock}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

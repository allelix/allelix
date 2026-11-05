"""Released ADRs are append-only (ADR-0001, amendment 2.6.0).

The lock in tests/data/adr_released.json was seeded from the v2.5.2 tag.
The guard tests run the same check on a synthetic ADR directory.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import adr_lock


def test_released_adr_lines_are_unchanged() -> None:
    lock = adr_lock.load(adr_lock.LOCK)
    assert "0001-record-architecture-decisions.md" in lock["files"]
    assert "0045-clinvar-aggregate-aware-scoring.md" in lock["files"]
    assert adr_lock.INDEX not in lock["files"]
    problems = adr_lock.check(adr_lock.ADR_DIR, lock)
    assert not problems, "\n".join(problems)


_ADR = """# ADR-0002: A decision

{status}
- **Date:** 2026-01-01

## Decision

We decide this.
It has two lines.
"""


@pytest.fixture
def adr_dir(tmp_path: Path) -> Path:
    d = tmp_path / "adr"
    d.mkdir()
    (d / "0002-a-decision.md").write_text(_ADR.format(status="- **Status:** Accepted"))
    (d / "0003-another.md").write_text("# ADR-0003: Another\n\n**Status:** Accepted\n\nText.\n")
    (d / "0004-plain.md").write_text("# ADR-0004: Plain\n\nStatus: Accepted\n\nText.\n")
    (d / "README.md").write_text("# Index\n")
    return d


def _lock(adr_dir: Path) -> dict:
    return adr_lock.refresh(adr_dir, {"files": {}})


def _edit(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new, 1))


def test_lock_covers_adrs_not_the_index(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    assert sorted(lock["files"]) == ["0002-a-decision.md", "0003-another.md", "0004-plain.md"]
    assert adr_lock.check(adr_dir, lock) == []
    assert adr_lock.uncovered(adr_dir, lock) == []


def test_append_passes(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    path = adr_dir / "0002-a-decision.md"
    path.write_text(path.read_text() + "\n## Amendment 2026-02-01\n\nMore.\n")
    assert adr_lock.check(adr_dir, lock) == []
    assert adr_lock.uncovered(adr_dir, lock) == ["0002-a-decision.md: locked 9 of 13 lines"]


def test_edited_line_fails_naming_file_and_line(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / "0002-a-decision.md", "We decide this.", "We decide that.")
    assert adr_lock.check(adr_dir, lock) == [
        "0002-a-decision.md: line 8 changed: 'We decide that.'"
    ]


def test_inserted_line_fails(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / "0002-a-decision.md", "## Decision\n", "## Decision\n\nInserted.\n")
    assert adr_lock.check(adr_dir, lock) == ["0002-a-decision.md: line 8 changed: 'Inserted.'"]


def test_deleted_line_fails(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / "0002-a-decision.md", "It has two lines.\n", "")
    assert adr_lock.check(adr_dir, lock) == [
        "0002-a-decision.md: line 9 deleted (released with 9 lines)"
    ]


def test_deleted_file_fails(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    (adr_dir / "0003-another.md").unlink()
    assert adr_lock.check(adr_dir, lock) == ["0003-another.md: released ADR deleted"]


@pytest.mark.parametrize(
    ("name", "old", "new"),
    [
        (
            "0002-a-decision.md",
            "- **Status:** Accepted",
            "- **Status:** Superseded by ADR-0003",
        ),
        (
            "0002-a-decision.md",
            "- **Status:** Accepted",
            "- **Status:** Superseded in part by ADR-0003 (the rule) and ADR-0004 (the scope)",
        ),
        ("0003-another.md", "**Status:** Accepted", "**Status:** Superseded by ADR-0002"),
        ("0004-plain.md", "Status: Accepted", "Status: Superseded in part by ADR-0002"),
    ],
)
def test_status_to_superseded_by_an_existing_adr_passes(
    adr_dir: Path, name: str, old: str, new: str
) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / name, old, new)
    assert adr_lock.check(adr_dir, lock) == []


@pytest.mark.parametrize(
    "new",
    [
        "- **Status:** Accepted; amended by ADR-0003",
        "- **Status:** Deprecated",
        "- **Status:** Rejected, superseded by ADR-0003",
        "**Status:** Superseded by ADR-0003",  # a different Status form
        "- **Status:** Superseded by ADR-0099",  # names a missing ADR
        "- **Status:** Superseded by ADR-0003 and ADR-0099",  # one of them missing
    ],
)
def test_other_status_changes_fail(adr_dir: Path, new: str) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / "0002-a-decision.md", "- **Status:** Accepted", new)
    assert adr_lock.check(adr_dir, lock) == [
        "0002-a-decision.md: line 3 changed: the Status line may change only to a "
        f"Superseded-by form naming an existing ADR: {new!r}"
    ]


def test_superseded_status_does_not_excuse_another_change(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    path = adr_dir / "0002-a-decision.md"
    _edit(path, "- **Status:** Accepted", "- **Status:** Superseded by ADR-0003")
    _edit(path, "We decide this.", "We decide that.")
    assert adr_lock.check(adr_dir, lock) == [
        "0002-a-decision.md: line 8 changed: 'We decide that.'"
    ]


def test_crlf_checkout_reads_as_released(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    path = adr_dir / "0002-a-decision.md"
    path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
    assert adr_lock.check(adr_dir, lock) == []


def test_refresh_adds_new_adrs_and_grows_entries(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    path = adr_dir / "0002-a-decision.md"
    path.write_text(path.read_text() + "\n## Amendment\n")
    (adr_dir / "0005-new.md").write_text("# ADR-0005: New\n\n**Status:** Accepted\n")
    grown = adr_lock.refresh(adr_dir, lock)
    assert grown["files"]["0002-a-decision.md"]["lines"] == 11
    assert grown["files"]["0005-new.md"]["status"] == "**Status:** Accepted"
    assert adr_lock.uncovered(adr_dir, grown) == []


def test_refresh_refuses_a_changed_prefix(adr_dir: Path) -> None:
    lock = _lock(adr_dir)
    _edit(adr_dir / "0003-another.md", "Text.", "Other text.")
    with pytest.raises(ValueError, match=r"0003-another\.md: line 5 changed"):
        adr_lock.refresh(adr_dir, lock)


def test_script_refuses_and_leaves_the_lock(adr_dir: Path, tmp_path: Path) -> None:
    lock_path = tmp_path / "lock.json"
    assert adr_lock.main(["--adr-dir", str(adr_dir), "--lock", str(lock_path)]) == 0
    before = lock_path.read_text()
    assert adr_lock.main(["--adr-dir", str(adr_dir), "--lock", str(lock_path), "--covers"]) == 0
    _edit(adr_dir / "0003-another.md", "Text.", "Other text.")
    assert adr_lock.main(["--adr-dir", str(adr_dir), "--lock", str(lock_path)]) == 1
    assert adr_lock.main(["--adr-dir", str(adr_dir), "--lock", str(lock_path), "--covers"]) == 1
    assert lock_path.read_text() == before

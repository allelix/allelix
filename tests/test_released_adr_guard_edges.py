"""A released-line guard identifies the first change and requires a real ADR file."""

import pytest
from test_adr_released_lock import adr_lock


@pytest.mark.parametrize("deleted_line", [1, 3, 5, 6])
def test_deleted_line_reports_the_first_changed_position(tmp_path, deleted_line):
    path = tmp_path / "0001-decision.md"
    lines = ["# ADR-0001", "", "- **Status:** Accepted", "", "First text.", "Last text."]
    path.write_text("\n".join(lines) + "\n")
    locked = adr_lock.refresh(tmp_path, {"files": {}})
    del lines[deleted_line - 1]
    path.write_text("\n".join(lines) + "\n")
    problems = adr_lock.check(tmp_path, locked)
    assert len(problems) == 1
    assert f"{path.name}: line {deleted_line} " in problems[0], problems


@pytest.mark.parametrize("prefix", ["- **Status:**", "**Status:**", "Status:"])
@pytest.mark.parametrize("target_is_file", [False, True], ids=["directory", "file"])
def test_superseding_target_must_be_an_adr_file(tmp_path, prefix, target_is_file):
    path = tmp_path / "0001-decision.md"
    path.write_text(f"# ADR-0001\n\n{prefix} Accepted\n\nDecision.\n")
    locked = adr_lock.refresh(tmp_path, {"files": {}})
    target = tmp_path / "0002-replacement.md"
    if target_is_file:
        target.write_text("# ADR-0002\n\nStatus: Accepted\n\nReplacement.\n")
    else:
        target.mkdir()
    path.write_text(path.read_text().replace("Accepted", "Superseded by ADR-0002"))
    problems = adr_lock.check(tmp_path, locked)
    if target_is_file:
        assert problems == []
    else:
        assert len(problems) == 1
        assert f"{path.name}: line 3 " in problems[0], problems


@pytest.mark.parametrize("other_line,first_changed", [(1, 1), (5, 3)])
def test_status_error_does_not_hide_an_earlier_changed_line(tmp_path, other_line, first_changed):
    path = tmp_path / "0001-decision.md"
    lines = ["# ADR-0001", "", "- **Status:** Accepted", "", "Decision."]
    path.write_text("\n".join(lines) + "\n")
    locked = adr_lock.refresh(tmp_path, {"files": {}})
    lines[2] = "- **Status:** Rejected"
    lines[other_line - 1] += " changed"
    path.write_text("\n".join(lines) + "\n")
    problems = adr_lock.check(tmp_path, locked)
    assert len(problems) == 1
    assert f"{path.name}: line {first_changed} " in problems[0], problems

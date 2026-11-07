"""A prose label is not a CSS declaration exemption."""

import pytest
from test_no_tracker_references import references


@pytest.mark.parametrize("label", ["Issue", "Fix", "# Regression"])
def test_prose_label_does_not_hide_tracker_reference(label: str) -> None:
    token = f"#{476}"
    assert references(f"{label}: {token}") == [token]


@pytest.mark.parametrize(
    "line",
    [
        "  --border: #333;",
        "color: #0066cc;",
        "&#8217; is an apostrophe",
        "rs1234 and chr1:#",
    ],
)
def test_non_tracker_controls(line: str) -> None:
    assert references(line) == []

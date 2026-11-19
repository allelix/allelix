# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Defined ClinVar status spellings must not take the unknown-status path."""

import pytest

from allelix.reports.html import _review_stars


@pytest.mark.parametrize(
    ("status", "expected"),
    [
        ("no classification for the individual variant", "☆☆☆☆"),
        ("criteria provided, multiple submitters", "★★☆☆"),
    ],
)
@pytest.mark.parametrize("spelling", ["original", "underscored", "capitalized"])
def test_documented_status_is_not_unknown(status: str, expected: str, spelling: str) -> None:
    """Recognize the published wording in each accepted spelling form."""
    if spelling == "underscored":
        status = status.replace(" ", "_")
    elif spelling == "capitalized":
        status = status.capitalize()
    assert _review_stars(status) == expected


def test_unrecognized_status_still_has_no_inferred_stars() -> None:
    """Keep the unknown-status fallback rather than infer stars from words."""
    assert _review_stars("criteria provided, multiple unknown reviewers") == ""

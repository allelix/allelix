# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Release-gate sanitize guard for the public CHANGELOG.

GH #118: CHANGELOG.md ships verbatim from the private allelix-dev
repo to the public allelix repo and PyPI. Bare `#NN` references
auto-link on GitHub against the surrounding repo's numbering — on
the public repo they ALL 404, because the issues/PRs only exist on
the private tracker. Same hazard for an `### Internal` section,
which is forensic provenance no public user needs.

These two assertions are the structural guard the issue requires:
   - any bare issue/PR reference (#NN) fails the gate
   - any `### Internal` heading fails the gate

The guard runs on every fast-tier pytest pass, so a PR that
re-introduces either form gets caught at PR time, not at tag time.

If you genuinely need internal traceability for an entry, put the
issue number in the git commit body (the private dev branch keeps
the squash commit + the private-only release commit body — both
retain forensic context for anyone tracing history). Public users
read the CHANGELOG; everything in there must stand on its own.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

_CHANGELOG = Path(__file__).resolve().parent.parent / "CHANGELOG.md"
_ISSUE_REF_PATTERN = re.compile(r"#\d+\b")
_INTERNAL_HEADING_PATTERN = re.compile(r"^###\s+Internal\b", re.MULTILINE)


@pytest.fixture(scope="module")
def changelog_text() -> str:
    return _CHANGELOG.read_text(encoding="utf-8")


def test_changelog_has_no_issue_or_pr_references(changelog_text: str) -> None:
    """GH #118: every `#NN` token in CHANGELOG.md auto-links to the
    surrounding repo on GitHub. Private-repo issue/PR numbers don't
    exist on the public repo, so every shipped `#NN` is a dead link.

    If you need to credit an issue, restructure the prose so it
    stands without the tracker number, and put the reference in the
    git commit body (forensic context the squash commit retains).
    """
    matches = _ISSUE_REF_PATTERN.findall(changelog_text)
    assert not matches, (
        f"CHANGELOG.md contains {len(matches)} bare issue/PR reference(s) "
        f"that would become dead links on the public repo: "
        f"{sorted(set(matches))[:10]}{'…' if len(set(matches)) > 10 else ''}. "
        f"Strip them (forensic context lives in the git commit body, not "
        f"in the user-facing CHANGELOG)."
    )


def test_changelog_has_no_internal_section(changelog_text: str) -> None:
    """GH #118: an `### Internal` heading carries forensic provenance
    (architectural notes, plumbing details, fixture changes) that no
    public user needs. The squash commit on main + the private
    release commit body both retain this context for anyone tracing
    history.
    """
    matches = _INTERNAL_HEADING_PATTERN.findall(changelog_text)
    assert not matches, (
        "CHANGELOG.md contains an '### Internal' heading. Move that "
        "content into the squash commit body on main (the load-bearing "
        "forensic surface) and drop the heading from the CHANGELOG."
    )


def test_changelog_keep_a_changelog_header_intact(changelog_text: str) -> None:
    """Light sanity check that the file structure isn't accidentally
    truncated by an over-aggressive sanitize pass — the Keep a
    Changelog header should always be present at the top."""
    first_lines = "\n".join(changelog_text.splitlines()[:3])
    assert "# Changelog" in first_lines
    assert "Keep a Changelog" in changelog_text

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


# Drift-vulnerable count patterns. These match the three shapes that have
# historically shipped in release content and silently rotted across DB
# refreshes: panel-coverage ratios ("3/20 found"), timing arrows
# ("~165s → ~133s"), and explicit before/after pairs ("went from X to Y").
# Each is an absolute measurement against a specific DB snapshot — the same
# code on a refreshed DB will produce different numbers, so any reader
# (or reviewer) re-running the analysis later sees what looks like a
# regression.
#
# The patterns are deliberately narrow: they flag comparison constructs
# specifically, not bare numbers in infrastructure descriptions ("10
# non-blank lines", "57M-row cache", "999-param convention"). The
# harness output format (``✓ <key>: N annotations, M unique keys, all
# invariants hold``) deliberately does not match, because that line IS
# the drift-tolerant signal — see ``test_data/check_ground_truth.py``.
_DRIFT_PATTERNS = (
    re.compile(r"\d+/\d+ found"),  # "3/20 found"
    re.compile(r"~\d+\s*s\b.*→.*~\d+\s*s\b"),  # "~165s → ~133s"
    re.compile(r"went from \d", re.IGNORECASE),  # "went from X to Y"
)

# Caveat keywords that opt a drift-pattern line back in: if a paragraph
# contains a comparison count AND one of these keywords nearby, the
# author has explicitly tagged the measurement as a snapshot. The check
# is per-section (not per-line) so a single calibration paragraph at
# the top of an [Unreleased] block covers all measurements inside.
_SNAPSHOT_CAVEAT_KEYWORDS = (
    "snapshot",
    "drift",
    "point-in-time",
    "harness floor",
    "harness invariant",
    "ground-truth harness",
)


def _extract_unreleased_section(text: str) -> str:
    """Return the ``[Unreleased]`` block (without surrounding sections).

    The guard scopes drift-pattern enforcement to ``[Unreleased]`` only.
    Already-shipped release sections are immutable history; rewriting
    them would create churn for no reader benefit (the corresponding
    PyPI sdist and squash commit body are forever pinned regardless).
    The discipline is forward-looking — every release shipped from this
    commit onward gets the check.
    """
    marker_open = "\n## [Unreleased]"
    if marker_open not in text:
        return ""
    after = text.split(marker_open, 1)[1]
    # The block runs until the next ``## [`` heading or EOF.
    next_section = re.search(r"\n## \[", after)
    if next_section:
        return after[: next_section.start()]
    return after


def test_changelog_unreleased_has_no_drift_vulnerable_counts(
    changelog_text: str,
) -> None:
    """GH #145-era discipline: raw drift-vulnerable counts in
    ``[Unreleased]`` must carry an explicit snapshot caveat, or the
    measurement must be quoted via the §19 ground-truth harness output
    format (which is drift-tolerant by design).

    See ``CONTRIBUTING.md`` § Release-content discipline for the
    rationale and the two acceptable shapes.
    """
    unreleased = _extract_unreleased_section(changelog_text)
    if not unreleased.strip():
        # Empty [Unreleased] section is fine — nothing to validate.
        return
    drift_hits: list[str] = []
    for pattern in _DRIFT_PATTERNS:
        for match in pattern.finditer(unreleased):
            drift_hits.append(match.group(0))
    if not drift_hits:
        return
    # Any drift hit requires a snapshot caveat keyword anywhere in the
    # [Unreleased] block. One calibration paragraph at the top of the
    # block covers all measurements inside.
    lowered = unreleased.lower()
    if any(kw in lowered for kw in _SNAPSHOT_CAVEAT_KEYWORDS):
        return
    raise AssertionError(
        f"CHANGELOG.md [Unreleased] contains {len(drift_hits)} drift-vulnerable "
        f"count(s) without a snapshot caveat: {drift_hits[:5]}"
        f"{'…' if len(drift_hits) > 5 else ''}. "
        f"Either remove the raw counts, wrap them with an explicit snapshot "
        f"caveat (any of: {', '.join(_SNAPSHOT_CAVEAT_KEYWORDS)}), or quote "
        f"the §19 ground-truth harness output format instead "
        f"(`✓ <key>: N annotations, M unique keys, all invariants hold`). "
        f"See CONTRIBUTING.md § Release-content discipline."
    )

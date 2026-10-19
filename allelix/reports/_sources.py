# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Every registered source's part in a run, and report-level warnings.

The CLI builds the sources list once, after the analysis, and the
renderers only read it. A source is ``used`` (it ran and is listed in
``annotators``), ``ready_not_used`` (it was ready and given to the run,
but nothing queried it), or ``skipped`` with a reason. Failure reasons
come from the cache's state (ADR-0048) and are shown as a notice; user
choices are shown in one quiet line.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

USED = "used"
READY_NOT_USED = "ready_not_used"
SKIPPED = "skipped"

# Cache states (CacheState values) that make a source unusable.
FAILURE_REASONS: tuple[str, ...] = (
    "missing",
    "outdated",
    "unreadable",
    "corrupt",
    "wal",
    "readonly_wal",
)
# The user's own choices, in the order the CLI applies them.
CHOICE_REASONS: tuple[str, ...] = (
    "disabled",
    "license",
    "excluded",
    "flag",
    "not_applicable_build",
)

_CHOICE_TEXT = {
    "disabled": "turned off in config",
    "license": "not licensed for this use",
    "excluded": "excluded",
    "flag": "turned off by a --no- flag",
    "not_applicable_build": "not used for this build",
}


@dataclass(frozen=True)
class SourceRecord:
    """One registered source and its part in the run.

    ``message`` is the console line for a failure (its reason and fix);
    it is shown in the HTML notice and is not part of the JSON.
    """

    name: str
    display_name: str
    state: str
    reason: str | None = None
    version: str | None = None
    message: str = ""

    @property
    def failed(self) -> bool:
        """Skipped because its cache cannot be used."""
        return self.state == SKIPPED and self.reason in FAILURE_REASONS

    @property
    def choice_text(self) -> str:
        """Short wording for a source the user chose not to run."""
        return _CHOICE_TEXT.get(self.reason or "", self.reason or "")

    def to_dict(self) -> dict[str, str | None]:
        """The JSON shape: name, display_name, state, reason, version."""
        return {
            "name": self.name,
            "display_name": self.display_name,
            "state": self.state,
            "reason": self.reason,
            "version": self.version,
        }


@dataclass(frozen=True)
class ReportWarning:
    """A report-level warning: a stable code, its text, and a count when one applies."""

    code: str
    message: str
    count: int | None = None

    def to_dict(self) -> dict[str, str | int | None]:
        """The JSON shape: code, message, count."""
        return {"code": self.code, "message": self.message, "count": self.count}


def build_sources(
    registry: Iterable[tuple[str, str]],
    *,
    choices: Mapping[str, str],
    failures: Mapping[str, tuple[str, str]],
    used: Mapping[str, str | None],
    not_applicable: Iterable[str] = (),
) -> list[SourceRecord]:
    """One record per registered source, in registry order.

    ``registry`` is ``(name, display_name)`` for every registered source.
    ``choices`` maps a source the user turned off to its reason;
    ``failures`` maps an unusable source to its cache state and console
    line; ``used`` maps each source that ran to its version. A user's
    choice wins over a failure, so a source turned off and not installed
    is reported as turned off. Any other source was ready and given to
    the run but not queried.
    """
    skipped_for_build = set(not_applicable)
    records: list[SourceRecord] = []
    for name, display_name in registry:
        if name in choices:
            records.append(SourceRecord(name, display_name, SKIPPED, choices[name]))
        elif name in failures:
            state, message = failures[name]
            records.append(SourceRecord(name, display_name, SKIPPED, state, message=message))
        elif name in used:
            records.append(SourceRecord(name, display_name, USED, version=used[name]))
        elif name in skipped_for_build:
            records.append(SourceRecord(name, display_name, SKIPPED, "not_applicable_build"))
        else:
            records.append(SourceRecord(name, display_name, READY_NOT_USED))
    return records


def clinvar_skipped_warning(sources: Iterable[SourceRecord]) -> ReportWarning | None:
    """The warning for a run in which ClinVar did not run, whatever the reason."""
    for record in sources:
        if record.name != "clinvar" or record.state != SKIPPED:
            continue
        why = record.message or f"ClinVar: {record.choice_text}."
        return ReportWarning(
            "clinvar_skipped",
            f"ClinVar was not used, so this report has no ClinVar findings. {why}",
        )
    return None

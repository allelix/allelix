# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Database and input text for a Rich markup string, printed exactly as written (2.6.2)."""

from __future__ import annotations

import re

# rich.markup.escape's own pattern: a "[" that would open or close a tag.
_TAG = re.compile(r"(\\*)(\[[a-z#/@][^[]*?])")
# A tag pair with no style: it ends a run of backslashes, which Rich halves
# before any tag, and leaves the style of what follows unchanged.
_END = "[none][/none]"


def literal(text: str) -> str:
    r"""``text`` for a Rich markup string, shown exactly as written wherever it is put.

    A ``[`` that would start a tag is escaped, as ``rich.markup.escape`` does.
    Rich also halves a run of backslashes right before a tag, so whether a
    trailing ``\`` survives depends on what follows: a closing tag, an
    optional suffix that may be empty, plain text or nothing.
    ``rich.markup.escape`` doubles one trailing backslash, which shows two
    where no tag follows and loses some of a longer run where one does. Here
    a trailing run is doubled and closed by an empty tag pair, so it always
    stands before a tag of its own and prints as written in every place.
    """
    escaped = _TAG.sub(lambda m: f"{m[1]}{m[1]}\\{m[2]}", text)
    run = len(escaped) - len(escaped.rstrip("\\"))
    return f"{escaped}{'\\' * run}{_END}" if run else escaped

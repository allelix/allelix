# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Controls for the D/I rule: what the fix must keep.

- A call written as a literal sequence (an insertion or deletion with its
  bases, such as ``C/CT``) is not a D/I call. It still matches a SNPedia
  page written in the same sequence, on the single-call and batch paths,
  on both builds. The rule names the letters ``D`` and ``I``; it does not
  exclude every allele that is not one base.
- A D/I call at a high-value rsID stays visible as not checked in the
  terminal and the HTML report, as in JSON, and has no source finding there.

The fixture is the one in ``test_di_calls_not_interpreted.py``. All rows
are synthetic.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from test_di_calls_not_interpreted import (
    _CHIP_CALLS,
    _FLOORS,
    _POS37,
    MIXED,
    _array,
    _both_paths,
    _snpedia,
    data,  # noqa: F401 - pytest fixture
    identity,  # noqa: F401 - pytest fixture
)

from allelix.cli import main
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path

# (allele1, allele2, summary): a plus page written in sequences, at MIXED.
_SEQUENCE_ROWS = (
    ("C", "CT", "sequence page C;CT"),
    ("CT", "CT", "sequence page CT;CT"),
)


def _add_sequence_rows(data_dir: Path) -> None:
    with contextlib.closing(sqlite3.connect(data_dir / "snpedia.sqlite")) as conn:
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, 3.0, 'Bad', ?, 'GENE', 'plus', "
            "'2026-10-10T00:00:00')",
            [(MIXED, a1, a2, summary) for a1, a2, summary in _SEQUENCE_ROWS],
        )
        conn.commit()


def _sequence_call(a1: str, a2: str, build: str) -> Variant:
    chrom, pos = _POS37[MIXED]
    if build == "GRCh38":
        pos += 100000
    return Variant(
        MIXED,
        chrom,
        pos,
        a1,
        a2,
        build=build,
        ref="C",
        input_alt="CT",
        ref_from_input=True,
    )


@pytest.mark.parametrize("with_identity", (True, False))
@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize(("a1", "a2", "summary"), _SEQUENCE_ROWS)
def test_a_sequence_call_still_matches_its_sequence_page(
    data,  # noqa: F811
    identity,  # noqa: F811
    a1,
    a2,
    summary,
    build,
    with_identity,
):
    _add_sequence_rows(data)
    single, batch = _both_paths(
        _snpedia(data, identity if with_identity else None), _sequence_call(a1, a2, build)
    )
    assert [a.description for a in single] == [a.description for a in batch]
    assert len(single) == 1
    assert summary in single[0].description


_NAMED = (
    "BRCA1 5382insC (rs80357906): called D/I, a deletion or insertion code that Allelix "
    "cannot check against its sources."
)


def _terminal_and_html(tmp_path: Path, data_dir: Path) -> tuple[str, str]:
    html = tmp_path / "r.html"
    result = CliRunner(env={"COLUMNS": "300"}).invoke(
        main,
        [
            "analyze",
            str(_array(tmp_path / "a.txt", _CHIP_CALLS)),
            "--build",
            "grch37",
            "--no-update",
            "--data-dir",
            str(data_dir),
            "--output",
            str(html),
            *_FLOORS,
        ],
    )
    assert result.exit_code == 0, result.output
    return " ".join(result.output.split()), " ".join(html.read_text().split())


def test_a_high_value_di_call_is_shown_not_checked_in_terminal_and_html(
    tmp_path,
    data,  # noqa: F811
):
    terminal, page = _terminal_and_html(tmp_path, data)
    assert _NAMED in terminal
    assert _NAMED in page


def test_a_high_value_di_call_has_no_source_finding_in_terminal_and_html(
    tmp_path,
    data,  # noqa: F811
):
    """Behavioural, beside the JSON check in the main file: fails until the fix."""
    terminal, page = _terminal_and_html(tmp_path, data)
    assert "high-value page DI" not in terminal
    assert "high-value page DI" not in page

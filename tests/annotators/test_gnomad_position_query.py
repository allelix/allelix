# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""gnomAD position lookup: one chromosome per statement, read-only safe.

The fixture is the loader's schema without ``idx_gnomad_position``, as on
a cache built before that index existed: the primary key ``(chrom, pos,
ref, alt)`` must serve the lookup on its own. All rows are synthetic.
"""

from __future__ import annotations

import contextlib
import os
import re
import sqlite3
import stat
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.gnomad import GNOMAD_DB_FILENAME, GnomadAnnotator
from allelix.databases.schema import GNOMAD_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

_SCHEMA = re.sub(
    r"CREATE INDEX IF NOT EXISTS idx_gnomad_position\s+ON gnomad_frequencies\(chrom, pos\);",
    "",
    GNOMAD_SCHEMA,
)

# (chrom, pos, ref, alt, rsid): two alleles at one position, a row with no
# rsID, and 2,500 positions on chromosome 1 to fill several statements.
_ROWS = [
    ("1", 100, "A", "G", "rs900780001"),
    ("1", 100, "A", "T", "rs900780002"),
    ("2", 200, "C", "T", None),
    ("X", 300, "G", "A", "rs900780003"),
    *(("1", 1_000 + i, "C", "G", f"rs9008{i:05d}") for i in range(2_500)),
]


def _cache(data_dir: Path) -> Path:
    assert "idx_gnomad_position" not in _SCHEMA
    db = data_dir / GNOMAD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.executescript(_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, 0.1)",
            _ROWS,
        )
        conn.commit()
        assert conn.execute("PRAGMA journal_mode").fetchone()[0] == "delete"
    return db


def _statements(gnomad: GnomadAnnotator, positions: set[tuple[str, int]]) -> list[str]:
    seen: list[str] = []
    gnomad._connection().set_trace_callback(seen.append)
    gnomad.bulk_resolve_rsids_from_positions(positions)
    gnomad._connection().set_trace_callback(None)
    return [s for s in seen if s.lstrip().upper().startswith("SELECT")]


_POSITIONS = {(c, p) for c, p, *_ in _ROWS} | {("2", 999), ("22", 5)}


class TestPositionQuery:
    def test_returns_every_row_with_an_rsid_and_nothing_else(self, tmp_path: Path) -> None:
        _cache(tmp_path)
        with GnomadAnnotator(tmp_path) as gnomad:
            found = gnomad.bulk_resolve_rsids_from_positions(_POSITIONS)
        assert found[("1", 100)] == [("A", "G", "rs900780001"), ("A", "T", "rs900780002")]
        assert found[("X", 300)] == [("G", "A", "rs900780003")]
        assert ("2", 200) not in found
        assert len(found) == 2 + 2_500

    def test_statement_names_one_chromosome(self, tmp_path: Path) -> None:
        """Shape, whatever SQLite plans: ``chrom = ? AND pos IN (...)``,
        never the row-value ``(chrom, pos) IN`` form."""
        _cache(tmp_path)
        with GnomadAnnotator(tmp_path) as gnomad:
            statements = _statements(gnomad, _POSITIONS)
        assert statements
        assert all("chrom = '" in s and "AND pos IN (" in s for s in statements)
        assert not any("(chrom, pos) IN" in s for s in statements)

    def test_statement_count_is_one_per_chromosome_per_chunk(self, tmp_path: Path) -> None:
        # chromosome 1: 2,501 positions (3 statements of up to 998);
        # 2, X and 22: one each.
        _cache(tmp_path)
        with GnomadAnnotator(tmp_path) as gnomad:
            assert len(_statements(gnomad, _POSITIONS)) == 3 + 3

    def test_plan_searches_on_chrom_and_pos(self, tmp_path: Path) -> None:
        _cache(tmp_path)
        with GnomadAnnotator(tmp_path) as gnomad:
            statements = _statements(gnomad, {("1", 100), ("1", 1_000)})
            conn = gnomad._connection()
            plans = [
                " ".join(row[3] for row in conn.execute(f"EXPLAIN QUERY PLAN {s}"))
                for s in statements
            ]
        print(f"SQLite {sqlite3.sqlite_version}: {plans}")
        assert plans
        assert all("(chrom=? AND pos=?)" in plan for plan in plans)


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root ignores file modes")
class TestReadOnlyCache:
    def test_resolver_runs_on_a_read_only_cache(self, tmp_path: Path) -> None:
        """No write on connect: a rollback-journal cache with file mode
        0444 and no ``idx_gnomad_position`` is read as it is."""
        db = _cache(tmp_path)
        db.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        try:
            with GnomadAnnotator(tmp_path) as gnomad:
                found = gnomad.bulk_resolve_rsids_from_positions({("X", 300)})
            assert found == {("X", 300): [("G", "A", "rs900780003")]}
            with contextlib.closing(sqlite3.connect(f"file:{db}?mode=ro", uri=True)) as conn:
                names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
            assert "idx_gnomad_position" not in names
        finally:
            db.chmod(stat.S_IRUSR | stat.S_IWUSR)

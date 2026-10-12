# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""AlphaMissense and CADD caches are ready only with the current schema stamp."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.databases._versions import ALPHAMISSENSE_SCHEMA_VERSION, CADD_SCHEMA_VERSION
from allelix.databases.alphamissense_loader import ALPHAMISSENSE_DB_FILENAME
from allelix.databases.cadd_loader import CADD_DB_FILENAME
from allelix.databases.schema import ALPHAMISSENSE_SCHEMA, CADD_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

SOURCES = {
    "alphamissense": (
        AlphaMissenseAnnotator,
        ALPHAMISSENSE_DB_FILENAME,
        ALPHAMISSENSE_SCHEMA,
        ALPHAMISSENSE_SCHEMA_VERSION,
    ),
    "cadd": (CaddAnnotator, CADD_DB_FILENAME, CADD_SCHEMA, CADD_SCHEMA_VERSION),
}


@pytest.mark.parametrize("source", sorted(SOURCES))
def test_ready_only_with_the_current_stamp(tmp_path: Path, source: str) -> None:
    """No stamp (schema version unknown) or an old stamp is not ready."""
    annotator_cls, filename, schema, version = SOURCES[source]
    db_path = tmp_path / filename

    def stamped(tag: str | None) -> bool:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.executescript(schema)
            conn.execute(
                "INSERT OR REPLACE INTO database_versions (name, source_url, version, "
                "downloaded_at, record_count, local_version_tag) VALUES (?, ?, ?, ?, ?, ?)",
                (source, "test://cache", "1", "2026-01-01T00:00:00Z", 1, tag),
            )
            conn.commit()
        with annotator_cls(tmp_path) as annotator:
            return annotator.is_ready()

    assert stamped(None) is False
    assert stamped(f"sv:{version - 1}") is False
    assert stamped(f"sv:{version}") is True

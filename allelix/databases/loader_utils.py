# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared utilities for pre-built SQLite cache loaders.

Extracted from gnomad_loader, alphamissense_loader, and snpedia_loader
which all follow the same download-decompress-stamp pattern for
HuggingFace-hosted SQLite caches.
"""

from __future__ import annotations

import contextlib
import gzip
import os
import shutil
import sqlite3
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

DISK_SPACE_MULTIPLIER = 6


def install_prebuilt_gz_cache(
    gz_path: Path,
    db_path: Path,
    record_name: str,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
    schema_version_tag: str | None = None,
) -> None:
    """Decompress a gzipped pre-built SQLite cache into place.

    Args:
        gz_path: Path to the downloaded .gz file.
        db_path: Destination path for the decompressed SQLite database.
        record_name: Database record name for version stamping.
        source_url: URL the file was downloaded from.
        remote_signal: Freshness signal to stamp in database_versions.
        schema_version_tag: If set, stamped as local_version_tag
            (e.g. ``"sv:2"``). Omit for caches that don't track
            schema versions (SNPedia).
    """
    gz_size = gz_path.stat().st_size
    free = shutil.disk_usage(db_path.parent).free
    needed = gz_size * DISK_SPACE_MULTIPLIER
    if free < needed:
        raise OSError(
            f"Not enough disk space to decompress {gz_path.name}: "
            f"{free / 1e9:.1f} GB free, need ~{needed / 1e9:.1f} GB. "
            "Free up space and retry."
        )

    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    with gzip.open(gz_path, "rb") as f_in, tmp_path.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
        if remote_signal:
            from allelix.databases.manager import stamp_remote_signal

            stamp_remote_signal(conn, record_name, remote_signal, source_url)

        if schema_version_tag:
            from allelix.databases.manager import _ensure_local_version_tag_column

            _ensure_local_version_tag_column(conn)
            conn.execute(
                "UPDATE database_versions SET local_version_tag = ? WHERE name = ?",
                (schema_version_tag, record_name),
            )

        conn.commit()

    os.replace(tmp_path, db_path)

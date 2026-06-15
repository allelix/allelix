# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""AlphaMissense pathogenicity cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update``. Contains 71M missense variant scores from AlphaMissense
with genomic coordinates and rsIDs (joined from gnomAD at build time).

The cache can also be built locally from the Zenodo TSV via
``scripts/build_alphamissense_cache.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases._versions import ALPHAMISSENSE_SCHEMA_VERSION
from allelix.databases.loader_utils import install_prebuilt_gz_cache

if TYPE_CHECKING:
    from pathlib import Path

ALPHAMISSENSE_DB_FILENAME = "alphamissense.sqlite"

ALPHAMISSENSE_CACHE_URL = (
    "https://huggingface.co/datasets/dial481/allelix-alphamissense"
    "/resolve/13a15e199536512b5e2d208d79c4f93c0a73f71f/alphamissense.sqlite.gz"
)

ALPHAMISSENSE_EXPECTED_SHA256 = "0cc1049d59b0aca61f397ad0650516555a271acffa65b7b8f23899bbd11c4386"


def install_prebuilt_cache(
    gz_path: Path,
    db_path: Path,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
) -> None:
    """Decompress a gzipped pre-built SQLite cache into place."""
    install_prebuilt_gz_cache(
        gz_path,
        db_path,
        "alphamissense",
        source_url=source_url,
        remote_signal=remote_signal,
        schema_version_tag=f"sv:{ALPHAMISSENSE_SCHEMA_VERSION}",
    )

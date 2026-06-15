# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""SNPedia pre-built cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update``. Contains ~216K raw wiki pages and ~105K parsed genotype
rows.

The cache can also be built locally via ``scripts/scrape_snpedia.py``
followed by ``scripts/parse_snpedia.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases.loader_utils import install_prebuilt_gz_cache

if TYPE_CHECKING:
    from pathlib import Path

SNPEDIA_CACHE_URL = (
    "https://huggingface.co/datasets/genomics-commons/snpedia"
    "/resolve/69a745401a0d63acb71fc759b9e79f6d5da79dd9/snpedia.sqlite.gz"
)

SNPEDIA_EXPECTED_SHA256 = "bd940b624143d03427baf9b2572da07257631bd6fb8b584b5ed0961f07cad104"


def install_prebuilt_cache(
    gz_path: Path,
    db_path: Path,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
) -> None:
    """Decompress a gzipped pre-built SNPedia SQLite cache into place."""
    install_prebuilt_gz_cache(
        gz_path,
        db_path,
        "snpedia",
        source_url=source_url,
        remote_signal=remote_signal,
    )

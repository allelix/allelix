# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""CADD variant deleteriousness cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update`` when CADD is enabled. Contains CADD PHRED scores for
SNV and indel positions present in Allelix's existing databases
(gnomAD, AlphaMissense, ClinVar GRCh38).

The cache can also be built locally from the full CADD prescored files
via ``scripts/build_cadd_cache.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases._versions import CADD_SCHEMA_VERSION
from allelix.databases.loader_utils import install_prebuilt_gz_cache

if TYPE_CHECKING:
    from pathlib import Path

CADD_DB_FILENAME = "cadd.sqlite"

CADD_CACHE_URL = (
    "https://huggingface.co/datasets/genomics-commons/cadd-scores"
    "/resolve/3157a4b5f65876eae2dbd7376dacfd15e2e2b08e/cadd.sqlite.gz"
)

CADD_EXPECTED_SHA256 = "ad8393f59eabd026c69ce2ef227b0d95eb9b428cd2960dc5fbddf142ca328b28"


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
        "cadd",
        source_url=source_url,
        remote_signal=remote_signal,
        schema_version_tag=f"sv:{CADD_SCHEMA_VERSION}",
    )

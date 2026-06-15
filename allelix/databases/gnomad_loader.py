# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""gnomAD exome frequency cache loader.

The pre-built SQLite cache is downloaded from HuggingFace during
``db update``. Contains all ~16M exome rsIDs from gnomAD v4.1 with
genomic coordinates (chrom/pos/ref/alt) for future AlphaMissense/CADD
integration.

The cache can also be built locally from gnomAD exome VCFs via
``scripts/build_gnomad_cache.py`` (streaming or local file mode).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.loader_utils import install_prebuilt_gz_cache

if TYPE_CHECKING:
    from pathlib import Path

GNOMAD_DB_FILENAME = "gnomad.sqlite"

GNOMAD_CACHE_URL = (
    "https://huggingface.co/datasets/allelix/allelix-gnomad"
    "/resolve/f0aadfb7940290c44930dc0d1b9b093bc089173f/gnomad.sqlite.gz"
)

GNOMAD_EXPECTED_SHA256 = "e001b6c472b89075f18c82a34ccfb1e8e5c524f8502b988db1a546d25b0c6fe4"


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
        "gnomad",
        source_url=source_url,
        remote_signal=remote_signal,
        schema_version_tag=f"sv:{GNOMAD_SCHEMA_VERSION}",
    )

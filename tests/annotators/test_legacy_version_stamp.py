# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GWAS and ClinPGx caches from before ``local_version_tag`` are stamped only
when the version baked into ``remote_signal`` is the current one."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.databases._versions import PHARMGKB_INTERPRETER_VERSION
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION, load_gwas_tsv
from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize("baked", ["current", "older", "none"])
@pytest.mark.parametrize("source", ["gwas", "pharmgkb"])
def test_legacy_stamp_only_for_the_current_baked_version(
    tmp_path: Path, mock_gwas_tsv: Path, mock_pharmgkb_dir: Path, source: str, baked: str
) -> None:
    """No cached TSV or ZIP is on disk, so a cache that cannot be stamped is not ready."""
    kind, current, annotator_cls = {
        "gwas": ("cv", _CATEGORIZER_VERSION, GWASCatalogAnnotator),
        "pharmgkb": ("iv", PHARMGKB_INTERPRETER_VERSION, PharmGKBAnnotator),
    }[source]
    db_path = tmp_path / f"{source}.sqlite"
    if source == "gwas":
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://legacy")
    else:
        load_pharmgkb_tsv(mock_pharmgkb_dir, db_path, remote_signal="etag:abc")
    suffix = {"current": f"|{kind}:{current}", "older": f"|{kind}:{current - 1}", "none": ""}
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.execute(
            "UPDATE database_versions SET local_version_tag = NULL, remote_signal = ? "
            "WHERE name = ?",
            (f"etag:abc{suffix[baked]}", source),
        )
        conn.commit()

    with annotator_cls(tmp_path) as annotator:
        assert annotator.is_ready() is (baked == "current")
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        tag, signal = conn.execute(
            "SELECT local_version_tag, remote_signal FROM database_versions WHERE name = ?",
            (source,),
        ).fetchone()
    if baked == "current":
        assert (tag, signal) == (f"{kind}:{current}", "etag:abc")
    else:
        assert tag is None

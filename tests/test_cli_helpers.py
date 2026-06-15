# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Unit tests for private CLI helpers."""

from __future__ import annotations

import os
import time
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

from allelix.cli._helpers import (
    _STALENESS_SECONDS,
    _chrom_sort_key,
    _maybe_refresh_databases,
    _percent,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestChromSortKey:
    def test_autosomes_then_sex_then_mt(self):
        keys = sorted(["MT", "X", "1", "22", "Y", "2"], key=_chrom_sort_key)
        assert keys == ["1", "2", "22", "X", "Y", "MT"]

    def test_unknown_chrom_falls_to_alphabetical(self):
        keys = sorted(["X", "1", "WEIRD", "2"], key=_chrom_sort_key)
        assert keys == ["1", "2", "X", "WEIRD"]


class TestPercent:
    def test_zero_total_does_not_divide(self):
        assert _percent(0, 0) == "0.00%"

    def test_basic(self):
        assert _percent(1, 2) == "50.00%"


def _mock_annotator(
    name: str = "clinvar",
    display_name: str = "ClinVar",
    requires_download: bool = True,
    is_ready: bool = True,
    fetch_remote: str | None = "md5:abc123",
    cached_remote: str | None = "md5:old999",
) -> MagicMock:
    ann = MagicMock()
    ann.name = name
    ann.display_name = display_name
    ann.requires_download = requires_download
    ann.is_ready.return_value = is_ready
    ann.fetch_remote_signal.return_value = fetch_remote
    ann.cached_remote_signal.return_value = cached_remote
    ann.version.return_value = "20260601"
    ann.__enter__ = MagicMock(return_value=ann)
    ann.__exit__ = MagicMock(return_value=False)
    return ann


class TestMaybeRefreshDatabases:
    """Tests for the pre-analysis database freshness gate."""

    def _write_stale_db(self, data_dir: Path, name: str = "clinvar") -> Path:
        db = data_dir / f"{name}.GRCh37.sqlite"
        db.write_text("")
        stale_time = time.time() - _STALENESS_SECONDS - 3600
        os.utime(db, (stale_time, stale_time))
        return db

    def _write_fresh_db(self, data_dir: Path, name: str = "clinvar") -> Path:
        db = data_dir / f"{name}.GRCh37.sqlite"
        db.write_text("")
        return db

    def test_stale_db_with_changed_remote_triggers_setup(self, tmp_path: Path) -> None:
        self._write_stale_db(tmp_path)
        ann = _mock_annotator(fetch_remote="md5:new", cached_remote="md5:old")
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        mock_setup.assert_called_once_with(ann)

    def test_fresh_db_skips_remote_check(self, tmp_path: Path) -> None:
        self._write_fresh_db(tmp_path)
        ann = _mock_annotator()
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        ann.fetch_remote_signal.assert_not_called()
        mock_setup.assert_not_called()

    def test_stale_db_same_remote_signal_skips_setup(self, tmp_path: Path) -> None:
        self._write_stale_db(tmp_path)
        ann = _mock_annotator(fetch_remote="md5:same", cached_remote="md5:same")
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        mock_setup.assert_not_called()

    def test_stale_db_offline_warns_no_setup(self, tmp_path: Path) -> None:
        self._write_stale_db(tmp_path)
        ann = _mock_annotator(fetch_remote=None)
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        mock_setup.assert_not_called()

    def test_annotator_not_ready_skips(self, tmp_path: Path) -> None:
        self._write_stale_db(tmp_path)
        ann = _mock_annotator(is_ready=False)
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        ann.fetch_remote_signal.assert_not_called()
        mock_setup.assert_not_called()

    def test_snpedia_excluded_via_requires_download(self, tmp_path: Path) -> None:
        self._write_stale_db(tmp_path, name="snpedia")
        ann = _mock_annotator(name="snpedia", requires_download=False)
        with (
            patch("allelix.cli._helpers.get_annotators", return_value=[ann]),
            patch("allelix.cli._helpers._run_setup") as mock_setup,
        ):
            _maybe_refresh_databases(tmp_path)
        ann.fetch_remote_signal.assert_not_called()
        mock_setup.assert_not_called()

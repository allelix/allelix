# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""R-4: ClinVar CLNSIG vocabulary drift detection.

Two tiers of coverage:

1. **Code consistency** (default CI): every CLNSIG value the annotator
   references in ``_CLNSIG_MAGNITUDE`` and ``_BENIGN_CLNSIGS`` is part
   of the canonical ClinVar vocabulary as captured in the snapshot
   file. Catches "we added a scoring rule for a CLNSIG ClinVar doesn't
   publish" — the kind of bug a refactor or test-data-driven change
   could introduce.

2. **Live drift** (``pytest -m slow``): every distinct
   ``clinical_significance`` value present in the installed ClinVar
   SQLite cache is in the snapshot. Fails when ClinVar ships a new
   CLNSIG term. Skipped when the cache isn't installed (CI runners
   don't have one).

Snapshot lives at ``allelix/data/clinvar_clnsig_snapshot.yaml``. When
it drifts, the remediation is documented in the YAML header.
"""

from __future__ import annotations

import contextlib
import sqlite3
from pathlib import Path

import pytest
import yaml

from allelix.annotators.clinvar import (
    _BENIGN_CLNSIGS,
    _CLNSIG_MAGNITUDE,
    _normalize_clnsig,
    clinvar_db_filename,
)
from allelix.databases import resolve_data_dir

_SNAPSHOT_PATH = Path("allelix/data/clinvar_clnsig_snapshot.yaml")


def _load_snapshot() -> set[str]:
    """Return the normalized vocabulary set from the snapshot file."""
    with open(_SNAPSHOT_PATH, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return {_normalize_clnsig(v) for v in data["values"]}


class TestCodeVocabularyMatchesSnapshot:
    """Code-side scoring rules must reference only snapshot vocabulary.

    Catches the case where someone adds a key to ``_CLNSIG_MAGNITUDE``
    that doesn't correspond to a real ClinVar term. Pure code check —
    runs in default CI without any external data.
    """

    def test_snapshot_file_exists_and_parses(self):
        """Pin the snapshot's location and structure."""
        assert _SNAPSHOT_PATH.exists(), f"Snapshot missing at {_SNAPSHOT_PATH}"
        snapshot = _load_snapshot()
        assert snapshot, "Snapshot is empty"
        assert "pathogenic" in snapshot, "Snapshot missing baseline 'pathogenic'"

    def test_every_clnsig_magnitude_key_is_in_snapshot(self):
        """Every scoring rule must reference a real ClinVar term."""
        snapshot = _load_snapshot()
        unknown = set(_CLNSIG_MAGNITUDE.keys()) - snapshot
        assert not unknown, (
            f"_CLNSIG_MAGNITUDE has keys not in snapshot: {sorted(unknown)}. "
            f"Either the term is real ClinVar vocabulary (add it to "
            f"{_SNAPSHOT_PATH}) or the scoring rule references a typo or "
            f"removed term (fix the code)."
        )

    def test_every_benign_clnsig_is_in_snapshot(self):
        """Benign suppression set must reference real ClinVar terms."""
        snapshot = _load_snapshot()
        unknown = _BENIGN_CLNSIGS - snapshot
        assert not unknown, (
            f"_BENIGN_CLNSIGS contains values not in snapshot: {sorted(unknown)}. "
            f"Either add the term to the snapshot or fix the code."
        )

    def test_normalized_snapshot_values_round_trip(self):
        """The snapshot stores raw ClinVar values; normalization is idempotent.

        Loading and normalizing again must produce the same set — pins
        the contract that ``_normalize_clnsig`` is a stable function
        callers can apply once.
        """
        snapshot = _load_snapshot()
        re_normalized = {_normalize_clnsig(v) for v in snapshot}
        assert snapshot == re_normalized


def _existing_clinvar_caches() -> list[Path]:
    """Return paths to installed ClinVar SQLite caches, if any."""
    data_dir = resolve_data_dir(None)
    candidates = [data_dir / clinvar_db_filename(b) for b in ("GRCh37", "GRCh38")]
    return [p for p in candidates if p.exists()]


@pytest.mark.slow
class TestLiveClnsigDrift:
    """Drift check against an installed ClinVar SQLite cache.

    Skipped when no cache is present (``allelix db update`` not run).
    Otherwise queries every distinct ``clinical_significance`` value,
    normalizes, and asserts it's in the snapshot.

    Fires on a fresh ClinVar release that introduces a new term — the
    sentinel that tells the dev to update the scoring rules and the
    snapshot.
    """

    def test_no_new_clnsig_values_in_live_cache(self):
        caches = _existing_clinvar_caches()
        if not caches:
            pytest.skip("No ClinVar cache installed; run `allelix db update` first")

        snapshot = _load_snapshot()
        unknown: set[str] = set()
        for cache_path in caches:
            with contextlib.closing(sqlite3.connect(cache_path)) as conn:
                rows = conn.execute(
                    "SELECT DISTINCT clinical_significance FROM clinvar_variants"
                ).fetchall()
            for (raw,) in rows:
                if raw is None:
                    continue
                # CLNSIG INFO may carry multiple pipe-separated values per row
                # (ADR-0011 — one per ALT). Split + normalize each.
                for value in raw.split("|"):
                    normalized = _normalize_clnsig(value)
                    if normalized not in snapshot:
                        unknown.add(normalized)

        assert not unknown, (
            f"ClinVar cache contains CLNSIG values not in snapshot:\n"
            f"  {sorted(unknown)}\n"
            f"Update {_SNAPSHOT_PATH} after verifying these are real "
            f"ClinVar terms (not malformed rows), and decide whether they "
            f"need scoring rules in _CLNSIG_MAGNITUDE."
        )

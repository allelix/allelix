# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar CLNSIG vocabulary drift detection.

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
import os
import sqlite3
from pathlib import Path

import pytest
import yaml

from allelix.annotators.clinvar import (
    _AGGREGATE_CAP,
    _BENIGN_CLNSIGS,
    _CLNSIG_MAGNITUDE,
    _normalize_clnsig,
    clinvar_db_filename,
)
from allelix.databases import resolve_data_dir
from allelix.databases.manager import _CLINVAR_PLACEHOLDER_CLNSIGS, aggregate_germline_term

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

    def test_every_snapshot_term_has_a_scoring_decision(self):
        """Lockstep: each known term is scored, filtered at load, or
        documented as appearing only in ClinVar's aggregate.

        A term added to the snapshot without a decision fails here. A term
        the loader meets that is in none of these is not shown, and the run
        counts it (it no longer scores 5.0).
        """
        snapshot = _load_snapshot()
        placeholders = {_normalize_clnsig(v) for v in _CLINVAR_PLACEHOLDER_CLNSIGS}
        undecided = snapshot - set(_CLNSIG_MAGNITUDE) - placeholders - _AGGREGATE_ONLY_TERMS
        assert not undecided, (
            f"Snapshot terms with no scoring decision: {sorted(undecided)}. Score them "
            f"in _CLNSIG_MAGNITUDE, filter them in _CLINVAR_PLACEHOLDER_CLNSIGS, or add "
            f"them to _AGGREGATE_ONLY_TERMS if no single submission carries them."
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


# Snapshot terms ClinVar uses only as a variant's aggregate, combining
# several submissions' classifications; no single submission carries them.
# Should one ever reach a cache row, it is hidden and counted.
_AGGREGATE_ONLY_TERMS: frozenset[str] = frozenset(
    _normalize_clnsig(v)
    for v in (
        "",
        "Pathogenic/Likely_pathogenic,_low_penetrance",
        "Pathogenic/Likely_pathogenic/established_risk_allele",
        "Pathogenic/Likely_pathogenic/Likely_pathogenic,_low_penetrance",
        "Pathogenic/Likely_pathogenic/likely_risk_allele",
        "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance",
        "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance/established_risk_allele",
        "Pathogenic/likely_risk_allele",
        "Pathogenic/Pathogenic,_low_penetrance",
        "Likely_pathogenic/Likely_pathogenic,_low_penetrance",
        "Likely_pathogenic/likely_risk_allele",
        "Likely_pathogenic/Pathogenic,_low_penetrance",
        "Uncertain_significance/uncertain_risk_allele",
        "Uncertain_significance/vus-high",
        "Uncertain_significance/vus-low",
        "Uncertain_significance/vus-mid",
        "no_classifications_from_unflagged_records",
    )
)


def _existing_clinvar_caches() -> list[Path]:
    """Return paths to installed ClinVar SQLite caches, if any."""
    data_dir = resolve_data_dir(None)
    candidates = [data_dir / clinvar_db_filename(b) for b in ("GRCh37", "GRCh38")]
    return [p for p in candidates if p.exists()]


@pytest.mark.slow
class TestClinvarClnsigDrift:
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
            if os.environ.get("ALLELIX_CI_REQUIRE_CACHE") == "1":
                pytest.fail("Required ClinVar cache is absent under ALLELIX_CI_REQUIRE_CACHE=1")
            pytest.skip("No ClinVar cache installed; run `allelix db update` first")

        snapshot = _load_snapshot()
        unknown: set[str] = set()
        uncapped: set[str] = set()
        for cache_path in caches:
            with contextlib.closing(sqlite3.connect(cache_path)) as conn:
                rows = conn.execute(
                    "SELECT DISTINCT clinical_significance FROM clinvar_variants"
                ).fetchall()
                aggregates = conn.execute(
                    "SELECT DISTINCT aggregate FROM clinvar_variants"
                ).fetchall()
            # ADR-0045: the cap keys on ClinVar's aggregate wording. A renamed
            # Conflicting or Benign aggregate would silently stop the cap.
            for (aggregate,) in aggregates:
                term = aggregate_germline_term(aggregate or "")
                if ("conflicting" in term or "benign" in term) and term not in _AGGREGATE_CAP:
                    uncapped.add(term)
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
        assert not uncapped, (
            f"ClinVar aggregates that look Conflicting or Benign but are not capped: "
            f"{sorted(uncapped)}. Add them to _AGGREGATE_CAP (ADR-0045) or explain why not."
        )


def test_ci_cache_requirement_converts_missing_cache_to_failure(monkeypatch):
    """The workflow flag must turn the local convenience skip into a hard failure."""
    monkeypatch.setenv("ALLELIX_CI_REQUIRE_CACHE", "1")
    monkeypatch.setitem(globals(), "_existing_clinvar_caches", lambda: [])

    with pytest.raises(pytest.fail.Exception, match="Required ClinVar cache is absent"):
        TestClinvarClnsigDrift().test_no_new_clnsig_values_in_live_cache()

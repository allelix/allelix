# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Bound metadata I/O per ClinVar instance/build without caching wrong answers.

Use real synthetic SQLite caches; count metadata reads, not wall time. Distinct
loci prevent the result cache from hiding per-call interpreter-tag reads.
Close/reuse and legacy-tag preparation must retain their existing semantics.
"""

import contextlib
import sqlite3
from collections import Counter

import pytest
from test_homref_ref_by_locus import _cv, _write_clinvar

from allelix.annotators import clinvar as module
from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

BUILDS = ("GRCh37", "GRCh38")
CURRENT = f"iv:{CLINVAR_INTERPRETER_VERSION}"
OLD = "iv:0"
N = 8


def _tag(data, build, tag, *, baked=False):
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.execute(
            "UPDATE database_versions SET local_version_tag=?, remote_signal=?",
            (tag, f"md5:synthetic|{CURRENT}" if baked else "md5:synthetic"),
        )
        conn.commit()


def _pair(data, build, *, carried):
    other = next(b for b in BUILDS if b != build)
    own, remote = [], []
    for index in range(N):
        rsid = f"rs399949{index:04d}"
        own.append(
            _cv(
                rsid,
                "1",
                1000 + index,
                "T" if carried else "A",
                "T" if carried else "G",
                allele_id=100 + index,
            )
        )
        remote.append(_cv(rsid, "1", 5000 + index, "C", "T", allele_id=100 + index))
    _write_clinvar(data, build, own)
    _write_clinvar(data, other, remote)
    return other


def _reads(monkeypatch):
    counts = Counter()
    original = module.get_database_info

    def counted(path, *args, **kwargs):
        counts[path.name] += 1
        return original(path, *args, **kwargs)

    monkeypatch.setattr(module, "get_database_info", counted)
    return counts


def _answer(method, carried, current):
    if method == "normal_allele_for":
        return "C" if carried and current else None
    return ("C" if current else "T") if carried else "A"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("method", ["reference_at", "normal_allele_for"])
@pytest.mark.parametrize("carried", [False, True], ids=["ordinary", "reference-carried"])
@pytest.mark.parametrize("current", [False, True], ids=["outdated", "current"])
def test_distinct_loci_do_not_reopen_metadata_per_call(
    tmp_path, monkeypatch, build, method, carried, current
):
    _pair(tmp_path, build, carried=carried)
    _tag(tmp_path, build, CURRENT if current else OLD)
    counts = _reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        for index in range(N):
            assert getattr(ann, method)(
                f"rs399949{index:04d}", build, ("1", 1000 + index)
            ) == _answer(method, carried, current)
        assert max(counts.values(), default=0) <= 1, dict(counts)


@pytest.mark.parametrize("current_build", BUILDS)
def test_build_tags_are_independent_when_lookups_interleave(tmp_path, monkeypatch, current_build):
    for build in BUILDS:
        # A flagged normal-allele REF==ALT row is deliberately excluded on a
        # current cache; an old interpreter retains its REF for compatibility.
        rows = [(*_cv(f"rs399950{i:04d}", "1", 1000 + i, "T", "T")[:-1], 1) for i in range(N)]
        _write_clinvar(tmp_path, build, rows)
        _tag(tmp_path, build, CURRENT if build == current_build else OLD)
    counts = _reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        for i in range(N):
            for build in BUILDS:
                assert ann.reference_at(f"rs399950{i:04d}", build, ("1", 1000 + i)) == (
                    None if build == current_build else "T"
                )
        assert max(counts.values(), default=0) <= 1, dict(counts)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("method", ["reference_at", "normal_allele_for"])
@pytest.mark.parametrize("initial_current", [False, True])
def test_close_reopens_changed_tag_and_clears_derived_answers(
    tmp_path, monkeypatch, build, method, initial_current
):
    _pair(tmp_path, build, carried=True)
    _tag(tmp_path, build, CURRENT if initial_current else OLD)
    counts = _reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        lookup = getattr(ann, method)
        assert lookup("rs3999490000", build, ("1", 1000)) == _answer(method, True, initial_current)
        ann.close()
        ann.close()  # close stays idempotent
        _tag(tmp_path, build, OLD if initial_current else CURRENT)
        counts.clear()
        assert lookup("rs3999490000", build, ("1", 1000)) == _answer(
            method, True, not initial_current
        )
        assert counts[clinvar_db_filename(build)] >= 1, "the changed tag must be read again"


@pytest.mark.parametrize("build", BUILDS)
def test_prepare_cache_does_not_leave_a_cached_outdated_decision(tmp_path, build):
    _pair(tmp_path, build, carried=True)
    _tag(tmp_path, build, None, baked=True)
    with ClinVarAnnotator(tmp_path) as ann:
        assert ann.reference_at("rs3999490000", build, ("1", 1000)) == "T"
        assert ann.prepare_cache() is True
        assert ann.reference_at("rs3999490000", build, ("1", 1000)) == "C"

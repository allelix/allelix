# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Bound local ClinVar reference-row SQL by batch, preserving locus answers.

Exercise ReferenceCarriedFallback.prefetch, the existing hook used by both
SNPedia and ClinPGx. Count ClinVar-table SELECTs on the input build, not wall time.
Other-build VariationID translation queries are counted separately and
excluded from the input-build chunk bound (a separate mechanism).
Synthetic SQLite only; no new API is needed to reproduce the failure.
"""

import contextlib
import math
import sqlite3
from collections import Counter

import pytest
from test_homref_ref_by_locus import _cv, _write_clinvar

from allelix.annotators import clinvar as module
from allelix.annotators.base import variant_locus
from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.annotators.reference_fallback import ReferenceCarriedFallback
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.models import Variant

BUILDS = ("GRCh37", "GRCh38")
CURRENT = f"iv:{CLINVAR_INTERPRETER_VERSION}"
OLD = "iv:0"


class _Unavailable:
    def is_ready(self):
        return False


def _fallback(ann):
    return ReferenceCarriedFallback(ann, _Unavailable(), _Unavailable())


def _tag(data, build, tag, *, baked=False):
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.execute(
            "UPDATE database_versions SET local_version_tag=?, remote_signal=?",
            (tag, f"md5:synthetic|{CURRENT}" if baked else "md5:synthetic"),
        )
        conn.commit()


def _fixture(data, build, n=21, *, current=True):
    other = next(b for b in BUILDS if b != build)
    rows, remote, calls, expected = [], [], [], []
    for i in range(n):
        rsid = f"rs399951{i:06d}"
        chrom, pos = str(1 + i % 2), 1000 + i
        call = Variant(rsid, chrom, pos, "A", "G", build=build)
        calls.append(call)
        vid = 10000 + i
        kind = i % 7
        if kind == 0:
            rows.extend(
                [
                    _cv(rsid, chrom, pos, "A", "G", allele_id=vid),
                    _cv(rsid, chrom, pos, "A", "C", allele_id=vid + 100000),
                ]
            )
            answer = ("A", None)
        elif kind == 1:
            rows.append(_cv(rsid, chrom, pos, "T", "T", allele_id=vid))
            remote.append(_cv(rsid, chrom, pos + 50000, "C", "T", allele_id=vid))
            answer = ("C", "C") if current else ("T", None)
        elif kind == 2:
            rows.append((*_cv(rsid, chrom, pos, "T", "T", allele_id=vid)[:-1], 1))
            answer = (None, None) if current else ("T", None)
        elif kind == 3:
            rows.append(_cv(rsid, chrom, pos, "G", "GA", allele_id=vid))
            answer = (None, None)  # insertion anchor is not SNV evidence
        elif kind == 4:
            # No own row; another rsID at the same position must not answer.
            rows.append(_cv(rsid + "9", chrom, pos, "C", "T", allele_id=vid))
            answer = (None, None)
        elif kind == 5:
            rows.extend(
                [
                    _cv(rsid, chrom, pos + 100000, "A", "G", allele_id=vid + 100000),
                    _cv(rsid, chrom, pos, "C", "T", allele_id=vid),
                ]
            )
            answer = ("C", None)
        else:
            rows.extend(
                [
                    _cv(rsid, "3", pos, "A", "G", allele_id=vid + 100000),
                    _cv(rsid, chrom, pos, "G", "T", allele_id=vid),
                ]
            )
            answer = ("G", None)
        expected.append(answer)
    _write_clinvar(data, build, rows)
    _write_clinvar(data, other, remote)
    _tag(data, build, CURRENT if current else OLD)
    return calls, expected


def _answers(provider, calls):
    return [
        (
            provider.reference_at(v.rsid, v.build, variant_locus(v)),
            provider.normal_allele_for(v.rsid, v.build, variant_locus(v)),
        )
        for v in calls
    ]


def _sql_reads(monkeypatch):
    """Observe real SQL on every opened/reopened annotator connection."""
    counts = Counter()
    original = module.connect_readonly

    def connected(path, *args, **kwargs):
        conn = original(path, *args, **kwargs)

        def trace(sql):
            compact = "".join(sql.upper().split())
            if "SELECT" in compact and "FROMCLINVAR_VARIANTS" in compact:
                counts[path.name] += 1

        conn.set_trace_callback(trace)
        return conn

    monkeypatch.setattr(module, "connect_readonly", connected)
    return counts


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("current", [True, False], ids=["current", "outdated"])
@pytest.mark.parametrize("n", [21, 501, 1001])
def test_reference_rows_sql_is_bounded_by_position_chunks(
    tmp_path, monkeypatch, build, current, n
):
    calls, expected = _fixture(tmp_path, build, n, current=current)
    with ClinVarAnnotator(tmp_path) as scalar:
        assert _answers(scalar, calls) == expected
    reads = _sql_reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(iter(calls))  # one-shot iterables must not lose calls
        assert _answers(provider, calls) == expected
        chunks = sum(
            math.ceil(size / 500) for size in Counter(v.chromosome for v in calls).values()
        )
        assert reads[clinvar_db_filename(build)] <= chunks, dict(reads)
        reads.clear()
        # The second consumer asks the same batch; subsets and duplicates
        # require no local-row SQL, including known-empty own-rsID answers.
        for batch in (calls, calls[::2], calls[:4] * 3):
            provider.prefetch(batch)
            assert _answers(provider, batch) == _answers(ann, batch)
        assert sum(reads.values()) == 0, dict(reads)


@pytest.mark.parametrize("build", BUILDS)
def test_prefetched_position_keeps_other_rsids_separate(tmp_path, monkeypatch, build):
    calls, _ = _fixture(tmp_path, build)
    reads = _sql_reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch([calls[4]])  # no own row at this locus
        reads.clear()
        locus = variant_locus(calls[4])
        assert provider.reference_at(calls[4].rsid, build, locus) is None
        assert provider.reference_at(calls[4].rsid + "9", build, locus) == "C"
        assert sum(reads.values()) == 0, dict(reads)


@pytest.mark.parametrize("build", BUILDS)
def test_unprefetched_locus_works_and_old_batch_is_not_retained(tmp_path, monkeypatch, build):
    calls, expected = _fixture(tmp_path, build)
    reads = _sql_reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls[:7])
        assert _answers(provider, [calls[14]]) == [expected[14]]  # outside current batch
        provider.prefetch(calls[7:14])
        reads.clear()
        assert provider.reference_at(calls[0].rsid, build, variant_locus(calls[0])) == "A"
        assert reads[clinvar_db_filename(build)] > 0, "old batches must not accumulate in memory"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("initial_current", [True, False])
def test_close_discards_prefetched_rows_and_tag(tmp_path, build, initial_current):
    calls, expected = _fixture(tmp_path, build, current=initial_current)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls)
        assert _answers(provider, calls) == expected
        ann.close()
        ann.close()
        with contextlib.closing(sqlite3.connect(tmp_path / clinvar_db_filename(build))) as conn:
            conn.execute(
                "UPDATE clinvar_variants SET ref='C',alt='G' WHERE rsid=?", (calls[0].rsid,)
            )
            conn.commit()
        _tag(tmp_path, build, OLD if initial_current else CURRENT)
        provider.prefetch(calls)
        assert provider.reference_at(calls[0].rsid, build, variant_locus(calls[0])) == "C"
        assert _answers(provider, [calls[2]]) == [
            (None, None) if not initial_current else ("T", None)
        ]


@pytest.mark.parametrize("build", BUILDS)
def test_prepare_cache_discards_prefetched_legacy_normal_flags(tmp_path, build):
    calls, _ = _fixture(tmp_path, build, current=False)
    _tag(tmp_path, build, None, baked=True)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls)
        assert _answers(provider, [calls[2]]) == [("T", None)]
        assert ann.prepare_cache() is True
        provider.prefetch(calls)
        assert _answers(provider, [calls[2]]) == [(None, None)]


@pytest.mark.parametrize("build", BUILDS)
def test_old_schema_without_normal_column_still_reads(tmp_path, build):
    calls, expected = _fixture(tmp_path, build, current=False)
    with contextlib.closing(sqlite3.connect(tmp_path / clinvar_db_filename(build))) as conn:
        conn.execute("ALTER TABLE clinvar_variants DROP COLUMN normal_allele")
        conn.commit()
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls)
        assert _answers(provider, calls) == expected


def test_build_and_chromosome_keys_do_not_alias(tmp_path, monkeypatch):
    rsid = "rs399951999999"
    _write_clinvar(
        tmp_path, "GRCh37", [_cv(rsid, "1", 100, "A", "G"), _cv(rsid, "2", 100, "C", "T")]
    )
    _write_clinvar(
        tmp_path, "GRCh38", [_cv(rsid, "1", 100, "G", "T"), _cv(rsid, "2", 100, "T", "C")]
    )
    calls = [
        Variant(rsid, f"chr{chrom}", 100, "A", "G", build=build)
        for build in BUILDS
        for chrom in (1, 2)
    ]
    reads = _sql_reads(monkeypatch)
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls)
        assert _answers(provider, calls) == [("A", None), ("C", None), ("G", None), ("T", None)]
        reads.clear()
        provider.prefetch(calls[:2])
        assert _answers(provider, calls[2:]) == [("G", None), ("T", None)]
        assert sum(reads.values()) == 0, dict(reads)


def test_missing_and_positionless_calls_do_not_borrow_a_locus(tmp_path):
    _write_clinvar(tmp_path, "GRCh37", [_cv("rs399951999999", "1", 100, "A", "G")])
    calls = [
        Variant("rs399951999999", "1", 100, "A", "G", build="GRCh38"),
        Variant("rs399951999999", "", 0, "A", "G", build="GRCh37"),
    ]
    with ClinVarAnnotator(tmp_path) as ann:
        provider = _fallback(ann)
        provider.prefetch(calls)
        assert _answers(provider, calls) == [(None, None), (None, None)]

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The direct analysis API closes supplied handles even before streaming."""

from __future__ import annotations

import contextlib
import sqlite3

import pytest

import allelix.reports._pipeline as pipeline
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize("phase", ["metadata", "build"])
@pytest.mark.parametrize("error_type", [ValueError, KeyboardInterrupt])
def test_early_failure_closes_actual_preopened_handle(tmp_path, monkeypatch, phase, error_type):
    with contextlib.closing(sqlite3.connect(tmp_path / "gnomad.sqlite")) as conn:
        conn.execute("CREATE TABLE synthetic (n INTEGER)")
    owner = GnomadAnnotator(tmp_path)
    connection = owner._connection()
    assert connection.execute("SELECT 1").fetchone() == (1,)
    parser = VcfParser()
    target, name = (
        (parser, "get_metadata") if phase == "metadata" else (pipeline, "decide_file_build")
    )
    monkeypatch.setattr(
        parser, "get_metadata", lambda path: {"sample_id": "synthetic", "build": "GRCh37"}
    )

    def refuse(*args, **kwargs):
        raise error_type("early refusal")

    monkeypatch.setattr(target, name, refuse)
    try:
        with pytest.raises(error_type, match="early refusal"):
            pipeline.run_analysis(
                tmp_path / "unused.vcf", parser, [], gnomad=owner, build_override="GRCh37"
            )
        # Observe actual closure before any test-owned cleanup or collection.
        with pytest.raises(sqlite3.ProgrammingError, match="closed"):
            connection.execute("SELECT 1")
    finally:
        owner.close()
        connection.close()


def _vcf(tmp_path):
    path = tmp_path / "watched.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=2,length=243199373>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "2\t100\trs990001\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("error_type", [RuntimeError, KeyboardInterrupt])
def test_analysis_partial_watch_registration_is_cleared(tmp_path, monkeypatch, error_type):
    parser = VcfParser()
    original = parser.watch
    error = error_type("failed after registration")
    events = []

    def watch(loci, rsids, **kwargs):
        original(loci, rsids, **kwargs)
        events.append(frozenset(rsids))
        if rsids:
            raise error

    monkeypatch.setattr(parser, "watch", watch)
    with pytest.raises(error_type) as raised:
        pipeline.run_analysis(
            _vcf(tmp_path), parser, [], build_override="GRCh37", high_value_rsids={"rs990001"}
        )
    assert raised.value is error
    assert parser._watch_rsids == frozenset()
    assert parser._watch_loci == {}
    assert events == [frozenset({"rs990001"}), frozenset()]


@pytest.mark.parametrize("reset_target", ["watch", "panel"])
@pytest.mark.parametrize("error_type", [ValueError, KeyboardInterrupt, None])
def test_analysis_reset_failure_preserves_primary_and_attempts_remaining_cleanup(
    tmp_path,
    monkeypatch,
    reset_target,
    error_type,
):
    with contextlib.closing(sqlite3.connect(tmp_path / "gnomad.sqlite")) as conn:
        conn.execute("CREATE TABLE synthetic (n INTEGER)")
    owner = GnomadAnnotator(tmp_path)
    connection = owner._connection()
    parser = VcfParser()
    watch, panel = parser.watch, parser.watch_panel
    primary = error_type("primary analysis failure") if error_type else None
    secondary = RuntimeError("watch reset failure")
    attempts = []

    def tracked_watch(loci, rsids, **kwargs):
        if not rsids:
            attempts.append("watch")
            if reset_target == "watch":
                raise secondary  # Deliberately fails BEFORE resetting its own state.
        watch(loci, rsids, **kwargs)

    def tracked_panel(rsids, *args, **kwargs):
        if not rsids:
            attempts.append("panel")
            if reset_target == "panel":
                raise secondary
        panel(rsids, *args, **kwargs)

    monkeypatch.setattr(parser, "watch", tracked_watch)
    monkeypatch.setattr(parser, "watch_panel", tracked_panel)
    if primary:

        def parse(path):
            raise primary
            yield  # Keep the parser iterator protocol.

        monkeypatch.setattr(parser, "parse", parse)
    try:
        with pytest.raises(error_type or RuntimeError) as raised:
            pipeline.run_analysis(
                _vcf(tmp_path),
                parser,
                [],
                gnomad=owner,
                build_override="GRCh37",
                high_value_rsids={"rs990001"},
                panel_rsids=frozenset({"rs990001"}),
            )
        assert raised.value is (primary if primary is not None else secondary)
        assert attempts == ["panel", "watch"]
        with pytest.raises(sqlite3.ProgrammingError, match="closed"):
            connection.execute("SELECT 1")
        # The deliberately failing reset was attempted, not claimed completed.
        assert bool(parser._watch_rsids) == (reset_target == "watch")
        assert bool(parser._panel_rsids) == (reset_target == "panel")
    finally:
        watch(set(), set())
        panel(frozenset())
        owner.close()
        connection.close()

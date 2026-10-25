# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Resolved-stream cleanup survives partial watches and failing owner cleanup."""

from __future__ import annotations

import contextlib
import sqlite3

import pytest

from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli._helpers import _resolved_variants
from allelix.parsers.vcf import VcfParser


def _vcf(tmp_path):
    path = tmp_path / "synthetic.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=2,length=243199373>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "2\t100\trs990001\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("error_type", [RuntimeError, KeyboardInterrupt])
def test_resolved_partial_watch_registration_is_cleared(tmp_path, monkeypatch, error_type):
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
    with (
        pytest.raises(error_type) as raised,
        _resolved_variants(
            parser,
            _vcf(tmp_path),
            data_dir=tmp_path,
            build_override="GRCh37",
            high_value={"rs990001"},
        ),
    ):
        pytest.fail("a failed registration cannot yield a stream")
    assert raised.value is error
    assert parser._watch_rsids == frozenset()
    assert parser._watch_loci == {}
    assert events == [frozenset({"rs990001"}), frozenset()]


@pytest.mark.parametrize("error_type", [ValueError, KeyboardInterrupt, None])
@pytest.mark.parametrize("failed_owner", [0, 1, 2])
@pytest.mark.parametrize("close_first", [False, True])
def test_resolved_attempts_all_owner_closes_and_preserves_primary(
    tmp_path,
    monkeypatch,
    error_type,
    failed_owner,
    close_first,
):
    with contextlib.closing(sqlite3.connect(tmp_path / "gnomad.sqlite")) as db:
        db.execute("CREATE TABLE synthetic (n INTEGER)")
    owners = [GnomadAnnotator(tmp_path) for _ in range(3)]
    connections = [owner._connection() for owner in owners]
    originals = [owner.close for owner in owners]
    events = []
    primary = error_type("primary resolved failure") if error_type else None
    secondary = RuntimeError("secondary provider cleanup")
    # Real acquired SQLite handles, with factory positions matching the CLI's
    # three owners. No stream query/biological behavior is replaced.
    for module, owner in zip(
        (
            "clinvar.ClinVarAnnotator",
            "gnomad.GnomadAnnotator",
            "rsid_identity.RsidIdentityAnnotator",
        ),
        owners,
        strict=True,
    ):
        monkeypatch.setattr("allelix.annotators." + module, lambda *args, owner=owner: owner)
    for i, owner in enumerate(owners):

        def close(index=i):
            events.append(index)
            if index != failed_owner or close_first:
                originals[index]()
            if index == failed_owner:
                raise secondary

        monkeypatch.setattr(owner, "close", close)
    try:
        with (
            pytest.raises(error_type or RuntimeError) as raised,
            _resolved_variants(
                VcfParser(),
                _vcf(tmp_path),
                data_dir=tmp_path,
                build_override="GRCh37",
            ),
        ):
            if primary:
                raise primary
        assert raised.value is (primary if primary is not None else secondary)
        assert events == [2, 1, 0]
        for i, conn in enumerate(connections):
            if i == failed_owner and not close_first:
                assert conn.execute("SELECT 1").fetchone() == (1,)
            else:
                with pytest.raises(sqlite3.ProgrammingError, match="closed"):
                    conn.execute("SELECT 1")
    finally:
        for close in originals:
            close()
        for conn in connections:
            conn.close()


@pytest.mark.parametrize("error_type", [ValueError, KeyboardInterrupt, None])
def test_resolved_watch_reset_failure_still_closes_every_owner(tmp_path, monkeypatch, error_type):
    with contextlib.closing(sqlite3.connect(tmp_path / "gnomad.sqlite")) as db:
        db.execute("CREATE TABLE synthetic (n INTEGER)")
    owners = [GnomadAnnotator(tmp_path) for _ in range(3)]
    connections = [owner._connection() for owner in owners]
    originals = [owner.close for owner in owners]
    events = []
    parser = VcfParser()
    watch = parser.watch
    primary = error_type("primary watched-stream failure") if error_type else None
    secondary = RuntimeError("secondary watch reset")
    for module, owner in zip(
        (
            "clinvar.ClinVarAnnotator",
            "gnomad.GnomadAnnotator",
            "rsid_identity.RsidIdentityAnnotator",
        ),
        owners,
        strict=True,
    ):
        monkeypatch.setattr("allelix.annotators." + module, lambda *args, owner=owner: owner)
    # Avoid any identity query: the fixture exercises acquired-owner cleanup only.
    monkeypatch.setattr(owners[-1], "is_ready", lambda: False)
    for i, owner in enumerate(owners):

        def close(index=i):
            events.append(index)
            originals[index]()

        monkeypatch.setattr(owner, "close", close)

    def reset(loci, rsids, **kwargs):
        if not rsids:
            events.append("watch")
            raise secondary
        watch(loci, rsids, **kwargs)

    monkeypatch.setattr(parser, "watch", reset)
    try:
        with (
            pytest.raises(error_type or RuntimeError) as raised,
            _resolved_variants(
                parser,
                _vcf(tmp_path),
                data_dir=tmp_path,
                build_override="GRCh37",
                high_value={"rs990001"},
            ),
        ):
            if primary:
                raise primary
        assert raised.value is (primary if primary is not None else secondary)
        assert events == ["watch", 2, 1, 0]
        assert parser._watch_rsids == frozenset({"rs990001"})  # Failed reset not completed.
        for conn in connections:
            with pytest.raises(sqlite3.ProgrammingError, match="closed"):
                conn.execute("SELECT 1")
    finally:
        watch(set(), set())
        for close in originals:
            close()
        for conn in connections:
            conn.close()


def _many_rows(tmp_path):
    path = _vcf(tmp_path)
    header = path.read_text(encoding="utf-8").split("2\t100\t", 1)[0]
    path.write_text(
        header
        + "".join(f"2\t{100 + i}\trs990001\tA\tG\t50\tPASS\t.\tGT\t0/1\n" for i in range(1001)),
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("mode", ["early", "exhausted", "error", "interrupt"])
def test_retained_resolved_stream_releases_actual_reader_on_context_exit(
    tmp_path, monkeypatch, mode
):
    import allelix.parsers.vcf as vcf

    path = _many_rows(tmp_path)
    opened = []
    original_open = vcf._open_vcf

    @contextlib.contextmanager
    def open_vcf(*args, **kwargs):
        with original_open(*args, **kwargs) as reader:
            opened.append(reader)
            yield reader

    monkeypatch.setattr(vcf, "_open_vcf", open_vcf)
    parser = VcfParser()
    original_parse = parser.parse
    consumed = []

    def parse(path):
        for variant in original_parse(path):
            consumed.append(variant.position)
            yield variant

    monkeypatch.setattr(parser, "parse", parse)
    stream = None
    primary = (
        ValueError("early consumer error")
        if mode == "error"
        else (KeyboardInterrupt("early consumer interrupt") if mode == "interrupt" else None)
    )
    try:
        expected = (
            pytest.raises(type(primary)) if primary is not None else contextlib.nullcontext()
        )
        with (
            expected,
            _resolved_variants(parser, path, data_dir=tmp_path, build_override="GRCh37") as stream,
        ):
            if mode == "exhausted":
                assert len(list(stream.variants)) == 1001
            else:
                assert next(stream.variants).position == 100
                if primary is not None:
                    raise primary
        assert stream is not None  # Retain every returned wrapper while observing the reader.
        assert len(opened) >= 2 and all(reader.closed for reader in opened)
        assert len(consumed) == (1001 if mode == "exhausted" else 1000)
    finally:
        if stream is not None:
            stream.variants.close()
        for reader in opened:
            reader.close()


@pytest.mark.parametrize("failed_iterator", ["output", "resolution", "build", "parser"])
@pytest.mark.parametrize("close_first", [False, True])
@pytest.mark.parametrize("error_type", [ValueError, KeyboardInterrupt, None])
def test_iterator_cleanup_failure_attempts_every_internal_owner(
    tmp_path,
    monkeypatch,
    failed_iterator,
    close_first,
    error_type,
):
    import allelix.annotators._resolution as resolution
    import allelix.cli._helpers as cli
    import allelix.parsers.vcf as vcf
    import allelix.reports._pipeline as pipeline

    path = _many_rows(tmp_path)
    with contextlib.closing(sqlite3.connect(tmp_path / "gnomad.sqlite")) as db:
        db.execute("CREATE TABLE synthetic (n INTEGER)")
    owner = GnomadAnnotator(tmp_path)
    connection = owner._connection()
    original_close = owner.close
    primary = error_type("primary iterator-consumer failure") if error_type else None
    secondary = RuntimeError("secondary iterator cleanup failure")
    opened, iterators, attempts = [], [], []
    original_open = vcf._open_vcf

    @contextlib.contextmanager
    def open_vcf(*args, **kwargs):
        with original_open(*args, **kwargs) as reader:
            opened.append(reader)
            yield reader

    monkeypatch.setattr(vcf, "_open_vcf", open_vcf)
    monkeypatch.setattr("allelix.annotators.gnomad.GnomadAnnotator", lambda *args: owner)

    def close_owner():
        attempts.append("provider")
        original_close()

    monkeypatch.setattr(owner, "close", close_owner)

    class TrackedIterator:
        def __init__(self, iterator, name):
            self.iterator = iterator
            self.name = name
            iterators.append(self)

        def __iter__(self):
            return self

        def __next__(self):
            return next(self.iterator)

        def close(self):
            attempts.append(self.name)
            if self.name != failed_iterator or close_first:
                self.iterator.close()
            if self.name == failed_iterator:
                raise secondary

    parser = VcfParser()
    for target, name, label in [
        (parser, "parse", "parser"),
        (pipeline, "with_effective_build", "build"),
        (resolution, "resolve_rsidless", "resolution"),
        (cli, "_stop_as_click", "output"),
    ]:
        original = getattr(target, name)

        def create(*args, original=original, label=label, **kwargs):
            return TrackedIterator(original(*args, **kwargs), label)

        monkeypatch.setattr(target, name, create)
    stream = None
    try:
        with (
            pytest.raises(error_type or RuntimeError) as raised,
            _resolved_variants(
                parser,
                path,
                data_dir=tmp_path,
                build_override="GRCh37",
            ) as stream,
        ):
            assert next(stream.variants).position == 100
            if primary is not None:
                raise primary
        assert raised.value is (primary if primary is not None else secondary)
        assert attempts == ["output", "resolution", "build", "parser", "provider"]
        with pytest.raises(sqlite3.ProgrammingError, match="closed"):
            connection.execute("SELECT 1")
        deliberately_open = failed_iterator == "parser" and not close_first
        assert [reader.closed for reader in opened].count(False) == int(deliberately_open)
    finally:
        # Only after observations: release the one deliberately failed owner,
        # plus any baseline leaks, without relying on GC.
        for iterator in reversed(iterators):
            iterator.iterator.close()
        original_close()
        connection.close()
        for reader in opened:
            reader.close()

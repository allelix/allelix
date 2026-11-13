"""Disabled caches cannot change startup writes or whether a report can run."""

import contextlib
import sqlite3

import pytest
from click.testing import CliRunner
from test_disabled_source_reports import CHOICES, _configure, _floors, _prepare

from allelix.annotators import get_annotators
from allelix.cli import main

COMMANDS = ("analyze", "methylation", "pharmacogenomics")


def _caches(data, source):
    return sorted(data.glob("clinvar.*.sqlite" if source == "clinvar" else source + ".sqlite"))


def _offline(monkeypatch, data):
    for annotator in get_annotators(data):
        monkeypatch.setattr(type(annotator), "fetch_remote_signal", lambda self: None)


def _args(command, path, data, output, options, *, no_update=False):
    return [
        command,
        str(path),
        "--data-dir",
        str(data),
        "--build",
        "grch38",
        "--output",
        str(output),
        *_floors(command),
        *options,
        *(["--no-update"] if no_update else []),
    ]


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("source,method", CHOICES)
def test_report_startup_leaves_switched_off_wal_caches_untouched(
    tmp_path, monkeypatch, command, source, method
):
    data, path = _prepare(tmp_path, "GRCh38")
    options = _configure(data, source, method, command=command)
    caches = _caches(data, source)
    assert caches
    for cache in caches:
        with contextlib.closing(sqlite3.connect(cache)) as conn:
            assert conn.execute("PRAGMA journal_mode=WAL").fetchone()[0] == "wal"
    before = {p.name: p.read_bytes() for p in caches}
    assert all(raw[18:20] == b"\x02\x02" for raw in before.values())
    _offline(monkeypatch, data)
    output = tmp_path / "report.json"
    result = CliRunner().invoke(main, _args(command, path, data, output, options))
    assert result.exit_code == 0, (result.output, result.exception)
    assert output.exists()
    assert {p.name: p.read_bytes() for p in caches} == before
    assert not any(p.name + ": converted" in result.output for p in caches)


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("source,method", CHOICES)
def test_only_a_switched_off_cache_cannot_make_the_report_runnable(
    tmp_path, command, source, method
):
    data, path = _prepare(tmp_path, "GRCh38")
    options = _configure(data, source, method, command=command)
    kept = _caches(data, source)
    assert kept
    for cache in data.glob("*.sqlite"):
        if cache not in kept:
            cache.rename(cache.with_suffix(".parked"))
    output = tmp_path / "report.json"
    args = _args(command, path, data, output, options, no_update=True)
    present = CliRunner().invoke(main, args)
    present_written = output.exists()
    if output.exists():
        output.unlink()
    for cache in kept:
        cache.rename(cache.with_suffix(".parked"))
    absent = CliRunner().invoke(main, args)
    assert absent.exit_code == 1, (absent.output, absent.exception)
    assert not output.exists()
    assert (present.exit_code, present_written) == (absent.exit_code, output.exists()), (
        present.output,
        present.exception,
    )
    assert present.output == absent.output


@pytest.mark.parametrize("command", COMMANDS)
def test_an_enabled_wal_cache_is_still_prepared(tmp_path, monkeypatch, command):
    data, path = _prepare(tmp_path, "GRCh38")
    _configure(data)
    cache = data / "gnomad.sqlite"
    with contextlib.closing(sqlite3.connect(cache)) as conn:
        assert conn.execute("PRAGMA journal_mode=WAL").fetchone()[0] == "wal"
    assert cache.read_bytes()[18:20] == b"\x02\x02"
    _offline(monkeypatch, data)
    result = CliRunner().invoke(main, _args(command, path, data, tmp_path / "report.json", ()))
    assert result.exit_code == 0, (result.output, result.exception)
    assert cache.read_bytes()[18:20] == b"\x01\x01"
    assert "gnomad.sqlite: converted" in result.output

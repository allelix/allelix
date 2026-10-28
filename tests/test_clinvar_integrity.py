# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Offline checksum, candidate validation and publication safety controls."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators import clinvar
from allelix.cli import main
from allelix.databases import manager
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.models import Variant
from tests.databases.test_clinvar_tsv_loader import _SS_HEADER, _VS_HEADER, _ss_row, _vs_row

if TYPE_CHECKING:
    from pathlib import Path

BUILDS = ("GRCh37", "GRCh38")
VS = "variant_summary.txt.gz"
SS = "submission_summary.txt.gz"


def _digest(payload: bytes) -> str:
    return hashlib.md5(payload, usedforsecurity=False).hexdigest()


def _payloads(counts: dict[str, int]) -> dict[str, bytes]:
    variants = []
    submissions = []
    for i in range(max(counts.values(), default=0)):
        for build, count in counts.items():
            if i < count:
                variants.append(
                    _vs_row(
                        assembly=build,
                        variation_id=str(100 + i),
                        allele_id=str(1000 + i),
                        rs=str(900940000 + i),
                        position_vcf=str(4000 + i),
                    )
                )
        submissions.append(_ss_row(variation_id=str(100 + i), scv=f"SCV{i + 1:09}.1"))
    return {
        VS: gzip.compress(("\n".join([_VS_HEADER, *variants]) + "\n").encode(), mtime=0),
        SS: gzip.compress(("\n".join([_SS_HEADER, *submissions]) + "\n").encode(), mtime=0),
    }


def _network(monkeypatch, counts, *, floors=True):
    payloads = _payloads(counts)
    sidecars = {name: f"{_digest(body)}  {name}\n" for name, body in payloads.items()}
    downloads = []

    def fetch(url, **kwargs):
        return sidecars[url.rsplit("/", 1)[-1].removesuffix(".md5")]

    def download(url, dest):
        name = url.rsplit("/", 1)[-1]
        downloads.append(name)
        dest.write_bytes(payloads[name])

    monkeypatch.setattr(clinvar, "fetch_remote_text", fetch)
    monkeypatch.setattr(clinvar, "download", download)
    if floors:
        monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, 1), raising=False)
    return payloads, sidecars, downloads


def _old_caches(data: Path, counts: dict[str, int], *, legacy=False) -> dict[str, bytes]:
    data.mkdir(exist_ok=True)
    sources = _payloads(counts)
    vs, ss = data.parent / "old-vs.gz", data.parent / "old-ss.gz"
    vs.write_bytes(sources[VS])
    ss.write_bytes(sources[SS])
    for build in counts:
        db = data / clinvar.clinvar_db_filename(build)
        manager.load_clinvar_tsv(vs, ss, db, build, remote_signal="md5:old")
        if not legacy:
            with contextlib.closing(sqlite3.connect(db)) as conn:
                conn.execute(
                    "UPDATE database_versions SET name = ? WHERE name = 'clinvar'",
                    (clinvar.clinvar_record_name(build),),
                )
                conn.commit()
    vs.unlink()
    ss.unlink()
    return {build: (data / clinvar.clinvar_db_filename(build)).read_bytes() for build in counts}


def _cache_entries(directory: Path) -> list[Path]:
    """The data directory's entries other than the ClinVar build log (2.5.2)."""
    return [p for p in directory.iterdir() if p.name != manager.CLINVAR_BUILD_LOG]


def _assert_old(data: Path, snapshots: dict[str, bytes]) -> None:
    for build, payload in snapshots.items():
        assert (data / clinvar.clinvar_db_filename(build)).read_bytes() == payload
    assert sorted(p.name for p in _cache_entries(data) if p.name != "config.toml") == sorted(
        clinvar.clinvar_db_filename(build) for build in snapshots
    )


@pytest.mark.parametrize("builds", [("GRCh37",), ("GRCh38",), BUILDS])
def test_selected_builds_publish_complete_metadata_once(tmp_path, monkeypatch, builds):
    payloads, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    with clinvar.ClinVarAnnotator(tmp_path, builds=builds) as ann:
        ann.setup()
        assert ann.is_ready()
        assert ann.cached_remote_signal() == f"md5:{_digest(payloads[SS])}"
    assert sorted(downloads) == sorted([VS, SS])
    assert sorted(p.name for p in _cache_entries(tmp_path)) == sorted(
        clinvar.clinvar_db_filename(build) for build in builds
    )
    for build in builds:
        info = manager.get_database_info(
            tmp_path / clinvar.clinvar_db_filename(build), clinvar.clinvar_record_name(build)
        )
        assert info is not None and info["record_count"] == 3
        assert info["version"] and info["downloaded_at"]
        assert info["local_version_tag"] == f"iv:{CLINVAR_INTERPRETER_VERSION}"


@pytest.mark.parametrize("source", [VS, SS])
@pytest.mark.parametrize("cut", ["changed", "member-boundary", "inside-member"])
def test_compressed_source_integrity_rejection_keeps_both_caches(
    tmp_path, monkeypatch, source, cut
):
    data = tmp_path / "data"
    snapshots = _old_caches(data, dict.fromkeys(BUILDS, 2))
    payloads, sidecars, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    if cut == "changed":
        sidecars[source] = f"{'0' * 32}  {source}"
    else:
        decoded = gzip.decompress(payloads[source])
        lines = decoded.splitlines(keepends=True)
        member = gzip.compress(b"".join(lines[:2]), mtime=0)
        complete = member + gzip.compress(b"".join(lines[2:]), mtime=0)
        sidecars[source] = f"{_digest(complete)}  {source}"
        payloads[source] = member if cut == "member-boundary" else complete[:-7]
        if cut == "member-boundary":
            assert gzip.decompress(payloads[source]) == b"".join(lines[:2])
    with (
        clinvar.ClinVarAnnotator(data) as ann,
        pytest.raises(OSError, match="Integrity check failed"),
    ):
        ann.setup()
    _assert_old(data, snapshots)


@pytest.mark.parametrize("source", [VS, SS])
@pytest.mark.parametrize(
    "body",
    [
        None,
        "",
        "<html>error</html>",
        "a" * 31,
        "z" * 32,
        "a" * 32 + " wrong.gz",
        "a" * 32 + " variant_summary.txt.gz\n" + "b" * 32 + " variant_summary.txt.gz",
        "a" * 32 + " variant_summary.txt.gz\nmalformed",
        "a" * 4097,
    ],
    ids=[
        "missing",
        "empty",
        "html",
        "short",
        "nonhex",
        "wrong-name",
        "conflicting",
        "malformed-line",
        "oversized",
    ],
)
def test_required_sidecar_fails_before_download(tmp_path, monkeypatch, body, source):
    _, sidecars, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    sidecars[source] = body.replace(VS, source) if isinstance(body, str) else body
    with (
        clinvar.ClinVarAnnotator(tmp_path) as ann,
        pytest.raises(
            RuntimeError,
            match=r"variant_summary.*checksum"
            if source == VS
            else "checksum for submission_summary is unavailable",
        ),
    ):
        ann.setup()
    assert downloads == []
    assert _cache_entries(tmp_path) == []


@pytest.mark.parametrize("source", [VS, SS])
@pytest.mark.parametrize("shape", ["uppercase", "absolute", "digest-only", "binary", "duplicate"])
def test_supported_checksum_shapes(tmp_path, monkeypatch, source, shape):
    payloads, sidecars, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    digest = _digest(payloads[source])
    sidecars[source] = {
        "uppercase": f"{digest.upper()}  {source}\r\n",
        "absolute": f"{digest}  /pub/clinvar/tab_delimited/{source}\r\n",
        "digest-only": digest + "\n",
        "binary": f"{digest} *{source}\n",
        "duplicate": f"{digest}  {source}\n{digest}  {source}\n",
    }[shape]
    with clinvar.ClinVarAnnotator(tmp_path) as ann:
        ann.setup()
        assert ann.is_ready()


def test_sidecar_http_read_is_bounded(monkeypatch, tmp_path):
    reads = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            pass

        def read(self, size=-1):
            reads.append(size)
            assert 0 < size <= 4097
            return b"a" * size

    monkeypatch.setattr(manager.urllib.request, "urlopen", lambda *_a, **_k: Response())
    assert clinvar.ClinVarAnnotator(tmp_path).fetch_remote_signal() is None
    assert reads == [4097]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("count,accepted", [(4, False), (5, True)])
def test_first_cache_fixed_floor_boundary(tmp_path, monkeypatch, build, count, accepted):
    _network(monkeypatch, {build: count})
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, 5))
    with clinvar.ClinVarAnnotator(tmp_path, builds=(build,)) as ann:
        if accepted:
            ann.setup()
            assert ann.record_count() == count
        else:
            with pytest.raises(manager.ClinVarFormatError, match=rf"{build}.*4.*5"):
                ann.setup()
            assert _cache_entries(tmp_path) == []


@pytest.mark.parametrize("build", BUILDS)
def test_production_setup_enforces_default_five_million_floor(tmp_path, monkeypatch, build):
    _network(monkeypatch, {build: 1}, floors=False)
    assert getattr(manager, "CLINVAR_MIN_ROWS", None) == dict.fromkeys(BUILDS, 5_000_000)
    with (
        clinvar.ClinVarAnnotator(tmp_path, builds=(build,)) as ann,
        pytest.raises(manager.ClinVarFormatError, match=rf"{build}.*1.*5,000,000"),
    ):
        ann.setup()
    assert _cache_entries(tmp_path) == []


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("new,accepted", [(9, False), (10, True)])
@pytest.mark.parametrize("legacy", [False, True])
def test_relative_floor_rounds_up_and_uses_each_old_build(
    tmp_path, monkeypatch, build, new, accepted, legacy
):
    data = tmp_path / "data"
    old = {"GRCh37": 20, "GRCh38": 20, build: 11}
    snapshots = _old_caches(data, old, legacy=legacy)
    counts = {"GRCh37": 19, "GRCh38": 19, build: new}
    _network(monkeypatch, counts)
    with clinvar.ClinVarAnnotator(data) as ann:
        if accepted:
            ann.setup()
            assert ann.is_ready()
        else:
            with pytest.raises(manager.ClinVarFormatError, match=rf"{build}.*9.*11.*10"):
                ann.setup()
            _assert_old(data, snapshots)


def test_small_existing_cache_does_not_weaken_fixed_floor(tmp_path, monkeypatch):
    data = tmp_path / "data"
    snapshots = _old_caches(data, dict.fromkeys(BUILDS, 2))
    _network(monkeypatch, dict.fromkeys(BUILDS, 4))
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, 5))
    with clinvar.ClinVarAnnotator(data) as ann, pytest.raises(manager.ClinVarFormatError):
        ann.setup()
    _assert_old(data, snapshots)


@pytest.mark.parametrize("bad", [None, "unknown", -1, 0, 1.5])
def test_bad_existing_metadata_is_not_treated_as_absent(tmp_path, monkeypatch, bad):
    data = tmp_path / "data"
    _old_caches(data, dict.fromkeys(BUILDS, 2))
    path = data / clinvar.clinvar_db_filename("GRCh38")
    with contextlib.closing(sqlite3.connect(path)) as conn:
        if bad is None:
            conn.execute("DELETE FROM database_versions")
        else:
            conn.execute("UPDATE database_versions SET record_count = ?", (bad,))
        conn.commit()
    snapshots = {
        build: (data / clinvar.clinvar_db_filename(build)).read_bytes() for build in BUILDS
    }
    _, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    with (
        clinvar.ClinVarAnnotator(data) as ann,
        pytest.raises(RuntimeError, match=r"GRCh38.*record count"),
    ):
        ann.setup()
    assert downloads == []
    _assert_old(data, snapshots)


@pytest.mark.parametrize(
    "failure", ["count", "empty", "parse", "metadata", "metadata-missing", "download"]
)
@pytest.mark.parametrize("existing", [False, True])
def test_second_candidate_prepublication_failure_preserves_all_selected_caches(
    tmp_path, monkeypatch, failure, existing
):
    data = tmp_path / "data"
    snapshots = _old_caches(data, dict.fromkeys(BUILDS, 2) if existing else {})
    payloads, sidecars, _ = _network(
        monkeypatch,
        {"GRCh37": 3, "GRCh38": 0 if failure == "empty" else 1 if failure == "count" else 3},
    )
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, 2))
    if failure == "parse":
        original = manager.iter_clinvar_tsv_records

        def records(vs, ss, build, **kwargs):
            if build == "GRCh38":
                raise manager.ClinVarFormatError("synthetic second-build parse failure")
            yield from original(vs, ss, build, **kwargs)

        monkeypatch.setattr(manager, "iter_clinvar_tsv_records", records)
    elif failure in ("metadata", "metadata-missing"):
        original_stamp = manager.stamp_remote_signal

        def stamp(conn, name, *args):
            if name.endswith("GRCh38"):
                if failure == "metadata":
                    raise sqlite3.OperationalError("synthetic metadata failure")
                return None
            return original_stamp(conn, name, *args)

        monkeypatch.setattr(manager, "stamp_remote_signal", stamp)
    elif failure == "download":

        def download(url, dest):
            if url.endswith(SS):
                raise OSError("synthetic second-download failure")
            dest.write_bytes(payloads[VS])

        monkeypatch.setattr(clinvar, "download", download)
    assert sidecars[SS]
    with (
        clinvar.ClinVarAnnotator(data) as ann,
        pytest.raises((RuntimeError, OSError, sqlite3.DatabaseError)),
    ):
        ann.setup()
    _assert_old(data, snapshots)


@pytest.mark.parametrize("failed_build", BUILDS)
def test_publication_failure_reports_actual_partial_result(tmp_path, monkeypatch, failed_build):
    data = tmp_path / "data"
    snapshots = _old_caches(data, dict.fromkeys(BUILDS, 2))
    payloads, _, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    original = manager.os.replace

    def replace(source, dest):
        if dest == data / clinvar.clinvar_db_filename(failed_build):
            raise OSError("synthetic replacement failure")
        original(source, dest)

    monkeypatch.setattr(manager.os, "replace", replace)
    with clinvar.ClinVarAnnotator(data) as ann, pytest.raises(RuntimeError) as error:
        ann.setup()
    assert failed_build in str(error.value)
    assert "old cache kept" in str(error.value)
    failed = data / clinvar.clinvar_db_filename(failed_build)
    assert failed.read_bytes() == snapshots[failed_build]
    first = data / clinvar.clinvar_db_filename("GRCh37")
    if failed_build == "GRCh37":
        _assert_old(data, snapshots)
        assert "none" in str(error.value).lower()
    else:
        assert "already updated: GRCh37" in str(error.value)
        info = manager.get_database_info(first, "clinvar.GRCh37")
        assert info is not None and info["record_count"] == 3
        assert info["remote_signal"] == f"md5:{_digest(payloads[SS])}"
        assert first.read_bytes() != snapshots["GRCh37"]
        assert sorted(p.name for p in _cache_entries(data)) == sorted(
            clinvar.clinvar_db_filename(b) for b in BUILDS
        )


def test_null_freshness_loader_still_stamps_complete_metadata(tmp_path):
    payloads = _payloads({"GRCh38": 2})
    vs, ss = tmp_path / VS, tmp_path / SS
    vs.write_bytes(payloads[VS])
    ss.write_bytes(payloads[SS])
    db = tmp_path / "cache.sqlite"
    manager.load_clinvar_tsv(
        vs, ss, db, "GRCh38", source_url="test://variant-summary", remote_signal=None
    )
    info = manager.get_database_info(db, "clinvar")
    assert info is not None and info["remote_signal"] is None
    assert info["record_count"] == 2 and info["version"] and info["downloaded_at"]
    assert info["source_url"] == "test://variant-summary"
    assert info["local_version_tag"] == f"iv:{CLINVAR_INTERPRETER_VERSION}"


@pytest.mark.parametrize("failure", ["checksum", "floor", "parse"])
def test_force_update_cannot_bypass_integrity(tmp_path, monkeypatch, failure):
    snapshots = _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    payloads, sidecars, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    if failure == "checksum":
        sidecars[VS] = None
    elif failure == "floor":
        monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, 3))
    else:
        payloads[VS] = gzip.compress(b"#changed header\n", mtime=0)
        sidecars[VS] = f"{_digest(payloads[VS])} {VS}"
    db_module = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(
        db_module, "get_annotators", lambda data, **kw: [clinvar.ClinVarAnnotator(data)]
    )
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(tmp_path), "--force"])
    assert result.exit_code == 1, result.output
    _assert_old(tmp_path, snapshots)


def test_healthy_unchanged_update_only_probes_freshness(tmp_path, monkeypatch):
    _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    for build in BUILDS:
        with contextlib.closing(
            sqlite3.connect(tmp_path / clinvar.clinvar_db_filename(build))
        ) as conn:
            conn.execute("UPDATE database_versions SET remote_signal = ?", ("md5:" + "a" * 32,))
            conn.commit()
    snapshots = {
        build: (tmp_path / clinvar.clinvar_db_filename(build)).read_bytes() for build in BUILDS
    }
    mtimes = {
        build: (tmp_path / clinvar.clinvar_db_filename(build)).stat().st_mtime_ns
        for build in BUILDS
    }
    calls = []

    def fetch(url, **kwargs):
        calls.append(url)
        return "a" * 32 + " " + SS

    def unexpected(*args, **kwargs):
        pytest.fail("unchanged update must not download, copy or rebuild")

    monkeypatch.setattr(clinvar, "fetch_remote_text", fetch)
    monkeypatch.setattr(clinvar, "download", unexpected)
    monkeypatch.setattr(clinvar, "load_clinvar_tsv", unexpected)
    db_module = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(
        db_module, "get_annotators", lambda data, **kw: [clinvar.ClinVarAnnotator(data)]
    )
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert "already current" in result.output
    assert calls == [manager.CLINVAR_SUBMISSION_SUMMARY_URL + ".md5"]
    _assert_old(tmp_path, snapshots)
    assert {
        build: (tmp_path / clinvar.clinvar_db_filename(build)).stat().st_mtime_ns
        for build in BUILDS
    } == mtimes


@pytest.mark.parametrize("source", [VS, SS])
def test_matching_checksum_does_not_bypass_damaged_gzip(tmp_path, monkeypatch, source):
    data = tmp_path / "data"
    snapshots = _old_caches(data, dict.fromkeys(BUILDS, 2))
    payloads, sidecars, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    payloads[source] = payloads[source][:-7]
    sidecars[source] = f"{_digest(payloads[source])} {source}"
    with clinvar.ClinVarAnnotator(data) as ann, pytest.raises((EOFError, OSError)):
        ann.setup()
    _assert_old(data, snapshots)


def _no_network(monkeypatch):
    def refuse(*args, **kwargs):
        raise AssertionError("no network access expected")

    monkeypatch.setattr(clinvar, "fetch_remote_text", refuse)
    monkeypatch.setattr(clinvar, "download", refuse)


def _rows(data: Path) -> list[tuple]:
    with clinvar.ClinVarAnnotator(data) as ann:
        assert ann.is_ready()
        return [
            (build, a.rsid, a.significance, a.gene, a.genotype_match, a.zygosity)
            for build in BUILDS
            for a in ann.batch_annotate(
                [
                    Variant(f"rs{900940000 + i}", "17", 4000 + i, "G", "A", build=build, ref="G")
                    for i in range(2)
                ]
            )
        ]


def test_cache_built_by_2_5_0_stays_ready_without_rebuild(tmp_path, monkeypatch):
    """The integrity checks change nothing stored, so iv:9 caches stay current.

    ``_old_caches`` stamps the way 2.5.0 did: a ``clinvar`` row renamed to the
    per-build name after the load.
    """
    snapshots = _old_caches(tmp_path / "old", dict.fromkeys(BUILDS, 2))
    for build in BUILDS:
        info = manager.get_database_info(
            tmp_path / "old" / clinvar.clinvar_db_filename(build),
            clinvar.clinvar_record_name(build),
        )
        assert info is not None and info["local_version_tag"] == "iv:9"
    assert CLINVAR_INTERPRETER_VERSION == 9
    with monkeypatch.context() as patcher:
        _no_network(patcher)
        with clinvar.ClinVarAnnotator(tmp_path / "old") as ann:
            ann.prepare_cache()
            assert ann.is_ready()
        old_rows = _rows(tmp_path / "old")
    _assert_old(tmp_path / "old", snapshots)

    _network(monkeypatch, dict.fromkeys(BUILDS, 2))
    with clinvar.ClinVarAnnotator(tmp_path / "fresh") as ann:
        ann.setup()
    assert old_rows == _rows(tmp_path / "fresh")
    assert [row[:3] for row in old_rows] == [
        (build, f"rs{900940000 + i}", "clinvar_pathogenic") for build in BUILDS for i in range(2)
    ]


def _stamp_signal(data: Path, signal: str) -> None:
    for build in BUILDS:
        with contextlib.closing(
            sqlite3.connect(data / clinvar.clinvar_db_filename(build))
        ) as conn:
            conn.execute("UPDATE database_versions SET remote_signal = ?", (signal,))
            conn.commit()


def _update(data: Path, monkeypatch):
    db_module = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(
        db_module, "get_annotators", lambda data, **kw: [clinvar.ClinVarAnnotator(data)]
    )
    return CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])


def test_ordinary_update_keeps_2_5_0_cache_with_equal_freshness(tmp_path, monkeypatch):
    _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    payloads, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    _stamp_signal(tmp_path, f"md5:{_digest(payloads[SS])}")
    snapshots = {
        build: (tmp_path / clinvar.clinvar_db_filename(build)).read_bytes() for build in BUILDS
    }
    result = _update(tmp_path, monkeypatch)
    assert result.exit_code == 0, result.output
    assert "already current" in result.output
    assert downloads == []
    _assert_old(tmp_path, snapshots)


def test_changed_freshness_rebuilds_2_5_0_cache_through_integrity_checks(tmp_path, monkeypatch):
    _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    payloads, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    signal = f"md5:{_digest(payloads[SS])}"
    verified = []
    original_verify = clinvar.verify_file_hash

    def verify(path, algorithm, expected):
        verified.append(path.name)
        return original_verify(path, algorithm, expected)

    monkeypatch.setattr(clinvar, "verify_file_hash", verify)
    first = _update(tmp_path, monkeypatch)
    assert first.exit_code == 0, first.output
    assert sorted(downloads) == sorted([VS, SS])
    assert sorted(verified) == sorted([VS, SS])
    for build in BUILDS:
        info = manager.get_database_info(
            tmp_path / clinvar.clinvar_db_filename(build), clinvar.clinvar_record_name(build)
        )
        assert info is not None and info["local_version_tag"] == "iv:9"
        assert info["record_count"] == 3 and info["remote_signal"] == signal
    second = _update(tmp_path, monkeypatch)
    assert second.exit_code == 0, second.output
    assert "already current" in second.output
    assert len(downloads) == 2


def _flaky_submission_checksum(monkeypatch, failures: int | None) -> list[str]:
    """Fail the submission_summary .md5 fetch ``failures`` times (None: always)."""
    sidecar_fetch = clinvar.fetch_remote_text
    calls: list[str] = []

    def fetch(url, **kwargs):
        if url.endswith(SS + ".md5"):
            calls.append(url)
            if failures is None or len(calls) <= failures:
                return None
        return sidecar_fetch(url, **kwargs)

    monkeypatch.setattr(clinvar, "fetch_remote_text", fetch)
    return calls


def _force(tmp_path):
    return CliRunner().invoke(
        main, ["db", "update", "--only", "clinvar", "--data-dir", str(tmp_path), "--force"]
    )


def test_force_retries_a_transient_checksum_failure_and_stores_it(tmp_path, monkeypatch):
    _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    payloads, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    calls = _flaky_submission_checksum(monkeypatch, failures=1)
    result = _force(tmp_path)
    assert result.exit_code == 0, result.output
    assert len(calls) == 2 and sorted(downloads) == sorted([VS, SS])
    signal = f"md5:{_digest(payloads[SS])}"
    for build in BUILDS:
        info = manager.get_database_info(
            tmp_path / clinvar.clinvar_db_filename(build), clinvar.clinvar_record_name(build)
        )
        assert info is not None and info["remote_signal"] == signal
        assert info["record_count"] == 3 and info["local_version_tag"] == "iv:9"


def test_force_never_refreshes_without_the_checksum(tmp_path, monkeypatch):
    """The reviewer's FORCE_NO_SIGNAL case: no unverified refresh, no --force loop."""
    snapshots = _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    _, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    _flaky_submission_checksum(monkeypatch, failures=None)
    result = _force(tmp_path)
    output = " ".join(result.output.split())
    assert result.exit_code == 1
    error = output.split("ClinVar not refreshed", 1)[1]
    assert error.startswith(": NCBI's checksum for submission_summary is unavailable")
    assert "Existing caches kept; try again later." in error
    assert "--force" not in error
    assert downloads == []
    _assert_old(tmp_path, snapshots)


def test_without_force_a_missing_checksum_keeps_the_cache(tmp_path, monkeypatch):
    snapshots = _old_caches(tmp_path, dict.fromkeys(BUILDS, 2))
    _, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 3))
    _flaky_submission_checksum(monkeypatch, failures=None)
    result = CliRunner().invoke(
        main, ["db", "update", "--only", "clinvar", "--data-dir", str(tmp_path)]
    )
    output = " ".join(result.output.split())
    assert result.exit_code == 0, result.output
    line = output.split("clinvar:", 1)[1]
    assert "cache kept. ClinVar is never refreshed unverified; try again later." in line
    assert "--force" not in line
    assert downloads == []
    _assert_old(tmp_path, snapshots)

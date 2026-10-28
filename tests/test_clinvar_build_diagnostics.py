# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Keep measured ClinVar ingestion evidence in the build log, including failed builds."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import logging
import re
import sqlite3

import pytest

from allelix.databases import manager


def _sources(tmp_path, build, compressed, *, count=2):
    """Two accepted variants, one joined group, one missing group, and rejected source rows."""
    vs_header = [manager._VS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(43)]
    vs_header[0] = "#AlleleID"
    vs_rows = []
    for identifier, assembly, pos in (
        (101, build, "1000"),
        (102, build, "2000"),
        (103, build, "not-a-position"),
        (104, "GRCh36", "4000"),
    ):
        row = [""] * 43
        for index, value in {
            0: str(identifier),
            1: "single nucleotide variant",
            2: "NM_SYNTHETIC:c.1A>G",
            4: "GENE",
            6: "Pathogenic",
            9: str(900252900 + identifier),
            16: assembly,
            18: "1",
            24: "criteria provided, single submitter",
            30: str(identifier),
            31: pos,
            32: "A",
            33: "G",
        }.items():
            row[index] = value
        vs_rows.append("\t".join(row))
    ss_header = [manager._SS_REQUIRED_COLUMNS.get(i, f"unused{i}") for i in range(16)]
    ss_header[0] = "#VariationID"
    ss_rows = []
    for number, (identifier, significance) in enumerate(
        [(101, "Pathogenic")] * count + [(101, "not specified"), (199, "Pathogenic")], 1
    ):
        row = [""] * 16
        for index, value in {
            0: str(identifier),
            1: significance,
            5: "na:Synthetic condition",
            6: "criteria provided, single submitter",
            10: f"SCV9002529{number:02}.1",
            15: "yes",
        }.items():
            row[index] = value
        ss_rows.append("\t".join(row))
    paths = []
    for name, header, rows in (
        ("variant_summary", vs_header, vs_rows),
        ("submission_summary", ss_header, ss_rows),
    ):
        payload = ("\t".join(header) + "\n" + "\n".join(rows) + "\n").encode()
        path = tmp_path / (name + (".txt.gz" if compressed else ".txt"))
        path.write_bytes(gzip.compress(payload, mtime=0) if compressed else payload)
        paths.append(path)
    return paths


def _old_cache(path):
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.execute("CREATE TABLE retained(value TEXT)")
        conn.execute("INSERT INTO retained VALUES ('old cache must survive')")
        conn.commit()
    return path.read_bytes()


def _log_build(log_path, call):
    logger = logging.getLogger("allelix")
    old_level = logger.level
    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s %(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    try:
        call()
    finally:
        handler.flush()
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(old_level)


def _metric(text, aliases, value):
    normalized = re.sub(r"[-_\s]+", " ", text.lower().replace(",", ""))
    names = "|".join(re.escape(alias) for alias in aliases)
    assert re.search(rf"\b(?:{names})\s*[:=]\s*{value}\b", normalized), (aliases, value, text)


def _assert_evidence(log_path, paths, build, *, yielded=2, joined=1, inserted=None):
    inserted = yielded if inserted is None else inserted
    text = log_path.read_text()
    assert build in text, text
    for path in paths:
        payload = path.read_bytes()
        digest = hashlib.md5(payload, usedforsecurity=False).hexdigest()
        # A digest/size must be attributable to its source, not an unlabeled pair of hashes.
        lines = [line for line in text.splitlines() if path.name in line]
        assert any(digest in line.lower() for line in lines), (path.name, digest, text)
        assert any(re.search(rf"\b{len(payload)}\b", line.replace(",", "")) for line in lines), (
            text
        )
    for aliases, value in (
        (
            ("variant rows accepted", "variant summary accepted", "variant summary rows accepted"),
            2,
        ),
        (
            (
                "submission rows accepted",
                "submission summary accepted",
                "submission summary rows accepted",
            ),
            yielded + 1,
        ),
        (("join hits", "joined groups"), joined),
        (("join misses", "unmatched groups"), 1),
        (("records yielded", "yielded records"), yielded),
        (("rows inserted", "inserted rows"), inserted),
        (("final count", "final sqlite count", "sqlite count"), inserted),
    ):
        _metric(text, aliases, value)


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize("outcome", ["success", "fixed-floor", "relative-floor", "rollback"])
def test_build_evidence_survives_every_exit(tmp_path, monkeypatch, build, compressed, outcome):
    vs, ss = _sources(tmp_path, build, compressed)
    cache = tmp_path / "clinvar.sqlite"
    before = _old_cache(cache)
    log = tmp_path / "build.log"
    kwargs = {"min_rows": 3} if outcome == "fixed-floor" else {}
    if outcome == "relative-floor":
        kwargs["previous_count"] = 3  # ceil(90% of 3) is 3, not 2.
    if outcome == "rollback":

        def fail_metadata(*args, **kwargs):
            raise RuntimeError("synthetic metadata failure after insertion")

        monkeypatch.setattr(manager, "stamp_remote_signal", fail_metadata)

    def build_cache():
        if outcome == "rollback":
            with pytest.raises(RuntimeError, match="synthetic metadata failure"):
                manager.load_clinvar_tsv(vs, ss, cache, build, **kwargs)
        elif outcome.endswith("floor"):
            with pytest.raises(manager.ClinVarFormatError, match="candidate has 2 rows"):
                manager.load_clinvar_tsv(vs, ss, cache, build, **kwargs)
        else:
            manager.load_clinvar_tsv(vs, ss, cache, build)

    _log_build(log, build_cache)
    assert not cache.with_name(cache.name + ".partial").exists()
    if outcome == "success":
        with contextlib.closing(sqlite3.connect(f"file:{cache}?mode=ro", uri=True)) as conn:
            assert conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0] == 2
    else:
        assert cache.read_bytes() == before
    # Read the persisted log after the transaction has closed and partial cleanup has run.
    _assert_evidence(log, (vs, ss), build)


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "count,minimum,previous,accepted",
    [
        (9, 0, 10, True),
        (8, 0, 10, False),
        (2, 2, None, True),
        (2, 3, None, False),
        (0, 0, None, False),
    ],
)
def test_diagnostics_do_not_weaken_the_fixed_relative_or_nonempty_floor(
    tmp_path, build, count, minimum, previous, accepted
):
    vs, ss = _sources(tmp_path, build, False, count=count)
    cache = tmp_path / "clinvar.sqlite"
    before = _old_cache(cache)
    if accepted:
        manager.load_clinvar_tsv(vs, ss, cache, build, min_rows=minimum, previous_count=previous)
        with contextlib.closing(sqlite3.connect(f"file:{cache}?mode=ro", uri=True)) as conn:
            assert conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0] == count
    else:
        with pytest.raises(manager.ClinVarFormatError):
            manager.load_clinvar_tsv(
                vs, ss, cache, build, min_rows=minimum, previous_count=previous
            )
        assert cache.read_bytes() == before


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("outcome", ["success", "floor", "publication-failure"])
def test_setup_keeps_evidence_after_owned_sources_are_cleaned(
    tmp_path, monkeypatch, build, outcome
):
    import shutil

    from allelix.annotators import clinvar as module
    from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename

    sources = tmp_path / "sources"
    sources.mkdir()
    vs, ss = _sources(sources, build, True)
    digests = {
        p.name: hashlib.md5(p.read_bytes(), usedforsecurity=False).hexdigest() for p in (vs, ss)
    }
    data = tmp_path / "data"
    data.mkdir()
    cache = data / clinvar_db_filename(build)
    owned = []

    def local_download(url, dest):
        source = vs if "variant_summary" in url else ss
        shutil.copyfile(source, dest)
        owned.append(dest)

    monkeypatch.setattr(module, "download", local_download)
    monkeypatch.setattr(module, "_fetch_checksum", lambda _url: digests[vs.name])
    monkeypatch.setattr(
        ClinVarAnnotator,
        "_fetch_remote_signal_for_tsv",
        staticmethod(lambda: "md5:" + digests[ss.name]),
    )
    # Tiny inputs exercise the unchanged production comparison with a small test-only floor.
    monkeypatch.setitem(manager.CLINVAR_MIN_ROWS, build, 3 if outcome == "floor" else 2)
    if outcome == "publication-failure":
        replace = module.os.replace

        def fail_publish(source, dest):
            if dest == cache:
                raise OSError("synthetic publication failure")
            return replace(source, dest)

        monkeypatch.setattr(module.os, "replace", fail_publish)
    log = tmp_path / "build.log"

    def build_cache():
        with ClinVarAnnotator(data, builds=(build,)) as annotator:
            if outcome == "floor":
                with pytest.raises(manager.ClinVarFormatError, match="candidate has 2 rows"):
                    annotator.setup()
            elif outcome == "publication-failure":
                with pytest.raises(RuntimeError, match="publication failed"):
                    annotator.setup()
            else:
                annotator.setup()
                assert annotator.is_ready()

    _log_build(log, build_cache)
    assert len(owned) == 2 and all(not p.exists() for p in owned)
    assert not list(data.glob(".clinvar-*"))
    if outcome != "success":
        assert not cache.exists()
    _assert_evidence(log, (vs, ss), build)


def test_production_fixed_floors_remain_five_million():
    assert manager.CLINVAR_MIN_ROWS == {"GRCh37": 5_000_000, "GRCh38": 5_000_000}


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_yielded_and_inserted_counts_are_measured_independently(tmp_path, monkeypatch, build):
    vs, ss = _sources(tmp_path, build, False)
    cache = tmp_path / "clinvar.sqlite"
    before = _old_cache(cache)
    # Simulate loss at the storage boundary, without changing the iterator's two yielded rows.
    # The floor must refuse one stored row, and the evidence must locate the loss after yielding.
    monkeypatch.setattr(
        manager,
        "CLINVAR_SCHEMA",
        manager.CLINVAR_SCHEMA
        + """
        CREATE TRIGGER synthetic_insert_loss BEFORE INSERT ON clinvar_variants
        WHEN EXISTS (SELECT 1 FROM clinvar_variants)
        BEGIN SELECT RAISE(IGNORE); END;
    """,
    )
    log = tmp_path / "build.log"

    def build_cache():
        with pytest.raises(manager.ClinVarFormatError, match="candidate has 1 rows"):
            manager.load_clinvar_tsv(vs, ss, cache, build, min_rows=2)

    _log_build(log, build_cache)
    assert cache.read_bytes() == before
    assert not cache.with_name(cache.name + ".partial").exists()
    _assert_evidence(log, (vs, ss), build, yielded=2, inserted=1)

"""O5: independent operations own their staging bytes, including on failure.

Scheduling only is injected; production writers, floors, stamps and publication
run. Nested publication is a deterministic interleaving, without timing sleeps.
The outside review's two-thread empty-cache publication is retained separately.
"""

import contextlib
import gzip
import io
import os
import sqlite3
from pathlib import Path

import pytest

from allelix.databases import gwas_loader, loader_utils, manager, pharmgkb_loader
from allelix.databases.rsid_identity_loader import write_identity_db
from allelix.databases.schema import GNOMAD_SCHEMA
from tests.test_clinvar_build_diagnostics import _sources

GWAS = dict(
    rsid="rs900000001",
    risk_allele="A",
    trait="Synthetic",
    p_value=1e-9,
    or_beta=1.1,
    ci_text="",
    gene="TEST",
    study_accession="GCST000000",
    pubmed_id="1",
    risk_allele_frequency=0.1,
    context="",
    mapped_trait_uri="",
    trait_category="other",
)
PHARM = dict(
    rsid="rs900000001",
    genotype="AA",
    gene="TEST",
    drugs="Drug",
    phenotype="Response",
    phenotype_category="Efficacy",
    annotation_text="Synthetic",
    level_of_evidence="1A",
    score=10,
    pgkb_annotation_id="1",
    allele_function="",
    function_class="decreased",
    is_nonfinding=False,
)


def harness(root, monkeypatch, kind):
    """Return (run(label), read_stamp), using real writers and tiny inputs."""
    root.mkdir()
    db = root / "cache.sqlite"
    if kind == "gwas":
        monkeypatch.setattr(
            gwas_loader, "iter_gwas_records", lambda p: iter([dict(GWAS, trait=p.stem)])
        )

        def run(label):
            return gwas_loader.load_gwas_tsv(
                root / f"{label}.tsv", db, min_rows=1, remote_signal=label
            )
    elif kind == "pharmgkb":
        monkeypatch.setattr(
            pharmgkb_loader,
            "iter_pharmgkb_records",
            lambda p, lookup: iter([dict(PHARM, annotation_text=p.stem)]),
        )

        def run(label):
            return pharmgkb_loader.load_pharmgkb_tsv(
                root / f"{label}.zip", db, min_rows=1, version=label
            )
    elif kind == "clinvar":
        paths = _sources(root, "GRCh38", False)

        def run(label):
            return manager.load_clinvar_tsv(*paths, db, "GRCh38", min_rows=1, remote_signal=label)
    elif kind == "prebuilt":
        for label in ("seed", "A", "B"):
            raw = root / f"{label}.sqlite"
            with contextlib.closing(sqlite3.connect(raw)) as conn, conn:
                conn.executescript(GNOMAD_SCHEMA)
                conn.execute(
                    "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
                    "VALUES ('1',1,'G','A','rs1',.1)"
                )
                conn.execute(
                    "INSERT INTO database_versions "
                    "(name,source_url,version,downloaded_at,record_count) "
                    "VALUES ('gnomad','test',?,'test',1)",
                    (label,),
                )
            (root / f"{label}.gz").write_bytes(gzip.compress(raw.read_bytes(), mtime=0))

        def run(label):
            return loader_utils.install_prebuilt_gz_cache(root / f"{label}.gz", db, "gnomad")
    elif kind == "download":

        class Response(io.BytesIO):
            @property
            def headers(self):
                return {"Content-Length": str(len(self.getvalue()))}

        monkeypatch.setattr(
            manager,
            "open_url",
            lambda req, **kw: Response(req.full_url.rsplit("/", 1)[-1].encode()),
        )

        def run(label):
            return manager.download("https://example.invalid/" + label, db)
    else:
        raise AssertionError(kind)

    def read():
        if kind == "download":
            return db.read_text()
        with contextlib.closing(sqlite3.connect(db)) as conn:
            assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)
            column = "remote_signal" if kind in ("clinvar", "gwas") else "version"
            stamp = conn.execute(f"SELECT {column},record_count FROM database_versions").fetchone()
            assert stamp is not None and stamp[1] > 0
            return stamp[0]

    return run, read, db


@pytest.mark.parametrize("kind", ["gwas", "pharmgkb", "clinvar", "prebuilt", "download"])
@pytest.mark.parametrize("outer", ["A", "B"])
@pytest.mark.parametrize(
    "inner_failure",
    [None, OSError, KeyboardInterrupt],
    ids=["both-success", "read-error", "interrupt"],
)
def test_one_writer_cannot_publish_or_delete_another_writers_staging(
    tmp_path, monkeypatch, kind, outer, inner_failure
):
    run, read, db = harness(tmp_path / "run", monkeypatch, kind)
    run("seed")
    real_replace = os.replace
    nested = False
    inner = "B" if outer == "A" else "A"
    inside = False
    inner_errors = []

    def replace(src, dst):
        nonlocal nested, inside
        if Path(dst) == db and inside and inner_failure:
            raise inner_failure("synthetic failed publication")
        if Path(dst) == db and not nested:
            nested = True
            inside = True
            try:
                run(inner)
            except BaseException as exc:
                inner_errors.append(exc)
            finally:
                inside = False
            assert read() == ("seed" if inner_failure else inner)
        return real_replace(src, dst)

    monkeypatch.setattr(os, "replace", replace)
    error = None
    try:
        run(outer)
    except BaseException as exc:
        error = exc
    assert nested, "must have reached the real publication boundary"
    assert error is None, f"outer writer lost its own completed bytes: {error!r}"
    assert len(inner_errors) == int(inner_failure is not None)
    if inner_failure:
        assert isinstance(inner_errors[0], inner_failure)
    assert read() == outer, "last successful publisher must publish its own complete cache"


@pytest.mark.parametrize("kind", ["gwas", "pharmgkb", "clinvar", "prebuilt", "download"])
def test_interrupt_preserves_destination_and_cleans_only_owned_staging(
    tmp_path, monkeypatch, kind
):
    run, read, db = harness(tmp_path / "run", monkeypatch, kind)
    run("seed")
    old = db.read_bytes()
    before = {p.relative_to(tmp_path) for p in tmp_path.rglob("*")}
    real_replace = os.replace

    def interrupt(src, dst):
        if Path(dst) == db:
            raise KeyboardInterrupt("synthetic SIGINT at publication boundary")
        return real_replace(src, dst)

    monkeypatch.setattr(os, "replace", interrupt)
    with pytest.raises(KeyboardInterrupt):
        run("A")
    assert db.read_bytes() == old and read() == "seed"
    assert {p.relative_to(tmp_path) for p in tmp_path.rglob("*")} == before


@pytest.mark.parametrize("failure", [False, True])
def test_offline_identity_builder_does_not_destroy_another_published_cache(tmp_path, failure):
    db = tmp_path / "identity.sqlite"
    kwargs = dict(version="synthetic", source_url="test://synthetic", downloaded_at="test")
    a = ("rs900000001", "1", 100, "G", "A", 200, "G", "A", 0)
    b = ("rs900000002", "1", 101, "G", "A", 201, "G", "A", 0)
    write_identity_db([b], db, **kwargs)

    def rows():
        write_identity_db([b], db, **kwargs)
        if failure:
            raise OSError("synthetic source failure")
        yield a

    error = None
    try:
        write_identity_db(rows(), db, **kwargs)
    except OSError as exc:
        error = exc
    except sqlite3.Error as exc:
        error = exc
    if failure:
        assert isinstance(error, OSError)
    else:
        assert error is None, f"builder lost its own output: {error!r}"
    with contextlib.closing(sqlite3.connect(db)) as conn:
        assert conn.execute("SELECT rsid FROM rsid_identity").fetchall() == [
            (b[0] if failure else a[0],)
        ]
        assert conn.execute("SELECT record_count FROM database_versions").fetchone() == (1,)

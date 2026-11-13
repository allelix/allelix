"""Ordinary update reports the retained release's skips, including local upgrades."""

import contextlib
import importlib
import sqlite3

import pytest
from click.testing import CliRunner
from test_gwas_risk_owner_disclosure import _assert_counts, _release

from allelix.annotators import gwas as module
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.cli import main
from allelix.databases.gwas_loader import gwas_source_hash, load_gwas_tsv, retained_source_stamp

SIGNAL = "etag:synthetic-retained-release"


def _cache(data, *, stale, other=2, unattributed=1):
    data.mkdir()
    source = _release(data / "gwas_catalog_associations.tsv", other, unattributed)
    digest = gwas_source_hash(source)
    target = data / "gwas.sqlite"
    load_gwas_tsv(source, target, remote_signal=SIGNAL, content_hash=digest)
    if stale:
        # A released cv:3 cache has no skip metadata and may hold the wrong
        # lowest-p winner. The retained table is this cache's own release.
        with contextlib.closing(sqlite3.connect(target)) as conn, conn:
            conn.execute("UPDATE database_versions SET local_version_tag='cv:3'")
            conn.execute("DROP TABLE gwas_source_skips")
            conn.execute("UPDATE gwas_associations SET risk_allele='G'")
    return digest


@pytest.fixture
def offline_gwas(monkeypatch):
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db.console, "_width", 500)
    monkeypatch.setattr(db, "get_annotators", lambda data, **kwargs: [GWASCatalogAnnotator(data)])
    monkeypatch.setattr(module, "GWAS_MIN_ROWS", 0)

    def no_download(*args, **kwargs):
        pytest.fail("The retained release needs no download")

    monkeypatch.setattr(module, "download", no_download)


@pytest.mark.parametrize("stale", (False, True), ids=("current", "cv3-local-upgrade"))
@pytest.mark.parametrize("remote", (SIGNAL, None), ids=("same-signal", "unavailable-signal"))
@pytest.mark.parametrize("other,unattributed", ((2, 1), (0, 0)))
def test_ordinary_update_discloses_the_retained_release_counts(
    tmp_path, monkeypatch, offline_gwas, stale, remote, other, unattributed
):
    data = tmp_path / "data"
    digest = _cache(data, stale=stale, other=other, unattributed=unattributed)
    monkeypatch.setattr(GWASCatalogAnnotator, "fetch_remote_signal", lambda self: remote)
    result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
    assert result.exit_code == 0, (result.output, result.exception)
    with GWASCatalogAnnotator(data) as ann:
        assert ann.is_ready(), "A retained cv:3 release must be re-ingested first"
    assert retained_source_stamp(data / "gwas.sqlite") == (SIGNAL, digest)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        assert conn.execute("SELECT risk_allele FROM gwas_associations").fetchall() == [("A",)]
    # The same persisted numbers are already available through status.
    status = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])
    assert status.exit_code == 0, (status.output, status.exception)
    _assert_counts(status.output, other, unattributed)
    _assert_counts(result.output, other, unattributed)


@pytest.mark.parametrize("stale", (False, True))
def test_only_another_source_does_not_disclose_or_rebuild_gwas(
    tmp_path, monkeypatch, offline_gwas, stale
):
    data = tmp_path / "data"
    _cache(data, stale=stale)
    before = (data / "gwas.sqlite").read_bytes()

    def no_probe(self):
        pytest.fail("--only clinvar must not probe GWAS")

    monkeypatch.setattr(GWASCatalogAnnotator, "fetch_remote_signal", no_probe)
    result = CliRunner().invoke(
        main, ["db", "update", "--only", "clinvar", "--data-dir", str(data)]
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert (data / "gwas.sqlite").read_bytes() == before
    assert "other_snp" not in result.output and "unattributed" not in result.output
    assert "GWAS Catalog" not in result.output

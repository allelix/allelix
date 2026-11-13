"""Risk-attribution skips belong to the source release, not the winning rows."""

import importlib
import re
import shutil
import zipfile

import pytest
from click.testing import CliRunner
from test_gwas_risk_owner import OTHER, RS, _row, _tsv

from allelix.annotators import gwas as module
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.cli import main
from allelix.databases.gwas_loader import load_gwas_tsv


def _release(path, other, unattributed):
    # Repeated invalid rows share the same dedup key; every rejected source
    # row is counted, not only the former lowest-p winner.
    return _tsv(
        path,
        [_row(OTHER + "-A", p="1e-50")] * other
        + [_row("unattributed-A", p="1e-40")] * unattributed
        + [_row(RS + "-A", p="1e-20")],
    )


def _assert_counts(text, other, unattributed):
    text = " ".join(text.replace(",", "").split())
    for label, count in ((r"other[ _-]+snp", other), (r"unattributed", unattributed)):
        pattern = rf"(?:\b{label}\s*[:=]?\s*{count}\b|\b{count}\s+{label}\b)"
        assert re.search(pattern, text, re.I), (
            f"Show source-release skip count {label}={count}, separate from loaded rows: {text}"
        )


@pytest.fixture
def only_gwas(monkeypatch):
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db.console, "_width", 500)
    monkeypatch.setattr(db, "get_annotators", lambda data, **kwargs: [GWASCatalogAnnotator(data)])
    monkeypatch.setattr(module, "GWAS_MIN_ROWS", 0)
    return db


@pytest.mark.parametrize("other,unattributed", ((2, 1), (1, 3)))
def test_db_status_reads_persisted_source_skip_counts_without_the_table(
    tmp_path, only_gwas, other, unattributed
):
    data = tmp_path / "data"
    data.mkdir()
    source = _release(tmp_path / "source.tsv", other, unattributed)
    load_gwas_tsv(source, data / "gwas.sqlite", remote_signal="etag:synthetic")
    source.unlink()
    before = (data / "gwas.sqlite").read_bytes()
    result = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])
    assert result.exit_code == 0, (result.output, result.exception)
    assert (data / "gwas.sqlite").read_bytes() == before, "Status must not rewrite the metadata"
    _assert_counts(result.output, other, unattributed)


def test_db_update_reports_the_new_release_counts_and_does_not_accumulate(
    tmp_path, monkeypatch, only_gwas
):
    data = tmp_path / "data"
    data.mkdir()
    active = [None]
    signal = [""]

    def download_local(url, destination, **kwargs):
        shutil.copyfile(active[0], destination)

    monkeypatch.setattr(module, "download", download_local)
    monkeypatch.setattr(GWASCatalogAnnotator, "fetch_remote_signal", lambda self: signal[0])
    observed = []
    for number, (other, unattributed) in enumerate(((2, 1), (1, 3)), start=1):
        source = _release(tmp_path / f"source{number}.tsv", other, unattributed)
        archive = tmp_path / f"source{number}.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.write(source, "gwas.tsv")
        active[0] = archive
        signal[0] = f"etag:synthetic{number}"
        result = CliRunner().invoke(main, ["db", "update", "--force", "--data-dir", str(data)])
        assert result.exit_code == 0, (result.output, result.exception)
        observed.append((result.output, other, unattributed))
        status = CliRunner().invoke(main, ["db", "status", "--data-dir", str(data)])
        assert status.exit_code == 0, (status.output, status.exception)
        observed.append((status.output, other, unattributed))
    for output, other, unattributed in observed:
        _assert_counts(output, other, unattributed)

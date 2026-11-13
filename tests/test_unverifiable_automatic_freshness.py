"""An unavailable freshness signal is not evidence that the user is offline."""

import os
import time
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _pharm

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.cli import main
from allelix.databases.gwas_loader import load_gwas_tsv

FIXTURES_DIR = Path(__file__).parent / "fixtures"

CLASSES = {
    "clinvar": ClinVarAnnotator,
    "gwas": GWASCatalogAnnotator,
    "pharmgkb": PharmGKBAnnotator,
}


def _data(path):
    path.mkdir()
    for build in ("GRCh37", "GRCh38"):
        _clinvar(
            path,
            build,
            [("rs1801133", "1", 1000, "G", "A", "SYN", "Pathogenic", "Pathogenic", 1, 1)],
        )
    _pharm(path, "rs1801133", "SYN", ["AG"])
    load_gwas_tsv(
        FIXTURES_DIR / "mock_gwas_catalog.tsv", path / "gwas.sqlite", source_url="test://synthetic"
    )
    (path / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=true\ngwas=true\nsnpedia=false\ngnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    return path


def _run(tmp_path, monkeypatch, command, source, *, stale=True, no_update=False):
    data = _data(tmp_path / "data")
    calls = []
    for name, cls in CLASSES.items():

        def signal(self, name=name):
            calls.append(name)
            return None

        monkeypatch.setattr(cls, "fetch_remote_signal", signal)

        def setup(self, **kwargs):
            pytest.fail("Unverifiable freshness must keep the cache without downloading")

        monkeypatch.setattr(cls, "setup", setup)
    if stale:
        old = time.time() - 30 * 86400
        for p in data.glob(source + "*.sqlite"):
            os.utime(p, (old, old))
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in data.glob("*.sqlite")}
    args = [
        command,
        str(FIXTURES_DIR / "mock_23andme.txt"),
        "--data-dir",
        str(data),
        "--build",
        "grch37",
        "--output",
        str(tmp_path / "report.json"),
    ]
    if no_update:
        args.append("--no-update")
    run = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    assert {
        p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in data.glob("*.sqlite")
    } == before
    return " ".join(run.output.lower().split()), calls


@pytest.mark.parametrize(
    "source,command",
    [
        (source, command)
        for source in ("clinvar", "gwas", "pharmgkb")
        for command in ("analyze", "methylation", "pharmacogenomics")
        if source != "gwas" or command == "analyze"
    ],
)
def test_each_auto_refresh_command_describes_unverifiable_source(
    tmp_path, monkeypatch, command, source
):
    output, calls = _run(tmp_path, monkeypatch, command, source)
    assert calls == [source], "The selected stale source must really be probed"
    assert "when online" not in output
    assert "verif" in output and "allelix db update" in output
    assert "cache" in output and ("retry" in output or "cached data is used" in output)


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("mode", ["fresh", "no_update"])
def test_fresh_or_offline_mode_does_not_probe_or_warn(tmp_path, monkeypatch, command, mode):
    output, calls = _run(
        tmp_path,
        monkeypatch,
        command,
        "gwas",
        stale=mode != "fresh",
        no_update=mode == "no_update",
    )
    assert not calls
    assert "when online" not in output and "freshness can't be verified" not in output


@pytest.mark.parametrize("command", ["methylation", "pharmacogenomics"])
def test_default_focused_command_does_not_probe_unused_gwas(tmp_path, monkeypatch, command):
    output, calls = _run(tmp_path, monkeypatch, command, "gwas")
    assert not calls, "A stale cache of an unused source must not trigger freshness traffic"
    assert "when online" not in output and "freshness can't be verified" not in output

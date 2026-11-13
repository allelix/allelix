"""New personal reports are owner-only; existing modes and nonreports stay intact."""

import contextlib
import os
import stat

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache, _write

from allelix.cli import main
from allelix.config import AllelixConfig, save_config
from allelix.exporters.plink import plink_paths
from allelix.reports import staged_output


def _mode(path):
    return stat.S_IMODE(path.stat().st_mode)


@contextlib.contextmanager
def _umask(value):
    previous = os.umask(value)
    try:
        yield
    finally:
        os.umask(previous)


@pytest.mark.parametrize("writer", ["json", "html", "atomic-text", "vcf", "plink"])
@pytest.mark.parametrize("mask", [0o000, 0o022, 0o077])
@pytest.mark.parametrize("prior", [None, 0o444, 0o600, 0o644, 0o666])
def test_report_modes_do_not_depend_on_extension_or_permissive_umask(
    tmp_path, writer, mask, prior
):
    source = tmp_path / "input.vcf"
    source.write_text(VCF)
    source.chmod(0o640)
    destination = tmp_path / "personal.output"
    outputs = list(plink_paths(destination).values()) if writer == "plink" else [destination]
    if prior is not None:
        for output in outputs:
            output.write_text("old report")
            output.chmod(prior)
    with _umask(mask):
        _write(writer, source, destination)
    expected = 0o600 if prior is None else prior
    assert [_mode(p) for p in outputs] == [expected] * len(outputs)
    assert all(p.read_bytes() != b"old report" for p in outputs)
    assert source.read_text() == VCF and _mode(source) == 0o640
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
@pytest.mark.parametrize("prior", [None, 0o640])
def test_report_and_vcf_out_use_the_same_privacy_rule(tmp_path, command, prior):
    source, data = tmp_path / "input.vcf", tmp_path / "data"
    source.write_text(VCF)
    _ready_cache(data)
    report, vcf = tmp_path / "report.json", tmp_path / "annotated.vcf"
    if prior is not None:
        for path in (report, vcf):
            path.write_text("old report")
            path.chmod(prior)
    args = [
        command,
        str(source),
        "--no-update",
        "--build",
        "grch37",
        "--data-dir",
        str(data),
        "--output",
        str(report),
        "--vcf-out",
        str(vcf),
    ]
    # --vcf-out belongs to analyze only; focused commands still own reports.
    if command != "analyze":
        args = args[:-2]
    with _umask(0o022):
        run = CliRunner().invoke(main, args)
    if command == "analyze" and prior is not None:
        assert run.exit_code != 0 and "refuses to overwrite" in run.output
        for path in (report, vcf):
            assert path.read_text() == "old report" and _mode(path) == prior
        return
    assert run.exit_code == 0, (run.output, run.exception)
    assert _mode(report) == (0o600 if prior is None else prior)
    if command == "analyze":
        assert _mode(vcf) == (0o600 if prior is None else prior)
        assert vcf.read_text().startswith("##fileformat=VCF")


@pytest.mark.parametrize("explicit_prefix", [False, True])
def test_cli_plink_all_three_new_files_are_private(tmp_path, explicit_prefix):
    source, data = tmp_path / "input.vcf", tmp_path / "data"
    source.write_text(VCF)
    _ready_cache(data)
    args = ["export", "plink", str(source), "--build", "grch37", "--data-dir", str(data)]
    if explicit_prefix:
        args += ["--output-prefix", str(tmp_path / "chosen")]
    with _umask(0o022):
        run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    outputs = [p for p in tmp_path.iterdir() if p.suffix in (".bed", ".bim", ".fam")]
    assert len(outputs) == 3
    assert {_mode(p) for p in outputs} == {0o600}


@pytest.mark.parametrize("prior", [None, 0o640])
def test_staging_is_private_and_failure_keeps_prior_bytes_and_mode(tmp_path, prior):
    dest = tmp_path / "report.txt"
    if prior is not None:
        dest.write_text("old report")
        dest.chmod(prior)
    with (
        _umask(0o022),
        pytest.raises(RuntimeError, match="stop before publish"),
        staged_output(dest) as stage,
    ):
        with stage.open_text() as stream:
            stream.write("private report data")
            stream.flush()
            assert _mode(stage.temporary) == 0o600
        raise RuntimeError("stop before publish")
    if prior is None:
        assert not dest.exists()
    else:
        assert dest.read_text() == "old report" and _mode(dest) == prior
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("mask", [0o022, 0o077])
@pytest.mark.parametrize("prior", [None, 0o640])
def test_config_is_not_a_report_and_keeps_its_existing_mode_contract(tmp_path, mask, prior):
    path = tmp_path / "config.toml"
    if prior is not None:
        path.write_text("old config")
        path.chmod(prior)
    with _umask(mask):
        save_config(tmp_path, AllelixConfig())
    assert _mode(path) == (prior if prior is not None else 0o666 & ~mask)
    assert "[sources]" in path.read_text()


@pytest.mark.parametrize("mask", [0o022, 0o077])
def test_cache_staging_keeps_its_nonreport_mode_policy(tmp_path, mask):
    from allelix.databases.staging import owned_staging

    target = tmp_path / "synthetic-cache.sqlite"
    with _umask(mask), owned_staging(target) as stage:
        stage.write_bytes(b"synthetic public cache bytes")
        os.replace(stage, target)
    assert _mode(target) == 0o666 & ~mask

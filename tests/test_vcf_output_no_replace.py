"""An annotated VCF never claims a destination another writer already owns."""

import os
import stat
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache

from allelix.cli import _helpers, main
from allelix.reports import StagedOutput, atomic_write_text
from allelix.reports.vcf import render_vcf

KINDS = ("file", "dangling", "symlink", "hardlink")
_LINK = os.link


def _occupy(path, kind):
    target = path.with_name("other-owner")
    if kind == "file":
        path.write_bytes(b"another writer's complete file\n")
        path.chmod(0o640)
    elif kind == "dangling":
        path.symlink_to(target.name)
    else:
        target.write_bytes(b"another writer's linked file\n")
        target.chmod(0o640)
        if kind == "symlink":
            path.symlink_to(target.name)
        else:
            _LINK(target, path)
    return _entry(path), _entry(target)


def _entry(path):
    if not os.path.lexists(path):
        return None
    info = path.lstat()
    content = os.readlink(path) if path.is_symlink() else path.read_bytes()
    return info.st_ino, stat.S_IMODE(info.st_mode), content


def _unchanged(path, before):
    assert (_entry(path), _entry(path.with_name("other-owner"))) == before
    assert not list(path.parent.glob(f".{path.name}.*.tmp")), "Owned staging must be cleaned"


def _write(source, output):
    return render_vcf(
        input_path=source,
        output_path=output,
        build="GRCh37",
        annotators_used=[],
        annotations_by_key={},
    )


def _race(monkeypatch, output, kind, boundary):
    saved = []

    def introduce():
        if not saved:
            assert not os.path.lexists(output)
            saved.append(_occupy(output, kind))

    if boundary == "publish":
        original = StagedOutput.publish

        def publish(self, *args, **kwargs):
            if self.destination == output:
                introduce()
            return original(self, *args, **kwargs)

        monkeypatch.setattr(StagedOutput, "publish", publish)
    else:
        # Insert the competing name at the final filesystem commit, after any
        # user-space existence check. A check followed by replace is not safe.
        for name in ("replace", "rename", "link"):
            original = getattr(os, name)

            def commit(src, dst, *args, _original=original, **kwargs):
                if Path(dst) == output:
                    introduce()
                return _original(src, dst, *args, **kwargs)

            monkeypatch.setattr(os, name, commit)
    return saved


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("boundary", ("publish", "filesystem-commit"))
def test_vcf_writer_does_not_replace_a_destination_created_during_writing(
    tmp_path, monkeypatch, kind, boundary
):
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(VCF)
    saved = _race(monkeypatch, output, kind, boundary)
    error = None
    try:
        _write(source, output)
    except (OSError, ValueError) as exc:
        error = exc
    assert saved, "The competing destination must be introduced during publication"
    _unchanged(output, saved[0])
    assert error is not None, "Refusal must reach the caller, not be reported as success"
    assert source.read_text() == VCF


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("explicit", (False, True), ids=("auto", "explicit"))
def test_cli_refuses_preexisting_names_including_dangling_links_before_analysis(
    tmp_path, monkeypatch, kind, explicit
):
    data = tmp_path / "data"
    _ready_cache(data)
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(
        VCF.replace(
            "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        )
    )
    before = _occupy(output, kind)
    calls = []
    original = _helpers.run_analysis

    def run(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)

    monkeypatch.setattr(_helpers, "run_analysis", run)
    args = [
        "analyze",
        str(source),
        "--no-update",
        "--data-dir",
        str(data),
        "--vcf-out",
        str(output),
    ]
    if explicit:
        args += ["--build", "grch37"]
    result = CliRunner().invoke(main, args)
    _unchanged(output, before)
    assert not calls, "lexists preflight must reject a dangling destination before analysis"
    assert result.exit_code != 0 and "Error" in result.output
    assert output.name in result.output


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("explicit", (False, True), ids=("auto", "explicit"))
def test_cli_reports_publication_race_as_failure_and_preserves_the_other_writer(
    tmp_path, monkeypatch, kind, explicit
):
    data = tmp_path / "data"
    _ready_cache(data)
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(
        VCF.replace(
            "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        )
    )
    saved = _race(monkeypatch, output, kind, "publish")
    args = [
        "analyze",
        str(source),
        "--no-update",
        "--data-dir",
        str(data),
        "--vcf-out",
        str(output),
    ]
    if explicit:
        args += ["--build", "grch37"]
    result = CliRunner().invoke(main, args)
    assert saved, (result.output, result.exception)
    _unchanged(output, saved[0])
    assert result.exit_code != 0, "A lost publication race must not report success"
    assert output.name in result.output and "Error" in result.output
    assert "Wrote annotated VCF" not in result.output


def test_unoccupied_vcf_is_complete_and_private(tmp_path):
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(VCF)
    assert _write(source, output) == 0
    text = output.read_text()
    assert "#CHROM\tPOS" in text and VCF.splitlines()[-1] in text
    assert stat.S_IMODE(output.stat().st_mode) == 0o600
    assert not list(tmp_path.glob(".result.vcf.*.tmp"))


@pytest.mark.parametrize("suffix", ("json", "html", "txt"))
def test_ordinary_report_replacement_keeps_its_existing_contract(tmp_path, suffix):
    output = tmp_path / f"report.{suffix}"
    output.write_text("old report")
    output.chmod(0o640)
    atomic_write_text(output, "new report")
    assert output.read_text() == "new report"
    assert stat.S_IMODE(output.stat().st_mode) == 0o640

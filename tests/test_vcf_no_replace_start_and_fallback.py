"""The entire CLI run and unsupported-hardlink paths retain no-replace safety."""

import errno
import os
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache
from test_vcf_output_no_replace import KINDS, _entry, _occupy, _unchanged, _write

from allelix.cli import _helpers, main


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("explicit", (False, True))
def test_name_created_during_analysis_is_not_an_existing_writer_api_destination(
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
    saved = []
    original = _helpers.run_analysis

    def analyze(*args, **kwargs):
        assert not os.path.lexists(output)
        saved.append(_occupy(output, kind))
        return original(*args, **kwargs)

    monkeypatch.setattr(_helpers, "run_analysis", analyze)
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--vcf-out",
        str(output),
    ]
    if explicit:
        args += ["--build", "grch37"]
    result = CliRunner().invoke(main, args)
    assert saved, (result.output, result.exception)
    _unchanged(output, saved[0])
    assert result.exit_code != 0 and "Error" in result.output
    assert output.name in result.output
    assert "Wrote annotated VCF" not in result.output


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("code", (errno.EPERM, errno.EOPNOTSUPP, errno.ENOSYS, errno.EMLINK))
def test_no_hardlinks_cannot_replace_a_competitor_after_exclusive_claim(
    tmp_path, monkeypatch, kind, code
):
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(VCF)
    original_link, original_open = os.link, os.open
    saved, attempted = [], []

    def link(src, dst, *args, **kwargs):
        if Path(dst) == output:
            attempted.append(True)
            raise OSError(code, "hard links unavailable", str(output))
        return original_link(src, dst, *args, **kwargs)

    def claim(path, flags, *args, **kwargs):
        fd = original_open(path, flags, *args, **kwargs)
        if Path(path) == output and flags & os.O_EXCL:
            # Another process removes the placeholder and installs its own
            # name after O_EXCL succeeds but before any later rename.
            output.unlink()
            saved.append(_occupy(output, kind))
        return fd

    monkeypatch.setattr(os, "link", link)
    monkeypatch.setattr(os, "open", claim)
    error = None
    try:
        _write(source, output)
    except (OSError, ValueError) as exc:
        error = exc
    assert attempted
    if saved:
        _unchanged(output, saved[0])
    else:
        # Refusing safely when atomic no-replace is unsupported is allowed.
        assert not os.path.lexists(output)
        assert not list(tmp_path.glob(".result.vcf.*.tmp"))
    assert error is not None, "Unsupported atomic publication must never claim unsafe success"
    assert output.name in str(error)


def test_failed_fallback_publication_leaves_no_empty_destination(tmp_path, monkeypatch):
    source, output = tmp_path / "input.vcf", tmp_path / "result.vcf"
    source.write_text(VCF)
    original_link, original_replace = os.link, os.replace

    def link(src, dst, *args, **kwargs):
        if Path(dst) == output:
            raise OSError(errno.EOPNOTSUPP, "hard links unavailable", str(output))
        return original_link(src, dst, *args, **kwargs)

    def replace(src, dst, *args, **kwargs):
        if Path(dst) == output:
            raise OSError(errno.EIO, "injected publication failure", str(output))
        return original_replace(src, dst, *args, **kwargs)

    monkeypatch.setattr(os, "link", link)
    monkeypatch.setattr(os, "replace", replace)
    with pytest.raises((OSError, ValueError)):
        _write(source, output)
    assert _entry(output) is None, "A failed write must not leave an empty final VCF"
    assert not list(tmp_path.glob(".result.vcf.*.tmp"))

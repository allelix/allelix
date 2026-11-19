"""An unsupported atomic VCF destination is refused before doing analysis work."""

import errno
import os
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache

from allelix.cli import _helpers, main


@pytest.mark.parametrize("code", [errno.EPERM, errno.EOPNOTSUPP, errno.ENOSYS, errno.EMLINK])
@pytest.mark.parametrize("auto", [False, True])
@pytest.mark.parametrize("no_update", [False, True])
def test_unsupported_hardlinks_refuse_before_cache_or_analysis(
    tmp_path, monkeypatch, code, auto, no_update
):
    data = tmp_path / "data"
    _ready_cache(data)
    source = tmp_path / "input.vcf"
    source.write_text(
        VCF.replace(
            "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        )
    )
    before = source.read_bytes()
    dest = tmp_path / "drive"
    dest.mkdir()
    output, report = dest / "annotated.vcf", tmp_path / "report.json"
    work = []
    attempted = []
    real_link = os.link

    def link(src, dst, *args, **kwargs):
        if Path(dst).parent == dest:
            attempted.append((Path(src), Path(dst)))
            raise OSError(code, "hard links unsupported", str(dst))
        return real_link(src, dst, *args, **kwargs)

    monkeypatch.setattr(os, "link", link)
    # Keep the real analysis, but make every entry into expensive work visible.
    for name in ("_prepare_caches", "_maybe_refresh_databases", "run_analysis"):
        original = getattr(_helpers, name)

        def track(*args, _name=name, _original=original, **kwargs):
            work.append(_name)
            if _name != "run_analysis":
                return None  # entirely offline, no cache/network refresh
            return _original(*args, **kwargs)

        monkeypatch.setattr(_helpers, name, track)
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--build",
        "auto" if auto else "GRCh37",
        "--output",
        str(report),
        "--vcf-out",
        str(output),
    ]
    if no_update:
        args += ["--no-update"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 1, (result.output, result.exception)
    assert attempted, "Probe the filesystem that will hold the VCF, not another directory"
    assert not work, f"Unsupported destination was discovered only after {work}"
    assert "Error" in result.output and str(output) in result.output
    assert "link" in result.output.lower(), result.output
    assert not report.exists() and not os.path.lexists(output)
    assert list(dest.iterdir()) == [], "Capability probes and staged outputs must be cleaned"
    assert source.read_bytes() == before


@pytest.mark.parametrize("kind", ["plain-report", "supported-vcf", "existing", "dangling", "gzip"])
def test_hardlink_preflight_keeps_other_output_contracts(tmp_path, monkeypatch, kind):
    data = tmp_path / "data"
    _ready_cache(data)
    source = tmp_path / "input.vcf"
    source.write_text(VCF)
    report = tmp_path / "report.json"
    output = tmp_path / ("out.vcf.gz" if kind == "gzip" else "out.vcf")
    if kind == "existing":
        output.write_bytes(b"keep me")
    if kind == "dangling":
        output.symlink_to(tmp_path / "absent")
    if kind == "plain-report":

        def unavailable(*args, **kwargs):
            raise AssertionError("Ordinary JSON output must not require hard links")

        monkeypatch.setattr(os, "link", unavailable)
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        "GRCh37",
        "--output",
        str(report),
    ]
    if kind != "plain-report":
        args += ["--vcf-out", str(output)]
    result = CliRunner().invoke(main, args)
    if kind in ("plain-report", "supported-vcf"):
        assert result.exit_code == 0, (result.output, result.exception)
        assert report.exists()
        if kind == "supported-vcf":
            assert output.is_file()
            assert not list(tmp_path.glob(".out.vcf.*"))
    else:
        assert result.exit_code == 1, (result.output, result.exception)
        assert not report.exists()
        if kind == "existing":
            assert output.read_bytes() == b"keep me"
        elif kind == "dangling":
            assert output.is_symlink() and os.readlink(output) == str(tmp_path / "absent")
        else:
            assert ".gz" in result.output and not output.exists()

"""A plain-text VCF writer must not label its bytes as gzip/BGZF."""

import gzip

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache

from allelix.cli import _helpers, main
from allelix.reports.vcf import render_vcf


def _source(tmp_path, compressed):
    source = tmp_path / ("input.vcf.gz" if compressed else "input.vcf")
    text = VCF.replace(
        "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
    )
    if compressed:
        with gzip.open(source, "wt") as out:
            out.write(text)
    else:
        source.write_text(text)
    return source


def _write(source, output):
    return render_vcf(
        input_path=source,
        output_path=output,
        build="GRCh37",
        annotators_used=[],
        annotations_by_key={},
    )


def _clear_error(text):
    text = text.lower()
    assert "vcf" in text
    assert "compress" in text or "gzip" in text or "bgzf" in text


@pytest.mark.parametrize("compressed", (False, True), ids=("plain-input", "gzip-input"))
@pytest.mark.parametrize("name", ("report.vcf.gz", "report.VCF.GZ"))
def test_writer_refuses_compressed_name_without_creating_any_output(tmp_path, compressed, name):
    source = _source(tmp_path, compressed)
    output = tmp_path / name
    before = source.read_bytes()
    error = None
    try:
        _write(source, output)
    except (ValueError, OSError) as exc:
        error = exc
    assert error is not None, "Plain VCF bytes must not be published under a gzip filename"
    _clear_error(str(error))
    assert not output.exists()
    assert not list(tmp_path.glob(f".{name}.*.tmp"))
    assert source.read_bytes() == before


@pytest.mark.parametrize("compressed", (False, True), ids=("plain-input", "gzip-input"))
@pytest.mark.parametrize("explicit", (False, True), ids=("auto", "explicit"))
@pytest.mark.parametrize("name", ("report.vcf.gz", "report.VCF.GZ"))
def test_cli_refuses_before_analysis_or_report_publication(
    tmp_path, monkeypatch, compressed, explicit, name
):
    source = _source(tmp_path, compressed)
    data = tmp_path / "data"
    _ready_cache(data)
    output, report = tmp_path / name, tmp_path / "report.json"
    calls = []
    original = _helpers.run_analysis

    def analyze(*args, **kwargs):
        calls.append(True)
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
        "--output",
        str(report),
    ]
    if explicit:
        args += ["--build", "grch37"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code != 0, result.output
    assert "Error" in result.output
    _clear_error(result.output)
    assert name in result.output
    assert not calls, "An invalid output format must be rejected before analysis"
    assert not output.exists() and not report.exists()


@pytest.mark.parametrize("compressed", (False, True), ids=("plain-input", "gzip-input"))
@pytest.mark.parametrize("name", ("report.vcf", "report.txt", "report"))
def test_plain_output_still_accepts_plain_or_compressed_input(tmp_path, compressed, name):
    source = _source(tmp_path, compressed)
    output = tmp_path / name
    _write(source, output)
    assert output.read_bytes().startswith(b"##fileformat=VCF")
    assert "rs900252401" in output.read_text()


def test_help_states_the_plain_output_and_compressed_name_refusal():
    result = CliRunner().invoke(main, ["analyze", "--help"])
    assert result.exit_code == 0
    text = " ".join(result.output.split()).lower()
    start = text.index("--vcf-out")
    end = text.find(" --", start + len("--vcf-out"))
    option = text[start:end] if end != -1 else text[start:]
    assert "uncompressed" in option or "plain-text" in option
    assert ".gz" in option and ("refus" in option or "reject" in option)

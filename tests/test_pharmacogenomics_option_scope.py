"""Only pharmacogenomics rejects options that cannot apply to its report."""

import re

import pytest
from click.testing import CliRunner

from allelix.cli import _helpers, main

FLAGS = [("--include-gwas",), ("--gwas-all",), ("--exclude-snpedia",)]
VALUES = [
    (name, value)
    for name, default in (("--gwas-min-magnitude", "9"), ("--snpedia-min-magnitude", "2"))
    for value in ("0", default, "4.5")
]
OPTIONS = FLAGS + VALUES + [(f"{name}={value}",) for name, value in VALUES]


@pytest.fixture
def dispatch(tmp_path, monkeypatch):
    source = tmp_path / "input.vcf"
    source.write_text("##fileformat=VCFv4.2\n")
    calls = []
    monkeypatch.setattr(_helpers, "_run_analysis_command", lambda **kw: calls.append(kw))
    return source, calls


def _assert_usage(run, option):
    assert run.exit_code == 2, (run.output, run.exception)
    assert option.split("=")[0] in run.output
    # A generic unknown-option error alone does not explain the breaking change.
    assert re.search(
        r"never applies?|does not apply|not applicable|cannot apply", run.output, re.I
    )


@pytest.mark.parametrize("option", OPTIONS)
@pytest.mark.parametrize("before_file", [False, True])
def test_pharma_refuses_explicit_options_even_when_the_value_equals_the_default(
    dispatch, tmp_path, option, before_file
):
    source, calls = dispatch
    output = tmp_path / "report.json"
    args = [*option, str(source)] if before_file else [str(source), *option]
    run = CliRunner().invoke(main, ["pharmacogenomics", *args, "--output", str(output)])
    _assert_usage(run, option[0])
    assert calls == [], "Refuse before cache refresh, input analysis or report writing"
    assert not output.exists()


@pytest.mark.parametrize(
    "options",
    [
        (),
        ("--no-update",),
        ("--min-magnitude", "0"),
        ("--include-benign", "--no-gnomad", "--show-unchecked"),
    ],
)
def test_pharma_implicit_defaults_and_unrelated_options_still_run(dispatch, options):
    source, calls = dispatch
    run = CliRunner().invoke(main, ["pharmacogenomics", str(source), *options])
    assert run.exit_code == 0, run.output
    assert len(calls) == 1 and calls[0]["category"] == "pharma"


@pytest.mark.parametrize("command", ["analyze", "methylation"])
@pytest.mark.parametrize("option", OPTIONS)
def test_other_report_commands_keep_the_released_option_interface(dispatch, command, option):
    source, calls = dispatch
    run = CliRunner().invoke(main, [command, str(source), *option])
    if command == "analyze" and option == ("--include-gwas",):
        # SUP clarification: this was never an analyze option; do not add it.
        assert run.exit_code == 2 and "--include-gwas" in run.output
        assert calls == []
        return
    assert run.exit_code == 0, run.output
    assert len(calls) == 1
    name = option[0].split("=")[0]
    passed = calls[0]
    if name == "--gwas-all":
        assert passed["gwas_all"] is True
    elif name == "--exclude-snpedia":
        assert "snpedia" in passed["exclude_sources"]
    elif name == "--include-gwas":
        assert "gwas" not in (passed["exclude_sources"] or ())
    else:
        value = option[1] if len(option) == 2 else option[0].split("=", 1)[1]
        assert passed[name[2:].replace("-", "_")] == float(value)


def test_analyze_still_includes_gwas_by_default(dispatch):
    source, calls = dispatch
    run = CliRunner().invoke(main, ["analyze", str(source)])
    assert run.exit_code == 0, run.output
    assert "gwas" not in (calls[0]["exclude_sources"] or ())
    assert calls[0]["gwas_min_magnitude"] == 9.0
    assert calls[0]["snpedia_min_magnitude"] == 2.0


@pytest.mark.parametrize("option", FLAGS + VALUES)
def test_rejection_also_happens_with_real_input_and_offline_ready_cache(tmp_path, option):
    from test_output_destinations_and_modes import VCF, _ready_cache

    source, data = tmp_path / "input.vcf", tmp_path / "data"
    source.write_text(VCF)
    _ready_cache(data)
    output = tmp_path / "report.json"
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in data.iterdir()}
    run = CliRunner().invoke(
        main,
        [
            "pharmacogenomics",
            str(source),
            "--no-update",
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--output",
            str(output),
            *option,
        ],
    )
    _assert_usage(run, option[0])
    assert not output.exists()
    assert {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in data.iterdir()} == before

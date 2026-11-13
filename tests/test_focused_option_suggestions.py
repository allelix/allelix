"""A refused option is not a useful spelling suggestion for that command."""

import re

import pytest
from click.testing import CliRunner
from test_publish_release_gates import ROOT

from allelix.cli import main

OPTIONS = (
    ("--include-gwas", "--include-gw", None),
    ("--gwas-all", "--gwas-al", None),
    ("--gwas-min-magnitude", "--gwas-min-magnitud", "0"),
    ("--exclude-snpedia", "--exclude-snp", None),
    ("--snpedia-min-magnitude", "--snpedia-min-magnitud", "0"),
)


@pytest.mark.parametrize("option,typo,value", OPTIONS)
def test_pharmacogenomics_does_not_recommend_a_refused_option(option, typo, value):
    result = CliRunner().invoke(main, ["pharmacogenomics", typo])
    assert result.exit_code == 2, result.output
    assert "No such option" in result.output and typo in result.output
    # The bad spelling is named; none of the five refused options may be suggested.
    suggestions = (
        result.output.split("Did you mean", 1)[-1] if "Did you mean" in result.output else ""
    )
    assert not any(re.search(re.escape(o) + r"(?![\w-])", suggestions) for o, _, _ in OPTIONS), (
        result.output
    )


@pytest.mark.parametrize("option,typo,value", OPTIONS)
def test_explicit_refused_option_still_names_it_before_input_read(option, typo, value):
    args = ["pharmacogenomics", option] + ([] if value is None else [value])
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 2 and option in result.output
    assert "never applies to pharmacogenomics" in result.output, result.output


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
def test_typo_for_a_supported_option_still_has_a_useful_hint(command):
    result = CliRunner().invoke(main, [command, "--min-magnitud"])
    assert result.exit_code == 2 and "Did you mean" in result.output
    assert "--min-magnitude" in result.output


@pytest.mark.parametrize("option,typo,value", OPTIONS)
def test_methylation_retains_its_existing_suggestions(option, typo, value):
    result = CliRunner().invoke(main, ["methylation", typo])
    assert result.exit_code == 2 and "Did you mean" in result.output
    assert option in result.output


@pytest.mark.parametrize("option,typo,value", OPTIONS[1:])
def test_analyze_retains_its_four_existing_suggestions(option, typo, value):
    result = CliRunner().invoke(main, ["analyze", typo])
    assert result.exit_code == 2 and "Did you mean" in result.output
    assert option in result.output


def test_adr_adds_a_current_focused_command_amendment_without_rewriting_history():
    text = (ROOT / "docs/adr/0024-gwas-magnitude-scoring.md").read_text()
    sections = re.split(r"^## ", text, flags=re.M)
    amendments = "\n".join(s for s in sections if s.startswith("Amendment 2.7.0"))
    assert amendments, (
        "Append a dated amendment; the old focused-GWAS paragraph is released history"
    )
    assert re.search(r"2026-\d{2}-\d{2}", amendments)
    assert "pharmacogenomics" in amendments and "--include-gwas" in amendments
    assert any(
        word in amendments.lower() for word in ("refus", "reject", "not accept", "unsupported")
    )
    assert "methylation" in amendments, "Do not imply that both focused commands refuse GWAS"

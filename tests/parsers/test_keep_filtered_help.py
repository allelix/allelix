"""--keep-filtered's help names both filters it overrides (ADR-0060)."""

from click.testing import CliRunner

from allelix.cli import main


def test_keep_filtered_help_names_site_filter_and_sample_ft():
    result = CliRunner(env={"COLUMNS": "200"}).invoke(main, ["analyze", "--help"])
    assert result.exit_code == 0
    text = " ".join(result.output.split())
    assert "keep rows whose FILTER, or whose selected sample's FT, is not PASS" in text

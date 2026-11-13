"""Filtering refused options must retain Click's context-added help options."""

import pytest
from click.testing import CliRunner

from allelix.cli import main


@pytest.mark.parametrize("command", ("analyze", "methylation", "pharmacogenomics"))
@pytest.mark.parametrize("typo", ("--hel", "--hepl"))
def test_supported_help_option_is_still_suggested(command, typo):
    result = CliRunner().invoke(main, [command, typo])
    assert result.exit_code == 2, result.output
    assert "No such option" in result.output and typo in result.output
    assert "Did you mean '--help'?" in result.output, result.output


@pytest.mark.parametrize("command", ("analyze", "methylation", "pharmacogenomics"))
def test_help_option_still_exits_before_reading_input(command):
    result = CliRunner().invoke(main, [command, "--help"])
    assert result.exit_code == 0, result.output
    assert "Usage:" in result.output and "Show this message and exit" in result.output

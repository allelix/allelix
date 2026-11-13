"""Report-command help describes the private new-file policy."""

import pytest
from click.testing import CliRunner

from allelix.cli import main


@pytest.mark.parametrize("command", ["analyze", "methylation", "pharmacogenomics"])
def test_output_help_says_new_report_files_are_private(command):
    result = CliRunner().invoke(main, [command, "--help"])
    assert result.exit_code == 0, result.output
    help_text = " ".join(result.output.split())
    assert "--output" in help_text
    assert "0600" in help_text, "--output help must state the new report-file mode"
    assert "a new file gets the umask's mode" not in help_text

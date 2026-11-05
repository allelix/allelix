# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Help text says what the commands do (2.6.0 documentation audit).

- ``--build`` on ``stats``, ``extract`` and ``export plink`` says an rsID and
  genotype list is refused there; the analysis commands say it is placed on
  GRCh38 (ADR-0069).
- ``--build auto``: confident known-SNP positions replace a contradicting label
  (ADR-0044); the help no longer says every contradiction stops the run, and
  it names the GRCh36 fail-safe and the vendor GRCh36 label (ADR-0025).
- ``--gwas-min-magnitude`` names the must-include exemption, ``--no-update``
  the whole update step, and ``config get``'s example the CADD default.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.config import _DEFAULT_SOURCES


def _flat_help(*args: str) -> str:
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, [*args, "--help"])
    assert result.exit_code == 0, result.output
    return " ".join(result.output.split())


PLACED = "it is placed on GRCh38, and any other build stops."
REFUSED = (
    "has no positions and is refused here; analyze, methylation and pharmacogenomics read it."
)


@pytest.mark.parametrize("command", [("analyze",), ("methylation",), ("pharmacogenomics",)])
def test_analysis_commands_say_an_rsid_list_is_placed(command: tuple[str, ...]) -> None:
    text = _flat_help(*command)
    assert PLACED in text
    assert REFUSED not in text


@pytest.mark.parametrize("command", [("stats",), ("extract",), ("export", "plink")])
def test_position_commands_say_an_rsid_list_is_refused(command: tuple[str, ...]) -> None:
    text = _flat_help(*command)
    assert REFUSED in text
    assert PLACED not in text


@pytest.mark.parametrize("command", [("analyze",), ("stats",), ("extract",), ("export", "plink")])
def test_auto_build_help_says_confident_positions_replace_a_label(
    command: tuple[str, ...],
) -> None:
    text = _flat_help(*command)
    assert "Confident known-SNP positions replace a file label they contradict" in text
    assert "and stops when that evidence contradicts the file's label" not in text


def test_gwas_floor_help_names_the_must_include_exemption() -> None:
    text = _flat_help("analyze")
    assert "except must-include rsIDs, which need only --min-magnitude" in text


def test_no_update_help_names_the_whole_update_step() -> None:
    text = _flat_help("analyze")
    assert "Skip the pre-analysis update step: no cache conversion, local rebuild" in text


def test_config_get_example_shows_the_cadd_default() -> None:
    assert _DEFAULT_SOURCES["cadd"] is False
    result = CliRunner(env={"COLUMNS": "100"}).invoke(main, ["config", "get", "--help"])
    assert "    allelix config get sources.cadd        # false" in result.output.splitlines()


@pytest.mark.parametrize("command", [("analyze",), ("stats",), ("extract",), ("export", "plink")])
def test_auto_build_help_names_the_grch36_exceptions(command: tuple[str, ...]) -> None:
    text = _flat_help(*command)
    assert (
        "Two GRCh36 exceptions: positions that point to GRCh36 without confidence, "
        "or a tie that includes GRCh36, take GRCh36 even against a modern label" in text
    )
    assert "an array file labelled GRCh36 runs on its label without positions" in text
    assert "(a VCF does not)" in text

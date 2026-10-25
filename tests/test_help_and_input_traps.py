# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Help text that renders as written, and two inputs that failed with no warning.

- ``config get``, ``config set`` and ``db path`` help keep their key lists and
  examples on separate lines (click's ``\\b`` marker works only in a plain
  docstring), and the ``db path`` example runs as shown.
- ``--category`` takes clinical, pharma or trait; a typo stops with the choices
  instead of an empty report with exit 0.
- ``--filter-file``: a ``#`` after an entry starts a comment, so
  ``rs1801133  # MTHFR`` keeps rs1801133.
- The GRCh36 notice says what still runs and what is not assessed.
- The Sequencing.com SV fixture carries a placeholder ``##fileDate``.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.cli._helpers import _parse_filter_file
from allelix.reports._pipeline import BuildDiagnostics
from allelix.reports.build_text import build_notices
from allelix.utils.build_detect import BUILD_GRCH36

FIXTURES = Path(__file__).parent / "fixtures"


def _help(*args: str) -> str:
    result = CliRunner(env={"COLUMNS": "100"}).invoke(main, [*args, "--help"])
    assert result.exit_code == 0, result.output
    return result.output


# ── Help rendering ────────────────────────────────────────────────────


@pytest.mark.parametrize(
    ("command", "lines"),
    [
        (
            ("config", "get"),
            [
                "  Keys:",
                "    sources.<name>       Show if a source is enabled",
                "    license.commercial   Show commercial mode",
                "  Examples:",
                "    allelix config get                     # dump entire config",
                "    allelix config get options.cadd_full   # false",
            ],
        ),
        (
            ("config", "set"),
            [
                "  Keys:",
                "    sources.<name>       Enable/disable a source (true/false)",
                "                         heterozygous_max_af (default 0.0)",
                "  Examples:",
                "    allelix config set sources.snpedia false",
                "    allelix config set plausibility.homozygous_max_af 0.00005",
            ],
        ),
        (
            ("db", "path"),
            ["      ALLELIX_DATA=$(allelix db path)", '      du -sh "$ALLELIX_DATA"'],
        ),
    ],
    ids=["config-get", "config-set", "db-path"],
)
def test_help_keeps_its_lines(command: tuple[str, ...], lines: list[str]) -> None:
    out = _help(*command).splitlines()
    for line in lines:
        assert line in out
    assert not any("\\b" in line for line in out)


def test_db_path_example_runs_as_shown(tmp_path: Path) -> None:
    out = _help("db", "path").splitlines()
    i = out.index("      ALLELIX_DATA=$(allelix db path)")
    script = "\n".join(line.strip() for line in out[i : i + 2])
    shown = script.replace("allelix db path", f"allelix db path --data-dir {tmp_path}")
    # The checkout's own entry point, else the one on PATH (CI installs there).
    venv_allelix = Path(__file__).resolve().parents[1] / ".venv" / "bin" / "allelix"
    allelix = str(venv_allelix) if venv_allelix.exists() else shutil.which("allelix")
    if allelix is None:
        pytest.skip("no allelix entry point in .venv or on PATH")
    shown = shown.replace("$(allelix ", f"$({allelix} ")
    result = subprocess.run(["sh", "-c", shown], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert str(tmp_path) in result.stdout


# ── --category ────────────────────────────────────────────────────────


def test_category_help_lists_the_three() -> None:
    out = " ".join(_help("analyze").split())
    assert "--category [clinical|pharma|trait]" in out
    assert "Show only one category: clinical, pharma or trait." in out


@pytest.mark.parametrize("typo", ["clinicals", "pharmacy", "traits", ""])
def test_category_typo_stops_with_the_choices(mock_mhg_path: Path, typo: str) -> None:
    result = CliRunner().invoke(main, ["analyze", str(mock_mhg_path), "--category", typo])
    assert result.exit_code == 2
    assert "Invalid value for '--category'" in result.output
    assert "'clinical', 'pharma', 'trait'" in result.output


@pytest.mark.parametrize(
    ("given", "passed"),
    [("clinical", "clinical"), ("Pharma", "pharma"), ("TRAIT", "trait")],
)
def test_category_choice_reaches_the_analysis(
    mock_mhg_path: Path, monkeypatch: pytest.MonkeyPatch, given: str, passed: str
) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        "allelix.cli._helpers._run_analysis_command", lambda **kw: captured.update(kw)
    )
    result = CliRunner().invoke(main, ["analyze", str(mock_mhg_path), "--category", given])
    assert result.exit_code == 0, result.output
    assert captured["category"] == passed


def test_no_category_still_means_all(mock_mhg_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(
        "allelix.cli._helpers._run_analysis_command", lambda **kw: captured.update(kw)
    )
    result = CliRunner().invoke(main, ["analyze", str(mock_mhg_path)])
    assert result.exit_code == 0, result.output
    assert captured["category"] is None


# ── --filter-file inline comments ─────────────────────────────────────


@pytest.mark.parametrize(
    "raw",
    [
        b"rs1801133  # MTHFR\nCOMT # gene\n",
        b"rs1801133\t#MTHFR\r\nCOMT#gene\r\n",
        b"\xef\xbb\xbfrs1801133  # MTHFR\nCOMT # gene\n",
        b"\xef\xbb\xbfrs1801133# MTHFR\r\nCOMT\t# gene\r\n",
        b"\xef\xbb\xbf# my panel\nrs1801133 # MTHFR\n\n   # indented note\nCOMT\n",
    ],
    ids=["lf", "crlf-tab", "bom-lf", "bom-crlf", "bom-comment-first"],
)
def test_inline_comment_keeps_the_entry(tmp_path: Path, raw: bytes) -> None:
    path = tmp_path / "filter.txt"
    path.write_bytes(raw)
    assert _parse_filter_file(path) == (frozenset({"COMT"}), frozenset({"rs1801133"}))


def test_comment_only_lines_give_nothing(tmp_path: Path) -> None:
    path = tmp_path / "filter.txt"
    path.write_bytes(b"\xef\xbb\xbf# header\n#rs4680\n  # rs1801133\n\n")
    assert _parse_filter_file(path) == (frozenset(), frozenset())


def test_filter_file_help_describes_inline_comments() -> None:
    out = " ".join(_help("analyze").split())
    assert "A # starts a comment, on its own line or after an entry" in out
    assert "Comments (#)" not in out


# ── GRCh36 notice ─────────────────────────────────────────────────────


def _grch36_notice() -> str:
    diag = BuildDiagnostics(
        effective_build=BUILD_GRCH36,
        source="file_label",
        detected_build=None,
        matched_count=0,
        inspected_count=0,
        position_confident=False,
        tied_builds=(),
        tied_votes=0,
        file_label=BUILD_GRCH36,
    )
    return next(n.text for n in build_notices(diag) if n.kind == "grch36")


def test_grch36_notice_says_what_runs_and_what_is_not_assessed() -> None:
    text = _grch36_notice()
    assert "ClinPGx matches by rsID and still runs" in text
    assert (
        "GWAS Catalog risk alleles and single-base SNPedia genotypes are not assessed "
        "on GRCh36 (counted below)"
    ) in text
    assert "GWAS Catalog and SNPedia match by rsID" not in text
    assert "ClinVar is skipped" in text


# ── Fixture header ────────────────────────────────────────────────────


def test_sv_fixture_file_date_is_a_placeholder() -> None:
    lines = (FIXTURES / "sequencing_com_sv.vcf").read_text().splitlines()
    assert lines[1] == "##fileDate=January-01-2000 at 00:00 UTC"

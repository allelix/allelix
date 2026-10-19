# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""--format is an expert override: allowed, warned and documented.

A file no reader recognises is rejected with a message that puts the
original file first. --format with a reader whose own detection does not
recognise the file reads it anyway, with a warning before anything else in
the terminal and on every report (HTML, JSON ``report_warnings`` code
``forced_format``). --format with a reader that does recognise the file is
not warned. Container checks run under --format as they do without it.

Synthetic files only: invented rows (rs3999990000 range), placeholder headers.
"""

from __future__ import annotations

import json
import zipfile
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import ParserNotFoundError, detect_parser

if TYPE_CHECKING:
    from pathlib import Path

REJECTION = (
    "Allelix doesn't recognise the format of 'headerless.txt'. Use the original file "
    "exactly as your testing company provided it, not edited or re-saved. If you're sure "
    "of its source, --format <name> reads it anyway, at your own risk: Allelix can't check "
    'a file it doesn\'t recognise, and the results may be wrong. See "Forcing a format" '
    "in docs/formats.md."
)
WARNING = (
    "'headerless.txt' doesn't look like a 23andMe file. Reading it as 23andMe because you "
    "chose --format 23andme. Allelix can't check its layout, strand or build the usual "
    "way, so the results may be wrong. Use the original file from your testing company if "
    "you can."
)
_ROWS = (
    "rs3999990701\t1\t1000101\tAG\nrs3999990702\t1\t1000202\tCC\nrs3999990703\t2\t2000303\tCT\n"
)


@pytest.fixture
def headerless(tmp_path: Path) -> Path:
    path = tmp_path / "headerless.txt"
    path.write_text(_ROWS)
    return path


def _run(*args: str):
    return CliRunner().invoke(main, list(args), env={"COLUMNS": "400"})


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_the_rejection_puts_the_original_file_first(headerless: Path) -> None:
    with pytest.raises(ParserNotFoundError) as exc:
        detect_parser(headerless)
    assert str(exc.value) == REJECTION
    result = _run("stats", str(headerless))
    assert result.exit_code != 0
    assert REJECTION in _flat(result.output)


class TestForcedMismatch:
    def test_analyze_warns_first_and_on_every_report(
        self, headerless: Path, clinvar_data_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "r.json"
        result = _run(
            "analyze", str(headerless), "--format", "23andme", "--build", "grch37", "--no-update",
            "--data-dir", str(clinvar_data_dir), "--output", str(out),
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        flat = _flat(result.output)
        assert flat.startswith(f"Warning: {WARNING}")
        warnings = json.loads(out.read_text())["report_warnings"]
        assert {"code": "forced_format", "message": WARNING, "count": None} in warnings

        html = tmp_path / "r.html"
        result = _run(
            "analyze", str(headerless), "--format", "23andme", "--build", "grch37", "--no-update",
            "--data-dir", str(clinvar_data_dir), "--output", str(html),
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        page = html.read_text()
        assert (
            '<div class="notice-warn"><strong>Warning.</strong> &#x27;headerless.txt&#x27; '
            "doesn&#x27;t look like a 23andMe file" in page
        )

    def test_the_terminal_report_carries_it(
        self, headerless: Path, clinvar_data_dir: Path
    ) -> None:
        result = _run(
            "analyze", str(headerless), "--format", "23andme", "--build", "grch37", "--no-update",
            "--data-dir", str(clinvar_data_dir),
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert _flat(result.output).count(f"Warning: {WARNING}") == 2

    @pytest.mark.parametrize(
        "command",
        [
            ["stats", "--build", "grch37"],
            ["extract", "--snps", "rs3999990701", "--build", "grch37"],
            ["export", "plink", "-o", "PLINK_OUT", "--build", "grch37"],
        ],
        ids=["stats", "extract", "export-plink"],
    )
    def test_every_command_with_format_warns(
        self, headerless: Path, tmp_path: Path, command: list[str]
    ) -> None:
        args = [a.replace("PLINK_OUT", str(tmp_path / "plink")) for a in command]
        result = _run(*args, "--format", "23andme", str(headerless))
        assert result.exit_code == 0, result.output
        assert _flat(result.output).startswith(f"Warning: {WARNING}")

    @pytest.mark.parametrize("command", ["methylation", "pharmacogenomics"])
    def test_focused_commands_warn(
        self, headerless: Path, clinvar_data_dir: Path, tmp_path: Path, command: str
    ) -> None:
        out = tmp_path / "r.json"
        result = _run(
            command, str(headerless), "--format", "23andme", "--build", "grch37", "--no-update",
            "--data-dir", str(clinvar_data_dir), "--output", str(out),
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert _flat(result.output).startswith(f"Warning: {WARNING}")
        codes = [w["code"] for w in json.loads(out.read_text())["report_warnings"]]
        assert "forced_format" in codes

    def test_compare_checks_each_file_on_its_own(
        self, headerless: Path, mock_23andme_path: Path
    ) -> None:
        result = _run(
            "compare", str(mock_23andme_path), str(headerless),
            "--format1", "23andme", "--format2", "23andme", "--build2", "grch37",
        )  # fmt: skip
        flat = _flat(result.output)
        assert flat.count("doesn't look like a 23andMe file") == 1
        assert "'headerless.txt' doesn't look like" in flat


class TestForcedMatch:
    def test_a_recognised_file_is_not_warned(
        self, mock_23andme_path: Path, clinvar_data_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "r.json"
        result = _run(
            "analyze", str(mock_23andme_path), "--format", "23andme", "--no-update",
            "--data-dir", str(clinvar_data_dir), "--output", str(out),
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert "doesn't look like" not in result.output
        codes = [w["code"] for w in json.loads(out.read_text())["report_warnings"]]
        assert "forced_format" not in codes


class TestContainersUnderFormat:
    @pytest.fixture
    def zips(self, tmp_path: Path, mock_23andme_path: Path) -> dict[str, Path]:
        two = tmp_path / "two.zip"
        with zipfile.ZipFile(two, "w") as z:
            z.write(mock_23andme_path, "a.txt")
            z.write(mock_23andme_path, "b.txt")
        one = tmp_path / "one.zip"
        with zipfile.ZipFile(one, "w") as z:
            z.write(mock_23andme_path, "a.txt")
        data = one.read_bytes()
        damaged = tmp_path / "damaged.zip"
        damaged.write_bytes(data[: len(data) // 2])
        flagged = bytearray(data)
        for signature, offset in ((b"PK\x03\x04", 6), (b"PK\x01\x02", 8)):
            flagged[flagged.find(signature) + offset] |= 1  # the encryption flag
        password = tmp_path / "password.zip"
        password.write_bytes(bytes(flagged))
        return {"two": two, "password": password, "damaged": damaged}

    @pytest.mark.parametrize(
        ("kind", "message"),
        [
            ("two", "is a zip file holding 2 files"),
            ("password", "is a password-protected zip file"),
            ("damaged", "looks truncated or damaged"),
        ],
    )
    @pytest.mark.parametrize("command", ["stats", "analyze"])
    def test_forced_runs_stop_as_detected_ones_do(
        self, zips: dict[str, Path], clinvar_data_dir: Path, kind: str, message: str, command: str
    ) -> None:
        extra = (
            ["--no-update", "--data-dir", str(clinvar_data_dir)] if command == "analyze" else []
        )
        for forced in ([], ["--format", "23andme"]):
            result = _run(command, *forced, *extra, str(zips[kind]))
            assert result.exit_code != 0
            assert message in _flat(result.output)
            assert "Traceback" not in result.output


@pytest.mark.parametrize(
    ("command", "option", "text"),
    [
        ("analyze", "--format", "at your own risk"),
        ("stats", "--format", "at your own risk"),
        ("compare", "--format1", "at your own risk"),
        ("compare", "--format2", "at your own risk"),
    ],
)
def test_help_says_forcing_is_at_your_own_risk(command: str, option: str, text: str) -> None:
    result = _run(command, "--help")
    flat = _flat(result.output)
    assert option in flat
    assert text in flat[flat.index(option) :]

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Utilities disclose imputation beside the input whose calls they display."""

from __future__ import annotations

import re

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import detect_parser
from tests.test_imputed_output_surfaces import _input, _method

KINDS = ("gfg-text", "gfg-vcf", "generic-vcf")


def _notices(output: str) -> list[str]:
    """Join wrapped method notes without borrowing labels from later tables."""
    lines = output.lower().splitlines()
    notes = []
    for index, line in enumerate(lines):
        if "input method:" not in line:
            continue
        paragraph = [line]
        for continuation in lines[index + 1 :]:
            if (
                not continuation.strip()
                or "input method:" in continuation
                or "build:" in continuation
                or any(char in continuation for char in "┏┓┃┌┐│")
                or "coverage summary" in continuation
            ):
                break
            paragraph.append(continuation)
        notes.append(" ".join(" ".join(paragraph).split()))
    return notes


def _assert_label(note: str, kind: str) -> None:
    fragment = (
        "vcf header declares imputation fields"
        if kind == "generic-vcf"
        else "genotyping array extended by imputation"
    )
    assert fragment in note, note
    assert "imputed, not measured" in note, note


@pytest.mark.parametrize("width", [80, 120])
@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("imputed", [False, True])
@pytest.mark.parametrize("command", ["stats", "extract", "plink"])
def test_single_input_utilities_disclose_the_method(
    tmp_path, clinvar_data_dir, monkeypatch, width, kind, imputed, command
):
    from allelix.cli import utility

    # Restore the raw override, including None; the computed width pins it on undo.
    monkeypatch.setattr(utility.console, "_width", width)
    path = _input(tmp_path, kind, imputed)
    assert detect_parser(path).get_metadata(path).get("method") == _method(kind, imputed)
    args = ["export", "plink"] if command == "plink" else [command]
    args += [str(path), "--build", "grch37", "--data-dir", str(clinvar_data_dir)]
    if command == "extract":
        args += ["--snps", "rs900260281"]
    prefix = tmp_path / "export"
    if command == "plink":
        args += ["--output-prefix", str(prefix)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    if command == "extract":
        assert "rs900260281" in result.output
        assert "A/G" in result.output
    elif command == "stats":
        assert "Genotype File Stats" in result.output
    else:
        assert all(prefix.with_suffix(ext).is_file() for ext in (".bed", ".bim", ".fam"))
    notes = _notices(result.output)
    assert len(notes) == int(imputed), result.output
    if imputed:
        _assert_label(notes[0], kind)


@pytest.mark.parametrize("width", [80, 120])
@pytest.mark.parametrize("kind1", KINDS)
@pytest.mark.parametrize("kind2", KINDS)
@pytest.mark.parametrize(
    "imputed1,imputed2", [(False, False), (True, False), (False, True), (True, True)]
)
def test_compare_attributes_each_method_to_its_input(
    tmp_path, monkeypatch, width, kind1, kind2, imputed1, imputed2
):
    from allelix.cli import utility

    # Restore the raw override, including None; the computed width pins it on undo.
    monkeypatch.setattr(utility.console, "_width", width)
    paths = []
    for number, kind, imputed in ((1, kind1, imputed1), (2, kind2, imputed2)):
        folder = tmp_path / f"side{number}"
        folder.mkdir()
        path = _input(folder, kind, imputed)
        assert detect_parser(path).get_metadata(path).get("method") == _method(kind, imputed)
        paths.append(path)
    result = CliRunner().invoke(
        main,
        ["compare", str(paths[0]), str(paths[1]), "--build1", "grch37", "--build2", "grch37"],
    )
    assert result.exit_code == 0, result.output
    assert "Genotype Concordance" in result.output
    notes = _notices(result.output)
    if not (imputed1 or imputed2):
        assert not notes, result.output
    for number, kind, imputed in ((1, kind1, imputed1), (2, kind2, imputed2)):
        associated = [
            note
            for note in notes
            if re.search(rf"\bfile\s*{number}\b", note) or re.search(r"\bboth\s+files\b", note)
        ]
        assert bool(associated) is imputed, result.output
        if imputed:
            assert len(associated) == 1, result.output
            _assert_label(associated[0], kind)
    assert all(re.search(r"\bfile\s*[12]\b|\bboth\s+files\b", note) for note in notes)

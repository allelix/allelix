# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Build notices on every command, from one shared helper (ADR-0044).

``analyze``, ``stats``, ``extract``, ``export plink`` and ``compare`` must
print the same notices for the same decision: the label mismatch, the
format-default note, the tie ("inconclusive") warning and the GRCh36
notice. The wording comes from ``allelix.reports.build_text``. Inputs are
small synthetic MyHappyGenes-format chips.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.reports._pipeline import BuildDiagnostics
from allelix.reports.build_text import build_notices, short_source, source_text
from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    KNOWN_SNP_POSITIONS,
)

if TYPE_CHECKING:
    from pathlib import Path


def _rows(votes: dict[str, int]) -> list[tuple[str, str, int]]:
    """Seed rows: ``votes[build]`` distinct seeds at that build's position."""
    rows = []
    pool = [r for r, e in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in e]
    for build, count in votes.items():
        for rsid in pool[:count]:
            rows.append((rsid, *KNOWN_SNP_POSITIONS[rsid][build]))
        pool = pool[count:]
    return sorted(rows, key=lambda r: (int(r[1]), r[2]))


def _chip(path: Path, label: str | None, votes: dict[str, int]) -> Path:
    header = f"# using human reference build {label} coordinates.\n" if label else ""
    body = "".join(f"{r}\t{c}\t{p}\tA\tG\n" for r, c, p in _rows(votes))
    path.write_text(
        "# MyHappyGenes [TEMPUS]\n"
        + header
        + "# Sample ID\tMHG1\n"
        + "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        + body
    )
    return path


# notice kind -> (label, seed votes, text every command must print)
CASES = {
    "mismatch": ("37.1", {BUILD_GRCH38: 4}, "Build mismatch: the file says GRCh37"),
    "default": (None, {BUILD_GRCH37: 1}, "Note: the file states no genome build"),
    "tie": (
        None,
        {BUILD_GRCH37: 2, BUILD_GRCH38: 2},
        "Build evidence inconclusive: the known-SNP positions are tied between GRCh37 and "
        "GRCh38 (2 each of 4)",
    ),
    "grch36": (None, {BUILD_GRCH36: 4}, "Warning: GRCh36 (hg18) detected"),
}


def _invoke(command: str, chip: Path, tmp_path: Path, data: Path) -> str:
    args = {
        "analyze": ["analyze", str(chip), "--data-dir", str(data), "--no-update"],
        "stats": ["stats", str(chip), "--data-dir", str(data)],
        "extract": ["extract", str(chip), "--snps", "rs1801133", "--data-dir", str(data)],
        "export plink": [
            "export",
            "plink",
            str(chip),
            "-o",
            str(tmp_path / "out"),
            "--data-dir",
            str(data),
        ],
        "compare": ["compare", str(chip), str(chip)],
    }[command]
    result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
    assert result.exit_code == 0, result.output
    return " ".join(result.output.split())


@pytest.mark.parametrize("kind", sorted(CASES))
@pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export plink", "compare"])
def test_every_command_prints_the_notice(
    command: str, kind: str, tmp_path: Path, clinvar_data_dir: Path
) -> None:
    label, votes, expected = CASES[kind]
    out = _invoke(command, _chip(tmp_path / "chip.txt", label, votes), tmp_path, clinvar_data_dir)
    assert expected in out
    if command == "compare":
        assert f"File 1 (chip.txt): {expected}" in out
        assert f"File 2 (chip.txt): {expected}" in out


def test_stats_names_a_grch36_tie(tmp_path: Path, clinvar_data_dir: Path) -> None:
    """A tie that takes the GRCh36 fail-safe is not "from known-SNP positions"."""
    chip = _chip(tmp_path / "chip.txt", None, {BUILD_GRCH36: 2, BUILD_GRCH37: 2})
    out = _invoke("stats", chip, tmp_path, clinvar_data_dir)
    assert "GRCh36 (GRCh36 safety choice; known-SNP positions tied)" in out
    assert "GRCh36 (from known-SNP positions)" not in out
    assert "Build evidence inconclusive" in out


# ── The wording itself: tie banner and tie notice ─────────────────


def _diag(**overrides) -> BuildDiagnostics:
    fields = {
        "effective_build": BUILD_GRCH37,
        "source": "default",
        "file_label": None,
        "detected_build": None,
        "matched_count": 0,
        "inspected_count": 4,
        "position_confident": False,
        "tied_builds": (BUILD_GRCH37, BUILD_GRCH38),
        "tied_votes": 2,
    }
    fields.update(overrides)
    return BuildDiagnostics(**fields)


def test_grch36_tie_banner_names_the_tie() -> None:
    diag = _diag(
        effective_build=BUILD_GRCH36,
        source="seeds",
        tied_builds=(BUILD_GRCH36, BUILD_GRCH37),
        tentative_grch36=True,
    )
    text = source_text(diag)
    assert (
        text
        == "GRCh36 safety choice: known-SNP positions tied between GRCh36 and GRCh37 (2 each of 4)"
    )
    assert "None" not in text
    assert "0/4" not in text
    assert short_source(diag) == "GRCh36 safety choice; known-SNP positions tied"


def test_tie_notice_does_not_splice_the_source_label() -> None:
    texts = [n.text for n in build_notices(_diag())]
    tie = next(t for t in texts if t.startswith("Build evidence inconclusive"))
    assert "format default; the file states no build" not in tie
    assert "Using GRCh37, this format's default, as the file states no build; it is one of" in tie


def test_three_way_tie_wording() -> None:
    diag = _diag(
        effective_build=BUILD_GRCH36,
        source="seeds",
        tied_builds=(BUILD_GRCH36, BUILD_GRCH37, BUILD_GRCH38),
        tied_votes=1,
        inspected_count=3,
        tentative_grch36=True,
    )
    texts = " ".join(n.text for n in build_notices(diag))
    assert "tied between GRCh36, GRCh37 and GRCh38 (1 each of 3)" in texts
    assert "the other build" not in texts
    assert "If the file uses a different build" in texts
    assert "GRCh36 and GRCh37 and GRCh38" not in source_text(diag)


# ── The GRCh36 notice ─────────────────────────────────────────────


def _grch36_text(**overrides) -> str:
    fields = {
        "effective_build": BUILD_GRCH36,
        "source": "seeds",
        "detected_build": BUILD_GRCH36,
        "matched_count": 4,
        "inspected_count": 4,
        "position_confident": True,
        "tied_builds": (),
        "tied_votes": 0,
    }
    fields.update(overrides)
    return next(n.text for n in build_notices(_diag(**fields)) if n.kind == "grch36")


def test_grch36_notice_says_detected_only_when_the_seeds_decided() -> None:
    assert "GRCh36 (hg18) detected from the known-SNP positions." in _grch36_text()
    tie = _grch36_text(
        detected_build=None,
        matched_count=0,
        tied_builds=(BUILD_GRCH36, BUILD_GRCH37),
        tied_votes=2,
        position_confident=False,
        tentative_grch36=True,
    )
    assert "chosen as the safe option because the known-SNP positions are tied" in tie
    assert "detected" not in tie
    label = _grch36_text(source="file_label", file_label=BUILD_GRCH36, detected_build=None)
    assert "the file states this build" in label
    assert "detected" not in label


def test_grch36_notice_claims_only_what_runs() -> None:
    """On GRCh36, ClinVar has no cache and gnomAD / AlphaMissense / CADD need
    GRCh38 coordinates; only the rsID-keyed sources run."""
    text = _grch36_text()
    assert "ClinPGx, GWAS Catalog and SNPedia match by rsID and still run" in text
    assert "gnomAD, AlphaMissense and CADD give no values" in text
    assert "ClinVar is skipped" in text
    assert "complete" not in text


def test_default_note_gives_the_counts() -> None:
    diag = _diag(
        detected_build=BUILD_GRCH37,
        matched_count=3,
        inspected_count=4,
        tied_builds=(),
        tied_votes=0,
    )
    note = next(n.text for n in build_notices(diag) if n.kind == "default")
    assert "3 of 4 known-SNP positions match GRCh37: not enough to be certain" in note
    assert "do not decide" not in note
    tie_note = next(n.text for n in build_notices(_diag()) if n.kind == "default")
    assert "its known-SNP positions are tied between GRCh37 and GRCh38 (2 each of 4)" in tie_note


def test_default_and_tie_notices_offer_grch36() -> None:
    """The same build choices as the stop message, GRCh36 included."""
    clause = (
        "pass --build grch37 or --build grch38 (for an old GRCh36 / hg18 file, "
        "--build grch36, see docs/grch36-liftover.md)."
    )
    texts = {n.kind: n.text for n in build_notices(_diag())}
    assert texts["default"].endswith(f"If the file uses another build, {clause}")
    assert texts["inconclusive"].endswith(f"If the file uses a different build, {clause}")

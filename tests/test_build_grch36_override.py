# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""``--build grch36`` and override validation (ADR-0025, ADR-0044).

A GRCh36 file that the build decision stops (for example an unlabelled
X-only export with no known-SNP position) must have a correct way through:
``--build grch36`` runs it on GRCh36 with position-keyed lookups off. Any
other override than GRCh36, GRCh37 or GRCh38 is rejected at the one choke
point, ``decide_build``. Inputs are synthetic.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import BuildDiagnostics, decide_file_build, run_analysis
from allelix.reports.build_text import build_notices
from allelix.utils.build_detect import (
    BUILD_GRCH36,
    NO_SEEDS,
    BuildDecisionError,
    decide_build,
)

if TYPE_CHECKING:
    from pathlib import Path


def _x_only_chip(path: Path) -> Path:
    """An unlabelled X-only chip: no seed can be present, so it stops by default."""
    path.write_text(
        "# MyHappyGenes [TEMPUS]\n# Sample ID\tMHG1\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        + "".join(f"rs99{n:05d}\tX\t{2_700_000 + n}\tA\tG\n" for n in range(20))
    )
    return path


# ── Validation at decide_build ────────────────────────────────────


@pytest.mark.parametrize("bad", ["hg19", "grch37", "GRCh39", "b37"])
def test_unsupported_override_is_rejected(bad: str) -> None:
    with pytest.raises(ValueError, match="Unsupported build override"):
        decide_build(
            file_name="x.txt",
            kind="array",
            override=bad,
            file_label=None,
            default="GRCh37",
            lengths=None,
            seeds=NO_SEEDS,
        )


def test_python_api_rejects_hg19(tmp_path: Path) -> None:
    chip = _x_only_chip(tmp_path / "x.txt")
    with pytest.raises(ValueError, match="Unsupported build override 'hg19'"):
        run_analysis(chip, MyHappyGenesParser(), [], build_override="hg19")
    parser = MyHappyGenesParser()
    with pytest.raises(ValueError, match="Unsupported build override"):
        decide_file_build(parser, chip, parser.get_metadata(chip), override="hg19")


def test_stop_message_offers_grch36(tmp_path: Path) -> None:
    chip = _x_only_chip(tmp_path / "x.txt")
    with pytest.raises(BuildDecisionError) as exc:
        run_analysis(chip, MyHappyGenesParser(), [])
    message = str(exc.value)
    assert "--build grch37 or --build grch38" in message
    assert "for an old GRCh36 / hg18 file, --build grch36, see docs/grch36-liftover.md" in message


# ── The GRCh36 notice for an override ─────────────────────────────


def test_grch36_override_notice(tmp_path: Path) -> None:
    result = run_analysis(
        _x_only_chip(tmp_path / "x.txt"), MyHappyGenesParser(), [], build_override="GRCh36"
    )
    assert result.build == BUILD_GRCH36
    diag = result.build_diagnostics
    assert diag is not None
    text = next(n.text for n in build_notices(diag) if n.kind == "grch36")
    assert text.startswith("Warning: GRCh36 (hg18) was forced with --build.")


def test_tentative_grch36_is_the_safe_option() -> None:
    diag = BuildDiagnostics(
        effective_build=BUILD_GRCH36,
        source="seeds",
        file_label=None,
        detected_build=BUILD_GRCH36,
        matched_count=2,
        inspected_count=3,
        position_confident=False,
        tentative_grch36=True,
    )
    text = next(n.text for n in build_notices(diag) if n.kind == "grch36")
    assert "was chosen as the safe option: the known-SNP positions lean to it (2/3)" in text
    assert "detected" not in text


# ── --build grch36 on every command ───────────────────────────────


@pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export plink", "compare"])
def test_build_grch36_on_every_command(
    command: str, tmp_path: Path, clinvar_data_dir: Path
) -> None:
    chip = _x_only_chip(tmp_path / "x.txt")
    data = ["--data-dir", str(clinvar_data_dir)]
    args = {
        "analyze": ["analyze", str(chip), *data, "--no-update", "--build", "grch36"],
        "stats": ["stats", str(chip), *data, "--build", "grch36"],
        "extract": ["extract", str(chip), "--snps", "rs9900001", *data, "--build", "grch36"],
        "export plink": [
            "export",
            "plink",
            str(chip),
            "-o",
            str(tmp_path / "out"),
            *data,
            "--build",
            "grch36",
        ],
        "compare": ["compare", str(chip), str(chip), "--build1", "grch36", "--build2", "grch36"],
    }[command]
    runner = CliRunner(env={"COLUMNS": "400"})
    stopped = runner.invoke(
        main, [a for a in args if a not in ("grch36",) and not a.startswith("--build")]
    )
    assert stopped.exit_code != 0  # without the override the file stops
    result = runner.invoke(main, args)
    assert result.exit_code == 0, result.output
    out = " ".join(result.output.split())
    assert "Warning: GRCh36 (hg18) was forced with --build." in out


# ── What runs on GRCh36 (the notice's claims) ─────────────────────


def test_grch36_runs_only_the_rsid_annotators(
    cm_stack, mock_mhg_path: Path, all_annotators_data_dir: Path
) -> None:
    """On GRCh36: no ClinVar rows, no gnomAD / AlphaMissense / CADD values;
    the rsID-keyed annotators (ClinPGx here) still produce rows. The same
    file on its own build gets ClinVar rows and the gnomAD value."""
    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.pharmgkb import PharmGKBAnnotator
    from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
    from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
    from allelix.databases.schema import GNOMAD_SCHEMA

    data = all_annotators_data_dir
    with contextlib.closing(sqlite3.connect(data / GNOMAD_DB_FILENAME)) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES ('1', 11796321, 'G', 'A', 'rs1801133', 0.35)"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('gnomad', 'test://mock', '4.1', "
            "'2026-01-01T00:00:00Z', 1, ?)",
            (f"sv:{GNOMAD_SCHEMA_VERSION}",),
        )
        conn.commit()

    def run(build: str | None):
        annotators = [
            cm_stack.enter_context(ClinVarAnnotator(data)),
            cm_stack.enter_context(PharmGKBAnnotator(data)),
        ]
        gnomad = cm_stack.enter_context(GnomadAnnotator(data))
        return run_analysis(
            mock_mhg_path, MyHappyGenesParser(), annotators, build_override=build, gnomad=gnomad
        )

    own = run(None)
    assert any(a.source == "clinvar" for a in own.annotations)
    assert any(a.allele_frequency is not None for a in own.annotations)

    grch36 = run(BUILD_GRCH36)
    assert grch36.build == BUILD_GRCH36
    assert not any(a.source == "clinvar" for a in grch36.annotations)
    assert any(a.source == "pharmgkb" for a in grch36.annotations)
    for a in grch36.annotations:
        assert a.allele_frequency is None
        assert a.am_pathogenicity is None
        assert a.cadd_phred is None

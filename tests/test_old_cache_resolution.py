# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An older ClinVar cache is not used to resolve unnamed VCF records.

A 2.4.x cache (``iv:5``) has no ``aggregate`` column. ``stats``, ``extract``
and ``export plink`` resolve records without an ID before any update step,
so they must skip such a cache as annotation does, not query it.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename
from allelix.cli import main
from tests.test_clinvar_matching_batch import _write_clinvar

if TYPE_CHECKING:
    from pathlib import Path

# An unnamed record at the synthetic ClinVar site 2:1000 G>A (rs900227001).
_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=2,length=243199373>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "2\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
)


def _data(tmp_path: Path, *, old: bool) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    _write_clinvar(data)
    if old:
        # What 2.4.2 left: the interpreter tag iv:5 and no aggregate column.
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
            conn.execute("ALTER TABLE clinvar_variants DROP COLUMN aggregate")
            conn.execute("UPDATE database_versions SET local_version_tag = 'iv:5'")
            conn.commit()
    return data


def _commands(tmp_path: Path, vcf: Path) -> dict[str, list[str]]:
    return {
        "stats": ["stats", str(vcf)],
        "extract": ["extract", str(vcf), "--snps", "rs900227001"],
        "export plink": ["export", "plink", str(vcf), "-o", str(tmp_path / "out")],
    }


@pytest.mark.parametrize("command", ["stats", "extract", "export plink"])
def test_old_clinvar_cache_is_skipped_not_queried(tmp_path: Path, command: str) -> None:
    data = _data(tmp_path, old=True)
    before = (data / clinvar_db_filename("GRCh37")).read_bytes()
    vcf = tmp_path / "unnamed.vcf"
    vcf.write_text(_VCF)
    args = [*_commands(tmp_path, vcf)[command], "--data-dir", str(data), "--build", "grch37"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert result.exception is None
    assert "Traceback" not in result.output and "no such column" not in result.output
    assert (data / clinvar_db_filename("GRCh37")).read_bytes() == before


def test_current_clinvar_cache_still_resolves_the_unnamed_record(tmp_path: Path) -> None:
    data = _data(tmp_path, old=False)
    vcf = tmp_path / "unnamed.vcf"
    vcf.write_text(_VCF)
    result = CliRunner().invoke(
        main,
        [
            "extract",
            str(vcf),
            "--snps",
            "rs900227001",
            "--data-dir",
            str(data),
            "--build",
            "grch37",
        ],
    )
    assert result.exit_code == 0, result.output
    row = next(line for line in result.output.splitlines() if "rs900227001" in line)
    assert "not in file" not in row and "G/A" in row

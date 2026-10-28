# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extract must complete required build evidence before printing an early hit."""

from __future__ import annotations

import re

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.utils.build_detect import KNOWN_SNP_POSITIONS


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("evidence", ["unmatched", "opposite", "matching"])
@pytest.mark.parametrize("padding", [0, 1100], ids=["short", "past-resolution-batch"])
def test_extract_early_hit_does_not_bypass_required_build_evidence(
    tmp_path, build, evidence, padding
):
    opposite = "GRCh38" if build == "GRCh37" else "GRCh37"
    seed_build = opposite if evidence == "opposite" else build
    rows = [("1", 1000, "rs900252699")]
    rows += [("1", 2000 + index, f"rs900260{index:04d}") for index in range(padding)]
    for rsid, placements in list(KNOWN_SNP_POSITIONS.items())[:4]:
        chrom, pos = placements[seed_build]
        if evidence == "unmatched":
            pos += 1
            assert (chrom, pos) not in placements.values()
        rows.append((chrom, pos, rsid))
    rows.sort(key=lambda row: (int(row[0]), row[1]))
    assert rows[0][2] == "rs900252699"
    length = 249250621 if build == "GRCh37" else 248956422
    source = tmp_path / "early-hit.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=1,length={length}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(
            f"{chrom}\t{pos}\t{rsid}\tA\tG\t.\tPASS\t.\tGT\t0/1\n" for chrom, pos, rsid in rows
        )
    )
    data = tmp_path / "data"
    data.mkdir()
    result = CliRunner(env={"COLUMNS": "400"}).invoke(
        main,
        [
            "extract",
            str(source),
            "--snps",
            "rs900252699",
            "--data-dir",
            str(data),
            "--build",
            build.lower(),
        ],
    )
    if evidence == "opposite":
        assert result.exit_code != 0, result.output
        assert opposite in result.output, result.output
        assert "Genotype" not in result.output, result.output
    else:
        assert result.exit_code == 0, (result.output, result.exception)
        assert "rs900252699" in result.output, result.output
        if evidence == "unmatched":
            text = " ".join(result.output.split())
            assert "Build warning" in text, text
            assert re.search(r"\b0\s*(?:/|of|out of)\s*4\b", text), text
        else:
            assert "Build warning" not in result.output, result.output

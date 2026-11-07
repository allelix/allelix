# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""compare's multi-sample error names --sample1 or --sample2, for the file concerned.

`compare` has no --sample; it has --sample1 and --sample2. Its error for a
multi-sample VCF named --sample, which `compare` then rejects. The first test
and its fixtures are reviewer C's (outside review of 2.6.0 RC4, finding
R4-08), unchanged. Synthetic data only.
"""

from __future__ import annotations

from click.testing import CliRunner

from allelix.cli import main

_VCF_ONE = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "1\t11856378\trs1801133\tG\tA\t50\tPASS\t.\tGT\t0/1\n"
)
_VCF_TWO = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=249250621>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
    "1\t11856378\trs1801133\tG\tA\t50\tPASS\t.\tGT\t0/1\t0/0\n"
)


def test_r08_compare_multisample_error_names_the_compare_option(tmp_path):
    one, two = tmp_path / "one.vcf", tmp_path / "two.vcf"
    one.write_text(_VCF_ONE)
    two.write_text(_VCF_TWO)
    result = CliRunner().invoke(main, ["compare", str(one), str(two)])
    assert result.exit_code != 0
    assert "Available samples" in result.output, result.output
    assert "--sample2" in result.output, (
        f"the advice must name compare's option for file 2:\n{result.output}"
    )


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_a_multi_sample_first_file_names_sample1(tmp_path):
    one, two = tmp_path / "one.vcf", tmp_path / "two.vcf"
    one.write_text(_VCF_TWO)
    two.write_text(_VCF_ONE)
    result = CliRunner().invoke(main, ["compare", str(one), str(two)])
    text = _flat(result.output)
    assert result.exit_code == 1, result.output
    assert "pass --sample1 <ID>" in text
    assert "--sample <ID>" not in text


def test_analyze_still_names_sample(tmp_path):
    """Control: analyze has --sample, and its error says so."""
    vcf = tmp_path / "two.vcf"
    vcf.write_text(_VCF_TWO)
    result = CliRunner().invoke(main, ["analyze", str(vcf), "--data-dir", str(tmp_path / "d")])
    assert result.exit_code == 1, result.output
    assert "pass --sample <ID>" in _flat(result.output)


def test_compare_with_sample2_runs(tmp_path):
    """Control: naming the column with --sample2 lets compare run."""
    one, two = tmp_path / "one.vcf", tmp_path / "two.vcf"
    one.write_text(_VCF_ONE)
    two.write_text(_VCF_TWO)
    result = CliRunner().invoke(main, ["compare", str(one), str(two), "--sample2", "A"])
    assert result.exit_code == 0, result.output

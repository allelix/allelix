# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Structural-variant / copy-number VCFs are rejected with the right advice.

A Sequencing.com 30x order delivers three VCFs. The SV and CNV files hold
only symbolic ALTs, which Allelix does not annotate, so analyzing them
used to "complete" with an empty report. The fixture is the provider's
header plus six synthetic SV records (one FILTER-failed, one breakend).
"""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SV_VCF = FIXTURES / "sequencing_com_sv.vcf"


class TestStructuralOnlyVcf:
    def test_counts_every_structural_record(self) -> None:
        parser = VcfParser()
        assert list(parser.parse(SV_VCF)) == []
        stats = parser.quality_stats
        assert stats.symbolic_rows == 6
        assert dict(stats.symbolic_kinds) == {"DEL": 2, "DUP": 1, "INV": 1, "INS": 1, "BND": 1}

    def test_pipeline_produces_no_annotations(self) -> None:
        result = run_analysis(SV_VCF, VcfParser(), [], build_override="GRCh38")
        assert result.annotations == []
        assert result.total_variants == 0
        assert result.vcf_quality is not None
        assert result.vcf_quality.symbolic_rows == 6

    def test_analyze_rejects_with_the_right_file_named(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        out = tmp_path / "report.json"
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main,
            [
                "analyze",
                str(SV_VCF),
                "--data-dir",
                str(clinvar_data_dir),
                "--no-update",
                "--build",
                "grch38",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code != 0
        message = " ".join(result.output.split())
        assert "only structural-variant / copy-number records" in message
        assert "snp-indel.genome.vcf.gz" in message
        assert "6 record(s)" in message
        assert not out.exists()


class TestOtherVcfsAreNotRejected:
    def test_gvcf_reference_blocks_are_not_structural(self, tmp_path: Path) -> None:
        vcf = tmp_path / "blocks.g.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            '##ALT=<ID=NON_REF,Description="Any other allele">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:GQ\t0/0:50\n"
            "chr1\t101\t.\tC\t<*>\t.\t.\tEND=200\tGT:GQ\t0/0:50\n"
        )
        parser = VcfParser()
        list(parser.parse(vcf))
        assert parser.quality_stats.symbolic_rows == 0

    def test_mixed_file_reports_snvs_and_counts_the_rest(self, tmp_path: Path) -> None:
        vcf = tmp_path / "mixed.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chr1\t1000\t.\tN\t<DEL>\t50\tPASS\tEND=2000\tGT\t0/1\n"
            "chr1\t3000\trs900760001\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        )
        parser = VcfParser()
        variants = list(parser.parse(vcf))
        assert [v.rsid for v in variants] == ["rs900760001"]
        assert parser.quality_stats.symbolic_rows == 1

    def test_breakend_alleles_never_become_genotypes(self, tmp_path: Path) -> None:
        vcf = tmp_path / "bnd.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chr1\t1000\t.\tG\tG]chr2:5000]\t50\tPASS\tSVTYPE=BND\tGT\t0/1\n"
            "chr1\t2000\t.\tT\t[chr3:9000[T\t50\tPASS\tSVTYPE=BND\tGT\t0/1\n"
            "chr1\t3000\t.\tA\tA.\t50\tPASS\tSVTYPE=BND\tGT\t0/1\n"
        )
        parser = VcfParser()
        assert list(parser.parse(vcf)) == []
        assert parser.quality_stats.symbolic_kinds["BND"] == 3

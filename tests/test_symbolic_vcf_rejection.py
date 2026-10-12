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

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SV_VCF = FIXTURES / "sequencing_com_sv.vcf"


# Header evidence that stops the build decision before any record is read.
HEADER_STOPS = {
    "mixed lengths": (
        "##contig=<ID=1,length=249250621>\n##contig=<ID=22,length=50818468>\n",
        "its contig lengths contradict each other",
    ),
    "unsupported lengths": ("##contig=<ID=1,length=123456789>\n", "unsupported assembly"),
    "lengths against assembly=": (
        "##contig=<ID=1,length=248956422,assembly=GRCh37>\n",
        "its header says GRCh37 but its contig lengths are GRCh38",
    ),
}


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


class TestWrongFileErrorsComeBeforeTheBuildStop:
    """A file Allelix cannot use whatever its build gets its own error,
    never the build stop that asks for --build."""

    @staticmethod
    def _run(args: list[str], data: Path) -> str:
        command, *rest = args
        extra = {
            "analyze": ["--data-dir", str(data), "--no-update"],
            "stats": ["--data-dir", str(data)],
            "extract": ["--snps", "rs1801133", "--data-dir", str(data)],
            "export": ["--data-dir", str(data)],
        }[command]
        result = CliRunner(env={"COLUMNS": "300"}).invoke(main, [command, *rest, *extra])
        assert result.exit_code != 0, result.output
        message = " ".join(result.output.split())
        assert "Cannot use" not in message, message
        assert "--build grch37 or --build grch38" not in message, message
        return message

    @staticmethod
    def _sv_with_lengths(tmp_path: Path) -> Path:
        text = SV_VCF.read_text()
        head, _, body = text.partition("#CHROM")
        path = tmp_path / "sv_lengths.vcf"
        path.write_text(head + "##contig=<ID=chr1,length=248956422>\n#CHROM" + body)
        return path

    def test_structural_only_without_build_or_lengths(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        for args in (
            ["analyze", str(SV_VCF)],
            ["stats", str(SV_VCF)],
            ["extract", str(SV_VCF)],
            ["export", "plink", str(SV_VCF), "-o", str(tmp_path / "out")],
        ):
            message = self._run(args, clinvar_data_dir)
            assert "only structural-variant / copy-number records" in message, args
            assert "upload the snp-indel.genome.vcf.gz file" in message, args

    def test_structural_only_with_contig_lengths(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        message = self._run(["analyze", str(self._sv_with_lengths(tmp_path))], clinvar_data_dir)
        assert "only structural-variant / copy-number records" in message

    @pytest.mark.parametrize("how", ["contig lengths", "--build"])
    def test_structural_only_on_every_command_without_a_prescan(
        self, tmp_path: Path, clinvar_data_dir: Path, how: str
    ) -> None:
        """With contig lengths or --build there is no pre-scan; stats, extract and
        export plink still give the structural-only error and write nothing."""
        vcf = self._sv_with_lengths(tmp_path) if how == "contig lengths" else SV_VCF
        forced = ["--build", "grch38"] if how == "--build" else []
        prefix = tmp_path / "out"
        for args in (
            ["stats", str(vcf), *forced],
            ["extract", str(vcf), *forced],
            ["export", "plink", str(vcf), "-o", str(prefix), *forced],
        ):
            message = self._run(args, clinvar_data_dir)
            assert "only structural-variant / copy-number records" in message, args
        assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("out")) == []

    def test_structural_only_in_compare_with_contig_lengths(self, tmp_path: Path) -> None:
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main, ["compare", str(self._sv_with_lengths(tmp_path)), str(FIXTURES / "mock_vcf.vcf")]
        )
        assert result.exit_code != 0
        assert "only structural-variant / copy-number records" in " ".join(result.output.split())

    @classmethod
    def _header_stop_vcf(cls, tmp_path: Path, stop: str, body: str) -> Path:
        contigs, _problem = HEADER_STOPS[stop]
        path = tmp_path / f"{body}.vcf"
        if body == "sv":
            head, _, records = SV_VCF.read_text().partition("#CHROM")
            path.write_text(head + contigs + "#CHROM" + records)
        else:
            rows = "1\t11856378\trs1801133\tG\tA\t50\tPASS\t.\tGT\t0/1\n" if body == "snv" else ""
            path.write_text(
                "##fileformat=VCFv4.2\n"
                + contigs
                + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + rows
            )
        return path

    @staticmethod
    def _command_args(command: str, vcf: Path, tmp_path: Path) -> list[str]:
        return {
            "analyze": ["analyze", str(vcf)],
            "stats": ["stats", str(vcf)],
            "extract": ["extract", str(vcf)],
            "export": ["export", "plink", str(vcf), "-o", str(tmp_path / "out")],
        }[command]

    @pytest.mark.parametrize("body", ["sv", "empty"])
    @pytest.mark.parametrize("stop", sorted(HEADER_STOPS))
    @pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export", "compare"])
    def test_wrong_file_before_a_header_stop(
        self, tmp_path: Path, clinvar_data_dir: Path, command: str, stop: str, body: str
    ) -> None:
        """Header evidence alone would stop the decision; the wrong-file error
        comes first, never the request for --build."""
        vcf = self._header_stop_vcf(tmp_path, stop, body)
        expected = (
            "only structural-variant / copy-number records"
            if body == "sv"
            else "empty.vcf contains no genotype records."
        )
        if command == "compare":
            result = CliRunner(env={"COLUMNS": "300"}).invoke(
                main, ["compare", str(vcf), str(FIXTURES / "mock_vcf.vcf")]
            )
            assert result.exit_code != 0
            message = " ".join(result.output.split())
            assert "Cannot use" not in message, message
        else:
            message = self._run(self._command_args(command, vcf, tmp_path), clinvar_data_dir)
        assert expected in message, message
        assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("out")) == []

    @pytest.mark.parametrize("stop", sorted(HEADER_STOPS))
    @pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export"])
    def test_a_usable_record_still_gets_the_header_stop(
        self, tmp_path: Path, clinvar_data_dir: Path, command: str, stop: str
    ) -> None:
        vcf = self._header_stop_vcf(tmp_path, stop, "snv")
        args = self._command_args(command, vcf, tmp_path)
        extra = (["--snps", "rs1801133"] if command == "extract" else []) + [
            "--data-dir",
            str(clinvar_data_dir),
        ]
        if command == "analyze":
            extra.append("--no-update")
        result = CliRunner(env={"COLUMNS": "300"}).invoke(main, [*args, *extra])
        assert result.exit_code != 0, result.output
        message = " ".join(result.output.split())
        assert f"Cannot use snv.vcf: {HEADER_STOPS[stop][1]}" in message, message

    @staticmethod
    def _vcf(path: Path, contig: str, rows: str = "") -> Path:
        path.write_text(
            "##fileformat=VCFv4.2\n"
            f"##contig=<{contig}>\n"
            '##FILTER=<ID=LowQual,Description="Low quality">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + rows
        )
        return path

    @pytest.mark.parametrize("how", ["contig lengths", "--build"])
    @pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export", "compare"])
    def test_empty_vcf_on_every_command(
        self, tmp_path: Path, clinvar_data_dir: Path, command: str, how: str
    ) -> None:
        """No genotype records: the same error whatever decided the build."""
        contig = "ID=1,length=248956422" if how == "contig lengths" else "ID=1"
        vcf = self._vcf(tmp_path / "empty.vcf", contig)
        forced = ["--build", "grch38"] if how == "--build" else []
        if command == "compare":
            forced = ["--build1", "grch38"] if forced else []
            args = ["compare", str(vcf), str(FIXTURES / "mock_vcf.vcf"), *forced]
            result = CliRunner(env={"COLUMNS": "300"}).invoke(main, args)
            assert result.exit_code != 0
            message = " ".join(result.output.split())
        else:
            rest = (
                ["plink", str(vcf), "-o", str(tmp_path / "out")]
                if command == "export"
                else [str(vcf)]
            )
            message = self._run([command, *rest, *forced], clinvar_data_dir)
        assert "empty.vcf contains no genotype records." in message
        assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("out")) == []

    def test_all_filtered_vcf_is_not_empty(self, tmp_path: Path, clinvar_data_dir: Path) -> None:
        """Records dropped by the quality gate are not "no genotype records"."""
        vcf = self._vcf(
            tmp_path / "filtered.vcf",
            "ID=1,length=248956422",
            "1\t1000000\t.\tA\tG\t5\tLowQual\t.\tGT\t0/1\n",
        )
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main, ["stats", str(vcf), "--data-dir", str(clinvar_data_dir)]
        )
        assert result.exit_code == 0, result.output
        assert "contains no genotype records" not in result.output

    def test_multi_sample_without_sample(self, tmp_path: Path, clinvar_data_dir: Path) -> None:
        vcf = tmp_path / "multi.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
            "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n"
        )
        for command in ("analyze", "stats"):
            message = self._run([command, str(vcf)], clinvar_data_dir)
            assert "Multi-sample VCF multi.vcf: pass --sample <ID>" in message

    def test_unrecognised_and_zero_byte_files(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        zero = tmp_path / "zero.txt"
        zero.write_text("")
        garbage = tmp_path / "garbage.txt"
        garbage.write_text("hello world\n")
        for path in (zero, garbage):
            for command in ("analyze", "stats"):
                message = self._run([command, str(path)], clinvar_data_dir)
                assert f"No parser recognized '{path.name}'" in message

    def test_header_only_file(self, tmp_path: Path, clinvar_data_dir: Path) -> None:
        chip = tmp_path / "header_only.txt"
        chip.write_text(
            "# MyHappyGenes [TEMPUS]\n# Sample ID\tMHG1\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        )
        for command in ("analyze", "stats"):
            message = self._run([command, str(chip)], clinvar_data_dir)
            assert "header_only.txt contains no genotype records." in message


class TestHg19ChrMIsNotEmpty:
    """MT records dropped for an hg19-length chrM are records: the file is not
    empty, and every command that drops them prints the chrM note."""

    @staticmethod
    def _vcf(path: Path, rows: str) -> Path:
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621>\n"
            "##contig=<ID=MT,length=16571>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + rows
        )
        return path

    MT_ROWS = (
        "MT\t73\trs3087742\tA\tG\t50\tPASS\t.\tGT\t1/1\n"
        "MT\t263\trs2853515\tA\tG\t50\tPASS\t.\tGT\t1/1\n"
    )
    NUCLEAR_ROW = "1\t11856378\trs1801133\tG\tA\t50\tPASS\t.\tGT\t0/1\n"

    @staticmethod
    def _invoke(command: str, vcf: Path, tmp_path: Path, data: Path) -> str:
        args = {
            "analyze": ["analyze", str(vcf), "--data-dir", str(data), "--no-update"],
            "stats": ["stats", str(vcf), "--data-dir", str(data)],
            "extract": ["extract", str(vcf), "--snps", "rs1801133", "--data-dir", str(data)],
            "export plink": [
                "export",
                "plink",
                str(vcf),
                "-o",
                str(tmp_path / "out"),
                "--data-dir",
                str(data),
            ],
        }[command]
        if command == "analyze":
            args += ["--output", str(tmp_path / "report.json")]
        result = CliRunner(env={"COLUMNS": "400"}).invoke(main, args)
        assert result.exit_code == 0, result.output
        return " ".join(result.output.split())

    @pytest.mark.parametrize("command", ["analyze", "stats", "extract"])
    def test_mt_only_file_gets_the_chrm_note(
        self, tmp_path: Path, clinvar_data_dir: Path, command: str
    ) -> None:
        vcf = self._vcf(tmp_path / "mt_only.vcf", self.MT_ROWS)
        out = self._invoke(command, vcf, tmp_path, clinvar_data_dir)
        assert "contains no genotype records" not in out
        assert "chrM is declared with the hg19 length (16,571 bp)" in out
        assert "2 mitochondrial record(s) were not annotated." in out

    def test_export_plink_refuses_an_empty_set(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """Every record was dropped, so there is nothing to write: say why, write nothing."""
        vcf = self._vcf(tmp_path / "mt_only.vcf", self.MT_ROWS)
        args = ["export", "plink", str(vcf), "-o", str(tmp_path / "out")]
        result = CliRunner(env={"COLUMNS": "400"}).invoke(
            main, [*args, "--data-dir", str(clinvar_data_dir)]
        )
        assert result.exit_code != 0
        out = " ".join(result.output.split())
        assert "chrM is declared with the hg19 length (16,571 bp)" in out
        assert (
            "Nothing to export from mt_only.vcf: no variant can be written to a PLINK set" in out
        )
        assert "2 mitochondrial records not exported (chrM has the hg19 length)" in out
        assert "No files were written." in out
        assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("out")) == []

    @pytest.mark.parametrize("command", ["analyze", "stats", "extract", "export plink"])
    def test_note_whenever_mt_records_were_dropped(
        self, tmp_path: Path, clinvar_data_dir: Path, command: str
    ) -> None:
        """A file with nuclear records too; extract finds its rsID before the MT
        records and still reads on (contig lengths decided), so it counts them."""
        vcf = self._vcf(tmp_path / "mixed.vcf", self.NUCLEAR_ROW + self.MT_ROWS)
        out = self._invoke(command, vcf, tmp_path, clinvar_data_dir)
        assert "2 mitochondrial record(s) were not annotated." in out

    def test_compare_keeps_the_mt_records(self, tmp_path: Path) -> None:
        """compare matches by rsID and drops nothing, so the file is not empty there."""
        vcf = self._vcf(tmp_path / "mt_only.vcf", self.MT_ROWS)
        result = CliRunner(env={"COLUMNS": "400"}).invoke(main, ["compare", str(vcf), str(vcf)])
        assert result.exit_code == 0, result.output
        assert "contains no genotype records" not in result.output


class TestCompareKeepsEachFileApart:
    """compare reads two files; each gets its own parser instance, so one
    file's checks never read the other file's parse state."""

    def test_resolve_parser_returns_a_new_instance_each_call(self) -> None:
        from allelix.cli._helpers import _resolve_parser

        vcf = FIXTURES / "mock_vcf.vcf"
        chip = FIXTURES / "mock_23andme.txt"
        assert _resolve_parser(vcf, None) is not _resolve_parser(vcf, None)
        assert _resolve_parser(chip, None) is not _resolve_parser(chip, None)

    def test_vcf_plus_vcf_with_a_structural_only_file(self) -> None:
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main, ["compare", str(SV_VCF), str(FIXTURES / "mock_vcf.vcf")]
        )
        assert result.exit_code != 0
        message = " ".join(result.output.split())
        assert "only structural-variant / copy-number records" in message
        assert "contains no genotype records" not in message

    @pytest.mark.parametrize(
        ("sample2", "concordant", "discordant"), [("SAMPLE_A", 2, 0), ("SAMPLE_B", 0, 1)]
    )
    def test_two_multi_sample_vcfs_keep_their_samples(
        self, sample2: str, concordant: int, discordant: int
    ) -> None:
        multi = FIXTURES / "mock_multisample.vcf"
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main,
            ["compare", str(multi), str(multi), "--sample1", "SAMPLE_A", "--sample2", sample2],
        )
        assert result.exit_code == 0, result.output
        rows = {
            parts[1]: int(parts[2])
            for line in result.output.splitlines()
            if (parts := [p.strip() for p in line.split("│")])
            and len(parts) > 3
            and parts[1] in ("Concordant", "Discordant")
        }
        assert rows == {"Concordant": concordant, "Discordant": discordant}

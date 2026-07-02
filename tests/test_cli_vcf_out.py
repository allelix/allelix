# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""End-to-end tests for the ``analyze --vcf-out`` code path (#150 PR 2).

Covers the CLI flag, the preflight validation (non-VCF input rejected,
existing output refused), and the pipeline hook that builds the
annotation dict + recovered-rsID map fed to the writer.

Uses the ``clinvar_data_dir`` fixture from ``tests/conftest.py`` — the
synthetic ClinVar cache with a handful of known-carrier rows so a real
analyze pipeline actually produces annotations for the test VCFs.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main

if TYPE_CHECKING:
    from _pytest.fixtures import FixtureRequest  # noqa: F401

FIXTURES = Path(__file__).resolve().parent / "fixtures"


class TestVcfOutHappyPath:
    def test_analyze_vcf_out_produces_valid_annotated_vcf(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """Feed a VCF with a known-pathogenic variant through analyze
        with ``--vcf-out``. Expect the annotated output to carry the
        Allelix provenance block, at least one ``ALLELIX_CLINVAR=``
        INFO stamp, and preserve every input header line."""
        out = tmp_path / "annotated.vcf"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        text = out.read_text(encoding="utf-8")
        # Provenance block present.
        assert "##ALLELIX_Version=" in text
        assert "##ALLELIX_VCF_SchemaVersion=0.1.0" in text
        assert "##ALLELIX_License=AGPL-3.0-or-later" in text
        assert "##ALLELIX_Build=" in text
        # ClinVar hit stamped on rs1801133 (MTHFR C677T carrier in mock).
        assert "ALLELIX_CLINVAR=" in text
        # Original input header line preserved.
        assert "##fileformat=VCFv4.2" in text
        # #CHROM row still present.
        assert "#CHROM\tPOS\tID\tREF\tALT" in text
        # User-facing confirmation printed.
        assert "Wrote annotated VCF" in result.output

    def test_pipeline_populates_vcf_write_inputs_only_when_requested(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """Default analyze (no --vcf-out) leaves ``vcf_write_inputs`` as
        None. --vcf-out sets both dicts non-empty."""
        from allelix.annotators import get_annotators
        from allelix.parsers import detect_parser
        from allelix.reports._pipeline import run_analysis

        parser = detect_parser(FIXTURES / "mock_vcf.vcf")
        annotators = [
            a for a in get_annotators(clinvar_data_dir) if a.is_ready() and a.name == "clinvar"
        ]

        # Default: collect flag off.
        default_result = run_analysis(FIXTURES / "mock_vcf.vcf", parser, annotators)
        assert default_result.vcf_write_inputs is None

        # Requested: collect flag on.
        requested_result = run_analysis(
            FIXTURES / "mock_vcf.vcf",
            parser,
            annotators,
            collect_vcf_write_inputs=True,
        )
        assert requested_result.vcf_write_inputs is not None
        assert requested_result.vcf_write_inputs.annotations_by_key
        # Keys are (chrom, pos, ref, alt) tuples with real VCF ALT values.
        for chrom, pos, ref, alt in requested_result.vcf_write_inputs.annotations_by_key:
            assert isinstance(chrom, str) and chrom
            assert isinstance(pos, int) and pos > 0
            assert isinstance(ref, str) and ref
            assert isinstance(alt, str) and alt


class TestVcfOutPreflight:
    def test_non_vcf_input_rejected_with_clear_error(
        self, tmp_path: Path, mock_mhg_path: Path, clinvar_data_dir: Path
    ) -> None:
        """--vcf-out against a 23andMe/MyHappyGenes file must error
        early with a message naming the incompatible format. The output
        file must NOT be created."""
        out = tmp_path / "should_not_exist.vcf"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code != 0
        assert "--vcf-out requires VCF" in result.output
        assert not out.exists()

    def test_existing_output_path_refused_no_clobber(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """--vcf-out at an existing path must refuse and leave the file
        untouched. No --force in v2.3.0."""
        existing = tmp_path / "existing.vcf"
        existing.write_text("pre-existing content\n", encoding="utf-8")
        mtime_before = existing.stat().st_mtime

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(existing),
            ],
        )
        assert result.exit_code != 0
        assert "refuses to overwrite" in result.output
        assert existing.read_text(encoding="utf-8") == "pre-existing content\n"
        assert existing.stat().st_mtime == mtime_before

    def test_vcf_out_composes_with_report_output(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """--vcf-out and --output are orthogonal; both should produce
        their respective files in one run."""
        vcf_out = tmp_path / "annotated.vcf"
        json_out = tmp_path / "report.json"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(vcf_out),
                "--output",
                str(json_out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert vcf_out.exists()
        assert json_out.exists()
        assert "Wrote annotated VCF" in result.output


class TestBcftoolsRoundTrip:
    """Verifies that ``bcftools view`` accepts the writer's output — the
    standard third-party validation for a VCF annotator.

    This is the CI gate the ADR-0036 promises: any regression in the
    ``##INFO=<>`` declaration syntax, per-allele ``Number=A`` value
    counts, provenance header format, or record-line structure would
    surface as a non-zero ``bcftools view`` exit before merge.

    Skips locally when ``bcftools`` is not on ``PATH``; CI installs it
    via the fast-tier apt step and runs the test on every push.
    """

    @pytest.fixture(autouse=True)
    def _require_bcftools(self) -> None:
        if shutil.which("bcftools") is None:
            pytest.skip("bcftools not installed; install to run the round-trip validation")

    def _run_bcftools(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bcftools", *args],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_bcftools_view_accepts_annotated_output(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """``bcftools view`` on the annotated output exits 0 with no
        errors on stderr and no warnings caused by the writer.

        Pre-existing fixture noise (mock_vcf.vcf declares only
        ``##contig=1,22`` but carries data on chroms 17 / 19) surfaces
        as ``Contig 'X' is not defined in the header`` warnings from
        bcftools — those are an input-header issue, not a writer
        regression, and are filtered here so this test catches only
        writer-introduced warnings.
        """
        out = tmp_path / "annotated.vcf"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        proc = self._run_bcftools("view", str(out))
        assert proc.returncode == 0, (
            f"bcftools view exited {proc.returncode}\n"
            f"stderr: {proc.stderr}\nstdout head: {proc.stdout[:400]}"
        )
        significant_lines = [
            line
            for line in proc.stderr.splitlines()
            if line.strip() and "is not defined in the header" not in line
        ]
        assert not significant_lines, (
            f"bcftools view emitted writer-attributable issues: {significant_lines}"
        )

    def test_bcftools_accepts_chr_prefixed_annotated_output(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """chr-prefix regression guard in the third-party CI gate.

        The unit tests in ``test_vcf_writer.py`` lock the
        ``normalize_chromosome`` boundary at the writer level; this
        pipes a chr-prefixed input through the full CLI and asks
        bcftools to accept the output, so the CI catches a boundary
        regression even if a future refactor breaks the unit tests'
        blast radius.

        Uses ``mock_vcf_chr_prefix.vcf`` (chr1 / chr22 with rs1801133,
        rs1801131, rs4680) — the same known-hit rows as ``mock_vcf.vcf``
        but chr-prefixed. Contigs declared in the header so bcftools
        doesn't emit the pre-existing fixture-drift warnings we filter
        for the bare-chrom test above.
        """
        out = tmp_path / "annotated.vcf"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf_chr_prefix.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        proc = self._run_bcftools("view", str(out))
        assert proc.returncode == 0, (
            f"bcftools view exited {proc.returncode}\n"
            f"stderr: {proc.stderr}\nstdout head: {proc.stdout[:400]}"
        )
        assert not proc.stderr.strip(), (
            f"bcftools view emitted issues on chr-prefixed output: {proc.stderr}"
        )
        # The rs1801133 row (chr1:11796321) is the load-bearing
        # assertion: this variant is in the synthetic ClinVar GRCh38
        # cache, so a working writer stamps it. A silent chr-prefix
        # regression would leave this row unstamped and drop the
        # test.
        stamped = next(
            line
            for line in out.read_text(encoding="utf-8").splitlines()
            if line.startswith("chr1\t11796321\t")
        )
        assert "ALLELIX_CLINVAR=" in stamped
        assert "rs1801133" in stamped

    def test_bcftools_header_lists_allelix_info_fields(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """``bcftools view -h`` prints every declared ``##INFO=<ID=ALLELIX_*>``
        line — sanity check that our header declarations survive
        bcftools's own parse."""
        out = tmp_path / "annotated.vcf"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        proc = self._run_bcftools("view", "-h", str(out))
        assert proc.returncode == 0, proc.stderr
        header = proc.stdout
        assert "##INFO=<ID=ALLELIX_CLINVAR" in header
        # Provenance survives bcftools too.
        assert "##ALLELIX_VCF_SchemaVersion=0.1.0" in header
        assert "##ALLELIX_License=AGPL-3.0-or-later" in header

    def test_bcftools_query_reads_allelix_info_field(
        self, tmp_path: Path, clinvar_data_dir: Path
    ) -> None:
        """``bcftools query`` on a known-stamped row extracts the same
        ``ALLELIX_CLINVAR`` value the writer emitted — proves bcftools
        parses our field declaration end-to-end, not just view-accepts
        it. rs1801133 (MTHFR C677T) is in the synthetic ClinVar cache
        and is present as a heterozygous carrier in ``mock_vcf.vcf``,
        so the writer definitely stamps this row.
        """
        out = tmp_path / "annotated.vcf"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_vcf.vcf"),
                "--data-dir",
                str(clinvar_data_dir),
                "--vcf-out",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        proc = self._run_bcftools(
            "query",
            "-f",
            "%CHROM\t%POS\t%ID\t%INFO/ALLELIX_CLINVAR\n",
            str(out),
        )
        assert proc.returncode == 0, proc.stderr
        rs1801133_line = next(
            (line for line in proc.stdout.splitlines() if "rs1801133" in line), None
        )
        assert rs1801133_line is not None, (
            f"rs1801133 row not present in bcftools query output: {proc.stdout!r}"
        )
        # bcftools returns the raw value from the INFO field verbatim.
        # ``Number=A`` on a single-ALT row is just one value (no comma).
        assert "clinvar_pathogenic|rs1801133|MTHFR" in rs1801133_line

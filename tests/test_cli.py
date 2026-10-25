# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""End-to-end tests for the CLI."""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli import main
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.cache_state import READY, CacheState, CacheStatus
from allelix.databases.schema import GNOMAD_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path


class TestStatsCommand:
    def test_stats_runs_on_fixture(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        assert "MyHappyGenes" in result.output
        assert "Total SNPs" in result.output
        assert "MHG000001" in result.output

    def test_stats_with_explicit_format(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(mock_mhg_path), "--format", "myhappygenes"])
        assert result.exit_code == 0, result.output

    def test_stats_counts_match_fixture_exactly(self, mock_mhg_path):
        """Pin specific numeric counts so a `+= 0` mutation fails the suite."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["stats", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        # Generated mock fixture (ADR-0021: rs104894490 NIPA1 added):
        #   2,016 SNPs / 103 no-call (5.11%) / 563 het (27.93%) / 1,350 hom (66.96%).
        assert "2,016" in result.output
        assert "5.11%" in result.output
        assert "27.93%" in result.output
        assert "66.96%" in result.output
        assert "1,350" in result.output

    def test_stats_count_invariant(self, mock_mhg_path):
        """Belt-and-braces: parser-level counts add up. Catches mutations directly."""
        from allelix.parsers import detect_parser

        parser = detect_parser(mock_mhg_path)
        total = no_calls = het = hom = 0
        for v in parser.parse(mock_mhg_path):
            total += 1
            if v.is_no_call:
                no_calls += 1
            elif v.is_heterozygous:
                het += 1
            else:
                hom += 1
        assert total == no_calls + het + hom
        assert (total, no_calls, het, hom) == (2016, 103, 563, 1350)

    def test_stats_unknown_format_errors(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(mock_mhg_path), "--format", "nonsense"])
        assert result.exit_code != 0
        assert "Unknown parser" in result.output

    def test_stats_build_comes_from_the_decision(self, tmp_path):
        """``stats`` shows the build ``analyze`` would use, and where it came
        from (ADR-0044): here hg19 contig lengths under chr-prefixed names."""
        runner = CliRunner(env={"COLUMNS": "200"})
        vcf = tmp_path / "minimal_chr.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=249250621>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "chr1\t11856378\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        )
        result = runner.invoke(main, ["stats", str(vcf)])
        assert result.exit_code == 0, result.output
        assert "GRCh37 (from the contig lengths)" in result.output

    def test_stats_build_override(self, tmp_path):
        runner = CliRunner(env={"COLUMNS": "200"})
        vcf = tmp_path / "bare.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        )
        result = runner.invoke(main, ["stats", str(vcf), "--build", "grch38"])
        assert result.exit_code == 0, result.output
        assert "GRCh38 (forced with --build)" in result.output

    def test_stats_stops_on_a_vcf_without_build_evidence(self, tmp_path):
        runner = CliRunner(env={"COLUMNS": "200"})
        vcf = tmp_path / "bare.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        )
        result = runner.invoke(main, ["stats", str(vcf)])
        assert result.exit_code != 0
        assert "Cannot use bare.vcf" in result.output
        assert "Total SNPs" not in result.output

    def test_stats_unrecognized_file_errors(self, tmp_path):
        f = tmp_path / "garbage.txt"
        f.write_text("hello world\n", encoding="utf-8")
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(f)])
        assert result.exit_code != 0
        assert "Allelix doesn't recognise the format of 'garbage.txt'" in result.output

    def test_stats_missing_file_errors(self, tmp_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(tmp_path / "missing.txt")])
        assert result.exit_code != 0


class TestSkippedLineReporting:
    """Malformed lines must surface to the user, not vanish silently."""

    def _dirty_file(self, tmp_path):
        f = tmp_path / "dirty.txt"
        f.write_text(
            "# MyHappyGenes [TEMPUS], human reference build 37.1 coordinates\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tA\tG\n"
            "junk\n"
            "rs2\t2\tBAD\tA\tG\n"
            # one known-SNP position, so the build line has evidence to check it
            "rs1801133\t1\t11856378\tA\tG\n",
            encoding="utf-8",
        )
        return f

    def test_stats_reports_skipped_lines(self, tmp_path):
        f = self._dirty_file(tmp_path)
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(f)])
        assert result.exit_code == 0, result.output
        assert "Skipped (malformed)" in result.output
        assert "warning:" in result.stderr
        assert "Line 5" in result.stderr  # the "junk" line
        assert "Line 6" in result.stderr  # the BAD position line

    def test_stats_no_skipped_row_when_clean(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        assert "Skipped (malformed)" not in result.output
        assert result.stderr == ""

    def test_stats_restores_logger_state(self, tmp_path):
        """Do not leak logger level or propagate setting across calls."""
        f = self._dirty_file(tmp_path)
        parser_logger = logging.getLogger("allelix.parsers")
        prev_level = parser_logger.level
        prev_propagate = parser_logger.propagate
        prev_handler_count = len(parser_logger.handlers)

        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(f)])
        assert result.exit_code == 0

        assert parser_logger.level == prev_level
        assert parser_logger.propagate == prev_propagate
        assert len(parser_logger.handlers) == prev_handler_count


class TestAnalyzeCommand:
    def test_analyze_stops_on_zero_vote_seeds(self, tmp_path, clinvar_data_dir: Path):
        """Through the CLI: a labelled chip whose seeds sit on no build stops."""
        chip = tmp_path / "nowhere.txt"
        chip.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# using human reference build 37.1 coordinates.\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1801133\t1\t7\tA\tG\n"
            "rs4680\t22\t7\tA\tG\n"
        )
        out = tmp_path / "report.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["analyze", str(chip), "--data-dir", str(clinvar_data_dir), "--no-update",
             "--output", str(out)],
        )  # fmt: skip
        assert result.exit_code != 0
        flat = " ".join(result.output.split())
        assert (
            "Cannot use nowhere.txt: its known-SNP positions match none of the supported" in flat
        )
        assert "2 known-SNP position(s) found, none at any supported build's coordinate" in flat
        assert not out.exists()

    def test_analyze_stops_on_an_unlabelled_chip_without_seeds(
        self, tmp_path, clinvar_data_dir: Path
    ):
        """An unlabelled chip with no seeds, through the CLI: exit non-zero, say why,
        write nothing."""
        chip = tmp_path / "x_only.txt"
        chip.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            + "".join(f"rs99{n:05d}\tX\t{2_700_000 + n}\tA\tG\n" for n in range(20))
        )
        out = tmp_path / "report.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        args = ["analyze", str(chip), "--data-dir", str(clinvar_data_dir), "--no-update"]
        result = runner.invoke(main, [*args, "--output", str(out)])
        assert result.exit_code != 0
        flat = " ".join(result.output.split())
        assert "Cannot use x_only.txt: the file states no genome build" in flat
        assert "no known-SNP positions in the file" in flat
        assert not out.exists()
        forced = runner.invoke(main, [*args, "--build", "grch37", "--output", str(out)])
        assert forced.exit_code == 0, forced.output
        assert out.exists()

    def test_errors_when_no_annotators_ready(self, mock_mhg_path, tmp_path: Path):
        runner = CliRunner()
        empty_dir = tmp_path / "empty"
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_path), "--data-dir", str(empty_dir)],
        )
        assert result.exit_code != 0
        assert "db update" in result.output

    def test_analyze_renders_known_pathogenic(self, mock_mhg_path, clinvar_data_dir: Path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_path), "--data-dir", str(clinvar_data_dir)],
        )
        assert result.exit_code == 0, result.output
        # MTHFR C677T is heterozygous in the mock fixture and Pathogenic in mock ClinVar
        assert "rs1801133" in result.output
        assert "MTHFR" in result.output
        assert "ClinVar" in result.output
        # terminal table strips the source prefix from significance
        # since the Source column already shows ClinVar; expect plain
        # ``pathogenic`` instead of ``clinvar_pathogenic``.
        assert "pathogenic" in result.output

    def test_analyze_does_not_flag_homozygous_reference(
        self, mock_mhg_path, clinvar_data_dir: Path
    ):
        # rs121918506 is REF=G ALT=T Pathogenic in mock ClinVar; mock MHG has G/G
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_path), "--data-dir", str(clinvar_data_dir)],
        )
        assert result.exit_code == 0, result.output
        assert "rs121918506" not in result.output

    def test_analyze_uses_env_var_data_dir(self, mock_mhg_path, clinvar_data_dir: Path):
        """Default data-dir resolution via $ALLELIX_DATA_DIR works in the CLI."""
        runner = CliRunner(env={"COLUMNS": "200", "ALLELIX_DATA_DIR": str(clinvar_data_dir)})
        result = runner.invoke(main, ["analyze", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        assert "rs1801133" in result.output

    def test_analyze_warns_on_build_header_data_mismatch(
        self, mock_mhg_mislabeled_path, clinvar_data_dir: Path
    ):
        """ADR-0021 end-to-end: file header claims GRCh37 but positions are
        GRCh38 (replicating the real MyHappyGenes/Tempus mislabel).
        `allelix analyze` must (a) detect GRCh38 from positions, (b) print
        the build mismatch warning, (c) annotate against the GRCh38 cache.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_mislabeled_path), "--data-dir", str(clinvar_data_dir)],
        )
        assert result.exit_code == 0, result.output
        # Width-independent — Rich wraps "Using GRCh38" mid-phrase at
        # COLUMNS=80 (the CliRunner default), even though this test sets
        # COLUMNS=200 in env. The warning text passes through Rich's
        # console.print which respects the runtime width.
        flat = " ".join(result.output.split())
        assert "Build mismatch" in flat
        assert "the file says GRCh37, but its known-SNP positions say GRCh38" in flat
        assert "Using GRCh38" in flat

    def test_analyze_no_warning_on_clean_grch37(
        self, mock_mhg_grch37_path, clinvar_data_dir: Path
    ):
        """ADR-0021: a file whose header and positions both agree on GRCh37
        triggers no warning. Detection silently confirms.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_grch37_path), "--data-dir", str(clinvar_data_dir)],
        )
        assert result.exit_code == 0, result.output
        assert "Build mismatch" not in result.output
        # The dim build banner still appears.
        assert "Build: GRCh37" in result.output

    def test_analyze_build_override_is_checked(
        self, mock_mhg_mislabeled_path, clinvar_data_dir: Path
    ):
        """On the deliberately-mislabeled fixture (label 37, positions 38),
        `--build grch38` forces the build and silences the mismatch warning;
        `--build grch37` stops, because the positions confidently say GRCh38.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        args = ["analyze", str(mock_mhg_mislabeled_path), "--data-dir", str(clinvar_data_dir)]
        result = runner.invoke(main, [*args, "--build", "grch38"])
        assert result.exit_code == 0, result.output
        assert "Build mismatch" not in result.output
        assert "Build: GRCh38 (forced with --build; the file says GRCh37)" in result.output
        result = runner.invoke(main, [*args, "--build", "grch37"])
        assert result.exit_code != 0
        output = " ".join(result.output.split())
        assert "--build grch37 was given but the file's known-SNP positions say GRCh38" in output

    def test_analyze_unknown_format_errors(self, mock_mhg_path, clinvar_data_dir: Path):
        """Cover the ParserNotFoundError branch in analyze."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--format",
                "nonsense",
            ],
        )
        assert result.exit_code != 0
        assert "Unknown parser" in result.output

    def test_analyze_reports_skipped_lines_on_dirty_input(
        self, tmp_path: Path, clinvar_data_dir: Path
    ):
        """Malformed-line warnings surface in the analyze command too."""
        f = tmp_path / "dirty.txt"
        f.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1801133\t1\t11856378\tG\tA\n"
            "this_is_garbage\n",
            encoding="utf-8",
        )
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["analyze", str(f), "--data-dir", str(clinvar_data_dir)])
        assert result.exit_code == 0, result.output
        assert "1 malformed line" in result.output
        assert "warning:" in result.stderr

    def test_analyze_with_both_annotators(self, mock_mhg_path, all_annotators_data_dir: Path):
        """End-to-end with ClinVar + ClinPGx both ready: rsids fire from each source."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_path), "--data-dir", str(all_annotators_data_dir)],
        )
        assert result.exit_code == 0, result.output
        # terminal strips source_prefix from significance display.
        # ``clinvar_pathogenic`` -> ``pathogenic``,
        # ``pharmgkb_loe_2a``   -> ``loe_2a`` next to the ClinPGx source col.
        assert "pathogenic" in result.output
        assert "loe_2a" in result.output
        # Both attribution labels present
        assert "ClinVar" in result.output
        assert "ClinPGx" in result.output
        # Summary line shows 3 databases
        assert "3 database" in result.output

    def test_min_magnitude_filter(self, mock_mhg_path, clinvar_data_dir: Path):
        # rs1801394 (MTRR) is Likely_benign (mag 2.0); should drop above threshold
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--min-magnitude",
                "5",
                "--include-benign",
            ],
        )
        assert result.exit_code == 0
        assert "rs1801394" not in result.output

    def test_benign_suppressed_by_default(self, mock_mhg_path, clinvar_data_dir: Path):
        """ADR-0008 amendment: Benign/Likely_benign are excluded by default."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "rs1801394" not in result.output

    def test_include_benign_flag(self, mock_mhg_path, clinvar_data_dir: Path):
        """--include-benign restores Benign/Likely_benign annotations."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--min-magnitude",
                "0",
                "--include-benign",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "rs1801394" in result.output

    def test_gwas_min_magnitude_default(self, mock_mhg_path, all_annotators_data_dir: Path):
        """Default --gwas-min-magnitude 9.0 filters all mock GWAS rows (max is 8.0)."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        # All mock GWAS rows are < 9.0 (max is Breast cancer at 8.0)
        assert "Pain sensitivity" not in result.output
        assert "Breast cancer" not in result.output
        assert "Homocysteine levels" not in result.output

    def test_gwas_min_magnitude_lowered(self, mock_mhg_path, all_annotators_data_dir: Path):
        """--gwas-min-magnitude 7.0 lets mag-8 GWAS rows through."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
                "--gwas-min-magnitude",
                "7.0",
            ],
        )
        assert result.exit_code == 0, result.output
        # Breast cancer (p=5e-30, OR=4.2 → mag 8.0) passes floor 7.0
        assert "Breast cancer" in result.output
        # Pain sensitivity (p=3e-9 → mag 6.0) still filtered
        assert "Pain sensitivity" not in result.output


class TestDbCommands:
    @pytest.fixture(autouse=True)
    def _no_identity_download(self, monkeypatch):
        """`db update` walks the whole registry; keep the pinned identity
        table's real download out of these tests."""
        from allelix.annotators.rsid_identity import RsidIdentityAnnotator

        monkeypatch.setattr(RsidIdentityAnnotator, "setup", lambda self: None)

    def test_status_with_empty_dir(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["db", "status", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "ClinVar" in result.output
        assert "no" in result.output  # not ready

    def test_status_with_populated_dir(self, clinvar_data_dir: Path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["db", "status", "--data-dir", str(clinvar_data_dir)])
        assert result.exit_code == 0
        assert "ClinVar" in result.output
        assert "yes" in result.output
        assert "20260101" in result.output  # version from VCF ##fileDate
        # ADR-0021: clinvar_data_dir populates both build caches.
        # 13 records per build x 2 builds = 26 composite total.
        assert "26" in result.output

    def test_db_update_only_clinvar_skips_every_other_source(self, tmp_path, monkeypatch):
        """The slow sentinel can refresh ClinVar without pulling enrichment caches."""

        class FakeAnnotator:
            requires_download = True
            server_driven_freshness = True

            def __init__(self, name):
                self.name = name
                self.setup_called = False

            def __enter__(self):
                return self

            def __exit__(self, *_exc):
                return None

            def is_ready(self):
                return False

            def cache_status(self):
                return READY if self.is_ready() else CacheStatus(CacheState.MISSING)

            def setup(self):
                self.setup_called = True

            def version(self):
                return "test"

            def cached_remote_signal(self):
                return None

        clinvar = FakeAnnotator("clinvar")
        gnomad = FakeAnnotator("gnomad")
        monkeypatch.setattr(
            "allelix.cli.db.get_annotators",
            lambda *_args, **_kwargs: [clinvar, gnomad],
        )

        result = CliRunner().invoke(
            main,
            ["db", "update", "--data-dir", str(tmp_path), "--only", "clinvar"],
        )

        assert result.exit_code == 0, result.output
        assert clinvar.setup_called is True
        assert gnomad.setup_called is False
        assert "clinvar ready" in result.output
        assert "gnomad" not in result.output

    @pytest.mark.parametrize("extra", [[], ["--force"]])
    def test_db_update_exits_nonzero_when_a_source_fails(self, tmp_path, monkeypatch, extra):
        """A failed install or refresh must not be masked by exit code 0;
        the other sources are still attempted."""

        class FakeAnnotator:
            requires_download = True
            server_driven_freshness = True

            def __init__(self, name, fails):
                self.name = name
                self.fails = fails
                self.setup_called = False

            def __enter__(self):
                return self

            def __exit__(self, *_exc):
                return None

            def is_ready(self):
                return bool(extra)

            def cache_status(self):
                return READY if self.is_ready() else CacheStatus(CacheState.MISSING)

            def setup(self):
                self.setup_called = True
                if self.fails:
                    raise RuntimeError("download failed")

            def version(self):
                return "test"

            def cached_remote_signal(self):
                return None

        broken = FakeAnnotator("gnomad", fails=True)
        working = FakeAnnotator("clinvar", fails=False)
        monkeypatch.setattr(
            "allelix.cli.db.get_annotators",
            lambda *_args, **_kwargs: [broken, working],
        )

        result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(tmp_path), *extra])

        assert result.exit_code == 1, result.output
        assert working.setup_called is True
        assert "gnomad: download failed" in result.output
        assert "Failed to install or refresh: gnomad" in result.output

    def test_db_update_with_file_url(
        self,
        tmp_path: Path,
        mock_pharmgkb_dir: Path,
        mock_cpic_lookup: dict[tuple[str, str], str],
        mock_gwas_tsv: Path,
        mock_gnomad_gz: Path,
        build_synthetic_clinvar_cache,
        monkeypatch,
    ):
        """Download+ingest path runs end-to-end against file:// URLs.

        All annotators must succeed -- `db update` iterates the whole registry,
        so each needs a local file:// URL or mock fixture. ClinVar setup is
        mocked to call build_synthetic_clinvar_cache; ClinPGx and GWAS use
        mock archives; gnomAD uses a mock gzipped SQLite cache.
        """
        import zipfile
        from urllib.request import pathname2url

        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases.gwas_loader import GWAS_DB_FILENAME
        from allelix.databases.manager import get_database_info

        pharmgkb_zip = tmp_path / "fixture_clinical_annotations.zip"
        with zipfile.ZipFile(pharmgkb_zip, "w") as zf:
            for f in mock_pharmgkb_dir.iterdir():
                zf.write(f, arcname=f.name)
        pharmgkb_url = f"file:{pathname2url(str(pharmgkb_zip.resolve()))}"
        monkeypatch.setattr(pharmgkb_module, "PHARMGKB_CLINICAL_URL", pharmgkb_url)

        # Never hit the real CPIC API from the test suite. The annotator
        # imports fetch_cpic_allele_functions by name, so patching it on the
        # pharmgkb module is what setup() actually resolves at call time.
        monkeypatch.setattr(
            pharmgkb_module,
            "fetch_cpic_allele_functions",
            lambda: dict(mock_cpic_lookup),
        )

        gwas_zip = tmp_path / "fixture_gwas_catalog.zip"
        with zipfile.ZipFile(gwas_zip, "w") as zf:
            zf.write(mock_gwas_tsv, arcname="gwas_catalog_associations.tsv")
        gwas_url = f"file:{pathname2url(str(gwas_zip.resolve()))}"
        monkeypatch.setattr(gwas_module, "GWAS_CATALOG_URL", gwas_url)

        from allelix.annotators import gnomad as gnomad_module

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)

        from allelix.annotators import clinvar as clinvar_mod
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.annotators.gwas import GWASCatalogAnnotator
        from allelix.annotators.pharmgkb import PharmGKBAnnotator

        def _fake_clinvar_setup(self):
            for build in self._builds:
                build_synthetic_clinvar_cache(
                    self._db_paths[build],
                    build,
                    source_url="test://clinvar-mock",
                    remote_signal="md5:test_signal",
                )

        monkeypatch.setattr(ClinVarAnnotator, "setup", _fake_clinvar_setup)
        # cached_remote_signal returns the shared TSV md5: both per-build
        # caches stamp the same value because they ingest one TSV pair.
        monkeypatch.setattr(clinvar_mod, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(
            PharmGKBAnnotator,
            "fetch_remote_signal",
            lambda _self: "pgkb:etag:test|cpic:test",
        )
        monkeypatch.setattr(
            GWASCatalogAnnotator,
            "fetch_remote_signal",
            lambda _self: "etag:test_signal",
        )
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        # The download() helper writes to data_dir/clinvar.vcf.gz and
        # data_dir/clinicalAnnotations.zip. Use a separate cache subdir so the
        # source ZIP fixture above doesn't collide.
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        runner = CliRunner()
        result = runner.invoke(main, ["db", "update", "--data-dir", str(cache_dir)])
        assert result.exit_code == 0, result.output

        clinvar_sqlite_grch37 = cache_dir / clinvar_db_filename("GRCh37")
        clinvar_sqlite_grch38 = cache_dir / clinvar_db_filename("GRCh38")
        pharmgkb_sqlite = cache_dir / "pharmgkb.sqlite"
        gwas_sqlite = cache_dir / GWAS_DB_FILENAME
        assert clinvar_sqlite_grch37.exists()
        assert clinvar_sqlite_grch38.exists()
        assert pharmgkb_sqlite.exists()
        assert gwas_sqlite.exists()
        # Staged downloads cleaned up after ingest; raw artifacts retained for
        # auto-reingest on interpreter version bumps (ZIP for ClinPGx, TSV for GWAS)
        assert not (cache_dir / "clinvar.vcf").exists()
        assert not (cache_dir / "clinvar.vcf.gz").exists()
        assert not (cache_dir / "gwas_catalog_associations.zip").exists()
        assert (cache_dir / "clinicalAnnotations.zip").exists()
        assert (cache_dir / "gwas_catalog_associations.tsv").exists()

        for build, db_path in (
            ("GRCh37", clinvar_sqlite_grch37),
            ("GRCh38", clinvar_sqlite_grch38),
        ):
            info = get_database_info(db_path, clinvar_record_name(build))
            assert info is not None, f"{build} cache missing version row"
            # grew the fixture set: 13 → 20 (added 3 identical-SCV
            # rows for rs999000111, 2 conflicting-SCV rows for rs999000222,
            # 1 sub-floor row for rs999000333, 1 unmapped-CLNSIG row for
            # rs999000444). The multi-allelic dedup fix added 2 more
            # rows for rs999000555 (multi-allelic distinct-alt dedup
            # guard) → 22 total. Cache row count stays at 22 even after
            # the dedup — dedup is annotator-side, not loader-side.
            assert info["record_count"] == 22
            assert info["version"] == "20260101"

        pharmgkb_info = get_database_info(pharmgkb_sqlite, "pharmgkb")
        assert pharmgkb_info is not None
        assert pharmgkb_info["record_count"] == 16

        gwas_info = get_database_info(gwas_sqlite, "gwas")
        assert gwas_info is not None
        assert gwas_info["record_count"] == 8

        gnomad_sqlite = cache_dir / "gnomad.sqlite"
        assert gnomad_sqlite.exists()
        gnomad_info = get_database_info(gnomad_sqlite, "gnomad")
        assert gnomad_info is not None
        assert gnomad_info["record_count"] == 3

    def test_db_update_bad_url_leaves_old_cache_intact(
        self, clinvar_data_dir: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """A failed update must not destroy a working cache.

        Uses --force because the cache is already ready; without it
        `db update` would correctly skip the download and the bad URL
        would never be exercised. ADR-0021: restrict to GRCh37 so the
        test only stress-tests one URL.
        """
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases import manager as manager_module
        from allelix.databases.manager import get_database_info

        sqlite_path = clinvar_data_dir / clinvar_db_filename("GRCh37")
        info_before = get_database_info(sqlite_path, clinvar_record_name("GRCh37"))
        assert info_before is not None

        # The annotator downloads variant_summary +
        # submission_summary instead of the per-build VCF. Pointing
        # either at a refused connection still exercises the
        # "failure preserves cache" invariant.
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_VARIANT_SUMMARY_URL",
            "http://127.0.0.1:1/missing-variant_summary.txt.gz",
        )
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_SUBMISSION_SUMMARY_URL",
            "http://127.0.0.1:1/missing-submission_summary.txt.gz",
        )
        # Need a non-None signal to get past the signal guard and into
        # the actual download attempt.
        from allelix.annotators.clinvar import ClinVarAnnotator

        monkeypatch.setattr(
            ClinVarAnnotator,
            "_fetch_remote_signal_for_tsv",
            staticmethod(lambda: "md5:" + "0" * 32),
        )

        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "db",
                "update",
                "--data-dir",
                str(clinvar_data_dir),
                "--force",
                "--build",
                "grch37",
            ],
        )
        assert result.exit_code == 1, result.output
        assert "clinvar:" in result.output

        info_after = get_database_info(sqlite_path, clinvar_record_name("GRCh37"))
        assert info_after == info_before

    def test_status_shows_pharmgkb_record_count(self, pharmgkb_data_dir: Path):
        """Status row for ClinPGx shows ready + version + record count (no ClinVar here)."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["db", "status", "--data-dir", str(pharmgkb_data_dir)])
        assert result.exit_code == 0
        assert "ClinPGx" in result.output
        assert "yes" in result.output
        assert "16" in result.output

    def test_db_update_skips_when_remote_signal_matches(
        self,
        tmp_path: Path,
        mock_gnomad_gz: Path,
        build_synthetic_clinvar_cache,
        monkeypatch,
    ):
        """Remote signal matches cached signal → skip without --force."""
        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename

        # Pre-populate the GRCh37 ClinVar cache with a known signal.
        build_synthetic_clinvar_cache(
            tmp_path / clinvar_db_filename("GRCh37"),
            "GRCh37",
            source_url="test",
            remote_signal="md5:cached_value",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "fetch_remote_signal",
            lambda self: "md5:cached_value",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "cached_remote_signal",
            lambda self: "md5:cached_value",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "is_ready",
            lambda self: True,
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "setup",
            lambda self: (_ for _ in ()).throw(
                AssertionError("setup() should not run when signal matches")
            ),
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main, ["db", "update", "--data-dir", str(tmp_path), "--build", "grch37"]
        )
        assert result.exit_code == 0, result.output
        assert "already current" in result.output

    def test_db_update_refreshes_when_remote_signal_differs(
        self,
        tmp_path: Path,
        mock_gnomad_gz: Path,
        build_synthetic_clinvar_cache,
        monkeypatch,
    ):
        """Remote signal differs from cached → trigger refresh without --force."""
        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename

        build_synthetic_clinvar_cache(
            tmp_path / clinvar_db_filename("GRCh37"),
            "GRCh37",
            source_url="test",
            remote_signal="md5:OLD",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "fetch_remote_signal",
            lambda self: "md5:NEW",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "cached_remote_signal",
            lambda self: "md5:OLD",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "is_ready",
            lambda self: True,
        )
        called = []
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "setup",
            lambda self: called.append(True),
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main, ["db", "update", "--data-dir", str(tmp_path), "--build", "grch37"]
        )
        assert result.exit_code == 0, result.output
        assert "remote signal changed" in result.output
        assert called  # setup() ran

    def test_db_update_skips_when_signal_unverifiable(
        self, clinvar_data_dir: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """Cache present + can't reach remote → skip with notice (don't crash)."""
        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator, "fetch_remote_signal", lambda self: None
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "setup",
            lambda self: (_ for _ in ()).throw(
                AssertionError("setup() should not run when signal is unverifiable")
            ),
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(main, ["db", "update", "--data-dir", str(clinvar_data_dir)])
        assert result.exit_code == 0, result.output
        assert "can't be verified" in result.output

    def test_db_update_legacy_cache_redownloads(
        self, tmp_path: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """A cache with no stored signal predates the signal mechanism
        (i.e., is old). It must be re-downloaded, NOT silently stamped as
        current — the previous behavior permanently marked stale caches as
        fresh."""
        import sqlite3

        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename

        # Hand-build a legacy database (no remote_signal column).
        legacy_db = tmp_path / clinvar_db_filename("GRCh37")
        with contextlib.closing(sqlite3.connect(legacy_db)) as conn:
            conn.executescript(
                """
                CREATE TABLE clinvar_variants (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    rsid TEXT NOT NULL, chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL, ref TEXT NOT NULL, alt TEXT NOT NULL,
                    clinical_significance TEXT, condition TEXT, gene TEXT,
                    review_status TEXT, allele_id INTEGER
                );
                CREATE TABLE database_versions (
                    name TEXT PRIMARY KEY, source_url TEXT NOT NULL,
                    version TEXT, downloaded_at TEXT NOT NULL,
                    record_count INTEGER NOT NULL
                );
                """
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?)",
                ("clinvar.GRCh37", "old", "old-version", "2024-01-01T00:00:00", 1),
            )
            conn.commit()

        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "fetch_remote_signal",
            lambda self: "md5:NEW",
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "cached_remote_signal",
            lambda self: None,
        )
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "is_ready",
            lambda self: True,
        )
        called = []
        monkeypatch.setattr(
            clinvar_module.ClinVarAnnotator,
            "setup",
            lambda self: called.append(True),
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main, ["db", "update", "--data-dir", str(tmp_path), "--build", "grch37"]
        )
        assert result.exit_code == 0, result.output
        assert "cache predates the freshness signal" in result.output
        assert called  # setup() MUST run — tagless cache is treated as stale

    def test_db_update_force_refreshes(
        self, clinvar_data_dir: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """--force must re-run setup() even when already ready."""
        from urllib.request import pathname2url

        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module

        called = []

        def fake_setup(self):
            called.append(True)

        monkeypatch.setattr(clinvar_module.ClinVarAnnotator, "setup", fake_setup)
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        from allelix.annotators.alphamissense import AlphaMissenseAnnotator

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["db", "update", "--data-dir", str(clinvar_data_dir), "--force"],
        )
        assert result.exit_code == 0, result.output
        assert called  # setup() ran

    def test_db_update_bad_url_shows_friendly_error(
        self, clinvar_data_dir: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """Failure path prints a friendly error — not propagate the raw error.

        Pinned facts:
          1. The error message appears in output with the annotator name prefix.
          2. The raw exception does NOT propagate to the caller.
        """
        import urllib.error
        from urllib.request import pathname2url

        from allelix.annotators import gnomad as gnomad_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator
        from allelix.databases import manager as manager_module

        # TSV URLs replace the per-build VCF map. Refused-connection
        # URLs on both files still exercise the friendly-error invariant.
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_VARIANT_SUMMARY_URL",
            "http://127.0.0.1:1/missing-variant_summary.txt.gz",
        )
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_SUBMISSION_SUMMARY_URL",
            "http://127.0.0.1:1/missing-submission_summary.txt.gz",
        )
        # Signal needs to succeed so the friendly-error path tests the
        # download failure, not the signal-guard failure.
        from allelix.annotators.clinvar import ClinVarAnnotator

        monkeypatch.setattr(
            ClinVarAnnotator,
            "_fetch_remote_signal_for_tsv",
            staticmethod(lambda: "md5:" + "0" * 32),
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
        monkeypatch.setattr(gnomad_module, "verify_file_hash", lambda *_a, **_kw: None)
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "db",
                "update",
                "--data-dir",
                str(clinvar_data_dir),
                "--force",
                "--build",
                "grch37",
            ],
        )
        assert result.exit_code == 1, result.output
        assert "clinvar:" in result.output
        assert not isinstance(result.exception, urllib.error.URLError)

    def test_db_path_prints_resolved_dir(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["db", "path", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert str(tmp_path) in result.output

    def test_db_path_no_styling_escape_codes(self, tmp_path: Path):
        """Output must be shell-safe for $(allelix db path) capture."""
        runner = CliRunner()
        result = runner.invoke(main, ["db", "path", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        # No ANSI escape sequences — Rich's color codes start with \x1b[
        assert "\x1b[" not in result.output

    def test_db_path_check_passes_on_writable_dir(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["db", "path", "--data-dir", str(tmp_path), "--check"])
        assert result.exit_code == 0

    def test_db_path_check_fails_on_readonly_dir(self, tmp_path: Path):
        import os

        from test_readonly_caches import _read_only

        # chmod 0 is bypassed by root — root sees a "readonly"
        # dir as writable (DAC override) and the test silently passes
        # for the wrong reason. CI runners and dev containers commonly
        # run as root; skip there with a clear message so the test
        # isn't a false-green on those environments. Hermetic
        # local-dev runs (non-root) keep the real assertion.
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            pytest.skip(
                "chmod-readonly is bypassed under root (DAC override); "
                "test is non-hermetic in that environment. Run as a "
                "non-root user to exercise the writability check."
            )
        with _read_only(tmp_path):
            runner = CliRunner()
            result = runner.invoke(main, ["db", "path", "--data-dir", str(tmp_path), "--check"])
            assert result.exit_code == 1
            # Error goes to stderr (plain print, not console.print); normalize
            # whitespace so the assertion survives any line-wrapping.
            stderr_flat = " ".join(result.stderr.split())
            assert "not writable" in stderr_flat

    def test_db_clean_empty_dir_is_noop(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
        assert result.exit_code == 0
        assert "Nothing to clean" in result.output

    def test_db_clean_missing_dir_resolves_to_empty_noop(self, tmp_path: Path):
        # resolve_data_dir creates the dir if missing; clean then finds
        # nothing to delete. Verify the directory is created AND the
        # command exits cleanly.
        runner = CliRunner()
        missing = tmp_path / "does-not-exist"
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(missing), "--yes"])
        assert result.exit_code == 0
        assert "Nothing to clean" in result.output
        assert missing.exists() and missing.is_dir()

    def test_db_clean_removes_caches_with_yes(self, tmp_path: Path):
        # Stage mock cache files + a sub-dir + a config.toml that MUST survive
        (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"x" * 1024)
        (tmp_path / "gnomad.sqlite").write_bytes(b"y" * 2048)
        (tmp_path / "cadd.sqlite").write_bytes(b"z" * 4096)
        nested = tmp_path / "snpedia"
        nested.mkdir()
        (nested / "page.txt").write_text("snpedia raw")
        (tmp_path / "config.toml").write_text("[sources]\nclinvar = true\n")

        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
        assert result.exit_code == 0
        assert "Removed" in result.output
        # Cache files gone
        assert not (tmp_path / "clinvar.GRCh38.sqlite").exists()
        assert not (tmp_path / "gnomad.sqlite").exists()
        assert not (tmp_path / "cadd.sqlite").exists()
        assert not nested.exists()
        # config.toml preserved
        assert (tmp_path / "config.toml").exists()
        assert (tmp_path / "config.toml").read_text() == "[sources]\nclinvar = true\n"

    def test_db_clean_keep_cadd_preserves_cadd_cache(self, tmp_path: Path):
        (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"x" * 1024)
        (tmp_path / "cadd.sqlite").write_bytes(b"z" * 4096)
        (tmp_path / "config.toml").write_text("")

        runner = CliRunner()
        result = runner.invoke(
            main, ["db", "clean", "--data-dir", str(tmp_path), "--keep-cadd", "--yes"]
        )
        assert result.exit_code == 0
        assert not (tmp_path / "clinvar.GRCh38.sqlite").exists()
        assert (tmp_path / "cadd.sqlite").exists()
        assert (tmp_path / "config.toml").exists()

    def test_db_clean_dry_run_deletes_nothing(self, tmp_path: Path):
        (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"x" * 1024)
        (tmp_path / "cadd.sqlite").write_bytes(b"z" * 4096)

        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "Dry run" in result.output
        assert "Would delete" in result.output
        assert (tmp_path / "clinvar.GRCh38.sqlite").exists()
        assert (tmp_path / "cadd.sqlite").exists()

    def test_db_clean_confirmation_prompt_blocks_on_no(self, tmp_path: Path):
        (tmp_path / "clinvar.GRCh38.sqlite").write_bytes(b"x" * 1024)

        runner = CliRunner()
        # Send "n\n" to the confirmation prompt
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path)], input="n\n")
        assert result.exit_code == 0
        assert "Aborted" in result.output
        assert (tmp_path / "clinvar.GRCh38.sqlite").exists()

    def test_db_clean_keep_cadd_with_only_cadd_present_is_noop(self, tmp_path: Path):
        (tmp_path / "cadd.sqlite").write_bytes(b"z" * 4096)

        runner = CliRunner()
        result = runner.invoke(
            main, ["db", "clean", "--data-dir", str(tmp_path), "--keep-cadd", "--yes"]
        )
        assert result.exit_code == 0
        # Rich wraps at 80 columns and the message holds the tmp path.
        output = " ".join(result.output.split())
        assert "Nothing to clean" in output
        assert "CADD preserved" in output
        assert (tmp_path / "cadd.sqlite").exists()

    def test_db_clean_refuses_unrecognized_dir_even_with_yes(self, tmp_path: Path):
        # Fat-fingered --data-dir pointing at a user's documents folder. No
        # recognized allelix artifacts. The --yes path MUST still refuse.
        (tmp_path / "resume.pdf").write_bytes(b"x" * 1024)
        (tmp_path / "tax_2025.xlsx").write_bytes(b"y" * 2048)
        subdir = tmp_path / "vacation_photos"
        subdir.mkdir()
        (subdir / "IMG_0001.jpg").write_bytes(b"z" * 4096)

        # Click 8.4 separates result.stdout / result.stderr automatically;
        # the guard error routes to stderr (matches `db path --check` pattern).
        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
        assert result.exit_code == 1
        # Normalize whitespace so the assertion survives any line-wrapping —
        # the resolved tmp_path can push the phrase across terminal-width
        # boundaries on narrow runners. (console.print() was
        # Rich-wrapping mid-phrase under CliRunner's 80-col default.)
        stderr_flat = " ".join(result.stderr.split())
        assert "doesn't look like an allelix cache" in stderr_flat
        assert "--force" in stderr_flat
        # Everything still present
        assert (tmp_path / "resume.pdf").exists()
        assert (tmp_path / "tax_2025.xlsx").exists()
        assert (subdir / "IMG_0001.jpg").exists()

    def test_db_clean_force_bypasses_guard(self, tmp_path: Path):
        # User explicitly opts in via --force.
        (tmp_path / "non_standard_cache.bin").write_bytes(b"x" * 1024)

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["db", "clean", "--data-dir", str(tmp_path), "--force", "--yes"],
        )
        assert result.exit_code == 0
        assert "Removed" in result.output
        assert not (tmp_path / "non_standard_cache.bin").exists()

    def test_db_clean_config_toml_alone_is_recognition_marker(self, tmp_path: Path):
        # A data dir with ONLY config.toml means a fresh install pre-update.
        # entries is empty (config.toml is preserved), so noop path applies.
        (tmp_path / "config.toml").write_text("[sources]\n")

        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
        assert result.exit_code == 0
        assert "Nothing to clean" in result.output
        assert (tmp_path / "config.toml").exists()

    def test_db_clean_dry_run_also_blocks_unrecognized_dir(self, tmp_path: Path):
        # The guard runs BEFORE --dry-run so even previewing against an
        # unrelated tree is refused — would print a misleading "would
        # delete" listing otherwise.
        (tmp_path / "random_file.txt").write_text("not allelix")

        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--dry-run"])
        assert result.exit_code == 1
        stderr_flat = " ".join(result.stderr.split())
        assert "doesn't look like an allelix cache" in stderr_flat

    def test_db_clean_refuses_dir_with_only_generic_config_toml(self, tmp_path: Path):
        """A `config.toml` alongside non-allelix content (the
        Rust / Hugo / Rocket / etc. project shape) must NOT pass the
        guard. Previously, the guard accepted any dir with a config.toml
        — so a typo'd `--data-dir ~/my-rust-project` could trigger
        `rm -rf` against a Cargo project root.

        Pinned: a generic config.toml + any non-allelix file →
        guard refuses with the "doesn't look like an allelix cache"
        error and exit code 1, regardless of what's in config.toml.
        """
        # A Rust-shaped config.toml (no allelix-specific sections).
        (tmp_path / "config.toml").write_text("[package]\nname = 'not-allelix'\n")
        (tmp_path / "Cargo.lock").write_text("# cargo lockfile\n")
        (tmp_path / "src").mkdir()
        (tmp_path / "src" / "main.rs").write_text("fn main() {}\n")

        runner = CliRunner()
        result = runner.invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
        assert result.exit_code == 1
        stderr_flat = " ".join(result.stderr.split())
        assert "doesn't look like an allelix cache" in stderr_flat
        # And the Cargo project survives the encounter unharmed.
        assert (tmp_path / "config.toml").exists()
        assert (tmp_path / "Cargo.lock").exists()
        assert (tmp_path / "src" / "main.rs").exists()


class TestAnalyzeOutputDispatch:
    def test_analyze_writes_json(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        out = tmp_path / "report.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        import json as _json

        payload = _json.loads(out.read_text())
        assert payload["input"]["sample_id"] == "MHG000001"
        assert any(a["attribution"] == "ClinVar" for a in payload["annotations"])

    def test_analyze_writes_html(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        out = tmp_path / "report.html"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        body = out.read_text()
        assert body.startswith("<!DOCTYPE html>")
        assert "rs1801133" in body
        assert "Informational only" in body

    def test_analyze_unknown_extension_errors(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        out = tmp_path / "report.xyz"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code != 0
        assert "Cannot infer report format" in result.output

    def test_analyze_report_format_override(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """--report-format json forces JSON regardless of file extension."""
        out = tmp_path / "report.weird"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(out),
                "--report-format",
                "json",
            ],
        )
        assert result.exit_code == 0, result.output
        import json as _json

        _json.loads(out.read_text())  # parses cleanly


class TestMethylationCommand:
    def test_filters_to_methylation_panel(self, mock_mhg_path, all_annotators_data_dir: Path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["methylation", str(mock_mhg_path), "--data-dir", str(all_annotators_data_dir)],
        )
        assert result.exit_code == 0, result.output
        # MTHFR is in the panel; BRCA1 is not.
        assert "MTHFR" in result.output
        assert "BRCA1" not in result.output

    def test_gwas_excluded_by_default(self, mock_mhg_path, all_annotators_data_dir: Path):
        """Methylation excludes GWAS by default — biology from ClinVar + ClinPGx only."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "methylation",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "GWAS" not in result.output
        assert "Homocysteine levels" not in result.output

    def test_include_gwas_flag(self, mock_mhg_path, all_annotators_data_dir: Path):
        """--include-gwas + --gwas-all restores all GWAS annotations in methylation report.

        --include-gwas re-enables the GWAS annotator (excluded by default
        in focused reports). --gwas-all disables trait-category filtering
        so measurement traits like "Homocysteine levels" appear.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "methylation",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
                "--include-gwas",
                "--gwas-all",
                "--gwas-min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        # Width-independent assertion — at narrow widths Rich wraps the
        # Condition column cell into "Homocysteine │ │ │ │ levels" (the
        # box-drawing chars and whitespace from neighboring columns
        # interleave between the two words). Strip everything except
        # letters and spaces, then collapse whitespace, then check.
        # Assertions on wrapped Rich output are width-fragile (see cli/db.py).
        import re

        letters_only = re.sub(r"[^A-Za-z ]+", " ", result.output)
        flat = " ".join(letters_only.split())
        assert "Homocysteine levels" in flat

    def test_default_output_below_threshold(self, mock_mhg_path, all_annotators_data_dir: Path):
        """Default methylation (no flags) on canonical mocks stays under 20 rows."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["methylation", str(mock_mhg_path), "--data-dir", str(all_annotators_data_dir)],
        )
        assert result.exit_code == 0, result.output
        lines = [
            line for line in result.output.splitlines() if "rs" in line.lower() and "│" in line
        ]
        assert len(lines) < 20, (
            f"Methylation default produced {len(lines)} annotation rows, "
            f"exceeding the 20-row sanity threshold."
        )

    def test_default_floor_surfaces_loe_3_hits(self, mock_mhg_path, all_annotators_data_dir: Path):
        """Methylation default --min-magnitude is 3.0, not analyze's 5.0.

        ClinPGx LoE 3 hits (magnitude 4.0) on COMT / MTHFR / MTR / CBS are
        exactly the signal a user running the focused subcommand wants to
        see. The old 5.0 default — inherited from the noise-prone general
        report — filtered them out.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["methylation", str(mock_mhg_path), "--data-dir", str(all_annotators_data_dir)],
        )
        assert result.exit_code == 0, result.output
        # loe_3 → "Reduced response" cell in the mock; magnitude 4.0.
        # Would not appear at the old 5.0 floor.
        assert "loe_3" in result.output

    def test_explicit_min_magnitude_5_still_filters_loe_3(
        self, mock_mhg_path, all_annotators_data_dir: Path
    ):
        """Users keep the escape hatch — passing --min-magnitude 5 restores the
        old behavior. The default changed, the flag didn't."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "methylation",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "5",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "loe_3" not in result.output


class TestPharmacogenomicsCommand:
    def test_filters_to_pharma_category(self, mock_mhg_path, all_annotators_data_dir: Path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "pharmacogenomics",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        # ClinPGx rows have category=pharma; ClinVar rows don't.
        # terminal strips the source_ prefix; check the stripped form.
        assert "ClinPGx" in result.output
        assert "loe_" in result.output
        # No ClinVar row should surface. ``ClinVar`` itself appears in the
        # database-status header banner regardless of filter, so check the
        # raw significance string instead — it would only appear in a
        # rendered table cell, which the pharma filter excludes.
        assert "clinvar_pathogenic" not in result.output

    def test_gwas_excluded_by_default(self, mock_mhg_path, all_annotators_data_dir: Path):
        """Pharmacogenomics excludes GWAS by default."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "pharmacogenomics",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "GWAS" not in result.output


class TestExtractCommand:
    """v0.5.2: spot-check diploid genotypes at specific rsids."""

    def test_extract_known_carrier(self, mock_mhg_path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["extract", str(mock_mhg_path), "--snps", "rs1801133"])
        assert result.exit_code == 0, result.output
        assert "rs1801133" in result.output
        # MHG fixture has rs1801133 as G/A heterozygous
        assert "G/A" in result.output
        assert "yes" in result.output  # heterozygous

    def test_extract_multiple_snps(self, mock_mhg_path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            ["extract", str(mock_mhg_path), "--snps", "rs1801133,rs4680,rs113993960"],
        )
        assert result.exit_code == 0, result.output
        assert "rs1801133" in result.output
        assert "rs4680" in result.output
        # rs113993960 is no-call in the v0.5.1-corrected fixture
        assert "rs113993960" in result.output
        assert "-/-" in result.output

    def test_extract_rsid_not_in_file(self, mock_mhg_path):
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["extract", str(mock_mhg_path), "--snps", "rs999000000"])
        assert result.exit_code == 0, result.output
        assert "not in file" in result.output

    def test_extract_empty_snps_errors(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(main, ["extract", str(mock_mhg_path), "--snps", ""])
        assert result.exit_code != 0
        assert "cannot be empty" in result.output


class TestDiffOption:
    """Integration tests for `--diff previous.json`."""

    def test_diff_no_changes_terminal(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """Same file analyzed twice with same filters → 'No changes'."""

        baseline = tmp_path / "baseline.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        assert baseline.exists()

        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "No changes" in result.output

    def test_diff_detects_removed_annotations(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """Raising min-magnitude should cause some annotations to disappear from current."""

        baseline = tmp_path / "baseline.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "0",
                "--output",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output

        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--min-magnitude",
                "9",
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "removed" in result.output

    def test_diff_json_output_includes_diff_key(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """JSON output with --diff includes a 'diff' key."""
        import json as _json

        baseline = tmp_path / "baseline.json"
        diffed = tmp_path / "diffed.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(baseline),
            ],
        )
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(diffed),
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        payload = _json.loads(diffed.read_text())
        assert "diff" in payload
        assert "summary" in payload["diff"]

    def test_diff_html_output_includes_diff_banner(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """HTML output with --diff includes the diff summary banner."""

        baseline = tmp_path / "baseline.json"
        diffed = tmp_path / "diffed.html"
        runner = CliRunner(env={"COLUMNS": "200"})
        runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(baseline),
            ],
        )
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(diffed),
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        body = diffed.read_text()
        assert "Diff:" in body
        assert "No changes" in body

    def test_diff_invalid_json_errors(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """--diff with an invalid JSON file → ClickException."""
        bad = tmp_path / "bad.json"
        bad.write_text("not json {{{")
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--diff",
                str(bad),
            ],
        )
        assert result.exit_code != 0
        assert "Cannot parse" in result.output

    def test_diff_wrong_schema_version_errors(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """--diff with a future schema version → ClickException."""
        import json as _json

        bad = tmp_path / "future.json"
        bad.write_text(_json.dumps({"schema_version": "99", "annotations": []}))
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--diff",
                str(bad),
            ],
        )
        assert result.exit_code != 0
        assert "schema version" in result.output

    def test_diff_on_methylation_command(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """--diff works on the methylation subcommand."""

        baseline = tmp_path / "baseline.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        runner.invoke(
            main,
            [
                "methylation",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(baseline),
            ],
        )
        result = runner.invoke(
            main,
            [
                "methylation",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "No changes" in result.output

    def test_diff_on_pharmacogenomics_command(
        self, mock_mhg_path, all_annotators_data_dir: Path, tmp_path: Path
    ):
        """--diff works on the pharmacogenomics subcommand."""

        baseline = tmp_path / "baseline.json"
        runner = CliRunner(env={"COLUMNS": "200"})
        runner.invoke(
            main,
            [
                "pharmacogenomics",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--output",
                str(baseline),
            ],
        )
        result = runner.invoke(
            main,
            [
                "pharmacogenomics",
                str(mock_mhg_path),
                "--data-dir",
                str(all_annotators_data_dir),
                "--diff",
                str(baseline),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "No changes" in result.output


class TestExcludeSnpedia:
    """--exclude-snpedia wires through to exclude_sources on all three commands."""

    def test_analyze_passes_exclude_sources(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--exclude-snpedia"])
        assert result.exit_code == 0, result.output
        assert captured["exclude_sources"] == frozenset({"snpedia"})

    def test_methylation_passes_exclude_sources(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["methylation", str(mock_mhg_path), "--exclude-snpedia"])
        assert result.exit_code == 0, result.output
        assert "snpedia" in captured["exclude_sources"]

    def test_pharmacogenomics_passes_exclude_sources(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["pharmacogenomics", str(mock_mhg_path), "--exclude-snpedia"])
        assert result.exit_code == 0, result.output
        assert "snpedia" in captured["exclude_sources"]

    def test_methylation_exclude_snpedia_and_gwas(self, mock_mhg_path, monkeypatch):
        """Both --exclude-snpedia and default GWAS exclusion coexist."""
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["methylation", str(mock_mhg_path), "--exclude-snpedia"])
        assert result.exit_code == 0, result.output
        assert captured["exclude_sources"] == frozenset({"snpedia", "gwas"})


class TestNoCaddFlag:
    """--no-cadd wires through to no_cadd on all three analysis commands."""

    def test_analyze_passes_no_cadd(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--no-cadd"])
        assert result.exit_code == 0, result.output
        assert captured["no_cadd"] is True

    def test_analyze_default_no_cadd_false(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        assert captured["no_cadd"] is False

    def test_methylation_passes_no_cadd(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["methylation", str(mock_mhg_path), "--no-cadd"])
        assert result.exit_code == 0, result.output
        assert captured["no_cadd"] is True

    def test_pharmacogenomics_passes_no_cadd(self, mock_mhg_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["pharmacogenomics", str(mock_mhg_path), "--no-cadd"])
        assert result.exit_code == 0, result.output
        assert captured["no_cadd"] is True


class TestParseFilterFile:
    """Unit tests for _parse_filter_file (parser classification)."""

    def test_rsid_lowercase(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("rs1801133\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset()
        assert rsids == frozenset({"rs1801133"})

    def test_rsid_uppercase_normalized_to_lowercase(self, tmp_path):
        """Input case-tolerant, output canonical: RS1801133 → rs1801133."""
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("RS1801133\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset()
        assert rsids == frozenset({"rs1801133"})

    def test_gene_lowercase_normalized_to_uppercase(self, tmp_path):
        """Input case-tolerant, output canonical: mthfr → MTHFR."""
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("mthfr\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"MTHFR"})
        assert rsids == frozenset()

    def test_mixed_messy_case_normalized(self, tmp_path):
        """End-to-end case-mixing across rsIDs and genes."""
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("Rs1801133\ncomt\nRSPO1\nRS4680\nmThFr\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"COMT", "RSPO1", "MTHFR"})
        assert rsids == frozenset({"rs1801133", "rs4680"})

    def test_gene_only(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("MTHFR\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"MTHFR"})
        assert rsids == frozenset()

    def test_gene_starting_with_rs_prefix_is_gene_not_rsid(self, tmp_path):
        """RSPO1, RSF1, RSC1A1 are real gene names — must not be classified as rsIDs."""
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("RSPO1\nRSF1\nRSC1A1\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"RSPO1", "RSF1", "RSC1A1"})
        assert rsids == frozenset()

    def test_mixed(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("rs1801133\nMTHFR\nrs4680\nCOMT\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"MTHFR", "COMT"})
        assert rsids == frozenset({"rs1801133", "rs4680"})

    def test_comments_and_blanks_ignored(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("# this is a comment\n\nMTHFR\n\n# another\nrs1801133\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset({"MTHFR"})
        assert rsids == frozenset({"rs1801133"})

    def test_empty_file_returns_empty_sets(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset()
        assert rsids == frozenset()

    def test_comments_only_returns_empty_sets(self, tmp_path):
        from allelix.cli._helpers import _parse_filter_file

        f = tmp_path / "filter.txt"
        f.write_text("# only a comment\n# another\n\n")
        genes, rsids = _parse_filter_file(f)
        assert genes == frozenset()
        assert rsids == frozenset()


class TestFilterFileOnAnalyze:
    """--filter-file is only on analyze; threads through _run_analysis_command."""

    def test_analyze_rsid_only(self, mock_mhg_path, tmp_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        f = tmp_path / "filter.txt"
        f.write_text("rs1801133\n")
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--filter-file", str(f)])
        assert result.exit_code == 0, result.output
        assert captured["genes"] == frozenset()
        assert captured["rsids"] == frozenset({"rs1801133"})

    def test_analyze_gene_only(self, mock_mhg_path, tmp_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        f = tmp_path / "filter.txt"
        f.write_text("MTHFR\n")
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--filter-file", str(f)])
        assert result.exit_code == 0, result.output
        assert captured["genes"] == frozenset({"MTHFR"})
        assert captured["rsids"] == frozenset()

    def test_analyze_mixed_or_combination(self, mock_mhg_path, tmp_path, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        f = tmp_path / "filter.txt"
        f.write_text("rs1801133\nCOMT\n")
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--filter-file", str(f)])
        assert result.exit_code == 0, result.output
        assert captured["genes"] == frozenset({"COMT"})
        assert captured["rsids"] == frozenset({"rs1801133"})

    def test_analyze_empty_filter_passes_empty_sets(self, mock_mhg_path, tmp_path, monkeypatch):
        """Empty filter file (only comments/blanks) threads empty frozensets through.

        The empty-set → match-nothing semantic on AnalysisResult.filter()
        is covered by a direct unit test in tests/reports/test_pipeline.py;
        here we verify only that the CLI layer forwards empty frozensets,
        not None.
        """
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        f = tmp_path / "filter.txt"
        f.write_text("# only comments\n\n")
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--filter-file", str(f)])
        assert result.exit_code == 0, result.output
        assert captured["genes"] == frozenset()
        assert captured["rsids"] == frozenset()

    def test_analyze_filter_file_nonexistent_path_errors(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["analyze", str(mock_mhg_path), "--filter-file", "/does/not/exist.txt"],
        )
        assert result.exit_code != 0

    def test_methylation_does_not_have_filter_file(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["methylation", str(mock_mhg_path), "--filter-file", "/tmp/x.txt"],
        )
        assert result.exit_code != 0
        assert "no such option" in result.output.lower()

    def test_pharmacogenomics_does_not_have_filter_file(self, mock_mhg_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["pharmacogenomics", str(mock_mhg_path), "--filter-file", "/tmp/x.txt"],
        )
        assert result.exit_code != 0
        assert "no such option" in result.output.lower()


class TestSampleOptionOnVcf:
    """``--sample`` is wired into analyze, methylation, and pharmacogenomics.

    Pre-flight checks happen in ``_resolve_parser``: a multi-sample VCF
    without ``--sample`` fails before the analysis pipeline kicks off,
    so the user sees a clear error instead of a mid-stream
    MultiSampleError exception.
    """

    @pytest.fixture
    def multisample_vcf(self):
        from pathlib import Path

        return Path("tests/fixtures/mock_multisample.vcf")

    @pytest.fixture
    def single_sample_vcf(self):
        from pathlib import Path

        return Path("tests/fixtures/mock_vcf.vcf")

    def test_analyze_multisample_without_sample_fails(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(multisample_vcf)])
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output
        assert "SAMPLE_A" in result.output
        assert "SAMPLE_B" in result.output
        assert "SAMPLE_C" in result.output

    def test_analyze_multisample_with_unknown_sample_fails(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(multisample_vcf), "--sample", "NOT_A_SAMPLE"])
        assert result.exit_code != 0
        assert "NOT_A_SAMPLE" in result.output
        assert "Available samples" in result.output

    def test_analyze_with_sample_passes_to_helper(self, multisample_vcf, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(multisample_vcf), "--sample", "SAMPLE_B"])
        assert result.exit_code == 0, result.output
        assert captured["sample"] == "SAMPLE_B"

    def test_methylation_has_sample_option(self, multisample_vcf, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["methylation", str(multisample_vcf), "--sample", "SAMPLE_A"])
        assert result.exit_code == 0, result.output
        assert captured["sample"] == "SAMPLE_A"

    def test_pharmacogenomics_has_sample_option(self, multisample_vcf, monkeypatch):
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(
            main, ["pharmacogenomics", str(multisample_vcf), "--sample", "SAMPLE_C"]
        )
        assert result.exit_code == 0, result.output
        assert captured["sample"] == "SAMPLE_C"

    def test_analyze_single_sample_vcf_no_sample_works(self, single_sample_vcf, monkeypatch):
        """Single-sample VCF doesn't need --sample (auto-picks the lone sample)."""
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(single_sample_vcf)])
        assert result.exit_code == 0, result.output
        assert captured["sample"] is None

    def test_analyze_sample_on_non_vcf_silently_ignored(self, mock_mhg_path, monkeypatch):
        """``--sample`` on an array file is silently ignored (no error).

        The option is VCF-specific by design but doesn't fail-loud when
        applied to a non-VCF format — the array parsers ignore it via
        ``_resolve_parser``.
        """
        captured: dict = {}

        def fake_run(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr("allelix.cli._helpers._run_analysis_command", fake_run)
        runner = CliRunner()
        result = runner.invoke(main, ["analyze", str(mock_mhg_path), "--sample", "ignored_value"])
        assert result.exit_code == 0, result.output


class TestExtractSampleOnVcf:
    """``allelix extract --snps`` accepts ``--sample`` for multi-sample VCFs."""

    @pytest.fixture
    def multisample_vcf(self):
        from pathlib import Path

        return Path("tests/fixtures/mock_multisample.vcf")

    def test_multisample_without_sample_fails(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(main, ["extract", str(multisample_vcf), "--snps", "rs1801133"])
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output

    def test_multisample_with_sample_works(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "extract",
                str(multisample_vcf),
                "--snps",
                "rs1801133",
                "--sample",
                "SAMPLE_A",
            ],
        )
        assert result.exit_code == 0, result.output
        # SAMPLE_A is 0/1 for rs1801133 (G/A) — heterozygous
        assert "rs1801133" in result.output

    def test_sequential_path_works_on_plain_vcf(self):
        """Without .tbi + pysam + gnomAD, sequential scan must still work."""
        from pathlib import Path

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "extract",
                str(Path("tests/fixtures/mock_vcf.vcf")),
                "--snps",
                "rs1801133,rs4680",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "rs1801133" in result.output
        assert "rs4680" in result.output


class TestExtractTabixFastPath:
    """The tabix fast-path takes over when pysam + .tbi + gnomAD are available.

    Tests use a controlled probe of the fast-path helper rather than
    invoking the CLI end-to-end — the helper is the integration point
    and its inputs are easy to mock.
    """

    def test_tabix_helper_returns_none_without_gz_extension(self, tmp_path):
        """Plain .vcf (no .gz) doesn't go through the tabix path."""
        from allelix.cli.utility import _try_tabix_extract

        f = tmp_path / "foo.vcf"
        f.write_text("##fileformat=VCFv4.2\n")
        result = _try_tabix_extract(f, {"rs1801133"}, None, None)
        assert result is None

    def test_tabix_helper_returns_none_without_tbi(self, tmp_path):
        """Missing .tbi index means no tabix fast-path."""
        from allelix.cli.utility import _try_tabix_extract

        f = tmp_path / "foo.vcf.gz"
        f.write_bytes(b"\x1f\x8b" + b"\x00" * 10)  # minimal gzip-shaped file
        # No .tbi alongside
        result = _try_tabix_extract(f, {"rs1801133"}, None, None)
        assert result is None


class TestSampleOnUtilityCommands:
    """``--sample`` threading on stats / compare / export plink.

    These three commands originally shipped without ``--sample`` even
    though they all hit ``_resolve_parser`` for VCF input. Multi-sample
    VCFs raised mid-pipeline. Pre-flight validation in
    ``_resolve_parser`` now catches the case clean.
    """

    @pytest.fixture
    def multisample_vcf(self):
        from pathlib import Path

        return Path("tests/fixtures/mock_multisample.vcf")

    def test_stats_without_sample_fails_on_multisample(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(multisample_vcf)])
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output

    def test_stats_with_sample_works(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(multisample_vcf), "--sample", "SAMPLE_A"])
        assert result.exit_code == 0, result.output

    def test_compare_without_sample1_fails_on_multisample_left(
        self, multisample_vcf, mock_mhg_path
    ):
        runner = CliRunner()
        result = runner.invoke(main, ["compare", str(multisample_vcf), str(mock_mhg_path)])
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output

    def test_compare_without_sample2_fails_on_multisample_right(
        self, multisample_vcf, mock_mhg_path
    ):
        runner = CliRunner()
        result = runner.invoke(main, ["compare", str(mock_mhg_path), str(multisample_vcf)])
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output

    def test_compare_with_both_samples_works(self, multisample_vcf):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "compare",
                str(multisample_vcf),
                str(multisample_vcf),
                "--sample1",
                "SAMPLE_A",
                "--sample2",
                "SAMPLE_B",
            ],
        )
        assert result.exit_code == 0, result.output

    def test_export_plink_without_sample_fails_on_multisample(self, multisample_vcf, tmp_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "export",
                "plink",
                str(multisample_vcf),
                "-o",
                str(tmp_path / "out"),
            ],
        )
        assert result.exit_code != 0
        assert "Multi-sample VCF" in result.output


class TestHighValueNoCalls:
    def test_stats_flags_dpyd_no_call(self, mock_mhg_path):
        """The MHG fixture has rs3918290 (DPYD) as a no-call; stats should flag it."""
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        assert "High-value no-calls" in result.output
        assert "rs3918290" in result.output
        assert "DPYD" in result.output

    def test_analyze_flags_dpyd_no_call(self, mock_mhg_path, clinvar_data_dir):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--min-magnitude",
                "0",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "high-value" in result.output.lower()
        assert "rs3918290" in result.output


class TestVersion:
    def test_version_flag(self):
        from allelix import __version__

        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        # Self-referential: bumping pyproject.toml shouldn't require touching
        # this test. The pyproject<->__version__ sync is pinned separately by
        # tests/test_version.py::test_pyproject_version_matches_metadata.
        assert __version__ in result.output


class TestConfigCommands:
    def test_config_show(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["config", "show", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "clinvar" in result.output
        assert "personal" in result.output

    def test_config_show_commercial_mode(self, tmp_path: Path):
        from allelix.config import AllelixConfig, save_config

        save_config(tmp_path, AllelixConfig(commercial=True))
        runner = CliRunner()
        result = runner.invoke(main, ["config", "show", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "commercial" in result.output

    def test_config_set_source(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "sources.gnomad", "false"],
        )
        assert result.exit_code == 0
        from allelix.config import load_config

        cfg = load_config(tmp_path)
        assert not cfg.sources["gnomad"]

    def test_config_set_commercial(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "license.commercial", "true"],
        )
        assert result.exit_code == 0
        from allelix.config import load_config

        cfg = load_config(tmp_path)
        assert cfg.commercial

    def test_config_set_invalid_key(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "bad.key", "true"],
        )
        assert result.exit_code != 0
        assert "Unknown key" in result.output

    def test_config_set_invalid_value(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "sources.clinvar", "yes"],
        )
        assert result.exit_code != 0
        assert "true" in result.output or "false" in result.output


class TestConfigGetCommand:
    def test_config_get_no_key_dumps_toml(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["config", "get", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "[sources]" in result.output
        assert "[license]" in result.output
        assert "[options]" in result.output

    def test_config_get_sources_cadd(self, tmp_path: Path):
        runner = CliRunner()
        runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "sources.cadd", "true"],
        )
        result = runner.invoke(
            main, ["config", "get", "--data-dir", str(tmp_path), "sources.cadd"]
        )
        assert result.exit_code == 0
        assert result.output.strip() == "true"

    def test_config_get_options_cadd_full(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main, ["config", "get", "--data-dir", str(tmp_path), "options.cadd_full"]
        )
        assert result.exit_code == 0
        assert result.output.strip() == "false"

    def test_config_get_license_commercial(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main, ["config", "get", "--data-dir", str(tmp_path), "license.commercial"]
        )
        assert result.exit_code == 0
        assert result.output.strip() == "false"

    def test_config_get_unknown_key(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(main, ["config", "get", "--data-dir", str(tmp_path), "foo.bar"])
        assert result.exit_code != 0
        assert "Unknown key" in result.output

    def test_config_get_unknown_source(self, tmp_path: Path):
        runner = CliRunner()
        result = runner.invoke(
            main, ["config", "get", "--data-dir", str(tmp_path), "sources.nonexistent"]
        )
        assert result.exit_code != 0
        assert "Unknown source" in result.output

    def test_config_get_license_source(self, tmp_path: Path):
        from allelix.config import AllelixConfig, save_config

        save_config(tmp_path, AllelixConfig(license_overrides={"cadd": True}))
        runner = CliRunner()
        result = runner.invoke(
            main, ["config", "get", "--data-dir", str(tmp_path), "license.cadd"]
        )
        assert result.exit_code == 0
        assert result.output.strip() == "true"


class TestLicensableGating:
    def test_block_purchasable_message_contains_url(self, tmp_path: Path):
        """Commercial mode, CADD not asserted → config show contains purchase URL."""
        from allelix.config import AllelixConfig, save_config

        save_config(tmp_path, AllelixConfig(commercial=True, sources={"cadd": True}))
        runner = CliRunner(env={"COLUMNS": "300"})
        result = runner.invoke(main, ["config", "show", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "els2.comotion.uw.edu" in result.output

    def test_consent_notice_branches_on_license_held(self):
        """license_held=True auto-accepts and does NOT show non-commercial affirmation."""
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers
        from allelix.cli import db as db_mod

        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf)
        db_mod.console = _helpers.console
        try:
            result = db_mod._confirm_cadd_license(license_held=True)
        finally:
            _helpers.console = original
            db_mod.console = original
        assert result is True
        output = buf.getvalue()
        assert "non-commercial" not in output.lower()
        assert "Commercial license asserted" in output

    def test_block_final_message_no_license_available(self, tmp_path: Path):
        """Commercial mode, SNPedia → config show says no commercial license available."""
        from allelix.config import AllelixConfig, save_config

        save_config(tmp_path, AllelixConfig(commercial=True))
        runner = CliRunner(env={"COLUMNS": "300"})
        result = runner.invoke(main, ["config", "show", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0
        assert "no commercial license is available" in result.output

    def test_config_set_license_non_licensable_rejected(self, tmp_path: Path):
        """Setting license.snpedia true is rejected — SNPedia is not licensable."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "license.snpedia", "true"],
        )
        assert result.exit_code != 0
        assert "not commercially licensable" in result.output

    def test_config_set_license_false_pops_key(self, tmp_path: Path):
        """Setting license.cadd false removes the key from the serialized config."""
        from allelix.config import AllelixConfig, load_config, save_config

        save_config(tmp_path, AllelixConfig(license_overrides={"cadd": True}))
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["config", "set", "--data-dir", str(tmp_path), "license.cadd", "false"],
        )
        assert result.exit_code == 0
        cfg = load_config(tmp_path)
        assert not cfg.license_held("cadd")
        assert "cadd" not in cfg.license_overrides


class TestExportPlinkCommand:
    def test_multi_allelic_prefers_forward_over_complement(self, tmp_path, monkeypatch):
        """Guard the two-pass coord selection loop in export_plink_cmd.

        At a multi-allelic site ref=G alts=[A,T], a donor carrying T/T
        must be exported against (G,T) — not (G,A) via complement.
        Reverting the CLI loop to single-pass would make this fail.
        """
        fixture = tmp_path / "donor.txt"
        fixture.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tTEST001\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tT\tT\n"
        )

        class _FakeGnomad(GnomadAnnotator):
            def is_ready(self):
                return True

            def bulk_resolve_coordinates(self, rsids):
                return {"rs1": [("1", 100, "G", "A"), ("1", 100, "G", "T")]}

        monkeypatch.setattr(
            "allelix.annotators.gnomad.GnomadAnnotator",
            _FakeGnomad,
        )

        runner = CliRunner()
        prefix = tmp_path / "out"
        result = runner.invoke(
            main,
            # gnomAD codes GRCh38 exports only.
            ["export", "plink", str(fixture), "-o", str(prefix), "--build", "grch38"],
        )
        assert result.exit_code == 0, result.output

        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[5] == "T", f"A2 should be T (forward match), got {parts[5]}"
        assert parts[4] == "G"

        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x03

    def test_multi_allelic_with_no_fitting_record_is_missing(self, tmp_path, monkeypatch):
        """C/C at a site gnomAD lists as G>A and G>T fits neither: missing, counted, said."""
        fixture = tmp_path / "donor.txt"
        fixture.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tTEST001\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tC\tC\n"
            "rs2\t1\t200\tG\tT\n"
        )

        class _FakeGnomad(GnomadAnnotator):
            def is_ready(self):
                return True

            def bulk_resolve_coordinates(self, rsids):
                return {
                    "rs1": [("1", 100, "G", "A"), ("1", 100, "G", "T")],
                    "rs2": [("1", 200, "G", "T")],
                }

        monkeypatch.setattr("allelix.annotators.gnomad.GnomadAnnotator", _FakeGnomad)

        result = CliRunner().invoke(
            main,
            ["export", "plink", str(fixture), "-o", str(tmp_path / "out"), "--build", "grch38"],
        )
        assert result.exit_code == 0, result.output
        assert "1 call(s) fit neither gnomAD allele" in result.output
        bim = [line.split("\t") for line in (tmp_path / "out.bim").read_text().splitlines()]
        assert [row[4:6] for row in bim] == [["G", "A"], ["G", "T"]]
        assert (tmp_path / "out.bed").read_bytes()[3:] == bytes([0x01, 0x02])

    def test_split_chromosome_sorted(self, tmp_path, monkeypatch):
        """CLI sorts variants so .bim has contiguous chromosome blocks.

        MHG files can have straggler variants after the main chromosome
        run (e.g. supplemental panel appended after chrY). PLINK1.9
        rejects split chromosomes. Guards the sort in cli.py.
        """
        fixture = tmp_path / "split.txt"
        fixture.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tTEST001\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tA\tA\n"
            "rs2\t2\t200\tG\tG\n"
            "rs3\t1\t300\tT\tT\n"
        )

        class _NoGnomad(GnomadAnnotator):
            def is_ready(self):
                return False

        monkeypatch.setattr("allelix.annotators.gnomad.GnomadAnnotator", _NoGnomad)

        runner = CliRunner()
        prefix = tmp_path / "out"
        result = runner.invoke(
            main,
            ["export", "plink", str(fixture), "-o", str(prefix), "--build", "grch37"],
        )
        assert result.exit_code == 0, result.output

        chroms = [
            line.split("\t")[0] for line in (tmp_path / "out.bim").read_text().strip().split("\n")
        ]
        assert chroms == ["1", "1", "2"], f"Expected contiguous chroms, got {chroms}"

    @pytest.mark.parametrize("build", ["grch38", "grch37"])
    def test_real_gnomad_cache_codes_alleles(self, tmp_path, build):
        """The real annotators the CLI builds run the export end to end.

        No test double: a synthetic gnomAD cache in the data directory gives
        the REF/ALT coding of a GRCh38 export, so the export path's resource
        handling is the production one. gnomAD is GRCh38 data: a GRCh37
        export with no rsID identity table uses its pair only where it is not
        palindromic and the call fits it as written, as both calls here do.
        """
        data = tmp_path / "data"
        data.mkdir()
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
            conn.executescript(GNOMAD_SCHEMA)
            conn.executemany(
                "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
                "VALUES (?,?,?,?,?,0.1)",
                [("1", 100, "G", "T", "rs1"), ("1", 200, "G", "A", "rs2")],
            )
            conn.execute(
                "INSERT INTO database_versions "
                "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
                "VALUES ('gnomad','test://','synthetic','2001-01-01',2,?)",
                (f"sv:{GNOMAD_SCHEMA_VERSION}",),
            )
            conn.commit()
        fixture = tmp_path / "donor.txt"
        fixture.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tTEST001\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tT\tT\n"
            "rs2\t1\t200\tG\tA\n"
        )

        result = CliRunner().invoke(
            main,
            [
                "export",
                "plink",
                str(fixture),
                "-o",
                str(tmp_path / "out"),
                "--build",
                build,
                "--data-dir",
                str(data),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "resolution failed" not in result.output
        bim = [line.split("\t") for line in (tmp_path / "out.bim").read_text().splitlines()]
        bed = (tmp_path / "out.bed").read_bytes()[3:]
        assert [row[4:6] for row in bim] == [["G", "T"], ["G", "A"]]
        assert bed == bytes([0x03, 0x02])
        if build == "grch37":
            assert "0 site(s) from the rsID identity table, 2 from gnomAD, 0 with fallback" in (
                " ".join(result.output.split())
            )


class TestBuildDiagnosticsMessages:
    """What the CLI prints about the build decision (ADR-0044)."""

    @staticmethod
    def _emit(**overrides) -> str:
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers
        from allelix.reports._pipeline import BuildDiagnostics

        fields = {
            "effective_build": "GRCh37",
            "source": "seeds",
            "file_label": None,
            "detected_build": "GRCh37",
            "matched_count": 10,
            "inspected_count": 10,
            "position_confident": True,
        }
        fields.update(overrides)

        class _R:
            build_diagnostics = BuildDiagnostics(**fields)

        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf, force_terminal=False, width=400)
        try:
            _helpers._emit_build_diagnostics(_R())
        finally:
            _helpers.console = original
        return " ".join(buf.getvalue().split())

    def test_seeds(self) -> None:
        out = self._emit()
        assert "Build: GRCh37 (from known-SNP positions: 10/10 match GRCh37)" in out
        assert "mismatch" not in out

    def test_override(self) -> None:
        out = self._emit(source="override", file_label="GRCh37", effective_build="GRCh38")
        assert "Build: GRCh38 (forced with --build; the file says GRCh37)" in out
        assert "mismatch" not in out

    def test_contig_lengths(self) -> None:
        out = self._emit(source="contig_lengths", lengths_detail="2 contig length(s) match GRCh37")
        assert "Build: GRCh37 (from the contig lengths: 2 contig length(s) match GRCh37)" in out

    def test_assembly_header_and_file_label(self) -> None:
        for source, wording in (
            ("assembly_header", "from the VCF assembly header"),
            ("file_label", "from the file's build line"),
        ):
            out = self._emit(source=source, file_label="GRCh37")
            assert wording in out
            # The label decided and nothing contradicted it: no notice.
            assert "mismatch" not in out
            assert "Note:" not in out

    def test_below_threshold_seeds_under_a_label_are_not_called_the_source(self) -> None:
        """3/4 agreeing with the label: the banner credits the label, with the counts."""
        out = self._emit(
            source="file_label",
            file_label="GRCh37",
            matched_count=3,
            inspected_count=4,
            position_confident=False,
        )
        assert (
            "Build: GRCh37 (from the file's build line; 3/4 known-SNP positions match GRCh37, "
            "tentative)" in out
        )
        assert "from known-SNP positions" not in out

    def test_tie_including_the_build_keeps_the_inconclusive_warning(self) -> None:
        """A tie that includes the deciding build: proceed on it, and warn."""
        out = self._emit(
            source="file_label",
            file_label="GRCh37",
            detected_build=None,
            matched_count=0,
            inspected_count=4,
            position_confident=False,
            tied_builds=("GRCh37", "GRCh38"),
            tied_votes=2,
        )
        assert "Build: GRCh37 (from the file's build line" in out
        assert (
            "Build evidence inconclusive: the known-SNP positions are tied between GRCh37 "
            "and GRCh38 (2 each of 4). Using GRCh37, which the file's build line states; it "
            "is one of the tied builds." in out
        )
        assert "the other build" not in out
        assert "--build grch37 or --build grch38" in out

    def test_grch36_tie_mismatch_notice(self) -> None:
        out = self._emit(
            effective_build="GRCh36",
            source="seeds",
            file_label="GRCh37",
            detected_build=None,
            matched_count=0,
            inspected_count=4,
            position_confident=False,
            label_overridden=True,
            tentative_grch36=True,
            tied_builds=("GRCh36", "GRCh38"),
            tied_votes=2,
        )
        assert (
            "Build mismatch: the file says GRCh37, but its known-SNP positions are tied "
            "between GRCh36 and GRCh38 (2 each of 4). Using GRCh36 as the safe choice" in out
        )
        assert "say None" not in out
        assert "Build evidence inconclusive" in out
        assert "GRCh36 (hg18) was chosen as the safe option" in out
        assert "detected" not in out

    def test_override_with_a_tie_does_not_warn(self) -> None:
        out = self._emit(source="override", tied_builds=("GRCh37", "GRCh38"), tied_votes=2)
        assert "inconclusive" not in out

    def test_default_prints_a_visible_note(self) -> None:
        out = self._emit(source="default", detected_build=None, matched_count=0, inspected_count=0)
        assert "this format's default; the file states no build" in out
        assert "Note: the file states no genome build" in out
        assert "--build grch37 or --build grch38" in out

    def test_confident_mismatch_notice(self) -> None:
        out = self._emit(
            effective_build="GRCh38",
            file_label="GRCh37",
            detected_build="GRCh38",
            matched_count=9,
            label_overridden=True,
        )
        assert (
            "Build mismatch: the file says GRCh37, but its known-SNP positions say GRCh38" in out
        )
        assert "Using GRCh38" in out

    def test_tentative_grch36_mismatch_notice(self) -> None:
        out = self._emit(
            effective_build="GRCh36",
            file_label="GRCh37",
            detected_build="GRCh36",
            matched_count=3,
            inspected_count=4,
            position_confident=False,
            label_overridden=True,
            tentative_grch36=True,
        )
        assert "positions tentatively say GRCh36 (3/4)" in out
        assert (
            "GRCh36 (hg18) was chosen as the safe option: the known-SNP positions lean to it (3/4)"
            in out
        )
        assert "detected" not in out

    def test_hg19_chrm_note(self) -> None:
        out = self._emit(source="contig_lengths", lengths_detail="x", mt_skipped=12)
        assert "chrM is declared with the hg19 length (16,571 bp)" in out
        assert "12 mitochondrial record(s) were not annotated" in out


class TestRuntimeNudges:
    """Post-pipeline diagnostic emission for
    strand-aware-inactive (array input + no ref) and GRCh38 rsID-less
    undercount cases. Both are once-per-run and side-effect-only."""

    def _capture(self, result, missing: tuple[str, ...] = ()) -> str:
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers

        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf, force_terminal=False, width=200)
        try:
            _helpers._emit_runtime_nudges(result, missing)
        finally:
            _helpers.console = original
        return buf.getvalue()

    def _result(self, **overrides):
        """Synthetic minimal AnalysisResult-shaped object. The CLI
        path only reads getattr() values; no need for the real dataclass."""
        from allelix.reports._pipeline import BuildDiagnostics

        defaults = dict(
            parser_name="23andme",
            total_variants=1000,
            no_ref_variant_count=0,
            rsidless_variant_count=0,
            build_diagnostics=BuildDiagnostics(
                file_label=None,
                detected_build="GRCh37",
                effective_build="GRCh37",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        defaults.update(overrides)

        class _R:
            pass

        r = _R()
        for k, v in defaults.items():
            setattr(r, k, v)
        return r

    # ---- strand-aware inactive ----

    def test_grch37_array_points_at_the_identity_table(self):
        """GRCh37 arrays gain reference context from the rsID identity table."""
        result = self._result(no_ref_variant_count=1000)
        out = self._capture(result)
        assert "strand-aware matching inactive" in out
        assert "run `allelix db update`" in out
        assert "rsID identity table" in out

    def test_grch36_array_does_not_recommend_a_download(self):
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            no_ref_variant_count=1000,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh36",
                detected_build="GRCh36",
                effective_build="GRCh36",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        out = self._capture(result)
        assert "no reference context is available for GRCh36 input" in out
        assert "db update" not in out

    def test_grch38_array_without_ref_keeps_actionable_download_hint(self):
        """A GRCh38 array can still gain reference context from gnomAD."""
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            no_ref_variant_count=1000,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh38",
                detected_build="GRCh38",
                effective_build="GRCh38",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        out = self._capture(result)
        assert "strand-aware matching inactive" in out
        assert "run `allelix db update`" in out

    def test_array_input_with_ref_populated_quiet(self):
        """Same input but gnomAD ref-resolution succeeded → no notice."""
        result = self._result(no_ref_variant_count=0)
        out = self._capture(result)
        assert "strand-aware" not in out

    def test_vcf_input_never_fires_strand_warning(self):
        """VCF inputs carry REF natively. Even if a few rows came in
        ref-less, the message is array-specific and should be silent."""
        result = self._result(parser_name="vcf", no_ref_variant_count=1000)
        out = self._capture(result)
        assert "strand-aware" not in out

    def test_partial_ref_population_does_not_fire(self):
        """Half the variants had ref populated (e.g. an array file that
        only partially overlaps gnomAD). Don't flap the warning on."""
        result = self._result(no_ref_variant_count=500, total_variants=1000)
        out = self._capture(result)
        assert "strand-aware" not in out

    # ---- GRCh38 rsID-less undercount ----

    def test_grch38_rsidless_fires_warning(self):
        """GRCh38 VCF with mostly ID=. rows → undercount warning."""
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            parser_name="vcf",
            total_variants=1000,
            rsidless_variant_count=900,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh38",
                detected_build="GRCh38",
                effective_build="GRCh38",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        out = self._capture(result, ("the rsID identity table",))
        assert "GRCh38 input without rsIDs" in out
        assert "the rsID identity table not ready" in out
        assert "Run `allelix db update`" in out
        # The measured coverage is equal once both are installed; no claim
        # that GRCh37 surfaces more.
        assert "GRCh37" not in out

    def test_grch38_rsidless_with_resolvers_ready_is_quiet(self):
        """With the identity table and gnomAD installed, coverage matches
        GRCh37, so the notice is not shown."""
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            parser_name="vcf",
            total_variants=1000,
            rsidless_variant_count=900,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh38",
                detected_build="GRCh38",
                effective_build="GRCh38",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        assert "GRCh38 input without rsIDs" not in self._capture(result, ())

    def test_grch38_with_rsids_no_warning(self):
        """GRCh38 VCF where most rows carry rsIDs (rare but real for
        GRCh38 files passed through a dbSNP annotator) → silent."""
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            parser_name="vcf",
            total_variants=1000,
            rsidless_variant_count=100,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh38",
                detected_build="GRCh38",
                effective_build="GRCh38",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        out = self._capture(result)
        assert "GRCh38 input without rsIDs" not in out

    def test_grch37_rsidless_no_warning(self):
        """The undercount story is GRCh38-specific. A GRCh37 input —
        even one without rsIDs — gets the equivalent coverage as
        always; this warning is interim, so don't fire."""
        result = self._result(
            parser_name="vcf",
            total_variants=1000,
            rsidless_variant_count=900,
        )
        out = self._capture(result)
        assert "GRCh38" not in out

    def test_array_input_no_grch38_warning(self):
        """Array inputs are always rsID-bearing — the GRCh38 case
        doesn't apply, even if the build happened to be GRCh38."""
        from allelix.reports._pipeline import BuildDiagnostics

        result = self._result(
            parser_name="23andme",
            total_variants=1000,
            rsidless_variant_count=0,
            build_diagnostics=BuildDiagnostics(
                file_label="GRCh38",
                detected_build="GRCh38",
                effective_build="GRCh38",
                source="seeds",
                matched_count=10,
                inspected_count=10,
                position_confident=True,
            ),
        )
        out = self._capture(result)
        assert "GRCh38 input without rsIDs" not in out

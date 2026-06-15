# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""End-to-end tests for the CLI."""

from __future__ import annotations

import contextlib
import logging
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main

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
        """C-3: pin specific numeric counts so a `+= 0` mutation fails the suite."""
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(main, ["stats", str(mock_mhg_path)])
        assert result.exit_code == 0, result.output
        # Generated mock fixture (Round 23 / ADR-0021: rs104894490 NIPA1 added):
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

    def test_stats_unrecognized_file_errors(self, tmp_path):
        f = tmp_path / "garbage.txt"
        f.write_text("hello world\n", encoding="utf-8")
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(f)])
        assert result.exit_code != 0
        assert "No parser recognized" in result.output

    def test_stats_missing_file_errors(self, tmp_path):
        runner = CliRunner()
        result = runner.invoke(main, ["stats", str(tmp_path / "missing.txt")])
        assert result.exit_code != 0


class TestSkippedLineReporting:
    """Pin C3: malformed lines must surface to the user, not vanish silently."""

    def _dirty_file(self, tmp_path):
        f = tmp_path / "dirty.txt"
        f.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1\t1\t100\tA\tG\n"
            "junk\n"
            "rs2\t2\tBAD\tA\tG\n",
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
        """m12/m13: do not leak logger level or propagate setting across calls."""
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
        assert "clinvar_pathogenic" in result.output

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
        """M-5: default data-dir resolution via $ALLELIX_DATA_DIR works in the CLI."""
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
        assert "Build mismatch" in result.output
        assert "header claims GRCh37" in result.output
        assert "position data is GRCh38" in result.output
        assert "Using GRCh38" in result.output

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

    def test_analyze_build_override_skips_detection(
        self, mock_mhg_mislabeled_path, clinvar_data_dir: Path
    ):
        """`--build grch37` forces the build and silences the mismatch
        warning, even on the deliberately-mislabeled fixture.
        """
        runner = CliRunner(env={"COLUMNS": "200"})
        result = runner.invoke(
            main,
            [
                "analyze",
                str(mock_mhg_mislabeled_path),
                "--data-dir",
                str(clinvar_data_dir),
                "--build",
                "grch37",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Build mismatch" not in result.output
        assert "Build: GRCh37 (override" in result.output

    def test_analyze_unknown_format_errors(self, mock_mhg_path, clinvar_data_dir: Path):
        """M-5: cover the ParserNotFoundError branch in analyze."""
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
        """M-5: malformed-line warnings surface in the analyze command too."""
        f = tmp_path / "dirty.txt"
        f.write_text(
            "# MyHappyGenes [TEMPUS]\n"
            "# Sample ID\tMHG_X\n"
            "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
            "rs1801133\t1\t11796321\tG\tA\n"
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
        # ClinVar pathogenic
        assert "clinvar_pathogenic" in result.output
        # ClinPGx pharmacogenomic (rs1801133 + AG triggers ClinPGx LoE 2A)
        assert "pharmgkb_loe_2a" in result.output
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
        assert "20260101" in result.output  # m-6: version from VCF ##fileDate
        # ADR-0021: clinvar_data_dir populates both build caches.
        # 13 records per build x 2 builds = 26 composite total.
        assert "26" in result.output

    def test_db_update_with_file_url(
        self,
        tmp_path: Path,
        mock_clinvar_vcf: Path,
        mock_pharmgkb_dir: Path,
        mock_cpic_lookup: dict[tuple[str, str], str],
        mock_gwas_tsv: Path,
        mock_gnomad_gz: Path,
        monkeypatch,
    ):
        """M-3: download+ingest path runs end-to-end against file:// URLs.

        All annotators must succeed -- `db update` iterates the whole registry,
        so each needs a local file:// URL or mock fixture. ClinVar, ClinPGx,
        and GWAS use mock archives; gnomAD uses a mock gzipped SQLite cache.
        """
        import zipfile
        from urllib.request import pathname2url

        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases import manager as manager_module
        from allelix.databases.gwas_loader import GWAS_DB_FILENAME
        from allelix.databases.manager import get_database_info

        # ADR-0021: per-build URL map. Point both builds at the same local
        # mock VCF so db update can populate either cache from one fixture.
        clinvar_url = f"file:{pathname2url(str(mock_clinvar_vcf.resolve()))}"
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_URL_BY_BUILD",
            {"GRCh37": clinvar_url, "GRCh38": clinvar_url},
        )

        pharmgkb_zip = tmp_path / "fixture_clinical_annotations.zip"
        with zipfile.ZipFile(pharmgkb_zip, "w") as zf:
            for f in mock_pharmgkb_dir.iterdir():
                zf.write(f, arcname=f.name)
        pharmgkb_url = f"file:{pathname2url(str(pharmgkb_zip.resolve()))}"
        monkeypatch.setattr(pharmgkb_module, "PHARMGKB_CLINICAL_URL", pharmgkb_url)

        # C-1: never hit the real CPIC API from the test suite. The annotator
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

        monkeypatch.setattr(
            ClinVarAnnotator,
            "_fetch_remote_signal_for",
            staticmethod(lambda _build: "md5:test_signal"),
        )
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
            assert info["record_count"] == 13
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
        """M-1+M-2: a failed update must not destroy a working cache.

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

        # Connection-refused URL on the GRCh37 build only.
        monkeypatch.setattr(
            manager_module,
            "CLINVAR_URL_BY_BUILD",
            {"GRCh37": "http://127.0.0.1:1/missing.vcf.gz"},
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
        assert result.exit_code == 0
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
        self, tmp_path: Path, mock_clinvar_vcf: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """Remote signal matches cached signal → skip without --force."""
        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases.manager import load_clinvar_vcf

        # Pre-populate the GRCh37 ClinVar cache with a known signal.
        load_clinvar_vcf(
            mock_clinvar_vcf,
            tmp_path / clinvar_db_filename("GRCh37"),
            source_url="test",
            remote_signal="md5:cached_value",
            record_name=clinvar_record_name("GRCh37"),
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
        self, tmp_path: Path, mock_clinvar_vcf: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """Remote signal differs from cached → trigger refresh without --force."""
        from allelix.annotators import clinvar as clinvar_module
        from allelix.annotators import gwas as gwas_module
        from allelix.annotators import pharmgkb as pharmgkb_module
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases.manager import load_clinvar_vcf

        load_clinvar_vcf(
            mock_clinvar_vcf,
            tmp_path / clinvar_db_filename("GRCh37"),
            source_url="test",
            remote_signal="md5:OLD",
            record_name=clinvar_record_name("GRCh37"),
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
        monkeypatch.setattr(AlphaMissenseAnnotator, "setup", lambda self: None)

        from allelix.annotators.snpedia import SNPediaAnnotator

        monkeypatch.setattr(SNPediaAnnotator, "setup", lambda self: None)

        from allelix.annotators.cadd import CaddAnnotator

        monkeypatch.setattr(CaddAnnotator, "setup", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(main, ["db", "update", "--data-dir", str(clinvar_data_dir)])
        assert result.exit_code == 0, result.output
        assert "can't be verified" in result.output

    def test_db_update_legacy_cache_stamps_signal(
        self, tmp_path: Path, mock_clinvar_vcf: Path, mock_gnomad_gz: Path, monkeypatch
    ):
        """Legacy caches (no stored signal) get signal stamped without re-download."""
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
        assert "stamped remote signal" in result.output
        assert not called  # setup() should NOT run — just stamp the signal

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
        """W-2: failure path prints a friendly error — not propagate the raw error.

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

        monkeypatch.setattr(
            manager_module,
            "CLINVAR_URL_BY_BUILD",
            {"GRCh37": "http://127.0.0.1:1/missing.vcf.gz"},
        )
        monkeypatch.setattr(pharmgkb_module.PharmGKBAnnotator, "requires_download", False)
        monkeypatch.setattr(gwas_module.GWASCatalogAnnotator, "requires_download", False)

        gnomad_gz_url = f"file:{pathname2url(str(mock_gnomad_gz.resolve()))}"
        monkeypatch.setattr(gnomad_module, "GNOMAD_CACHE_URL", gnomad_gz_url)
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
        assert result.exit_code == 0
        assert "clinvar:" in result.output
        assert not isinstance(result.exception, urllib.error.URLError)


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
        assert "Homocysteine levels" in result.output

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
        assert "ClinPGx" in result.output
        assert "pharmgkb_loe_" in result.output
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
        is covered by a direct unit test in tests/test_pipeline_filter.py;
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

        class _FakeGnomad:
            def __init__(self, *a, **kw):
                pass

            def is_ready(self):
                return True

            def bulk_resolve_coordinates(self, rsids):
                return {"rs1": [("1", 100, "G", "A"), ("1", 100, "G", "T")]}

            def close(self):
                pass

        monkeypatch.setattr(
            "allelix.annotators.gnomad.GnomadAnnotator",
            _FakeGnomad,
        )

        runner = CliRunner()
        prefix = tmp_path / "out"
        result = runner.invoke(
            main,
            ["export", "plink", str(fixture), "-o", str(prefix), "--build", "grch37"],
        )
        assert result.exit_code == 0, result.output

        bim = (tmp_path / "out.bim").read_text().strip()
        parts = bim.split("\t")
        assert parts[5] == "T", f"A2 should be T (forward match), got {parts[5]}"
        assert parts[4] == "G"

        bed = (tmp_path / "out.bed").read_bytes()
        assert bed[3] == 0x03

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

        monkeypatch.setattr(
            "allelix.annotators.gnomad.GnomadAnnotator",
            type(
                "_NoGnomad",
                (),
                {
                    "__init__": lambda self, *a, **kw: None,
                    "is_ready": lambda self: False,
                    "close": lambda self: None,
                },
            ),
        )

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


class TestBuildDiagnosticsFallback:
    """The CLI loudly recommends --build when auto-detect has no signal.

    Real VCFs from variant callers (GATK / DeepVariant) emit ID=. and skip
    the ##contig assembly tag, so both auto-detect paths come up empty.
    Defaulting silently makes annotation positions diverge from reality —
    the user must be told to pass --build explicitly.
    """

    def test_no_signal_fires_yellow_warning(self) -> None:
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers
        from allelix.reports._pipeline import BuildDiagnostics

        diag = BuildDiagnostics(
            header_build=None,
            detected_build=None,
            effective_build="GRCh37",
            override=False,
            matched_count=0,
            inspected_count=0,
        )

        class _R:
            build_diagnostics = diag

        # Swap the module-level console for a recording one so we can read output.
        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf, force_terminal=False, width=200)
        try:
            _helpers._emit_build_diagnostics(_R())
        finally:
            _helpers.console = original

        output = buf.getvalue()
        assert "Could not auto-detect" in output
        assert "--build grch37" in output and "--build grch38" in output

    def test_override_does_not_warn(self) -> None:
        """When the user explicitly passed --build, no warning should fire."""
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers
        from allelix.reports._pipeline import BuildDiagnostics

        diag = BuildDiagnostics(
            header_build=None,
            detected_build=None,
            effective_build="GRCh38",
            override=True,
            matched_count=0,
            inspected_count=0,
        )

        class _R:
            build_diagnostics = diag

        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf, force_terminal=False, width=200)
        try:
            _helpers._emit_build_diagnostics(_R())
        finally:
            _helpers.console = original

        assert "Could not auto-detect" not in buf.getvalue()

    def test_header_only_does_not_warn(self) -> None:
        """When the header has an assembly tag, the dim status line is enough."""
        from io import StringIO

        from rich.console import Console

        from allelix.cli import _helpers
        from allelix.reports._pipeline import BuildDiagnostics

        diag = BuildDiagnostics(
            header_build="GRCh38",
            detected_build=None,
            effective_build="GRCh38",
            override=False,
            matched_count=0,
            inspected_count=0,
        )

        class _R:
            build_diagnostics = diag

        buf = StringIO()
        original = _helpers.console
        _helpers.console = Console(file=buf, force_terminal=False, width=200)
        try:
            _helpers._emit_build_diagnostics(_R())
        finally:
            _helpers.console = original

        assert "Could not auto-detect" not in buf.getvalue()

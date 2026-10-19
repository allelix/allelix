# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Plausibility flag: zygosity x gnomAD allele frequency (ADR-0041).

Covers the three buckets (plausible, suspicious, implausible), the
not-assessed cases (no gnomAD value, Y / MT, a homozygous X call), the
config-tunable thresholds, the pipeline end to end on GRCh38 and on
GRCh37 through the identity table, and the HTML / JSON / terminal
surfaces. All genotypes and frequencies are synthetic.
"""

from __future__ import annotations

import contextlib
import io
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.config import AllelixConfig, load_config, save_config
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION, GNOMAD_SCHEMA_VERSION
from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA
from allelix.models import Annotation
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import run_analysis
from allelix.reports.html import render_html
from allelix.reports.json_report import SCHEMA_VERSION, render_json
from allelix.reports.terminal import render_terminal
from allelix.utils.plausibility import (
    CHROMOSOME,
    GNOMAD_NOT_RUN,
    IMPLAUSIBLE,
    NO_FREQUENCY,
    NOT_ASSESSED,
    PLAUSIBLE,
    SUSPICIOUS,
    PlausibilityThresholds,
    classify,
    most_severe,
)

if TYPE_CHECKING:
    from pathlib import Path


class TestClassify:
    @pytest.mark.parametrize(
        ("af", "copies", "chrom", "expected"),
        [
            # The three buckets.
            (0.00001, 2, "3", (IMPLAUSIBLE, "")),
            (0.0, 2, "3", (IMPLAUSIBLE, "")),
            (0.0, 1, "3", (SUSPICIOUS, "")),
            (0.2, 2, "3", (PLAUSIBLE, "")),
            (0.00001, 1, "3", (PLAUSIBLE, "")),
            # Boundaries: homozygous is strict "<", heterozygous is "<=".
            (0.0001, 2, "3", (PLAUSIBLE, "")),
            (0.0, 1, "22", (SUSPICIOUS, "")),
            # Not assessed, with the reason (ADR-0050).
            (None, 2, "3", (NOT_ASSESSED, NO_FREQUENCY)),
            (None, 1, "3", (NOT_ASSESSED, NO_FREQUENCY)),
            (0.00001, 0, "3", (NOT_ASSESSED, CHROMOSOME)),
            # X: a homozygous call may be a male's hemizygous call.
            (0.00001, 2, "X", (NOT_ASSESSED, CHROMOSOME)),
            (0.0, 1, "X", (SUSPICIOUS, "")),
            # Y and MT are haploid.
            (0.0, 2, "Y", (NOT_ASSESSED, CHROMOSOME)),
            (0.0, 1, "Y", (NOT_ASSESSED, CHROMOSOME)),
            (0.0, 2, "MT", (NOT_ASSESSED, CHROMOSOME)),
            # No frequency on a chromosome with no rule: the chromosome decides.
            (None, 2, "Y", (NOT_ASSESSED, CHROMOSOME)),
        ],
    )
    def test_buckets(
        self, af: float | None, copies: int, chrom: str, expected: tuple[str, str]
    ) -> None:
        assert classify(af, copies, chrom) == expected

    def test_thresholds_are_tunable(self) -> None:
        strict = PlausibilityThresholds(homozygous_max_af=0.000001, heterozygous_max_af=0.001)
        assert classify(0.00001, 2, "3", strict) == (PLAUSIBLE, "")
        assert classify(0.0005, 1, "3", strict) == (SUSPICIOUS, "")

    def test_most_severe(self) -> None:
        assert most_severe([PLAUSIBLE, IMPLAUSIBLE, SUSPICIOUS]) == IMPLAUSIBLE
        assert most_severe([NOT_ASSESSED, PLAUSIBLE]) == PLAUSIBLE
        assert most_severe([NOT_ASSESSED, NOT_ASSESSED]) == NOT_ASSESSED
        assert most_severe([]) == ""


# ---------------------------------------------------------------------------
# Pipeline, end to end. Synthetic ClinVar rows (pathogenic, so none are
# suppressed), a synthetic gnomAD cache keyed by GRCh38 position, and a
# synthetic chip. rsIDs are in an unused 9005xxxxx range.
# ---------------------------------------------------------------------------

_SINGLE = "criteria_provided,_single_submitter"

# (rsid, chrom, pos38, ref, alt, user alleles, gnomAD AF or None, expected)
_SITES: tuple[tuple[str, str, int, str, str, str, float | None, str], ...] = (
    ("rs900580001", "3", 3_800_100, "G", "A", "AA", 0.00001, IMPLAUSIBLE),
    ("rs900580002", "3", 3_800_200, "C", "T", "CT", 0.0, SUSPICIOUS),
    ("rs900580003", "3", 3_800_300, "T", "C", "CC", 0.2, PLAUSIBLE),
    ("rs900580004", "3", 3_800_400, "A", "G", "AG", 0.00002, PLAUSIBLE),
    ("rs900580005", "3", 3_800_500, "C", "G", "GG", None, NOT_ASSESSED),
    ("rs900580006", "X", 3_800_600, "G", "A", "AA", 0.00001, NOT_ASSESSED),
    ("rs900580007", "X", 3_800_700, "C", "T", "CT", 0.0, SUSPICIOUS),
)
_POS37_OFFSET = 100_000  # synthetic GRCh37 position = GRCh38 position - offset


def _write_clinvar(data: Path, build: str, offset: int) -> None:
    rows = [
        (rsid, chrom, pos - offset, ref, alt, "Pathogenic", f"Cond_{rsid}", "GENE", _SINGLE, i)
        for i, (rsid, chrom, pos, ref, alt, *_rest) in enumerate(_SITES, 1)
    ]
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://plausibility', '1', '2026-10-10', ?, NULL, ?)",
            (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _write_gnomad(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / GNOMAD_DB_FILENAME)) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            [
                (chrom, pos, ref, alt, rsid, af)
                for rsid, chrom, pos, ref, alt, _gt, af, _exp in _SITES
                if af is not None
            ],
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) "
            "VALUES ('gnomad', 'test://plausibility', '4.1', '2026-10-10', ?, ?)",
            (len(_SITES), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()


def _write_chip(path: Path, offset: int) -> Path:
    lines = [
        "# MyHappyGenes [TEMPUS]",
        "# Sample ID\tPLAUS000001",
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward",
    ]
    lines += [
        f"{rsid}\t{chrom}\t{pos - offset}\t{gt[0]}\t{gt[1]}"
        for rsid, chrom, pos, _ref, _alt, gt, _af, _exp in _SITES
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def grch38_data(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    _write_clinvar(data, "GRCh38", 0)
    _write_gnomad(data)
    return data


def _plausibility_by_rsid(result) -> dict[str, str]:
    return {a.rsid: a.plausibility for a in result.annotations}


class TestPipeline:
    def test_grch38_buckets(
        self, tmp_path: Path, grch38_data: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        result = run_analysis(
            _write_chip(tmp_path / "chip.txt", 0),
            MyHappyGenesParser(),
            [cm_stack.enter_context(ClinVarAnnotator(grch38_data, builds=("GRCh38",)))],
            build_override="GRCh38",
            gnomad=cm_stack.enter_context(GnomadAnnotator(grch38_data)),
        )
        assert _plausibility_by_rsid(result) == {s[0]: s[7] for s in _SITES}

    def test_flagged_rows_are_not_suppressed(
        self, tmp_path: Path, grch38_data: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        result = run_analysis(
            _write_chip(tmp_path / "chip.txt", 0),
            MyHappyGenesParser(),
            [cm_stack.enter_context(ClinVarAnnotator(grch38_data, builds=("GRCh38",)))],
            build_override="GRCh38",
            gnomad=cm_stack.enter_context(GnomadAnnotator(grch38_data)),
        )
        assert len(result.annotations) == len(_SITES)

    def test_thresholds_reach_the_pipeline(
        self, tmp_path: Path, grch38_data: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        result = run_analysis(
            _write_chip(tmp_path / "chip.txt", 0),
            MyHappyGenesParser(),
            [cm_stack.enter_context(ClinVarAnnotator(grch38_data, builds=("GRCh38",)))],
            build_override="GRCh38",
            gnomad=cm_stack.enter_context(GnomadAnnotator(grch38_data)),
            plausibility=PlausibilityThresholds(
                homozygous_max_af=0.000001, heterozygous_max_af=0.00005
            ),
        )
        got = _plausibility_by_rsid(result)
        assert got["rs900580001"] == PLAUSIBLE  # 1e-5 is no longer below the cut-off
        assert got["rs900580004"] == SUSPICIOUS  # 2e-5 is now at or below it

    def test_without_gnomad_nothing_is_assessed(
        self, tmp_path: Path, grch38_data: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        result = run_analysis(
            _write_chip(tmp_path / "chip.txt", 0),
            MyHappyGenesParser(),
            [cm_stack.enter_context(ClinVarAnnotator(grch38_data, builds=("GRCh38",)))],
            build_override="GRCh38",
        )
        assert set(_plausibility_by_rsid(result).values()) == {NOT_ASSESSED}
        assert {a.plausibility_reason for a in result.annotations} == {GNOMAD_NOT_RUN}

    def test_grch37_chip_through_identity_table(
        self, tmp_path: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        """The case the flag was designed for: a GRCh37 consumer chip."""
        data = tmp_path / "data"
        data.mkdir()
        _write_clinvar(data, "GRCh37", _POS37_OFFSET)
        _write_gnomad(data)
        write_identity_db(
            [
                (rsid, chrom, pos - _POS37_OFFSET, ref, alt, pos, ref, alt, 0)
                for rsid, chrom, pos, ref, alt, *_rest in _SITES
            ],
            data / RSID_IDENTITY_DB_FILENAME,
            version="dbSNP test",
            source_url="test://plausibility",
            downloaded_at="2026-10-10T00:00:00Z",
        )
        result = run_analysis(
            _write_chip(tmp_path / "chip37.txt", _POS37_OFFSET),
            MyHappyGenesParser(),
            [cm_stack.enter_context(ClinVarAnnotator(data, builds=("GRCh37",)))],
            build_override="GRCh37",
            gnomad=cm_stack.enter_context(GnomadAnnotator(data)),
            identity=cm_stack.enter_context(RsidIdentityAnnotator(data)),
        )
        assert _plausibility_by_rsid(result) == {s[0]: s[7] for s in _SITES}


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


class TestConfig:
    def test_defaults(self, tmp_path: Path) -> None:
        cfg = load_config(tmp_path)
        assert cfg.plausibility == PlausibilityThresholds(0.0001, 0.0)
        assert "[plausibility]" in (tmp_path / "config.toml").read_text()

    def test_roundtrip(self, tmp_path: Path) -> None:
        save_config(tmp_path, AllelixConfig(plausibility=PlausibilityThresholds(5e-05, 1e-06)))
        assert load_config(tmp_path).plausibility == PlausibilityThresholds(5e-05, 1e-06)

    def test_older_config_without_section_gets_defaults(self, tmp_path: Path) -> None:
        (tmp_path / "config.toml").write_text("[sources]\nclinvar = true\n")
        assert load_config(tmp_path).plausibility == PlausibilityThresholds()

    @pytest.mark.parametrize("bad", ["1.5", "-0.1", '"abc"', "true", "nan"])
    def test_invalid_value_falls_back_to_default(self, tmp_path: Path, bad: str) -> None:
        (tmp_path / "config.toml").write_text(
            f"[plausibility]\nhomozygous_max_af = {bad}\nheterozygous_max_af = 0.001\n"
        )
        cfg = load_config(tmp_path)
        assert cfg.plausibility == PlausibilityThresholds(0.0001, 0.001)

    def test_cli_set_and_get(self, tmp_path: Path) -> None:
        runner = CliRunner()
        args = ["--data-dir", str(tmp_path)]
        result = runner.invoke(
            main, ["config", "set", *args, "plausibility.homozygous_max_af", "0.00005"]
        )
        assert result.exit_code == 0, result.output
        result = runner.invoke(main, ["config", "get", *args, "plausibility.homozygous_max_af"])
        assert result.output.strip() == "5e-05"
        assert load_config(tmp_path).plausibility.homozygous_max_af == 0.00005

    @pytest.mark.parametrize(
        ("key", "value", "message"),
        [
            ("plausibility.homozygous_max_af", "2", "between 0 and 1"),
            ("plausibility.homozygous_max_af", "often", "between 0 and 1"),
            ("plausibility.max_af", "0.1", "Unknown plausibility key"),
        ],
    )
    def test_cli_rejects_bad_input(self, tmp_path: Path, key: str, value: str, message: str):
        result = CliRunner().invoke(
            main, ["config", "set", "--data-dir", str(tmp_path), key, value]
        )
        assert result.exit_code != 0
        assert message in result.output


# ---------------------------------------------------------------------------
# Report surfaces
# ---------------------------------------------------------------------------


def _ann(rsid: str, plausibility: str) -> Annotation:
    return Annotation(
        source="clinvar",
        rsid=rsid,
        significance="Pathogenic",
        category="disease",
        magnitude=8.0,
        description="",
        attribution="ClinVar",
        genotype_match="AA",
        condition=f"Cond_{rsid}",
        gene="GENE",
        alt="A",
        allele_frequency=0.00001,
        plausibility=plausibility,
    )


_FLAGGED = [
    _ann("rs900580001", IMPLAUSIBLE),
    _ann("rs900580002", SUSPICIOUS),
    _ann("rs900580003", PLAUSIBLE),
]


def _result(tmp_path: Path, annotations: list[Annotation]):
    from allelix.reports._pipeline import AnalysisResult

    return AnalysisResult(
        file_path=tmp_path / "chip.txt",
        parser_name="myhappygenes",
        parser_display_name="MyHappyGenes",
        sample_id="PLAUS000001",
        build="GRCh38",
        total_variants=len(annotations),
        skipped_count=0,
        annotators_used=[("clinvar", "1")],
        annotations=annotations,
    )


class TestReports:
    def test_json_carries_the_field_and_bumps_schema(self, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        render_json(_result(tmp_path, _FLAGGED), output_path=out)
        payload = json.loads(out.read_text())
        assert payload["schema_version"] == SCHEMA_VERSION
        assert int(SCHEMA_VERSION) >= 7
        got = {a["rsid"]: a["plausibility"] for a in payload["annotations"]}
        assert got == {
            "rs900580001": IMPLAUSIBLE,
            "rs900580002": SUSPICIOUS,
            "rs900580003": PLAUSIBLE,
        }

    def test_html_warning_icon_on_flagged_rows_only(self, tmp_path: Path) -> None:
        out = tmp_path / "report.html"
        render_html(_result(tmp_path, _FLAGGED), output_path=out)
        body = out.read_text()
        rows = {
            rsid: body.split(f'data-search-text="{rsid}', 1)[0].rsplit("<tr ", 1)[1]
            for rsid in ("rs900580001", "rs900580002", "rs900580003")
        }
        assert 'data-plausibility="implausible"' in rows["rs900580001"]
        assert 'data-plausibility="suspicious"' in rows["rs900580002"]
        assert 'data-plausibility="plausible"' in rows["rs900580003"]
        assert body.count('class="plaus plaus-implausible"') == 1
        assert body.count('class="plaus plaus-suspicious"') == 1
        assert "plaus-plausible" not in body

    def test_terminal_marks_flagged_genotypes(self, tmp_path: Path) -> None:
        buf = io.StringIO()
        console = Console(file=buf, width=200, color_system=None)
        render_terminal(_result(tmp_path, _FLAGGED), console)
        text = buf.getvalue()
        lines = {rsid: next(ln for ln in text.splitlines() if rsid in ln) for rsid in _RSIDS}
        assert "AA ⚠" in lines["rs900580001"]
        assert "AA ?" in lines["rs900580002"]
        assert "⚠" not in lines["rs900580003"] and "?" not in lines["rs900580003"]
        assert "implausible: homozygous for an allele very rare in gnomAD" in text

    def test_terminal_has_no_legend_without_flags(self, tmp_path: Path) -> None:
        buf = io.StringIO()
        console = Console(file=buf, width=200, color_system=None)
        render_terminal(_result(tmp_path, [_ann("rs900580003", PLAUSIBLE)]), console)
        assert "implausible:" not in buf.getvalue()


_RSIDS = ("rs900580001", "rs900580002", "rs900580003")

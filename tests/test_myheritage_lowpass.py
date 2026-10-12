# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""MyHeritage low-pass WGS exports: provenance and header build.

Kits processed from January 2026 are sequenced at about 2x and imputed.
The export is the usual MyHeritage CSV with a ``##`` preamble that
declares the method and the build. The fixture is that five-line header
plus synthetic rows.
"""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers import detect_parser
from allelix.parsers.myheritage import MyHeritageParser
from allelix.reports._pipeline import run_analysis

FIXTURES = Path(__file__).resolve().parent / "fixtures"
LOWPASS = FIXTURES / "myheritage_lowpass_wgs.csv"


class TestMetadata:
    def test_method_and_build(self) -> None:
        meta = MyHeritageParser().get_metadata(LOWPASS)
        assert meta["method"] == "low-pass-wgs-imputed"
        assert meta["build"] == "GRCh37"

    def test_array_export_declares_neither(self) -> None:
        meta = MyHeritageParser().get_metadata(FIXTURES / "mock_myheritage.csv")
        assert "method" not in meta
        assert meta["build"] is None  # no default here (ADR-0044)

    def test_a_build38_reference_is_a_build38_header(self, tmp_path: Path) -> None:
        path = tmp_path / "x.csv"
        path.write_text(LOWPASS.read_text().replace("##reference=build37", "##reference=build38"))
        meta = MyHeritageParser().get_metadata(path)
        assert meta["build"] == "GRCh38"

    def test_unknown_method_is_kept_verbatim(self, tmp_path: Path) -> None:
        path = tmp_path / "x.csv"
        path.write_text(
            LOWPASS.read_text().replace(
                "##method=Low-pass Whole Genome Sequencing", "##method=Something New"
            )
        )
        assert MyHeritageParser().get_metadata(path)["method"] == "Something New"


class TestParsing:
    def test_detected_and_parsed(self) -> None:
        assert detect_parser(LOWPASS).name == "myheritage"
        calls = [(v.rsid, v.allele1, v.allele2) for v in MyHeritageParser().parse(LOWPASS)]
        assert calls == [
            ("rs4477212", "A", "A"),
            ("rs3094315", "A", "G"),
            ("rs9001001", "-", "-"),
            ("rs9001002", "C", "T"),
        ]

    def test_header_build_reaches_the_build_decision(self) -> None:
        result = run_analysis(LOWPASS, MyHeritageParser(), [])
        assert result.build_diagnostics is not None
        assert result.build_diagnostics.file_label == "GRCh37"
        assert result.build_diagnostics.source == "file_label"
        assert result.input_method == "low-pass-wgs-imputed"


class TestReports:
    def _analyze(self, tmp_path: Path, name: str) -> tuple[str, Path]:
        out = tmp_path / name
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main,
            [
                "analyze",
                str(LOWPASS),
                "--data-dir",
                str(_clinvar_data_dir(tmp_path)),
                "--no-update",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        return result.output, out

    def test_json_states_the_method(self, tmp_path: Path) -> None:
        output, out = self._analyze(tmp_path, "r.json")
        assert json.loads(out.read_text())["input"]["method"] == "low-pass-wgs-imputed"
        assert "genotypes are imputed" in " ".join(output.split())

    def test_html_states_provenance_and_swaps_the_caveat(self, tmp_path: Path) -> None:
        _output, out = self._analyze(tmp_path, "r.html")
        body = out.read_text()
        assert "<strong>Input method.</strong>" in body
        assert "genotypes are imputed" in body
        assert "<h3>Imputed genotypes</h3>" in body
        assert "Pseudogene cross-hybridization" not in body

    def test_array_input_keeps_the_array_caveat(self, tmp_path: Path) -> None:
        out = tmp_path / "r.html"
        result = CliRunner().invoke(
            main,
            [
                "analyze",
                str(FIXTURES / "mock_myheritage.csv"),
                "--data-dir",
                str(_clinvar_data_dir(tmp_path)),
                "--no-update",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        body = out.read_text()
        assert "Pseudogene cross-hybridization" in body
        assert "<strong>Input method.</strong>" not in body


def _clinvar_data_dir(tmp_path: Path) -> Path:
    """A minimal ready ClinVar cache (both builds) so `analyze` runs."""
    import contextlib
    import sqlite3

    from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
    from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
    from allelix.databases.schema import CLINVAR_SCHEMA

    data = tmp_path / "data"
    if data.exists():
        return data
    data.mkdir()
    for build in ("GRCh37", "GRCh38"):
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
                "record_count, remote_signal, local_version_tag) "
                "VALUES (?, 'test://lowpass', '1', '2026-10-10', 0, NULL, ?)",
                (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    return data

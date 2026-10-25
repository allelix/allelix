"""Every registered source's part in a run, and report-level warnings.

The CLI builds ``sources`` once from the cache states (ADR-0048) and the
user's choices; JSON, HTML and the annotated VCF only read it.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from pathlib import Path

import pytest
from click.testing import CliRunner
from test_readonly_caches import _data_dir, _denied_permissions, _flat, _read_only

from allelix.cli import main
from allelix.reports._sources import (
    READY_NOT_USED,
    SKIPPED,
    USED,
    build_sources,
    clinvar_skipped_warning,
)

_FIXTURES = Path(__file__).parent / "fixtures"
_REGISTRY = [
    "clinvar",
    "pharmgkb",
    "gwas",
    "snpedia",
    "gnomad",
    "alphamissense",
    "cadd",
    "rsid_identity",
]


@pytest.fixture
def data(
    tmp_path: Path, build_synthetic_clinvar_cache, mock_pharmgkb_dir: Path, mock_gwas_tsv: Path
) -> Path:
    """Every cache type; AlphaMissense in WAL mode, CADD off by default config."""
    return _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)


def _run(
    data: Path,
    out: Path,
    *extra: str,
    input_file: str = "mock_23andme.txt",
    build: str = "grch37",
) -> tuple[int, str]:
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(_FIXTURES / input_file),
            "--build",
            build,
            "--no-update",
            "--data-dir",
            str(data),
            "--output",
            str(out),
            *extra,
        ],
    )
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return result.exit_code, result.output


def _sources(report: Path) -> dict[str, tuple[str, str | None]]:
    payload = json.loads(report.read_text())
    return {s["name"]: (s["state"], s["reason"]) for s in payload["sources"]}


def _write_config(data: Path, text: str) -> None:
    (data / "config.toml").write_text(text)


class TestSourcesBlock:
    def test_every_registered_source_has_one_record(self, data: Path, tmp_path: Path) -> None:
        out = tmp_path / "r.json"
        code, output = _run(data, out)
        assert code == 0, output
        payload = json.loads(out.read_text())
        assert [s["name"] for s in payload["sources"]] == _REGISTRY
        assert set(payload["sources"][0]) == {"name", "display_name", "state", "reason", "version"}
        sources = _sources(out)
        assert sources["clinvar"] == (USED, None)
        assert sources["alphamissense"] == (SKIPPED, "wal")
        assert sources["cadd"] == (SKIPPED, "disabled")
        # ``annotators`` keeps meaning the sources that ran.
        ran = {a["name"] for a in payload["annotators"]}
        assert ran == {name for name, (state, _) in sources.items() if state == USED}
        assert payload["report_warnings"] == []

    def test_stale_clinvar_is_named_and_warned(self, data: Path, tmp_path: Path) -> None:
        for build in ("GRCh37", "GRCh38"):
            with contextlib.closing(sqlite3.connect(data / f"clinvar.{build}.sqlite")) as conn:
                conn.execute("UPDATE database_versions SET local_version_tag = 'iv:4'")
                conn.commit()
        out = tmp_path / "r.json"
        code, output = _run(data, out)
        assert code == 0, output
        assert _sources(out)["clinvar"] == (SKIPPED, "outdated")
        warnings = json.loads(out.read_text())["report_warnings"]
        assert [w["code"] for w in warnings] == ["clinvar_skipped"]
        assert "ClinVar was not used" in warnings[0]["message"]
        assert "out of date (built for iv:4, needs iv:9)" in warnings[0]["message"]
        assert _flat("Warning: ClinVar was not used") in _flat(output)

        html = tmp_path / "r.html"
        code, output = _run(data, html)
        assert code == 0, output
        page = html.read_text()
        assert '<div class="notice-warn"><strong>Warning.</strong> ClinVar was not used' in page
        assert "<strong>Sources not used.</strong>" in page
        assert "ClinVar: out of date (built for iv:4, needs iv:9)" in page
        assert "Not run by choice: CADD (turned off in config)" in page

    @pytest.mark.parametrize(
        ("name", "reason"),
        [("gwas", "missing"), ("pharmgkb", "unreadable"), ("pharmgkb", "corrupt")],
    )
    def test_failure_reasons(self, data: Path, tmp_path: Path, name: str, reason: str) -> None:
        path = data / f"{name}.sqlite"
        if reason == "missing":
            path.unlink()
        elif reason == "corrupt":
            path.write_bytes(b"x" * 4096)
        out = tmp_path / "r.json"
        denial = (
            _denied_permissions(path, read=True)
            if reason == "unreadable"
            else contextlib.nullcontext()
        )
        with denial:
            code, output = _run(data, out)
        assert code == 0, output
        assert _sources(out)[name] == (SKIPPED, reason)

    def test_readonly_wal(self, data: Path, tmp_path: Path) -> None:
        out = tmp_path / "r.json"
        with _read_only(data):
            code, output = _run(data, out)
        assert code == 0, output
        assert _sources(out)["alphamissense"] == (SKIPPED, "readonly_wal")

    def test_user_choices(self, data: Path, tmp_path: Path) -> None:
        _write_config(data, "[license]\ncommercial = true\n")
        out = tmp_path / "r.json"
        code, output = _run(data, out, "--no-gnomad")
        assert code == 0, output
        sources = _sources(out)
        assert sources["snpedia"] == (SKIPPED, "license")
        assert sources["gnomad"] == (SKIPPED, "flag")

        _write_config(data, "")
        code, output = _run(data, out, "--exclude-snpedia")
        assert code == 0, output
        assert _sources(out)["snpedia"] == (SKIPPED, "excluded")

    def test_a_choice_wins_over_a_missing_cache(self, data: Path, tmp_path: Path) -> None:
        (data / "cadd.sqlite").unlink()
        out = tmp_path / "r.json"
        code, output = _run(data, out)
        assert code == 0, output
        assert _sources(out)["cadd"] == (SKIPPED, "disabled")
        assert "Skipping CADD" not in output

    def test_annotated_vcf_names_the_skipped_sources(self, data: Path, tmp_path: Path) -> None:
        vcf_out = tmp_path / "out.vcf"
        code, output = _run(
            data,
            tmp_path / "r.json",
            "--vcf-out",
            str(vcf_out),
            input_file="mock_vcf.vcf",
            build="grch38",
        )
        assert code == 0, output
        header = [ln for ln in vcf_out.read_text().splitlines() if ln.startswith("##ALLELIX_")]
        assert "##ALLELIX_SkippedSources=alphamissense:wal,cadd:disabled" in header
        assert "##ALLELIX_VCF_SchemaVersion=0.5.0" in header

    def test_threshold_note_names_only_sources_that_ran(self, data: Path, tmp_path: Path) -> None:
        html = tmp_path / "r.html"
        extra = ("--min-magnitude", "3", "--snpedia-min-magnitude", "1")
        code, output = _run(data, html, *extra)
        assert code == 0, output
        assert "SNPEDIA (mag &ge; 1.0)" in html.read_text()
        code, output = _run(data, html, *extra, "--exclude-snpedia")
        assert code == 0, output
        assert "SNPEDIA" not in html.read_text()


class TestClinVarOneBuild:
    """``db update --build grch37`` installs one ClinVar build; analyze uses it."""

    def test_the_installed_build_is_used(self, data: Path, tmp_path: Path) -> None:
        (data / "clinvar.GRCh38.sqlite").unlink()
        out = tmp_path / "r.json"
        code, output = _run(data, out)
        assert code == 0, output
        assert _sources(out)["clinvar"] == (USED, None)
        assert json.loads(out.read_text())["report_warnings"] == []
        assert "Skipping ClinVar" not in output

    def test_input_on_the_other_build_skips_clinvar(self, data: Path, tmp_path: Path) -> None:
        (data / "clinvar.GRCh38.sqlite").unlink()
        out = tmp_path / "r.json"
        code, output = _run(data, out, input_file="mock_vcf.vcf", build="grch38")
        assert code == 0, output
        assert _sources(out)["clinvar"] == (SKIPPED, "missing")
        payload = json.loads(out.read_text())
        assert "clinvar" not in {a["name"] for a in payload["annotators"]}
        assert [w["code"] for w in payload["report_warnings"]] == ["clinvar_skipped"]
        assert _flat("Skipping ClinVar: not installed.") in _flat(output)
        assert _flat("(input is GRCh38)") in _flat(output)


class TestBuildSources:
    def test_states(self) -> None:
        registry = [(n, n.upper()) for n in ("clinvar", "gnomad", "rsid_identity", "cadd")]
        records = build_sources(
            registry,
            choices={"cadd": "disabled"},
            failures={},
            used={"clinvar": "v1"},
            not_applicable=["rsid_identity"],
        )
        assert [(r.name, r.state, r.reason, r.version) for r in records] == [
            ("clinvar", USED, None, "v1"),
            ("gnomad", READY_NOT_USED, None, None),
            ("rsid_identity", SKIPPED, "not_applicable_build", None),
            ("cadd", SKIPPED, "disabled", None),
        ]
        assert clinvar_skipped_warning(records) is None

    def test_a_choice_wins_over_a_failure(self) -> None:
        records = build_sources(
            [("clinvar", "ClinVar")],
            choices={"clinvar": "excluded"},
            failures={"clinvar": ("missing", "ClinVar: not installed.")},
            used={},
        )
        assert (records[0].state, records[0].reason) == (SKIPPED, "excluded")
        warning = clinvar_skipped_warning(records)
        assert warning is not None
        assert warning.code == "clinvar_skipped"
        assert warning.message.endswith("ClinVar: excluded.")

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Database and input text is printed as it is, not read as Rich markup (2.6.2).

Text in square brackets that looks like a Rich tag was read as one: an opening
tag such as ``[bold]`` vanished from the line, and a closing tag that matched
nothing (``[/x]``) raised MarkupError and ended ``analyze`` with exit 1 after
the analysis had run. A ``:dna:`` code became an emoji. The high-value no-call
lines lost their bracketed reason the same way. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename
from allelix.cli import main
from tests.test_annotation_parity_and_named_disclosures import _mthfr, _vcf

_RUNNER_ENV = {"COLUMNS": "500"}


def _invoke(*args: str):
    result = CliRunner(env=_RUNNER_ENV).invoke(main, list(args))
    assert result.exception is None or isinstance(result.exception, SystemExit), repr(
        result.exception
    )
    assert result.exit_code == 0, result.output
    return result.output


def _data(tmp_path, *, trait="Synthetic condition", condition="Synthetic condition"):
    data = tmp_path / "data"
    pos = _mthfr(data, "GRCh37")
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute("UPDATE gwas_associations SET trait = ?", (trait,))
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
        conn.execute("UPDATE clinvar_variants SET condition = ?", (condition,))
        conn.commit()
    return data, pos


def _analyze(tmp_path, data, pos, *extra: str, genotype="0/1"):
    vcf = _vcf(tmp_path / "input.vcf", [("rs1801133", "1", pos, "G", "A", genotype, ".")])
    return _invoke(
        "analyze",
        str(vcf),
        "--build",
        "grch37",
        "--data-dir",
        str(data),
        "--no-update",
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        *extra,
    )


def test_a_closing_tag_in_a_gwas_trait_does_not_crash_analyze(tmp_path):
    data, pos = _data(tmp_path, trait="Synthetic [/x] trait")

    output = _analyze(tmp_path, data, pos)

    assert "Synthetic [/x] trait" in output


def test_bracketed_database_text_is_shown_not_dropped(tmp_path):
    data, pos = _data(
        tmp_path,
        trait="Trait [bold] kept",
        condition="Condition [red]kept[/red] too",
    )

    output = _analyze(tmp_path, data, pos)

    assert "Trait [bold] kept" in output
    assert "Condition [red]kept[/red] too" in output


def test_an_emoji_code_in_database_text_is_shown_as_written(tmp_path):
    data, pos = _data(tmp_path, trait="Trait :dna: code")

    output = _analyze(tmp_path, data, pos)

    assert "Trait :dna: code" in output


def test_a_high_value_no_call_keeps_its_bracketed_reason(tmp_path):
    data, pos = _data(tmp_path)

    output = _analyze(tmp_path, data, pos, genotype="./.")

    assert "high-value SNP(s) returned no-call" in output
    assert "[no genotype called]" in output


def test_a_removed_row_from_the_baseline_report_is_shown_as_written(tmp_path):
    data, pos = _data(tmp_path)
    baseline = tmp_path / "baseline.json"
    _analyze(tmp_path, data, pos, "--output", str(baseline))
    report = json.loads(baseline.read_text())
    gone = dict(report["annotations"][0])
    gone.update(rsid="rs3999442001", condition="Gone [/x] condition", gene="GENE[/y]")
    report["annotations"].append(gone)
    baseline.write_text(json.dumps(report))

    output = _analyze(tmp_path, data, pos, "--diff", str(baseline))

    assert "Removed Annotations (1)" in output
    assert "Gone [/x] condition" in output
    assert "GENE[/y]" in output


def test_db_status_shows_a_bracketed_version_as_written(tmp_path):
    data, _ = _data(tmp_path)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute("UPDATE database_versions SET version = 'release [/x] 1'")
        conn.commit()

    output = _invoke("db", "status", "--data-dir", str(data))

    assert "release [/x] 1" in output


def test_stats_shows_a_bracketed_file_name_and_sample_id_as_written(tmp_path):
    vcf = tmp_path / "genome [bold] one.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample[/x]\n"
        "1\t11856378\trs1801133\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
    )

    output = _invoke("stats", str(vcf), "--build", "grch37")

    assert "genome [bold] one.vcf" in output
    assert "Sample[/x]" in output


# Text ending in a backslash is printed with one backslash, as written: escaping
# for Rich must not double it where no tag follows (2.6.2 RC2).


def _no_doubled_backslash(output: str) -> None:
    assert "\\\\" not in output, [line for line in output.splitlines() if "\\\\" in line]


def test_db_status_and_db_update_show_a_version_ending_in_a_backslash_as_written(tmp_path):
    data, pos = _data(tmp_path)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute("UPDATE database_versions SET version = 'release 1\\'")
        conn.commit()

    status = _invoke("db", "status", "--data-dir", str(data))
    analyzed = _analyze(tmp_path, data, pos)

    assert "release 1\\" in status
    _no_doubled_backslash(status)
    assert "release 1\\)" in analyzed  # "Analyzing against: ... (release 1\)"
    _no_doubled_backslash(analyzed)


def test_stats_shows_a_file_name_and_sample_id_ending_in_a_backslash_as_written(tmp_path):
    vcf = tmp_path / "genome\\"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\\\n"
        "1\t11856378\trs1801133\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
    )

    output = _invoke("stats", str(vcf), "--format", "vcf", "--build", "grch37")

    assert "Genotype File Stats: genome\\" in output
    assert "Sample\\" in output
    _no_doubled_backslash(output)


def test_db_clean_and_config_show_a_data_dir_ending_in_a_backslash_as_written(tmp_path):
    data = tmp_path / "data\\"
    data.mkdir()
    (data / "gwas.sqlite").write_text("synthetic")

    clean = _invoke("db", "clean", "--dry-run", "--data-dir", str(data))
    config = _invoke("config", "show", "--data-dir", str(data))

    # Join wrapped lines without Rich's alignment padding before both checks.
    # This also catches a doubled backslash run split across two title lines.
    clean = "".join(line.strip() for line in clean.splitlines())
    config = "".join(line.strip() for line in config.splitlines())
    assert "data\\)" in clean
    _no_doubled_backslash(clean)
    _no_doubled_backslash(config)


@pytest.mark.parametrize(
    "text",
    ["plain", "Trait [bold] kept", "Gone [/x] condition", "a\\[b]", "", "x[/y]\\"]
    + [f"data[bold]{'\\' * n}" for n in range(1, 5)],
)
@pytest.mark.parametrize(
    "template",
    [
        "{}",
        "[dim]{}[/dim]",
        "[dim]at {}{}[/dim]",
        "[cyan]{}[/cyan] tail",
        "x {} y",
        "{}[bold]B[/bold]",
    ],
)
def test_literal_prints_as_written_before_a_tag_an_empty_suffix_text_or_nothing(text, template):
    from rich.console import Console

    from allelix.utils.markup import literal

    console = Console(width=300, color_system=None, record=True, emoji=False)
    console.print(template.format(literal(text), ""), end="")
    plain = template.replace("[dim]", "").replace("[/dim]", "").replace("[cyan]", "")
    plain = plain.replace("[/cyan]", "").replace("[bold]", "").replace("[/bold]", "")
    assert console.export_text() == plain.format(text, "")

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Long observed alleles must not hide the other answers in terminal tables."""

from __future__ import annotations

import io
import re
from dataclasses import asdict
from pathlib import Path

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.cli import main, utility
from allelix.models import Annotation, Variant
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.diff import ChangedAnnotation, DiffResult
from allelix.reports.terminal import render_terminal, render_terminal_diff


def _console(width: int) -> tuple[Console, io.StringIO]:
    stream = io.StringIO()
    return Console(file=stream, width=width, color_system=None, force_terminal=False), stream


def _cells(output: str) -> list[str]:
    """Join wrapped lines within each displayed column, not across columns."""
    rows = [
        re.split("[│┃]", line)[1:-1] for line in output.splitlines() if "│" in line or "┃" in line
    ]
    assert rows, output
    return ["".join("".join(row[i].split()) for row in rows) for i in range(len(rows[0]))]


def _answer(genotype: str, text: str) -> None:
    """Keep every allele, or show an explicit prefix and length for each copy."""
    compact = "".join(text.split())
    if genotype in compact:
        return
    alleles = genotype.split("/")
    for allele in set(alleles):
        copies = alleles.count(allele)
        markers = re.findall(rf"(?<![0-9]){len(allele)}(?:bp|nt|bases)\b", compact, re.IGNORECASE)
        assert compact.count(allele[:6]) >= copies and len(markers) >= copies, (
            f"Genotype was silently truncated: {text}"
        )


def _ref(length: int, padded: bool) -> str:
    core = ("ACGT" * (length // 4 + 1))[:length]
    return "T" + core + "C" if padded else core


@pytest.mark.parametrize("width", [80, 250])
@pytest.mark.parametrize("length", [60, 300])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("haploid", [False, True])
def test_extract_long_reference_keeps_every_answer(monkeypatch, width, length, padded, haploid):
    ref = _ref(length, padded)
    long_call = Variant("rs900260101", "22", 20_010_383, ref, ref, haploid=haploid)
    short_call = Variant("rs900260102", "22", 20_010_401, "A", "T")
    no_call = Variant("rs900260103", "22", 20_010_402, "-", "-")
    console, stream = _console(width)
    monkeypatch.setattr(utility, "console", console)
    found = {v.rsid: v for v in [long_call, short_call, no_call]}
    utility._render_extract_table(Path("mock.vcf"), set(found), found)
    output = stream.getvalue()
    cells = _cells(output)
    assert len(cells) == 6, output
    for index, expected in enumerate(["rs900260101", "22", "20,010,383", None, "no", "no"]):
        if expected is not None:
            assert expected in cells[index], output
    for index, expected in [(0, "rs900260102"), (2, "20,010,401"), (4, "yes")]:
        assert expected in cells[index], output
    assert "Position" in cells[2] and "Het?" in cells[4] and "No-call?" in cells[5], output
    assert "rs900260103" in cells[0] and "yes" in cells[5], output
    _answer(long_call.genotype, cells[3])
    assert "…" not in "".join(cells[i] for i in (0, 1, 2, 4, 5)), output


def _vcf(path: Path, ref: str, gt: str, *, alt: str = "A", pos: int = 20_010_383) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n##reference=GRCh38\n"
        "##contig=<ID=22,length=50818468>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMOCK\n"
        f"22\t{pos}\trs900260101\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{gt}\n"
    )


@pytest.mark.parametrize("width", [80, 250])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("gt", ["0", "0/0"])
def test_extract_cli_preserves_long_reference_answer(tmp_path, monkeypatch, width, padded, gt):
    ref = _ref(300, padded)
    vcf = tmp_path / "mock.vcf"
    position = 20_010_382 if padded else 20_010_383
    _vcf(vcf, ref, gt, alt="TAC" if padded else "A", pos=position)
    console, stream = _console(width)
    monkeypatch.setattr(utility, "console", console)
    result = CliRunner().invoke(
        main,
        ["extract", str(vcf), "--snps", "rs900260101", "--data-dir", str(tmp_path / "db")],
    )
    assert result.exit_code == 0, result.output
    cells = _cells(stream.getvalue())
    assert f"{position:,}" in cells[2], stream.getvalue()
    assert "Het?no" in cells[4], stream.getvalue()
    assert "No-call?no" in cells[5], stream.getvalue()
    _answer(f"{ref}/{ref}", cells[3])


def _annotation(genotype: str, rsid: str = "rs900260101") -> Annotation:
    return Annotation(
        source="clinvar",
        rsid=rsid,
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=9.0,
        description="Synthetic finding",
        attribution="ClinVar",
        genotype_match=genotype,
        gene="MOCK",
        condition="Mock trait",
    )


@pytest.mark.parametrize("width", [80, 250])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("haploid", [False, True])
@pytest.mark.parametrize("length", [3, 300])
@pytest.mark.parametrize("view", ["analyze", "new", "changed", "removed"])
def test_report_tables_keep_source_and_classification(width, padded, haploid, length, view):
    ref = _ref(length, padded)
    genotype = ref if haploid else f"{ref}/{ref}"
    annotation = _annotation(genotype)
    short = _annotation("A/G", "rs900260102")
    console, stream = _console(width)
    if view == "analyze":
        result = AnalysisResult(
            file_path=Path("mock.vcf"),
            parser_name="vcf",
            parser_display_name="VCF",
            sample_id="MOCK",
            build="GRCh38",
            total_variants=2,
            skipped_count=0,
            annotators_used=[],
            annotations=[annotation, short],
        )
        assert render_terminal(result, console) == 2
    else:
        diff = DiffResult()
        if view == "new":
            diff.new = [annotation, short]
        elif view == "changed":
            diff.changed = [
                ChangedAnnotation(
                    current=a,
                    previous_significance="clinvar_benign",
                    previous_magnitude=1.0,
                    previous_genotype_match="A/G",
                )
                for a in (annotation, short)
            ]
        else:
            diff.removed = [asdict(annotation), asdict(short)]
        assert render_terminal_diff(diff, console) == 2
    output = stream.getvalue()
    cells = _cells(output)
    # Source, classification and the neighboring short row must survive a long GT.
    for value in ("rs900260101", "rs900260102", "MOCK", "ClinVar", "pathogenic", "9.0"):
        assert any(value in cell for cell in cells), output
    assert any("Mocktrait" in cell for cell in cells), output
    if view == "changed":
        assert any("benign" in cell for cell in cells), output
        assert any("1.0" in cell for cell in cells), output
    gt_column = 7 if view == "changed" else 5
    _answer(genotype, cells[gt_column])
    assert "A/G" in cells[gt_column], output


@pytest.mark.parametrize("width", [80, 250])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("haploid", [False, True])
@pytest.mark.parametrize("command", ["stats", "compare"])
def test_summary_tables_do_not_depend_on_allele_length(
    tmp_path, monkeypatch, width, padded, haploid, command
):
    """These tables show counts rather than calls; long alleles cannot change them."""
    ref = _ref(300, padded)
    vcf = tmp_path / "mock.vcf"
    # Carry the long allele at a non-reference record so the stats stream includes it.
    _vcf(
        vcf,
        "TAC" if padded else "A",
        "1" if haploid else "1/1",
        alt=ref,
        pos=20_010_382 if padded else 20_010_383,
    )
    console, stream = _console(width)
    monkeypatch.setattr(utility, "console", console)
    args = [command, str(vcf)] + ([str(vcf)] if command == "compare" else [])
    if command == "stats":
        args += ["--data-dir", str(tmp_path / "db")]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    compact = "".join(stream.getvalue().split())
    if command == "stats":
        assert "TotalSNPs│1" in compact, stream.getvalue()
        assert "Homozygous│1(100.00%)" in compact, stream.getvalue()
    else:
        assert "Concordant│1│100.00%" in compact, stream.getvalue()
        assert "Discordant│0│0.00%" in compact, stream.getvalue()

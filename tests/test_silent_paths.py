# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Paths that could read a file wrongly with no warning, now fixed.

- No-call codes: ``00`` and ``NN`` (``0`` and ``N`` per allele column, or
  alone on a one-letter row) are no-calls in every array reader, as ``--``
  is, so ``stats`` and the high-value no-call warning count them.
- A ``--filter`` list saved with a byte order mark keeps its first entry.
- Genes for Good's imputed VCF (its ``##Genes_for_Good_version`` line and
  phased calls) is marked imputed, as its text download is.
- A GRCh37 VCF that declares no length for the mitochondrion and names it
  ``chrM`` or ``M`` (UCSC hg19's names, the 16,571 bp sequence) has its MT
  records left out with a note, instead of being read as the rCRS.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from test_sources import data  # noqa: F401  (the data-folder fixture)

from allelix.cli import main
from allelix.cli._helpers import _parse_filter_file
from allelix.parsers.ancestrydna import AncestryDNAParser
from allelix.parsers.decodeme import DeCODEmeParser
from allelix.parsers.ftdna import FTDNAParser
from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser
from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser
from allelix.parsers.genesforgood import GenesForGoodParser
from allelix.parsers.livingdna import LivingDNAParser
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.myheritage import MyHeritageParser
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import (
    BuildDiagnostics,
    StreamBuild,
    decide_file_build,
    with_effective_build,
)
from allelix.reports.build_text import CHRM, build_notices
from allelix.utils.build_detect import _CONTIG_LENGTHS, KNOWN_SNP_POSITIONS

if TYPE_CHECKING:
    from collections.abc import Callable

    from allelix.models import Variant
    from allelix.parsers.base import GenotypeParser
    from allelix.utils.build_detect import BuildDecision

_FIXTURES = Path(__file__).parent / "fixtures"

# ── No-call codes ─────────────────────────────────────────────────────

# Synthetic rows: (rsID, chromosome, position, genotype as the file writes it).
_DIPLOID_NO_CALLS = [
    ("rs3999990101", "2", 3999101, "--"),
    ("rs3999990102", "2", 3999102, "00"),
    ("rs3999990103", "2", 3999103, "NN"),
    ("rs3999990104", "2", 3999104, "nn"),
]
_HAPLOID_NO_CALLS = [
    ("rs3999990105", "MT", 3999105, "-"),
    ("rs3999990106", "MT", 3999106, "0"),
    ("rs3999990107", "MT", 3999107, "N"),
]
_CALLED = [("rs3999990108", "2", 3999108, "AG"), ("rs3999990109", "MT", 3999109, "A")]


def _tab(r: str, c: str, p: int, g: str) -> str:
    return f"{r}\t{c}\t{p}\t{g}"


def _csv(r: str, c: str, p: int, g: str) -> str:
    return f'"{r}","{c}","{p}","{g}"'


def _columns(r: str, c: str, p: int, g: str) -> str:
    """Two allele columns; a one-letter row writes the letter in both."""
    a1, a2 = (g[0], g[1]) if len(g) == 2 else (g, g)
    return f"{r}\t{c}\t{p}\t{a1}\t{a2}"


def _decodeme(r: str, c: str, p: int, g: str) -> str:
    return f"{r},A/G,{c},{p},+,{g}"


# Reader, its fixture (the envelope), the row writer, and whether it reads
# one-letter rows (the two-column readers write both columns instead).
_READERS: list[tuple[type[GenotypeParser], str, Callable[[str, str, int, str], str], bool]] = [
    (TwentyThreeAndMeParser, "mock_23andme.txt", _tab, True),
    (GenesForGoodParser, "mock_genesforgood.txt", _tab, True),
    (AncestryDNAParser, "mock_ancestrydna.txt", _columns, False),
    (FTDNAParser, "mock_ftdna.csv", _csv, True),
    (FTDNAIlluminaParser, "mock_ftdna_illumina.txt", _tab, True),
    (FTDNAFamFinderParser, "mock_ftdna_famfinder.txt", _columns, False),
    (LivingDNAParser, "mock_livingdna.csv", _tab, True),
    (MyHeritageParser, "mock_myheritage.csv", _csv, True),
    (MyHappyGenesParser, "mock_myhappygenes.txt", _columns, False),
    (DeCODEmeParser, "mock_decodeme.csv", _decodeme, True),
]


def _with_rows(tmp_path: Path, fixture: str, rows: list[str]) -> Path:
    """The fixture with ``rows`` appended, in the fixture's own line ending."""
    raw = (_FIXTURES / fixture).read_bytes()
    eol = b"\r\n" if b"\r\n" in raw else b"\n"
    if not raw.endswith(eol):
        raw += eol
    path = tmp_path / fixture
    path.write_bytes(raw + b"".join(r.encode() + eol for r in rows))
    return path


@pytest.mark.parametrize(
    ("reader", "fixture", "row", "one_letter"), _READERS, ids=lambda v: getattr(v, "name", "")
)
def test_no_call_codes_in_every_array_reader(
    tmp_path: Path,
    reader: type[GenotypeParser],
    fixture: str,
    row: Callable[[str, str, int, str], str],
    one_letter: bool,
) -> None:
    rows = _DIPLOID_NO_CALLS + (_HAPLOID_NO_CALLS if one_letter else []) + _CALLED
    path = _with_rows(tmp_path, fixture, [row(*r) for r in rows])
    calls = {v.rsid: v for v in reader().parse(path)}
    for rsid, *_ in _DIPLOID_NO_CALLS + (_HAPLOID_NO_CALLS if one_letter else []):
        assert calls[rsid].is_no_call, rsid
        assert (calls[rsid].allele1, calls[rsid].allele2) == ("-", "-"), rsid
        assert not calls[rsid].haploid, rsid
    assert not calls["rs3999990108"].is_no_call
    assert (calls["rs3999990108"].allele1, calls["rs3999990108"].allele2) == ("A", "G")
    assert not calls["rs3999990109"].is_no_call


def test_one_called_allele_beside_a_no_call_code(tmp_path: Path) -> None:
    path = _with_rows(tmp_path, "mock_ancestrydna.txt", ["rs3999990110\t2\t3999110\tA\t0"])
    (call,) = [v for v in AncestryDNAParser().parse(path) if v.rsid == "rs3999990110"]
    assert (call.allele1, call.allele2) == ("A", "-")
    assert call.is_no_call


def _stats(path: Path) -> str:
    result = CliRunner().invoke(main, ["stats", str(path)])
    assert result.exit_code == 0, result.output
    return " ".join(result.output.split())


@pytest.mark.parametrize("code", ["00", "NN"])
def test_stats_counts_the_code_as_a_no_call(tmp_path: Path, code: str) -> None:
    """``stats`` and the high-value warning see ``00``/``NN`` as ``--`` (APOE rs429358)."""
    base = (_FIXTURES / "mock_23andme.txt").read_text()
    assert "rs429358\t19\t45411941\t" in base
    dashes = tmp_path / "dashes.txt"
    coded = tmp_path / "coded.txt"
    for path, value in ((dashes, "--"), (coded, code)):
        path.write_text(
            "".join(
                f"rs429358\t19\t45411941\t{value}\n" if line.startswith("rs429358\t") else line
                for line in base.splitlines(keepends=True)
            )
        )
    expected = _stats(dashes)
    assert "rs429358" in expected
    assert _stats(coded) == expected.replace(dashes.name, coded.name)


# ── A --filter list saved with a byte order mark ──────────────────────


@pytest.mark.parametrize(
    "raw",
    [
        b"\xef\xbb\xbfrs1801133\nMTHFR\n",
        b"\xef\xbb\xbfrs1801133\r\nMTHFR\r\n",
        b"\xef\xbb\xbfrs1801133\rMTHFR\r",
        b"rs1801133\nMTHFR\n",
        b"rs1801133\r\nMTHFR\r\n",
    ],
    ids=["bom-lf", "bom-crlf", "bom-cr", "lf", "crlf"],
)
def test_filter_keeps_its_first_rsid(tmp_path: Path, raw: bytes) -> None:
    path = tmp_path / "filter.txt"
    path.write_bytes(raw)
    assert _parse_filter_file(path) == (frozenset({"MTHFR"}), frozenset({"rs1801133"}))


def test_filter_with_a_bom_first_gene(tmp_path: Path) -> None:
    path = tmp_path / "filter.txt"
    path.write_bytes(b"\xef\xbb\xbfMTHFR\nrs4680\n")
    assert _parse_filter_file(path) == (frozenset({"MTHFR"}), frozenset({"rs4680"}))


# ── Genes for Good's imputed VCF ──────────────────────────────────────

_GFG_LINE = "##Genes_for_Good_version=v0.0 placeholder\n"


def _gfg_vcf(path: Path, gts: list[str], *, gfg_line: bool = True) -> Path:
    body = "".join(
        f"1\t{3999200 + i}\trs{3999990200 + i}\tA\tG\t.\tPASS\t.\tGT\t{gt}\n"
        for i, gt in enumerate(gts)
    )
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + (_GFG_LINE if gfg_line else "")
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + body
    )
    return path


@pytest.mark.parametrize(
    ("gts", "gfg_line", "method"),
    [
        (["0|1", "1|1", "0|0"], True, "array-imputed"),
        (["0|1", "./.", ".|1", "1", "0|0"], True, "array-imputed"),
        (["0|1", "1|1", "0/0"], True, None),
        (["0/1", "1/1"], True, None),
        (["./.", "1"], True, None),
        (["0|1", "1|1", "0|0"], False, None),
    ],
    ids=[
        "phased",
        "phased-with-missing-and-haploid",
        "one-unphased",
        "unphased",
        "none-called",
        "no-gfg-line",
    ],
)
def test_genes_for_good_vcf_is_marked_imputed(
    tmp_path: Path, gts: list[str], gfg_line: bool, method: str | None
) -> None:
    path = _gfg_vcf(tmp_path / "gfg.vcf", gts, gfg_line=gfg_line)
    assert VcfParser().get_metadata(path).get("method") == method


def test_phase_check_reads_a_bounded_number_of_calls(tmp_path: Path) -> None:
    """Past the first 1,000 phased calls an unphased one is not looked at."""
    gts = ["0|1"] * 1000 + ["0/1"]
    path = _gfg_vcf(tmp_path / "gfg.vcf", gts)
    assert VcfParser().get_metadata(path).get("method") == "array-imputed"
    path = _gfg_vcf(tmp_path / "gfg2.vcf", ["0|1"] * 999 + ["0/1"])
    assert VcfParser().get_metadata(path).get("method") is None


def test_genes_for_good_vcf_report_says_imputed(tmp_path: Path, data: Path) -> None:  # noqa: F811
    path = _gfg_vcf(tmp_path / "gfg.vcf", ["0|1", "1|1"])
    out = tmp_path / "r.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(path),
            "--build",
            "grch37",
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(out.read_text())["input"]["method"] == "array-imputed"
    assert "genotypes are imputed" in " ".join(result.output.split())


# ── The mitochondrion on GRCh37: decided per record name ─────────────
#
# A record named chrM or M (UCSC hg19's names; hg19's chrM is a 16,571 bp
# sequence, not the rCRS ClinVar uses) is read only when that same name is
# declared with the rCRS length (16,569); otherwise it is left out, counted,
# with the note. MT and chrMT are GRCh37's own names and are read. A
# declared hg19 length (16,571) leaves out every MT record, as before.

_SEEDS37 = [v["GRCh37"] for v in KNOWN_SNP_POSITIONS.values()][:4]
_MT_POSITIONS = (73, 263, 750)

# The header's declaration for the mitochondrion: (name, length) or None.
_DECLARATIONS = [
    None,
    ("MT", None),
    ("MT", 16569),
    ("MT", 16571),
    ("chrMT", 16569),
    ("chrM", None),
    ("chrM", 16569),
    ("chrM", 16571),
    ("M", 16569),
]
_RECORD_NAMES = ["MT", "chrMT", "chrM", "M"]
READ, UCSC, HG19 = "read", "left out: UCSC name", "left out: hg19 length"


def _expected(declared: tuple[str, int | None] | None, record: str) -> str:
    """The rule, written out once."""
    if declared is not None and declared[1] == 16571:
        return HG19
    if record.upper() in ("CHRM", "M"):
        return READ if declared == (record, 16569) else UCSC
    return READ


def _mt_vcf(
    path: Path,
    *,
    nuclear: str,
    mito: str,
    contigs: str = "",
    seeds: list[tuple[str, int]] = _SEEDS37,
) -> Path:
    rows = [f"{nuclear}{c}\t{p}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n" for c, p in seeds]
    rows += [f"{mito}\t{p}\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n" for p in _MT_POSITIONS]
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + contigs
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "".join(rows)
    )
    return path


def _contig_lines(
    prefix: str, declared: tuple[str, int | None] | None, *, nuclear_lengths: bool
) -> str:
    """Nuclear ``##contig`` lines (with GRCh37 lengths, or none) and the mito declaration."""
    lines = ""
    if nuclear_lengths:
        for c in sorted({c for c, _ in _SEEDS37}):
            lines += f"##contig=<ID={prefix}{c},length={_CONTIG_LENGTHS['GRCh37'][c]}>\n"
    if declared is not None:
        name, length = declared
        lines += f"##contig=<ID={name}" + (f",length={length}" if length else "") + ">\n"
    return lines


def _run(path: Path) -> tuple[BuildDecision, list[Variant], BuildDiagnostics]:
    parser = VcfParser()
    decision = decide_file_build(parser, path, parser.get_metadata(path), override=None)
    record = StreamBuild(decision, parser=parser)
    kept = list(
        with_effective_build(
            parser.parse(path), decision=decision, file_name=path.name, record=record
        )
    )
    return decision, kept, record.diagnostics()


@pytest.mark.parametrize("lengths_decide", [True, False], ids=["lengths-decide", "pre-scan"])
@pytest.mark.parametrize("prefix", ["chr", ""], ids=["chr-nuclear", "plain-nuclear"])
@pytest.mark.parametrize("record", _RECORD_NAMES)
@pytest.mark.parametrize(
    "declared",
    _DECLARATIONS,
    ids=lambda d: "none" if d is None else f"{d[0]}-{d[1] or 'nolength'}",
)
def test_mito_matrix(
    tmp_path: Path,
    declared: tuple[str, int | None] | None,
    record: str,
    prefix: str,
    lengths_decide: bool,
) -> None:
    path = _mt_vcf(
        tmp_path / "x.vcf",
        nuclear=prefix,
        mito=record,
        contigs=_contig_lines(prefix, declared, nuclear_lengths=lengths_decide),
    )
    decision, kept, diag = _run(path)
    assert decision.build == "GRCh37"
    assert (decision.source == "contig_lengths") is lengths_decide
    expected = _expected(declared, record)
    mt_kept = sum(v.chromosome == "MT" for v in kept)
    notes = [n.text for n in build_notices(diag) if n.kind == CHRM]
    if expected == READ:
        assert (mt_kept, diag.mt_skipped, notes) == (3, 0, [])
    else:
        assert (mt_kept, diag.mt_skipped) == (0, 3)
        (note,) = notes
        if expected == UCSC:
            assert note.startswith("Note: 3 mitochondrial record(s) named chrM or M")
        else:
            assert note.startswith("Note: chrM is declared with the hg19 length")
    assert len(kept) - mt_kept == len(_SEEDS37)


def test_matrix_expectations_as_written() -> None:
    """The cells the PR body lists, spelt out (guards the rule helper itself)."""
    assert _expected(None, "MT") == READ
    assert _expected(None, "chrM") == UCSC
    assert _expected(("MT", 16569), "chrM") == UCSC
    assert _expected(("chrM", 16569), "chrM") == READ
    assert _expected(("chrM", 16569), "M") == UCSC
    assert _expected(("M", 16569), "M") == READ
    assert _expected(("chrMT", 16569), "chrM") == UCSC
    assert _expected(("chrM", None), "chrM") == UCSC
    assert _expected(("MT", 16571), "MT") == HG19
    assert _expected(("chrM", 16571), "chrMT") == HG19


@pytest.mark.parametrize("record", ["chrM", "MT"])
def test_another_declared_length_is_not_the_rcrs(tmp_path: Path, record: str) -> None:
    """chrM declared with a length that is neither 16,569 nor 16,571: not shown to be the rCRS."""
    path = _mt_vcf(
        tmp_path / "x.vcf",
        nuclear="chr",
        mito=record,
        contigs=_contig_lines("chr", ("chrM", 16570), nuclear_lengths=True),
    )
    _decision, kept, diag = _run(path)
    mt_kept = sum(v.chromosome == "MT" for v in kept)
    assert (mt_kept, diag.mt_skipped) == ((0, 3) if record == "chrM" else (3, 0))


def test_grch38_chrm_is_read(tmp_path: Path) -> None:
    seeds38 = [v["GRCh38"] for v in KNOWN_SNP_POSITIONS.values()][:4]
    path = _mt_vcf(tmp_path / "x.vcf", nuclear="chr", mito="chrM", seeds=seeds38)
    decision, kept, diag = _run(path)
    assert decision.build == "GRCh38"
    assert sum(v.chromosome == "MT" for v in kept) == 3
    assert diag.mt_skipped == 0


def test_flag_is_set_again_on_every_decision(tmp_path: Path) -> None:
    parser = VcfParser()
    seeds38 = [v["GRCh38"] for v in KNOWN_SNP_POSITIONS.values()][:4]
    grch37 = _mt_vcf(tmp_path / "a.vcf", nuclear="chr", mito="chrM")
    grch38 = _mt_vcf(tmp_path / "b.vcf", nuclear="chr", mito="chrM", seeds=seeds38)
    hg19 = _mt_vcf(
        tmp_path / "c.vcf",
        nuclear="chr",
        mito="chrM",
        contigs=_contig_lines("chr", ("chrM", 16571), nuclear_lengths=True),
    )
    decide_file_build(parser, grch37, parser.get_metadata(grch37), override=None)
    assert parser.drop_ucsc_mito is True
    decide_file_build(parser, grch38, parser.get_metadata(grch38), override=None)
    assert parser.drop_ucsc_mito is False
    decide_file_build(parser, grch37, parser.get_metadata(grch37), override=None)
    decide_file_build(parser, hg19, parser.get_metadata(hg19), override=None)
    assert parser.drop_ucsc_mito is False  # the declared hg19 length already leaves MT out


def test_prescan_does_not_drop(tmp_path: Path) -> None:
    """A parser left flagged by an earlier file reads every record in the next pre-scan."""
    parser = VcfParser()
    parser.drop_ucsc_mito = True
    path = _mt_vcf(tmp_path / "x.vcf", nuclear="", mito="chrM")
    seen: list[int] = []
    original = parser.parse

    def spy(file_path: Path):
        seen.append(int(parser.drop_ucsc_mito))
        yield from original(file_path)

    parser.parse = spy  # type: ignore[method-assign]
    decide_file_build(parser, path, parser.get_metadata(path), override=None)
    assert seen == [0]


@pytest.mark.parametrize("command", ["analyze", "stats"])
def test_ucsc_drop_reaches_the_report(tmp_path: Path, data: Path, command: str) -> None:  # noqa: F811
    path = _mt_vcf(
        tmp_path / "x.vcf",
        nuclear="",
        mito="chrM",
        contigs=_contig_lines("", ("MT", 16569), nuclear_lengths=True),
    )
    out = tmp_path / "r.json"
    args = [command, str(path)]
    if command == "analyze":
        args += ["--data-dir", str(data), "--no-update", "--output", str(out)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    flat = " ".join(result.output.split())
    assert "3 mitochondrial record(s) named chrM or M were not annotated" in flat
    if command == "analyze":
        payload = json.loads(out.read_text())
        assert payload["input"]["build_evidence"]["mt_records_skipped"] == 3

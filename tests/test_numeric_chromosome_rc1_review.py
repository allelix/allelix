"""Outside reviewer A's ruled A-02/A-03/A-04 CLI regressions, synthetic only.

Adapted for the repository import layout; behavior assertions preserved.
A-01 is a documentation correction, not a demand to change VCF parsing.
"""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity

from allelix.cli import main

SOURCES = (
    "clinvar",
    "pharmgkb",
    "snpedia",
    "gwas",
    "gnomad",
    "cadd",
    "alphamissense",
    "rsid_identity",
)
L37 = {"1": 249250621, "X": 155270560, "Y": 59373566}
RS_X, POS_X = "rs900263001", 3_000_001


def _config(data: Path, on: set[str]) -> None:
    data.mkdir(parents=True, exist_ok=True)
    (data / "config.toml").write_text(
        "[sources]\n" + "".join(f"{n} = {str(n in on).lower()}\n" for n in SOURCES)
    )


def _data(tmp_path: Path) -> Path:
    """ClinVar (both builds) and the identity table: one X site and one autosome site."""
    data = tmp_path / "data"
    _config(data, {"clinvar", "rsid_identity"})
    rows = [(RS_X, "X", POS_X, "G", "A"), ("rs900263005", "1", 1_000_005, "G", "A")]
    for build in ("GRCh37", "GRCh38"):
        _clinvar(
            data,
            build,
            [(r, c, p, a, b, "SYN", "Pathogenic", "Pathogenic", 1, 1) for r, c, p, a, b in rows],
        )
    _identity(data, [(r, c, p, a, b, p, a, b, 0) for r, c, p, a, b in rows])
    return data


def _vcf(path: Path, contigs: list[tuple[str, int | None]], records: list[tuple]) -> Path:
    head = "##fileformat=VCFv4.2\n" + "".join(
        f"##contig=<ID={c}" + (f",length={n}" if n else "") + ">\n" for c, n in contigs
    )
    body = "".join("\t".join(map(str, r)) + "\n" for r in records)
    path.write_text(head + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + body)
    return path


def _invoke(args: list[str]):
    return CliRunner(env={"COLUMNS": "250"}).invoke(main, args)


def _export_plink(src: Path, data: Path, prefix: Path, build: str) -> list[list[str]]:
    result = _invoke(
        [
            "export",
            "plink",
            str(src),
            "--build",
            build,
            "--data-dir",
            str(data),
            "--output-prefix",
            str(prefix),
        ]
    )
    assert result.exit_code == 0, (result.output, result.exception)
    return [line.split("\t") for line in Path(str(prefix) + ".bim").read_text().splitlines()]


# A-02 -------------------------------------------------------------------------
# New in 2.6.3: a PAR contig (25 or XY) that declares its own length is compared
# with X's length, so a file whose 1, X and Y lengths are GRCh37's is refused as
# "unsupported assembly". 2.6.2 ignored the PAR contig and decided GRCh37.
def test_a02_par_contig_length_does_not_refuse_a_grch37_file(tmp_path):
    data = _data(tmp_path)
    src = _vcf(
        tmp_path / "par.vcf",
        [("1", L37["1"]), ("X", L37["X"]), ("Y", L37["Y"]), ("25", 2_699_520)],
        [
            ("X", POS_X, RS_X, "G", "A", ".", "PASS", ".", "GT", "0/1"),
            ("25", 100_001, "rs900263002", "C", "T", ".", "PASS", ".", "GT", "0/1"),
            ("1", 1_000_005, "rs900263005", "G", "A", ".", "PASS", ".", "GT", "0/1"),
        ],
    )
    out = tmp_path / "r.json"
    result = _invoke(
        ["analyze", str(src), "--data-dir", str(data), "--no-update", "--output", str(out)]
    )
    assert result.exit_code == 0, result.output
    report = json.loads(out.read_text())
    assert report["input"]["build"] == "GRCh37"


# A-03 -------------------------------------------------------------------------
# New in 2.6.3, not in the CHANGELOG: the ID `export plink` generates for an
# rsID-less record on a numeric contig changed from `23:pos:ref:alt` (2.6.2) to
# `X:pos:ref:alt`, and `25:...` to `X:...` on a row whose CHROM column is 25.
def test_a03_generated_plink_id_keeps_the_files_chromosome(tmp_path):
    data = tmp_path / "empty"
    _config(data, set())
    src = _vcf(
        tmp_path / "unnamed.vcf",
        [("23", None), ("25", None)],
        [
            ("23", POS_X, ".", "G", "A", ".", "PASS", ".", "GT", "0/1"),
            ("25", 100_001, ".", "C", "T", ".", "PASS", ".", "GT", "0/1"),
        ],
    )
    bim = _export_plink(src, data, tmp_path / "out", "grch38")
    ids = sorted(row[1] for row in bim)
    assert ids == ["23:3000001:G:A", "25:100001:C:T"], bim


# A-04 -------------------------------------------------------------------------
# Pre-existing: `extract` prefers a call at the rsID's placement over one elsewhere
# (ADR-0054), but compares the reader's raw chromosome with the placement's. For a
# reader that keeps the raw code (MyHappyGenes, 23andMe, MyHeritage, Living DNA,
# FTDNA CSV), a call written `23` at the placement loses to a call elsewhere; the
# same file written `X` (or a deCODEme file written `23`, in rc1) shows the placed call.
def _mhg(path: Path, chrom: str) -> Path:
    path.write_text(
        "# MyHappyGenes [TEMPUS]\n# Sample ID\tMHG1\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        f"{RS_X}\t1\t5000\tA\tA\n"
        f"{RS_X}\t{chrom}\t{POS_X}\tG\tA\n"
    )
    return path


def _extract_row(src: Path, data: Path) -> list[str]:
    result = _invoke(
        ["extract", str(src), "--snps", RS_X, "--build", "grch37", "--data-dir", str(data)]
    )
    assert result.exit_code == 0, (result.output, result.exception)
    rows = [
        [c.strip() for c in line.split("│")[1:-1]]
        for line in result.output.splitlines()
        if RS_X in line
    ]
    assert len(rows) == 1, result.output
    return rows[0]


def test_a04_extract_prefers_the_placed_call_for_a_numeric_code(tmp_path):
    data = _data(tmp_path)
    control = _extract_row(_mhg(tmp_path / "letters.txt", "X"), data)
    assert control[2:4] == ["3,000,001", "G/A"], control
    actual = _extract_row(_mhg(tmp_path / "numeric.txt", "23"), data)
    assert actual[2:4] == control[2:4], actual

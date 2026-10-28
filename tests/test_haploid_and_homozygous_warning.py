"""File-declared haploid calls and the multi-homozygote warning (ADR-0050).

A one-letter array call or a one-index VCF GT is one copy: the parsers
mark it, matching still sees the allele twice, and the reports show it
once as hemizygous. Several rare homozygous ClinVar Pathogenic / Likely
pathogenic calls on autosomes give the report a warning. All genotypes
and ClinVar rows are synthetic.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import is_pathogenic_aggregate
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA
from allelix.parsers.ftdna import FTDNAParser
from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser
from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser
from allelix.parsers.livingdna import LivingDNAParser
from allelix.parsers.myheritage import MyHeritageParser
from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser

_FIXTURES = Path(__file__).parent / "fixtures"


# (parser, mock fixture, one-letter X row, two-letter row)
_ARRAY_FORMATS = [
    (
        TwentyThreeAndMeParser,
        "mock_23andme.txt",
        "rs900950001\tX\t5000\tT",
        "rs900950002\t1\t6000\tCT",
    ),
    (LivingDNAParser, "mock_livingdna.csv", "rs900950001\tX\t5000\tT", "rs900950002\t1\t6000\tCT"),
    (
        MyHeritageParser,
        "mock_myheritage.csv",
        '"rs900950001","X","5000","T"',
        '"rs900950002","1","6000","CT"',
    ),
    (
        FTDNAParser,
        "mock_ftdna.csv",
        '"rs900950001","X","5000","T"',
        '"rs900950002","1","6000","CT"',
    ),
    (
        FTDNAIlluminaParser,
        "mock_ftdna_illumina.txt",
        "rs900950001\tX\t5000\tT",
        "rs900950002\t1\t6000\tCT",
    ),
    (
        FTDNAFamFinderParser,
        "mock_ftdna_famfinder.txt",
        "rs900950001\tMT\t5000\tT",
        "rs900950002\t1\t6000\tC\tT",
    ),
]


class TestParsersMarkHaploid:
    @pytest.mark.parametrize(("parser_cls", "fixture", "haploid_row", "pair_row"), _ARRAY_FORMATS)
    def test_one_letter_call_is_haploid(
        self, tmp_path: Path, parser_cls, fixture: str, haploid_row: str, pair_row: str
    ) -> None:
        text = (_FIXTURES / fixture).read_text().rstrip("\n")
        path = tmp_path / fixture
        path.write_text(f"{text}\n{haploid_row}\n{pair_row}\n")
        variants = {v.rsid: v for v in parser_cls().parse(path)}
        one, two = variants["rs900950001"], variants["rs900950002"]
        assert (one.allele1, one.allele2, one.haploid) == ("T", "T", True)
        assert (two.allele1, two.allele2, two.haploid) == ("C", "T", False)
        # A no-call is never haploid.
        assert not any(v.haploid for v in variants.values() if v.is_no_call)

    def test_vcf_one_index_gt_is_haploid(self, tmp_path: Path) -> None:
        path = tmp_path / "x.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "X\t5000\trs900950001\tC\tT\t50\tPASS\t.\tGT\t1\n"
            "1\t6000\trs900950002\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
            "X\t7000\trs900950003\tC\tT\t50\tPASS\t.\tGT\t.\n"
        )
        variants = {v.rsid: v for v in VcfParser().parse(path)}
        assert variants["rs900950001"].haploid
        assert (variants["rs900950001"].allele1, variants["rs900950001"].allele2) == ("T", "T")
        assert not variants["rs900950002"].haploid
        assert "rs900950003" not in variants or not variants["rs900950003"].haploid


# ---------------------------------------------------------------------------
# End to end. ClinVar rows in the unused 900950xxx range.
# ---------------------------------------------------------------------------

_SINGLE = "criteria_provided,_single_submitter"


def _clinvar(data: Path, rows: list[tuple[str, str, int, str, str, str]]) -> None:
    """Both builds' caches with the same rows: (rsid, chrom, pos, ref, alt, aggregate)."""
    for build in ("GRCh37", "GRCh38"):
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            conn.executemany(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
                "submitters_same, submitters_total) "
                "VALUES (?, ?, ?, ?, ?, 'Pathogenic', 'Test condition', 'GENE', ?, ?, ?, 1, 1)",
                [
                    (rsid, chrom, pos, ref, alt, _SINGLE, i, aggregate)
                    for i, (rsid, chrom, pos, ref, alt, aggregate) in enumerate(rows, 1)
                ],
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, 'test://', '1', '2026-10-10', ?, "
                "NULL, ?)",
                (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()


def _rare_in_gnomad(data: Path, rows: tuple[tuple[str, str, int, str, str, str], ...]) -> None:
    """A ready gnomAD giving each row's ALT a rare frequency, reached on GRCh37.

    The identity table places each row at the same position on both builds,
    same strand, so the GRCh37 key translates to the gnomAD key.
    """
    write_identity_db(
        [(rsid, chrom, pos, ref, alt, pos, ref, alt, 0) for rsid, chrom, pos, ref, alt, _ in rows],
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2001-01-01T00:00:00Z",
    )
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, 0.000001)",
            [(chrom, pos, ref, alt, rsid) for rsid, chrom, pos, ref, alt, _ in rows],
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('gnomad', 'test://', '1', '2001-01-01', ?, "
            "'sv:1')",
            (len(rows),),
        )
        conn.commit()


def _array(path: Path, calls: list[tuple[str, str, int, str]]) -> Path:
    # A 23andMe banner, so --format 23andme is a choice, not a forced mismatch.
    lines = [
        "# This data file generated by 23andMe at: Mon Jan 01 00:00:00 2001",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    lines += [f"{rsid}\t{chrom}\t{pos}\t{gt}" for rsid, chrom, pos, gt in calls]
    path.write_text("\n".join(lines) + "\n")
    return path


def _report(tmp_path: Path, data: Path, array: Path, *extra: str) -> dict:
    out = tmp_path / "r.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(array),
            "--format",
            "23andme",
            "--build",
            "grch37",
            "--no-update",
            "--data-dir",
            str(data),
            "--output",
            str(out),
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output
    return json.loads(out.read_text())


@pytest.fixture
def data(tmp_path: Path) -> Path:
    d = tmp_path / "data"
    d.mkdir()
    return d


class TestHemizygousInReports:
    def test_one_letter_x_call_is_hemizygous(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, [("rs900950010", "X", 4000, "C", "T", "Pathogenic")])
        array = _array(tmp_path / "a.txt", [("rs900950010", "X", 4000, "T")])
        payload = _report(tmp_path, data, array)
        (row,) = payload["annotations"]
        assert row["genotype_match"] == "T"
        assert row["zygosity"] == "hemizygous"
        assert row["plausibility"] == "not_assessed"
        assert row["plausibility_reason"] == "gnomad_not_run"

    def test_haploid_vcf_call_matched_after_trimming(self, tmp_path: Path, data: Path) -> None:
        """A one-index GT whose allele ClinVar matches only after trimming
        (``GAT>GA`` as ``AT>A`` one base on, ADR-0047) is hemizygous too."""
        _clinvar(
            data,
            [
                ("rs900950031", "X", 5001, "AT", "A", "Pathogenic"),
                ("rs900950032", "X", 6000, "G", "A", "Pathogenic"),
            ],
        )
        vcf = tmp_path / "x.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "X\t5000\t.\tGAT\tGA\t50\tPASS\t.\tGT\t1\n"
            "X\t6000\t.\tG\tA\t50\tPASS\t.\tGT\t1\n"
        )
        out = tmp_path / "r.json"
        result = CliRunner().invoke(
            main,
            [
                "analyze",
                str(vcf),
                "--build",
                "grch37",
                "--no-update",
                "--data-dir",
                str(data),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        rows = {a["position"]: a for a in json.loads(out.read_text())["annotations"]}
        assert (rows[6000]["genotype_match"], rows[6000]["zygosity"]) == ("A", "hemizygous")
        # The trimmed match: ClinVar's row is at 5001; the record is at 5000.
        assert (rows[5001]["genotype_match"], rows[5001]["zygosity"]) == ("GA", "hemizygous")

    def test_two_letter_x_call_stays_homozygous(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, [("rs900950010", "X", 4000, "C", "T", "Pathogenic")])
        array = _array(tmp_path / "a.txt", [("rs900950010", "X", 4000, "TT")])
        (row,) = _report(tmp_path, data, array)["annotations"]
        assert (row["genotype_match"], row["zygosity"]) == ("TT", "homozygous")


_HOM_ROWS = (
    ("rs900950021", "2", 1000, "G", "A", "Pathogenic"),
    ("rs900950022", "3", 1000, "G", "A", "Likely pathogenic"),
    ("rs900950023", "4", 1000, "G", "A", "Pathogenic/Likely pathogenic; risk factor"),
    ("rs900950024", "5", 1000, "G", "A", "Conflicting classifications of pathogenicity"),
    ("rs900950025", "X", 1000, "G", "A", "Pathogenic"),
)


class TestMultipleHomozygousWarning:
    def test_three_rare_homozygous_plp_calls_warn(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, _HOM_ROWS)
        _rare_in_gnomad(data, _HOM_ROWS)
        calls = [(rsid, chrom, pos, "AA") for rsid, chrom, pos, *_ in _HOM_ROWS]
        payload = _report(tmp_path, data, _array(tmp_path / "a.txt", calls))
        warnings = [
            w for w in payload["report_warnings"] if w["code"] == "multiple_homozygous_plp"
        ]
        # The Conflicting aggregate and the X call do not count.
        assert [w["count"] for w in warnings] == [3]
        assert warnings[0]["message"].startswith("3 variants are called homozygous")

    def test_the_console_prints_it(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, _HOM_ROWS)
        _rare_in_gnomad(data, _HOM_ROWS)
        calls = [(rsid, chrom, pos, "AA") for rsid, chrom, pos, *_ in _HOM_ROWS]
        array = _array(tmp_path / "a.txt", calls)
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main,
            [
                "analyze",
                str(array),
                "--format",
                "23andme",
                "--build",
                "grch37",
                "--no-update",
                "--data-dir",
                str(data),
                "--output",
                str(tmp_path / "r.json"),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "Warning: 3 variants are called homozygous" in result.output

    def test_the_display_filter_cannot_hide_it(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, _HOM_ROWS)
        _rare_in_gnomad(data, _HOM_ROWS)
        calls = [(rsid, chrom, pos, "AA") for rsid, chrom, pos, *_ in _HOM_ROWS]
        payload = _report(
            tmp_path, data, _array(tmp_path / "a.txt", calls), "--min-magnitude", "10"
        )
        assert payload["annotations"] == []
        assert [w["count"] for w in payload["report_warnings"]] == [3]

    def test_missing_frequencies_do_not_warn(self, tmp_path: Path, data: Path) -> None:
        """Without gnomAD no call was assessed, so none counts as rare (2.5.2)."""
        _clinvar(data, _HOM_ROWS)
        calls = [(rsid, chrom, pos, "AA") for rsid, chrom, pos, *_ in _HOM_ROWS]
        payload = _report(tmp_path, data, _array(tmp_path / "a.txt", calls))
        assert payload["report_warnings"] == []

    def test_two_do_not_warn(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, _HOM_ROWS)
        calls = [(rsid, chrom, pos, "AA") for rsid, chrom, pos, *_ in _HOM_ROWS[:2]]
        calls += [(rsid, chrom, pos, "GA") for rsid, chrom, pos, *_ in _HOM_ROWS[2:3]]
        payload = _report(tmp_path, data, _array(tmp_path / "a.txt", calls))
        assert payload["report_warnings"] == []


@pytest.mark.parametrize(
    ("aggregate", "expected"),
    [
        ("Pathogenic", True),
        ("Likely pathogenic", True),
        ("Pathogenic/Likely pathogenic", True),
        ("Pathogenic/Likely pathogenic; risk factor", True),
        ("Pathogenic, low penetrance", True),
        ("Likely_pathogenic", True),
        ("Pathogenic/Likely risk allele", False),
        ("Conflicting classifications of pathogenicity", False),
        ("Uncertain significance", False),
        ("", False),
    ],
)
def test_pathogenic_aggregate(aggregate: str, expected: bool) -> None:
    assert is_pathogenic_aggregate(aggregate) is expected


class TestNotAssessedShown:
    def test_html_and_terminal(self, tmp_path: Path, data: Path) -> None:
        _clinvar(data, [("rs900950010", "2", 4000, "C", "T", "Pathogenic")])
        array = _array(tmp_path / "a.txt", [("rs900950010", "2", 4000, "TT")])
        html = tmp_path / "r.html"
        args = ["analyze", str(array), "--format", "23andme", "--build", "grch37"]
        args += ["--no-update", "--data-dir", str(data)]
        result = CliRunner().invoke(main, [*args, "--output", str(html)])
        assert result.exit_code == 0, result.output
        page = html.read_text()
        assert 'class="plaus plaus-not_assessed"' in page
        assert "Plausibility not assessed: gnomAD was not used in this run" in page
        assert '"plausibility_reason": "gnomAD was not used in this run"' in page
        result = CliRunner(env={"COLUMNS": "200"}).invoke(main, args)
        assert result.exit_code == 0, result.output
        assert (
            "plausibility not assessed; the JSON and HTML reports give the reason" in result.output
        )


def test_wrong_allele_checker_reads_half_and_hemizygous_calls() -> None:
    import sys

    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
    from check_wrong_allele import _carried

    assert _carried({"genotype_match": "T/."}) == ["T"]
    assert _carried({"genotype_match": "T"}) == ["T"]
    assert _carried({"genotype_match": "AG"}) == ["A", "G"]
    assert _carried({"genotype_match": "AT/A"}) == ["AT", "A"]

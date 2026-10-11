# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""VCF quality gate: spanning deletions, caller filters, genotype quality.

All records are synthetic. Positions are arbitrary and carry no real
genotype data.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from allelix.exporters.plink import export_plink
from allelix.models import Variant
from allelix.parsers.vcf import VcfParser

if TYPE_CHECKING:
    from pathlib import Path

_HEADER = (
    "##fileformat=VCFv4.3\n"
    "##contig=<ID=1,length=249250621,assembly=GRCh37>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
)


def _vcf(tmp_path: Path, rows: list[str], name: str = "in.vcf") -> Path:
    path = tmp_path / name
    path.write_text(_HEADER + "".join(r + "\n" for r in rows))
    return path


# ── Spanning deletion ``*`` (VCF 4.3) ─────────────────────────────


_STAR_ROWS = [
    # 1/2 against ALT=G,*: G on one haplotype, a spanning deletion on the other.
    "1\t1000\trs900750001\tA\tG,*\t50\tPASS\t.\tGT\t1/2",
    # 1/1 against ALT=*: both haplotypes deleted.
    "1\t2000\trs900750002\tC\t*\t50\tPASS\t.\tGT\t1/1",
    # 0/1 against ALT=T,*: the * is listed but not called.
    "1\t3000\trs900750003\tA\tT,*\t50\tPASS\t.\tGT\t0/1",
    # 0/2 against ALT=T,*: REF plus a spanning deletion.
    "1\t4000\trs900750004\tG\tT,*\t50\tPASS\t.\tGT\t0/2",
]


class TestSpanningDeletion:
    def test_star_never_reaches_a_genotype(self, tmp_path: Path) -> None:
        variants = list(VcfParser().parse(_vcf(tmp_path, _STAR_ROWS)))
        assert all("*" not in (v.allele1, v.allele2) for v in variants)

    def test_star_only_row_is_dropped(self, tmp_path: Path) -> None:
        rsids = {v.rsid for v in VcfParser().parse(_vcf(tmp_path, _STAR_ROWS))}
        assert "rs900750002" not in rsids

    def test_call_with_a_star_haplotype_is_dropped(self, tmp_path: Path) -> None:
        """Same as a call on <NON_REF>: no Variant shape says "one allele,
        other haplotype deleted", so nothing is claimed."""
        rsids = {v.rsid for v in VcfParser().parse(_vcf(tmp_path, _STAR_ROWS))}
        assert "rs900750001" not in rsids
        assert "rs900750004" not in rsids

    def test_uncalled_star_in_alt_list_is_harmless(self, tmp_path: Path) -> None:
        by_rsid = {v.rsid: v for v in VcfParser().parse(_vcf(tmp_path, _STAR_ROWS))}
        v = by_rsid["rs900750003"]
        assert (v.allele1, v.allele2) == ("A", "T")

    def test_plink_bim_never_carries_star(self, tmp_path: Path) -> None:
        prefix = tmp_path / "out"
        export_plink(VcfParser().parse(_vcf(tmp_path, _STAR_ROWS)), prefix, "GRCh37")
        bim_alleles = {
            allele
            for line in prefix.with_suffix(".bim").read_text().splitlines()
            for allele in line.split("\t")[4:6]
        }
        assert "*" not in bim_alleles
        assert bim_alleles <= {"A", "C", "G", "T", "0"}


class TestPlinkAlleleGuard:
    def test_non_nucleotide_single_char_alleles_are_skipped(self, tmp_path: Path) -> None:
        """Array I/D indel codes and a stray * are single characters but not
        nucleotides; they are counted with the indels, never written."""
        variants = [
            Variant("rs900750010", "1", 100, "A", "G"),
            Variant("rs900750011", "1", 200, "I", "D"),
            Variant("rs900750012", "1", 300, "A", "*"),
        ]
        prefix = tmp_path / "out"
        written, no_calls, non_snv, _mono = export_plink(iter(variants), prefix, "GRCh37")
        assert (written, no_calls, non_snv) == (1, 0, 2)
        assert prefix.with_suffix(".bim").read_text().split("\t")[1] == "rs900750010"


# ── Caller filters and genotype quality ───────────────────────────

_QUAL_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=3,length=198295559,assembly=GRCh38>\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=LowQual,Description="Low quality">\n'
    '##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="VQSR tranche">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
)

# The issue's fixture: four hets at ClinVar positions with PASS, ".",
# LowQual and a VQSR tranche, plus a PASS het with GQ 5.
_GATE_SITES = (
    ("rs900740001", 3_900_100, "G", "A", "PASS", "0/1:99:30"),
    ("rs900740002", 3_900_200, "C", "T", ".", "0/1:99:30"),
    ("rs900740003", 3_900_300, "T", "C", "LowQual", "0/1:99:30"),
    ("rs900740004", 3_900_400, "A", "G", "VQSRTrancheSNP99.90to100.00", "0/1:99:30"),
    ("rs900740005", 3_900_500, "G", "T", "PASS", "0/1:5:30"),
)


def _gate_vcf(tmp_path: Path, extra: tuple[str, ...] = ()) -> Path:
    path = tmp_path / "gate.vcf"
    rows = [
        f"3\t{pos}\t{rsid}\t{ref}\t{alt}\t50\t{filt}\t.\tGT:GQ:DP\t{sample}"
        for rsid, pos, ref, alt, filt, sample in _GATE_SITES
    ]
    path.write_text(_QUAL_HEADER + "\n".join([*rows, *extra]) + "\n")
    return path


def _calls(variants: list[Variant]) -> dict[str, tuple[str, str]]:
    return {v.rsid: (v.allele1, v.allele2) for v in variants}


class TestCallerFilterGate:
    def test_default_skips_filtered_rows_and_no_calls_low_gq(self, tmp_path: Path) -> None:
        parser = VcfParser()
        calls = _calls(list(parser.parse(_gate_vcf(tmp_path))))
        assert calls == {
            "rs900740001": ("G", "A"),
            "rs900740002": ("C", "T"),
            "rs900740005": ("-", "-"),
        }
        stats = parser.quality_stats
        assert stats.filtered_rows == 2
        assert dict(stats.filtered_by_label) == {
            "LowQual": 1,
            "VQSRTrancheSNP99.90to100.00": 1,
        }
        assert stats.low_quality_calls == 1
        assert not stats.kept_filtered

    def test_keep_filtered_keeps_and_counts(self, tmp_path: Path) -> None:
        parser = VcfParser(keep_filtered=True)
        calls = _calls(list(parser.parse(_gate_vcf(tmp_path))))
        assert calls["rs900740003"] == ("T", "C")
        assert calls["rs900740004"] == ("A", "G")
        assert parser.quality_stats.filtered_rows == 2
        assert parser.quality_stats.kept_filtered

    def test_min_gq_is_tunable_and_zero_disables(self, tmp_path: Path) -> None:
        assert _calls(list(VcfParser(min_gq=0).parse(_gate_vcf(tmp_path))))["rs900740005"] == (
            "G",
            "T",
        )
        strict = _calls(list(VcfParser(min_gq=100).parse(_gate_vcf(tmp_path))))
        assert strict["rs900740001"] == ("-", "-")

    def test_zero_depth_is_a_no_call_even_without_a_gq_floor(self, tmp_path: Path) -> None:
        row = "3\t3900600\trs900740006\tC\tA\t50\tPASS\t.\tGT:GQ:DP\t1/1:99:0"
        parser = VcfParser(min_gq=0)
        calls = _calls(list(parser.parse(_gate_vcf(tmp_path, (row,)))))
        assert calls["rs900740006"] == ("-", "-")
        assert parser.quality_stats.low_quality_calls == 1

    def test_missing_gq_or_dp_does_not_gate(self, tmp_path: Path) -> None:
        rows = (
            "3\t3900700\trs900740007\tC\tA\t50\tPASS\t.\tGT:GQ:DP\t0/1:.:.",
            "3\t3900800\trs900740008\tC\tA\t50\tPASS\t.\tGT\t0/1",
        )
        calls = _calls(list(VcfParser().parse(_gate_vcf(tmp_path, rows))))
        assert calls["rs900740007"] == ("C", "A")
        assert calls["rs900740008"] == ("C", "A")

    def test_multiple_labels_count_once_each(self, tmp_path: Path) -> None:
        row = "3\t3900900\trs900740009\tC\tA\t50\tLowQual;q10\t.\tGT:GQ:DP\t0/1:99:30"
        parser = VcfParser()
        list(parser.parse(_gate_vcf(tmp_path, (row,))))
        assert parser.quality_stats.filtered_rows == 3
        assert parser.quality_stats.filtered_by_label["LowQual"] == 2
        assert parser.quality_stats.filtered_by_label["q10"] == 1

    def test_stats_reset_between_runs(self, tmp_path: Path) -> None:
        parser = VcfParser()
        path = _gate_vcf(tmp_path)
        list(parser.parse(path))
        list(parser.parse(path))
        assert parser.quality_stats.filtered_rows == 2


def _gate_data_dir(tmp_path: Path) -> Path:
    """Synthetic ClinVar caches with a pathogenic row at every gate site."""
    import contextlib
    import sqlite3

    from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
    from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
    from allelix.databases.schema import CLINVAR_SCHEMA

    data = tmp_path / "data"
    data.mkdir()
    rows = [
        (rsid, "3", pos, ref, alt, "Pathogenic", f"Cond_{rsid}", "GENE",
         "criteria_provided,_single_submitter", i)
        for i, (rsid, pos, ref, alt, _filt, _sample) in enumerate(_GATE_SITES, 1)
    ]  # fmt: skip
    for build in ("GRCh37", "GRCh38"):
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
                "VALUES (?, 'test://gate', '1', '2026-10-10', ?, NULL, ?)",
                (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    return data


def _analyze(tmp_path: Path, *extra: str):
    import json

    from click.testing import CliRunner

    from allelix.cli import main

    out_json = tmp_path / "report.json"
    out_vcf = tmp_path / "annotated.vcf"
    result = CliRunner(env={"COLUMNS": "200"}).invoke(
        main,
        [
            "analyze",
            str(_gate_vcf(tmp_path)),
            "--data-dir",
            str(_gate_data_dir(tmp_path)),
            "--build",
            "grch38",
            "--no-update",
            "--output",
            str(out_json),
            "--vcf-out",
            str(out_vcf),
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output
    return result.output, json.loads(out_json.read_text()), out_vcf.read_text()


def _stamped_rsids(vcf_text: str) -> set[str]:
    return {
        line.split("\t")[2]
        for line in vcf_text.splitlines()
        if not line.startswith("#") and "ALLELIX_CLINVAR=" in line
    }


class TestCallerFilterEndToEnd:
    def test_default_reports_two_and_counts_three(self, tmp_path: Path) -> None:
        output, report, vcf_text = _analyze(tmp_path)
        assert {a["rsid"] for a in report["annotations"]} == {"rs900740001", "rs900740002"}
        assert report["input"]["vcf_quality"] == {
            "filtered_rows": 2,
            "filtered_by_label": {"LowQual": 1, "VQSRTrancheSNP99.90to100.00": 1},
            "kept_filtered": False,
            "low_quality_calls": 1,
            "min_gq": 20,
            "spanning_deletion_calls": 0,
            "symbolic_rows": 0,
        }
        assert _stamped_rsids(vcf_text) == {"rs900740001", "rs900740002"}
        assert "ALLELIX_FILTER" not in vcf_text
        assert "2 row(s) excluded by the variant caller's filters" in output
        assert "1 call(s) reported as no-calls" in output

    def test_filtered_rows_pass_through_the_annotated_vcf_unchanged(self, tmp_path: Path) -> None:
        _output, _report, vcf_text = _analyze(tmp_path)
        source = _gate_vcf(tmp_path).read_text().splitlines()
        for rsid in ("rs900740003", "rs900740004"):
            original = next(ln for ln in source if f"\t{rsid}\t" in ln)
            assert original in vcf_text.splitlines()

    def test_keep_filtered_reports_and_marks_them(self, tmp_path: Path) -> None:
        output, report, vcf_text = _analyze(tmp_path, "--keep-filtered")
        assert {a["rsid"] for a in report["annotations"]} == {
            "rs900740001",
            "rs900740002",
            "rs900740003",
            "rs900740004",
        }
        assert report["input"]["vcf_quality"]["kept_filtered"] is True
        assert _stamped_rsids(vcf_text) == {
            "rs900740001",
            "rs900740002",
            "rs900740003",
            "rs900740004",
        }
        lowqual = next(ln for ln in vcf_text.splitlines() if "\trs900740003\t" in ln)
        assert "ALLELIX_FILTER=LowQual" in lowqual
        passing = next(ln for ln in vcf_text.splitlines() if "\trs900740001\t" in ln)
        assert "ALLELIX_FILTER" not in passing
        assert "##INFO=<ID=ALLELIX_FILTER," in vcf_text
        assert "kept because of --keep-filtered" in output

    def test_min_gq_flag(self, tmp_path: Path) -> None:
        _output, report, _vcf = _analyze(tmp_path, "--min-gq", "0")
        assert "rs900740005" in {a["rsid"] for a in report["annotations"]}
        assert report["input"]["vcf_quality"]["min_gq"] == 0

    def test_html_input_quality_banner(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from allelix.cli import main

        out = tmp_path / "report.html"
        result = CliRunner().invoke(
            main,
            [
                "analyze",
                str(_gate_vcf(tmp_path)),
                "--data-dir",
                str(_gate_data_dir(tmp_path)),
                "--build",
                "grch38",
                "--no-update",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        body = out.read_text()
        assert "<strong>Input quality.</strong>" in body
        assert "excluded by the variant caller&#x27;s filters" in body or (
            "excluded by the variant caller's filters" in body
        )

    def test_non_vcf_input_has_null_vcf_quality(self, tmp_path: Path) -> None:
        from allelix.reports._pipeline import AnalysisResult
        from allelix.reports.json_report import render_json

        out = tmp_path / "r.json"
        render_json(
            AnalysisResult(
                file_path=tmp_path / "chip.txt",
                parser_name="myhappygenes",
                parser_display_name="MyHappyGenes",
                sample_id="S",
                build="GRCh37",
                total_variants=0,
                skipped_count=0,
                annotators_used=[],
            ),
            output_path=out,
        )
        import json

        assert json.loads(out.read_text())["input"]["vcf_quality"] is None


# ── Reviewer additions: GQ boundary, every command honours the gate ──


class TestGateBoundaryAndCommands:
    def test_gq_19_is_a_no_call_and_gq_20_passes(self, tmp_path: Path) -> None:
        rows = (
            "3\t3901000\trs900740010\tC\tA\t50\tPASS\t.\tGT:GQ:DP\t0/1:19:30",
            "3\t3901100\trs900740011\tC\tA\t50\tPASS\t.\tGT:GQ:DP\t0/1:20:30",
        )
        calls = _calls(list(VcfParser().parse(_gate_vcf(tmp_path, rows))))
        assert calls["rs900740010"] == ("-", "-")
        assert calls["rs900740011"] == ("C", "A")

    def test_export_plink_honours_the_gate(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from allelix.cli import main

        prefix = tmp_path / "out"
        result = CliRunner().invoke(
            main,
            [
                "export",
                "plink",
                str(_gate_vcf(tmp_path)),
                "-o",
                str(prefix),
                "--build",
                "grch38",
                "--data-dir",
                str(tmp_path / "empty"),
            ],
        )
        assert result.exit_code == 0, result.output
        bim_ids = {ln.split("\t")[1] for ln in prefix.with_suffix(".bim").read_text().splitlines()}
        assert bim_ids == {"rs900740001", "rs900740002"}

    def test_stats_honours_the_gate(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from allelix.cli import main

        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main, ["stats", str(_gate_vcf(tmp_path))]
        )
        assert result.exit_code == 0, result.output
        # Two calls plus the GQ-5 no-call; the two filtered rows are not counted.
        chrom_row = next(
            ln for ln in result.output.splitlines() if ln.startswith("│ 3 ") and ln.count("│") == 3
        )
        assert chrom_row.split("│")[2].strip() == "3"


class TestSpanningDeletionReporting:
    def test_one_haplotype_star_calls_are_counted(self, tmp_path: Path) -> None:
        from allelix.reports import vcf_quality_lines

        parser = VcfParser()
        list(parser.parse(_vcf(tmp_path, _STAR_ROWS)))
        # 1/2 G,* and 0/2 T,* have '*' on one haplotype; 1/1 * has it on both.
        assert parser.quality_stats.spanning_deletion_calls == 2
        assert (
            "2 record(s) with a spanning deletion on one haplotype were not assessed."
            in vcf_quality_lines(parser.quality_stats)
        )

    def test_star_never_reaches_report_or_annotated_vcf_fields(self, tmp_path: Path) -> None:
        """A '*' listed in ALT but not called: the G call is annotated, the
        '*' slot of every Allelix INFO field stays empty."""
        import json

        from click.testing import CliRunner

        from allelix.cli import main

        data = _gate_data_dir(tmp_path)
        vcf = tmp_path / "star.vcf"
        vcf.write_text(
            _QUAL_HEADER + "3\t3900100\trs900740001\tG\tA,*\t50\tPASS\t.\tGT:GQ:DP\t0/1:99:30\n"
        )
        out_json, out_vcf = tmp_path / "r.json", tmp_path / "a.vcf"
        result = CliRunner().invoke(
            main,
            [
                "analyze",
                str(vcf),
                "--data-dir",
                str(data),
                "--build",
                "grch38",
                "--no-update",
                "--output",
                str(out_json),
                "--vcf-out",
                str(out_vcf),
            ],
        )
        assert result.exit_code == 0, result.output
        report = json.loads(out_json.read_text())
        assert [a["genotype_match"] for a in report["annotations"]] == ["AG"]
        row = next(ln for ln in out_vcf.read_text().splitlines() if "\trs900740001\t" in ln)
        allelix_fields = [f for f in row.split("\t")[7].split(";") if f.startswith("ALLELIX_")]
        assert allelix_fields
        for entry in allelix_fields:
            values = entry.split("=", 1)[1].split(",")
            assert values[1] == "."
            assert "*" not in values[0]


# ── gVCF no-call semantics and high-value positions ────────

_GVCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##ALT=<ID=NON_REF,Description="Any other allele">\n'
    "##contig=<ID=chr22,length=50818468,assembly=GRCh38>\n"
    '##FILTER=<ID=LowQual,Description="Low quality">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
)

_GVCF_ROWS = (
    # GQ 0 / DP 0 block, the unambiguous case.
    "chr22\t1000\t.\tT\t<NON_REF>\t.\t.\tEND=1100\tGT:DP:GQ\t0/0:0:0",
    # Confident reference block.
    "chr22\t2000\t.\tA\t<NON_REF>\t.\t.\tEND=2100\tGT:DP:GQ\t0/0:30:50",
    # Boundary: GQ 19 demoted, GQ 20 kept.
    "chr22\t3000\t.\tC\t<NON_REF>\t.\t.\tEND=3100\tGT:DP:GQ\t0/0:10:19",
    "chr22\t4000\t.\tG\t<NON_REF>\t.\t.\tEND=4100\tGT:DP:GQ\t0/0:10:20",
    # A block without GQ is not demoted.
    "chr22\t5000\t.\tT\t<NON_REF>\t.\t.\tEND=5100\tGT:DP\t0/0:12",
    # Variant records: low GQ, caller no-call, confident call, filtered.
    "chr22\t6000\t.\tC\tT,<NON_REF>\t50\t.\t.\tGT:DP:GQ\t0/1:20:5",
    "chr22\t7000\t.\tG\tA,<NON_REF>\t.\t.\t.\tGT:DP:GQ\t./.:0:0",
    "chr22\t8000\t.\tA\tG,<NON_REF>\t50\t.\t.\tGT:DP:GQ\t0/1:30:99",
    "chr22\t9000\t.\tC\tA,<NON_REF>\t50\tLowQual\t.\tGT:DP:GQ\t0/1:30:99",
)

# Synthetic high-value rsIDs placed by the identity table (GRCh38 position).
_HV_SITES = {
    "rs900720001": 1050,  # inside the GQ 0 / DP 0 block
    "rs900720002": 1000,  # first base of that block
    "rs900720003": 1100,  # END of that block
    "rs900720004": 1101,  # one past it: no record covers it
    "rs900720005": 2050,  # confident reference
    "rs900720006": 3050,  # GQ 19 block
    "rs900720007": 4050,  # GQ 20 block
    "rs900720008": 5050,  # block without GQ
    "rs900720009": 6000,  # variant call with GQ 5
    "rs900720010": 7000,  # caller no-call
    "rs900720011": 8000,  # confident variant call
    "rs900720012": 9000,  # filtered row
    "rs900720013": 9500,  # no record
}

_EXPECTED_GVCF = {
    "rs900720001": ("no_call", "low_gq_block"),
    "rs900720002": ("no_call", "low_gq_block"),
    "rs900720003": ("no_call", "low_gq_block"),
    "rs900720004": ("no_call", "no_coverage"),
    "rs900720005": ("reference", ""),
    "rs900720006": ("no_call", "low_gq_block"),
    "rs900720007": ("reference", ""),
    "rs900720008": ("reference", ""),
    "rs900720009": ("no_call", "low_gq_call"),
    "rs900720010": ("no_call", "no_call"),
    "rs900720011": ("called", ""),
    "rs900720012": ("no_call", "filtered"),
    "rs900720013": ("no_call", "no_coverage"),
}


def _hv_data_dir(tmp_path: Path) -> Path:
    from allelix.databases.rsid_identity_loader import (
        RSID_IDENTITY_DB_FILENAME,
        write_identity_db,
    )

    data = tmp_path / "hvdata"
    data.mkdir(exist_ok=True)
    write_identity_db(
        [
            (rsid, "22", pos - 50_000, "N", "N", pos, "N", "N", 0)
            for rsid, pos in _HV_SITES.items()
        ],
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://hv",
        downloaded_at="2026-10-10T00:00:00Z",
    )
    return data


def _hv_states(tmp_path: Path, vcf_text: str, cm_stack) -> dict[str, tuple[str, str]] | None:
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.reports._pipeline import run_analysis

    vcf = tmp_path / "hv.vcf"
    vcf.write_text(vcf_text)
    result = run_analysis(
        vcf,
        VcfParser(),
        [],
        build_override="GRCh38",
        identity=cm_stack.enter_context(RsidIdentityAnnotator(_hv_data_dir(tmp_path))),
        high_value_rsids=set(_HV_SITES),
    )
    return result.high_value_vcf_states


class TestGvcfHighValueStates:
    def test_every_state_and_reason(self, tmp_path: Path, cm_stack) -> None:
        states = _hv_states(tmp_path, _GVCF_HEADER + "\n".join(_GVCF_ROWS) + "\n", cm_stack)
        assert states == _EXPECTED_GVCF

    def test_plain_vcf_leaves_uncovered_positions_unknown(self, tmp_path: Path, cm_stack) -> None:
        plain_header = _GVCF_HEADER.replace(
            '##ALT=<ID=NON_REF,Description="Any other allele">\n', ""
        )
        variant_rows = [r for r in _GVCF_ROWS if "END=" not in r]
        states = _hv_states(tmp_path, plain_header + "\n".join(variant_rows) + "\n", cm_stack)
        assert states == {
            "rs900720009": ("no_call", "low_gq_call"),
            "rs900720010": ("no_call", "no_call"),
            "rs900720011": ("called", ""),
            "rs900720012": ("no_call", "filtered"),
        }

    def test_scan_lists_no_calls_with_reasons(self, tmp_path: Path, cm_stack) -> None:
        from allelix.reports.high_value import HighValueSNP, format_warnings, scan_no_calls

        states = _hv_states(tmp_path, _GVCF_HEADER + "\n".join(_GVCF_ROWS) + "\n", cm_stack)
        high_value = {
            rsid: HighValueSNP(rsid=rsid, gene="GENE", cluster="", note="note")
            for rsid in _HV_SITES
        }
        warnings = scan_no_calls([], high_value, states)
        assert {w.snp.rsid: w.reason for w in warnings} == {
            rsid: reason for rsid, (state, reason) in _EXPECTED_GVCF.items() if state == "no_call"
        }
        line = next(ln for ln in format_warnings(warnings) if ln.startswith("rs900720001"))
        assert "[reference call with low genotype quality or no reads]" in line

    def test_a_variant_call_beats_a_no_call(self) -> None:
        from allelix.reports.high_value import HighValueSNP, scan_no_calls

        hv = {"rs900720011": HighValueSNP("rs900720011", "GENE", "", "note")}
        called = Variant("rs900720011", "22", 8000, "A", "G")
        assert scan_no_calls([called], hv, {"rs900720011": ("no_call", "low_gq_block")}) == []

    def test_cli_json_reports_gvcf_and_reasons(self, tmp_path: Path, monkeypatch) -> None:
        import json

        from click.testing import CliRunner

        from allelix.cli import main
        from allelix.reports.high_value import HighValueSNP

        monkeypatch.setattr(
            "allelix.cli._helpers.load_high_value_snps",
            lambda: {r: HighValueSNP(r, "GENE", "", "note") for r in _HV_SITES},
        )
        vcf = tmp_path / "hv.g.vcf"
        vcf.write_text(_GVCF_HEADER + "\n".join(_GVCF_ROWS) + "\n")
        data = _hv_data_dir(tmp_path)
        out = tmp_path / "r.json"
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main,
            [
                "analyze",
                str(vcf),
                "--data-dir",
                str(data),
                "--build",
                "grch38",
                "--no-update",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        report = json.loads(out.read_text())
        assert report["input"]["gvcf"] is True
        got = {e["rsid"]: e["reason"] for e in report["high_value_no_calls"]}
        assert got["rs900720001"] == "low_gq_block"
        assert got["rs900720004"] == "no_coverage"
        assert "rs900720005" not in got


class TestLocusWatcherAgreesWithBruteForce:
    """The streaming watcher must classify watched positions exactly as a
    brute-force scan of every record does: sorted, shuffled (fallback),
    and with a chromosome that comes back."""

    @staticmethod
    def _random_gvcf(seed: int) -> list[str]:
        import random

        rng = random.Random(seed)
        rows: list[str] = []
        for chrom in ("chr21", "chr22"):
            pos = 1
            while pos < 3000:
                if rng.random() < 0.25:
                    gq = rng.choice([0, 5, 19, 20, 60, 99])
                    gt = rng.choice(["0/1", "1/1", "./.", "0/0"])
                    filt = rng.choice([".", ".", "PASS", "LowQual"])
                    ref = rng.choice(["A", "CA", "GTT"])
                    rows.append(
                        f"{chrom}\t{pos}\t.\t{ref}\tT,<NON_REF>\t50\t{filt}\t.\t"
                        f"GT:DP:GQ\t{gt}:{rng.choice([0, 12])}:{gq}"
                    )
                    # Sometimes the next record starts inside this one's REF,
                    # as gVCFs do after a deletion.
                    pos += 1 if rng.random() < 0.3 else len(ref)
                else:
                    length = rng.randint(1, 120)
                    gq = rng.choice([0, 10, 19, 20, 45, 99])
                    rows.append(
                        f"{chrom}\t{pos}\t.\tA\t<NON_REF>\t.\t.\tEND={pos + length - 1}\t"
                        f"GT:DP:GQ\t0/0:{rng.choice([0, 25])}:{gq}"
                    )
                    pos += length
                if rng.random() < 0.05:
                    pos += rng.randint(1, 40)  # an uncovered gap
        return rows

    @staticmethod
    def _brute_force(rows: list[str], loci: set[tuple[str, int]]) -> dict:
        from allelix.parsers._helpers import normalize_chromosome
        from allelix.parsers.vcf import (
            VcfQualityStats,
            _parse_record,
            _row_stop,
            strongest_state,
        )

        out: dict = {}
        for row in rows:
            cols = row.split("\t")
            state = _parse_record(cols, row, 0, VcfQualityStats(), want_state=True)[1]
            if state is None:
                continue
            chrom, start = normalize_chromosome(cols[0]), int(cols[1])
            stop = _row_stop(cols, start)
            for c, p in loci:
                if c == chrom and start <= p <= stop:
                    out[(c, p)] = strongest_state(out.get((c, p)), state)
        return out

    @pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
    @pytest.mark.parametrize("layout", ["sorted", "shuffled", "chromosome_returns"])
    def test_matches(self, tmp_path: Path, seed: int, layout: str) -> None:
        import random

        rows = self._random_gvcf(seed)
        if layout == "shuffled":
            random.Random(seed).shuffle(rows)
        elif layout == "chromosome_returns":
            first = [r for r in rows if r.startswith("chr21\t")]
            half = len(first) // 2
            rows = first[:half] + [r for r in rows if r.startswith("chr22\t")] + first[half:]
        rng = random.Random(seed * 7)
        loci = {(c, rng.randint(1, 3100)) for c in ("21", "22") for _ in range(60)}
        vcf = tmp_path / "w.g.vcf"
        vcf.write_text(_GVCF_HEADER + "\n".join(rows) + "\n")
        parser = VcfParser()
        parser.watch(loci, set())
        list(parser.parse(vcf))
        assert parser.quality_stats.locus_states == self._brute_force(rows, loci)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ADR-0047: VCF calls that used to be dropped without a word.

Untrimmed alleles, half calls (``./1``), homozygous deletions to one base,
and spanning-deletion sites. All records and cache rows are synthetic;
positions are arbitrary and carry no real genotype data.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.models import HALF_CALL, Annotation, Variant
from allelix.parsers.vcf import VcfParser
from allelix.reports import vcf_quality_lines
from allelix.reports._pipeline import _flag_plausibility
from allelix.utils.allele import trim_allele_pair
from allelix.utils.plausibility import DEFAULT_THRESHOLDS

if TYPE_CHECKING:
    from pathlib import Path

BUILD = "GRCh38"

# (position, REF, ALT, classification): the ClinVar rows, written the way
# ClinVar writes them (trimmed, one anchor base for indels).
_ROWS = (
    (1_000, "GT", "G", "Pathogenic"),  # BRCA2-style one-base deletion
    (2_000, "C", "T", "Pathogenic"),  # BRCA1-style SNV
    (3_001, "C", "T", "Pathogenic"),  # reached only after a prefix trim
    (4_000, "A", "ATG", "Pathogenic"),  # TG insertion of a compound call
    (4_000, "ATG", "A", "Likely pathogenic"),  # TG deletion of the same call
    (5_000, "AC", "A", "Pathogenic"),  # homozygous deletion to one base
    (6_000, "C", "T", "Pathogenic"),  # half calls
)

# (position, REF, ALT column, GT): the VCF records, written the way a caller
# or a merge can write them.
_RECORDS = (
    (1_000, "GTT", "GT,G", "0/1"),
    (2_000, "CA", "TA", "0/1"),
    (3_000, "ACG", "ATG", "0/1"),
    (4_000, "ATG", "ATGTG,A", "2/1"),
    (5_000, "AC", "A", "1/1"),
    (6_000, "C", "T", "./1"),
    (6_000, "C", "T", "1/."),
    (7_000, "A", "G,*", "1/2"),
    (8_000, "C", "*", "1/1"),
)

_HEADER = (
    "##fileformat=VCFv4.3\n"
    "##contig=<ID=1,length=248956422,assembly=GRCh38>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
)


def _vcf(tmp_path: Path, records=_RECORDS) -> Path:
    path = tmp_path / "in.vcf"
    path.write_text(
        _HEADER
        + "".join(
            f"1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n" for pos, ref, alt, gt in records
        )
    )
    return path


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    rows = [
        (f"rs90076{i:04d}", "1", pos, ref, alt, sig, f"Condition {i}", f"GENE{i}",
         "criteria provided, single submitter", i)
        for i, (pos, ref, alt, sig) in enumerate(_ROWS, 1)
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
                "VALUES (?, 'test://drops', '1', '2026-10-10', ?, NULL, ?)",
                (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    return data


def _parse(tmp_path: Path, records=_RECORDS) -> tuple[list[Variant], VcfParser]:
    parser = VcfParser()
    variants = list(parser.parse(_vcf(tmp_path, records)))
    for v in variants:
        v.build = BUILD
    return variants, parser


def _annotate(data_dir: Path, variants: list[Variant]) -> list[Annotation]:
    with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
        return list(clinvar.batch_annotate(variants))


def _at(tmp_path: Path, data_dir: Path, pos: int) -> list[Annotation]:
    variants, _parser = _parse(tmp_path, [r for r in _RECORDS if r[0] == pos])
    return _annotate(data_dir, variants)


class TestTrimAllelePair:
    @pytest.mark.parametrize(
        ("pair", "trimmed"),
        [
            ((1_000, "GTT", "GT"), (1_000, "GT", "G")),
            ((2_000, "CA", "TA"), (2_000, "C", "T")),
            ((3_000, "ACG", "ATG"), (3_001, "C", "T")),
            ((4_000, "ATG", "ATGTG"), (4_000, "A", "ATG")),
            ((5_000, "AC", "A"), (5_000, "AC", "A")),
            ((6_000, "C", "T"), (6_000, "C", "T")),
        ],
    )
    def test_trims_suffix_then_prefix_keeping_an_anchor(self, pair, trimmed) -> None:
        assert trim_allele_pair(*pair) == trimmed


class TestUntrimmedAlleles:
    def test_deletion_written_with_a_shared_suffix_matches(self, tmp_path, data_dir) -> None:
        [a] = _at(tmp_path, data_dir, 1_000)
        assert (a.locus, a.locus_ref, a.alt) == (("1", 1_000), "GT", "G")
        assert a.genotype_match == "GTT/GT"
        assert a.zygosity == "Heterozygous"

    def test_snv_written_with_a_shared_suffix_matches(self, tmp_path, data_dir) -> None:
        [a] = _at(tmp_path, data_dir, 2_000)
        assert (a.locus, a.locus_ref, a.alt) == (("1", 2_000), "C", "T")

    def test_prefix_trim_moves_the_match_one_base_on(self, tmp_path, data_dir) -> None:
        [a] = _at(tmp_path, data_dir, 3_000)
        assert (a.locus, a.locus_ref, a.alt) == (("1", 3_001), "C", "T")
        assert a.record_key == ("1", 3_000, "ACG", "ATG")

    def test_each_allele_of_a_compound_call_is_trimmed_on_its_own(
        self, tmp_path, data_dir
    ) -> None:
        found = {(a.locus_ref, a.alt) for a in _at(tmp_path, data_dir, 4_000)}
        assert found == {("A", "ATG"), ("ATG", "A")}

    def test_two_variants_one_call_are_not_collapsed(self, tmp_path, data_dir) -> None:
        """TTCTC>TTC,T 2/1 carries two deletions: TTC>T (trimmed) and
        TTCTC>T. ClinVar can give both one rsID, one classification and one
        condition; they are still two variants, each with its own row."""
        for build in ("GRCh37", "GRCh38"):
            with contextlib.closing(
                sqlite3.connect(data_dir / clinvar_db_filename(build))
            ) as conn:
                conn.executemany(
                    "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                    "clinical_significance, condition, gene, review_status, allele_id, "
                    "submitters_same, submitters_total) VALUES "
                    "('rs900769000', '1', 9000, ?, 'T', 'Pathogenic', 'Condition 9', "
                    "'GENE9', 'criteria provided, single submitter', ?, ?, ?)",
                    [("TTC", 90, 1, 2), ("TTCTC", 91, 2, 3)],
                )
                conn.commit()
        variants, _parser = _parse(tmp_path, [(9_000, "TTCTC", "TTC,T", "2/1")])
        found = _annotate(data_dir, variants)
        assert sorted(
            (a.locus_ref, a.alt, a.submitters_same, a.submitters_total) for a in found
        ) == [
            ("TTC", "T", 1, 2),
            ("TTCTC", "T", 2, 3),
        ]

    def test_array_call_is_never_trimmed(self, data_dir) -> None:
        """A REF from a reference table is not the record's: no trimming."""
        call = Variant("rs900760002", "1", 2_000, "CA", "TA", build=BUILD, ref="CA")
        assert _annotate(data_dir, [call]) == []

    def test_count_of_calls_matched_after_trimming(self, tmp_path, data_dir) -> None:
        variants, _parser = _parse(tmp_path)
        with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
            list(clinvar.batch_annotate(variants))
            assert clinvar.run_counts()["trimmed_matches"] == 4


class TestHomozygousDeletionToOneBase:
    def test_vcf_call_matches(self, tmp_path, data_dir) -> None:
        [a] = _at(tmp_path, data_dir, 5_000)
        assert (a.locus_ref, a.alt, a.zygosity) == ("AC", "A", "Homozygous")

    def test_array_call_keeps_the_anchor_rule(self, data_dir) -> None:
        """ADR-0011: an array's A/A must not match the AC>A deletion, even
        with a REF filled in from the identity table."""
        call = Variant("rs900760006", "1", 5_000, "A", "A", build=BUILD, ref="AC")
        assert _annotate(data_dir, [call]) == []


class TestHalfCalls:
    def test_parser_yields_the_called_allele(self, tmp_path) -> None:
        variants, parser = _parse(tmp_path)
        halves = [v for v in variants if v.position == 6_000]
        assert [(v.allele1, v.allele2) for v in halves] == [("-", "T"), ("T", "-")]
        assert all(v.is_half_call for v in halves)
        assert parser.quality_stats.half_calls == 2

    def test_clinvar_assesses_the_called_allele(self, tmp_path, data_dir) -> None:
        found = _at(tmp_path, data_dir, 6_000)
        assert [(a.alt, a.genotype_match, a.zygosity) for a in found] == [
            ("T", "T/.", HALF_CALL),
            ("T", "T/.", HALF_CALL),
        ]

    def test_array_half_call_stays_a_no_call(self, data_dir) -> None:
        call = Variant("rs900760007", "1", 6_000, "-", "T", build=BUILD, ref="C")
        assert not call.is_half_call
        assert _annotate(data_dir, [call]) == []

    def test_low_quality_half_call_is_a_no_call(self, tmp_path) -> None:
        path = tmp_path / "gq.vcf"
        path.write_text(_HEADER + "1\t6000\t.\tC\tT\t50\tPASS\t.\tGT:GQ\t./1:5\n")
        parser = VcfParser()
        [v] = list(parser.parse(path))
        assert (v.allele1, v.allele2) == ("-", "-")
        assert (parser.quality_stats.half_calls, parser.quality_stats.low_quality_calls) == (0, 1)

    def test_gwas_assesses_the_called_allele(self, gwas_data_dir: Path) -> None:
        call = Variant(
            "rs1801133", "1", 11_796_321, "-", "A", build=BUILD, ref="G", ref_from_input=True
        )
        with GWASCatalogAnnotator(gwas_data_dir, filter_traits=False) as gwas:
            found = gwas.annotate(call)
        assert found
        assert {a.genotype_match for a in found} == {"A/."}

    def test_clinvar_resolver_names_a_half_call(self, tmp_path, data_dir) -> None:
        variants, _parser = _parse(tmp_path, [r for r in _RECORDS if r[0] == 6_000])
        with ClinVarAnnotator(data_dir, builds=(BUILD,)) as clinvar:
            clinvar.bulk_resolve_rsids(variants)
        assert {v.rsid for v in variants} == {"rs900760007"}

    def test_identity_resolver_names_a_half_call(self) -> None:
        from allelix.annotators._resolution import _resolve_from_identity

        class _Identity:
            def alleles_at_loci(self, build, loci):
                return {("1", 6_000): [("C", "T", "rs900760007")]}

        call = Variant("", "1", 6_000, "-", "T", build=BUILD, ref="C", ref_from_input=True)
        coords: dict = {}
        _resolve_from_identity(_Identity(), [call], coords)  # type: ignore[arg-type]
        assert call.rsid == "rs900760007"
        assert coords == {"rs900760007": ("1", 6_000, "C", "T")}

    def test_plausibility_is_not_assessed(self) -> None:
        a = Annotation(
            source="clinvar",
            rsid="rs900760007",
            significance="clinvar_pathogenic",
            category="clinical",
            magnitude=9.0,
            description="",
            attribution="ClinVar",
            genotype_match="T/.",
            alt="T",
            locus_ref="C",
            carried_alleles=("-", "T"),
            locus=("1", 6_000),
            allele_frequency=0.9,
        )
        _flag_plausibility([a], [(a, ("1", 6_000, "C", "T"))], DEFAULT_THRESHOLDS, gnomad_ran=True)
        assert (a.plausibility, a.plausibility_reason) == ("not_assessed", "half_call")


class TestDropsAreCounted:
    def test_spanning_deletion_sites_are_listed(self, tmp_path) -> None:
        _variants, parser = _parse(tmp_path)
        stats = parser.quality_stats
        assert stats.spanning_deletion_calls == 2
        assert stats.spanning_deletion_sites == ["1:7000 A>G,* 1/2", "1:8000 C>* 1/1"]

    def test_rows_without_a_usable_genotype_are_counted(self, tmp_path) -> None:
        path = tmp_path / "odd.vcf"
        path.write_text(
            _HEADER
            + "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/3\n"
            + "1\t200\t.\tA\tG\t50\tPASS\t.\tDP\t30\n"
            + "1\tx\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "1\t300\t.\tA\n"
        )
        parser = VcfParser()
        assert list(parser.parse(path)) == []
        stats = parser.quality_stats
        assert (stats.unusable_genotype_calls, stats.no_genotype_rows, stats.malformed_rows) == (
            1,
            1,
            2,
        )

    @pytest.mark.parametrize(
        ("row", "field"),
        [
            ("1\t400\t.\tA\tG,<NON_REF>\t50\tPASS\t.\tGT\t0/2", "symbolic_allele_calls"),
            ("1\t500\t.\tA\tG,<*>\t50\tPASS\t.\tGT\t1/2", "symbolic_allele_calls"),
            ("1\t600\t.\tG\t.\t50\tPASS\t.\tGT\t0/1", "unusable_genotype_calls"),
            ("1\t610\t.\tG\t.\t50\tPASS\t.\tGT\t1", "unusable_genotype_calls"),
        ],
    )
    def test_calls_on_alleles_that_name_no_sequence_are_counted(
        self, tmp_path, row, field
    ) -> None:
        path = tmp_path / "odd.vcf"
        path.write_text(_HEADER + row + "\n")
        parser = VcfParser()
        assert list(parser.parse(path)) == []
        assert getattr(parser.quality_stats, field) == 1
        assert parser.quality_stats.symbolic_rows == 0

    def test_a_watched_record_is_counted_once(self, tmp_path) -> None:
        """A record at a watched rsID takes the full path, not the early exit."""
        path = tmp_path / "watched.vcf"
        path.write_text(_HEADER + "1\t620\trs900760099\tG\t.\t50\tPASS\t.\tGT\t0/1\n")
        parser = VcfParser()
        parser.watch(set(), {"rs900760099"})
        list(parser.parse(path))
        assert parser.quality_stats.unusable_genotype_calls == 1

    def test_no_calls_reference_blocks_and_svs_are_not_counted_as_drops(self, tmp_path) -> None:
        path = tmp_path / "plain.vcf"
        path.write_text(
            _HEADER
            + "1\t700\t.\tA\tG\t50\tPASS\t.\tGT\t./.\n"
            + "1\t800\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/0\n"
            + "1\t900\t.\tG\t.\t50\tPASS\t.\tGT\t0/0\n"
            + "1\t950\t.\tA\t<DEL>\t50\tPASS\t.\tGT\t0/1\n"
        )
        parser = VcfParser()
        list(parser.parse(path))
        stats = parser.quality_stats
        assert (stats.symbolic_allele_calls, stats.unusable_genotype_calls) == (0, 0)
        assert stats.symbolic_rows == 1

    def test_console_lines_name_each_reason(self, tmp_path) -> None:
        _variants, parser = _parse(tmp_path)
        lines = vcf_quality_lines(parser.quality_stats)
        assert any(line.startswith("2 call(s) with a spanning deletion (*)") for line in lines)
        assert any(line.startswith("2 half call(s)") for line in lines)


def _analyze(tmp_path: Path, data_dir: Path) -> tuple[str, dict, str]:
    out_json = tmp_path / "report.json"
    out_vcf = tmp_path / "annotated.vcf"
    result = CliRunner(env={"COLUMNS": "200"}).invoke(
        main,
        [
            "analyze",
            str(_vcf(tmp_path)),
            "--data-dir",
            str(data_dir),
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
    return result.output, json.loads(out_json.read_text()), out_vcf.read_text()


class TestEndToEnd:
    def test_json_input_block_counts_every_reason(self, tmp_path, data_dir) -> None:
        _output, report, _vcf_text = _analyze(tmp_path, data_dir)
        quality = report["input"]["vcf_quality"]
        assert quality["half_calls"] == 2
        assert quality["trimmed_matches"] == 4
        assert quality["spanning_deletion_sites"] == ["1:7000 A>G,* 1/2", "1:8000 C>* 1/1"]

    def test_json_rows_use_clinvar_position_and_alleles(self, tmp_path, data_dir) -> None:
        _output, report, _vcf_text = _analyze(tmp_path, data_dir)
        found = {(a["position"], a["ref"], a["alt"]) for a in report["annotations"]}
        assert {(1_000, "GT", "G"), (2_000, "C", "T"), (3_001, "C", "T")} <= found
        assert (5_000, "AC", "A") in found

    def test_vcf_out_stamps_the_allele_as_written(self, tmp_path, data_dir) -> None:
        """GTT>GT,G 0/1: the deletion matched as GT>G, so the stamp goes on
        the first ALT (GT), not on the second (G)."""
        _output, _report, vcf_text = _analyze(tmp_path, data_dir)
        line = next(ln for ln in vcf_text.splitlines() if ln.startswith("1\t1000\t"))
        stamp = next(f for f in line.split("\t")[7].split(";") if f.startswith("ALLELIX_CLINVAR="))
        first, second = stamp.removeprefix("ALLELIX_CLINVAR=").split(",")
        assert first.startswith("clinvar_pathogenic")
        assert second == "."

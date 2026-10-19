# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar matching: REF identity, hom-ref suppression by locus, reference-carried sites.

- A VCF is forward strand: a REF that is the complement of ClinVar's is not
  a match, and is counted as a REF mismatch.
- ClinPGx and SNPedia suppress a hom-ref call by ClinVar's REF only where
  ClinVar places the rsID at the call's locus.
- A hom-ref VCF call at a site where the reference carries a ClinVar
  variant allele is that variant, homozygous.
- A match carries its row's own rsID; an indel record's call keeps a slash.

All rsIDs, loci and rows here are synthetic.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.base import is_clinvar_homref
from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.models import Variant
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.reports._pipeline import AnalysisResult

GRCH37 = "GRCh37"

# (rsid, chromosome, position, ref, alt, significance, allele_id)
_ROWS = (
    # An ordinary SNV.
    ("rs900227001", "2", 1000, "G", "A", "Pathogenic", 1),
    # The reference carries the variant allele here (written REF == ALT).
    ("rs900230001", "2", 2000, "T", "T", "Pathogenic", 2),
    # Two ALTs of one site, each with its own rsID.
    ("rs900244001", "2", 3000, "A", "G", "Pathogenic", 3),
    ("rs900244002", "2", 3000, "A", "T", "Likely_pathogenic", 4),
    # A one-base deletion.
    ("rs900244003", "2", 4000, "AC", "A", "Pathogenic", 5),
)

_AGGREGATE = {"Pathogenic": "Pathogenic", "Likely_pathogenic": "Likely pathogenic"}


def _write_clinvar(data: Path) -> None:
    path = data / clinvar_db_filename(GRCH37)
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total) "
            "VALUES (?, ?, ?, ?, ?, ?, 'Synthetic condition', 'GENE', "
            "'criteria_provided,_single_submitter', ?, ?, 1, 1)",
            [(*row, _AGGREGATE[row[5]]) for row in _ROWS],
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://clinvar', '20261016', '2026-10-16', ?, NULL, ?)",
            (clinvar_record_name(GRCH37), len(_ROWS), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


@pytest.fixture
def data(tmp_path: Path) -> Path:
    out = tmp_path / "data"
    out.mkdir()
    _write_clinvar(out)
    return out


_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=2>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
)
_GVCF_HEADER = '##ALT=<ID=NON_REF,Description="Any other allele">\n'
_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def _vcf(path: Path, records: list[str], *, gvcf: bool = False) -> Path:
    body = "".join(f"{r}\n" for r in records)
    path.write_text(_HEADER + (_GVCF_HEADER if gvcf else "") + _COLUMNS + body)
    return path


def _run(data: Path, vcf: Path, out: Path | None = None) -> AnalysisResult:
    with ClinVarAnnotator(data) as clinvar:
        result = run_analysis(vcf, VcfParser(), [clinvar], build_override=GRCH37)
    if out is not None:
        from allelix.reports.json_report import render_json

        render_json(result, output_path=out)
    return result


def _clinvar(result: AnalysisResult) -> list[tuple[str, str, str]]:
    return sorted((a.rsid, a.significance, a.genotype_match) for a in result.annotations)


# ── A complemented REF is a REF mismatch ────────────────────────


def test_vcf_record_with_a_complemented_ref_is_not_matched(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t1000\trs900227001\tC\tT\t50\tPASS\t.\tGT:GQ\t0/1:50"])
    out = tmp_path / "r.json"
    result = _run(data, vcf, out)
    assert _clinvar(result) == []
    quality = json.loads(out.read_text())["input"]["vcf_quality"]
    assert quality["clinvar_ref_mismatches"] == 1


def test_vcf_record_with_clinvars_ref_matches(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t1000\trs900227001\tG\tA\t50\tPASS\t.\tGT:GQ\t0/1:50"])
    result = _run(data, vcf)
    assert _clinvar(result) == [("rs900227001", "clinvar_pathogenic", "AG")]
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 0


# ── Hom-ref suppression only at ClinVar's own locus ─────────────


def test_reference_at_needs_the_locus(data: Path) -> None:
    with ClinVarAnnotator(data) as clinvar:
        assert clinvar.reference_for("rs900227001", GRCH37) == "G"
        assert clinvar.reference_at("rs900227001", GRCH37, ("2", 1000)) == "G"
        assert clinvar.reference_at("rs900227001", GRCH37, ("2", 1001)) is None
        assert clinvar.reference_at("rs900227001", GRCH37, ("3", 1000)) is None
        here = Variant("rs900227001", "2", 1000, "G", "G", build=GRCH37)
        elsewhere = Variant("rs900227001", "2", 5000, "G", "G", build=GRCH37)
        unplaced = Variant("rs900227001", "2", 0, "G", "G", build=GRCH37)
        assert is_clinvar_homref(here, clinvar.reference_at)
        assert not is_clinvar_homref(elsewhere, clinvar.reference_at)
        assert not is_clinvar_homref(unplaced, clinvar.reference_at)


def test_snpedia_is_not_suppressed_by_clinvar_elsewhere(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.execute(
            "CREATE TABLE snpedia_genotypes (rsid TEXT, allele1 TEXT, allele2 TEXT, "
            "magnitude REAL, repute TEXT, summary TEXT, gene TEXT, orientation TEXT)"
        )
        conn.execute(
            "INSERT INTO snpedia_genotypes VALUES "
            "('rs900227001', 'G', 'G', 2.0, 'Good', 'common genotype', 'GENE', 'plus')"
        )
        conn.commit()
    with ClinVarAnnotator(data) as clinvar:
        snpedia = SNPediaAnnotator(data, clinvar_ref_provider=clinvar.reference_at)
        with snpedia:
            here = Variant("rs900227001", "2", 1000, "G", "G", build=GRCH37)
            elsewhere = Variant("rs900227001", "2", 5000, "G", "G", build=GRCH37)
            assert snpedia.annotate(here) == []
            assert [a.genotype_match for a in snpedia.annotate(elsewhere)] == ["GG"]


# ── Hom-ref calls at reference-carried sites ────────────────────


def test_hom_ref_call_at_a_reference_carried_site_is_the_variant(
    data: Path, tmp_path: Path
) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t2000\t.\tT\tC\t50\tPASS\t.\tGT:GQ\t0/0:50"])
    result = _run(data, vcf)
    assert _clinvar(result) == [("rs900230001", "clinvar_pathogenic", "TT")]
    (row,) = result.annotations
    assert row.zygosity_code == "homozygous"
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_carried_calls == 1


def test_reference_block_over_a_reference_carried_site(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(
        tmp_path / "x.g.vcf",
        ["2\t1990\t.\tA\t<NON_REF>\t.\tPASS\tEND=2010\tGT:GQ\t0/0:50"],
        gvcf=True,
    )
    result = _run(data, vcf)
    assert _clinvar(result) == [("rs900230001", "clinvar_pathogenic", "TT")]


def test_low_quality_hom_ref_is_not_assessed(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t2000\t.\tT\tC\t50\tPASS\t.\tGT:GQ\t0/0:5"])
    result = _run(data, vcf)
    assert _clinvar(result) == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_carried_calls == 0


def test_plain_vcf_without_a_record_counts_the_site(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t1000\trs900227001\tG\tA\t50\tPASS\t.\tGT:GQ\t0/1:50"])
    out = tmp_path / "r.json"
    result = _run(data, vcf, out)
    assert ("rs900230001", "clinvar_pathogenic", "TT") not in _clinvar(result)
    quality = json.loads(out.read_text())["input"]["vcf_quality"]
    assert quality["reference_carried_not_in_file"] == 1
    assert quality["reference_carried_calls"] == 0


def test_gvcf_without_coverage_is_not_counted_as_not_in_file(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(
        tmp_path / "x.g.vcf",
        ["2\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=100\tGT:GQ\t0/0:50"],
        gvcf=True,
    )
    result = _run(data, vcf)
    assert _clinvar(result) == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_carried_not_in_file == 0


def test_a_carrier_record_is_matched_as_before(data: Path, tmp_path: Path) -> None:
    """A person without the variant has a 1/1 record T>C there: no ClinVar row fits."""
    vcf = _vcf(tmp_path / "x.vcf", ["2\t2000\t.\tT\tC\t50\tPASS\t.\tGT:GQ\t1/1:50"])
    result = _run(data, vcf)
    assert _clinvar(result) == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_carried_calls == 0


# ── The row's own rsID; an indel record's call keeps a slash ────


def test_multi_allelic_record_rows_carry_their_own_rsids(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t3000\trs900244001\tA\tG,T\t50\tPASS\t.\tGT:GQ\t1/2:50"])
    result = _run(data, vcf)
    assert _clinvar(result) == [
        ("rs900244001", "clinvar_pathogenic", "GT"),
        ("rs900244002", "clinvar_likely_pathogenic", "GT"),
    ]


def test_homozygous_one_base_deletion_keeps_a_slash(data: Path, tmp_path: Path) -> None:
    vcf = _vcf(tmp_path / "x.vcf", ["2\t4000\trs900244003\tAC\tA\t50\tPASS\t.\tGT:GQ\t1/1:50"])
    result = _run(data, vcf)
    assert _clinvar(result) == [("rs900244003", "clinvar_pathogenic", "A/A")]


def test_gwas_row_on_a_homozygous_one_base_deletion_keeps_a_slash(tmp_path: Path) -> None:
    from allelix.annotators.gwas import GWASCatalogAnnotator
    from allelix.databases.gwas_loader import _CATEGORIZER_VERSION, GWAS_DB_FILENAME
    from allelix.databases.schema import GWAS_SCHEMA

    with contextlib.closing(sqlite3.connect(tmp_path / GWAS_DB_FILENAME)) as conn:
        conn.executescript(GWAS_SCHEMA)
        conn.execute(
            "INSERT INTO gwas_associations (rsid, risk_allele, trait, p_value, or_beta, "
            "ci_text, gene, study_accession, pubmed_id, trait_category) "
            "VALUES ('rs900244003', 'A', 'Synthetic trait', 1e-30, 1.5, '[1.2-1.8]', "
            "'GENE', 'GCST0', '1', 'disease')"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES ('gwas', 'test://gwas', 'test', '2026-10-16T00:00:00Z', 1, '', ?)",
            (f"cv:{_CATEGORIZER_VERSION}",),
        )
        conn.commit()
    deletion = Variant(
        "rs900244003", "2", 4000, "A", "A", build=GRCH37, ref="AC", ref_from_input=True
    )
    with GWASCatalogAnnotator(tmp_path) as gwas:
        assert [a.genotype_match for a in gwas.annotate(deletion)] == ["A/A"]


def test_report_lines_for_the_new_vcf_counts() -> None:
    from allelix.parsers.vcf import VcfQualityStats
    from allelix.reports import vcf_quality_lines

    stats = VcfQualityStats()
    stats.clinvar_ref_mismatches = 2
    stats.reference_carried_calls = 3
    stats.reference_carried_not_in_file = 4
    text = " ".join(vcf_quality_lines(stats))
    assert "2 record(s) at a ClinVar position have a REF that is not ClinVar's" in text
    assert "3 homozygous-reference call(s)" in text
    assert "4 site(s) where the reference genome carries a ClinVar variant allele" in text

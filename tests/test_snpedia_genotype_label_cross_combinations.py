"""SNPedia's genotype label (``diploid_label``) crossed with input shapes not in the contract.

Lower-case records, a multi-ALT sequence call, a chip call (no REF from the
input) and a haploid call matched through a minus page.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.models import Variant
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from allelix.utils.allele import complement

if TYPE_CHECKING:
    from pathlib import Path

RS = "rs900000001"


def _cache(path: Path, a1: str, a2: str, orientation: str = "plus") -> None:
    if orientation == "minus":
        a1, a2 = complement(a1), complement(a2)
    a1, a2 = sorted((a1, a2))
    with contextlib.closing(sqlite3.connect(path / "snpedia.sqlite")) as conn, conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        conn.execute(
            "INSERT INTO snpedia_genotypes "
            "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
            "VALUES (?,?,?,6,'Bad','Synthetic finding','TEST',?)",
            (RS, a1, a2, orientation),
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
            "VALUES ('snpedia','test://synthetic','synthetic','2026-11-02',1,?)",
            (f"pv:{_PARSER_VERSION}",),
        )


def _analyze(path: Path, ref: str, alt: str, gt: str) -> list:
    src = path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t1000\t{RS}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
    )
    return run_analysis(
        src, VcfParser(), [SNPediaAnnotator(path)], build_override="GRCh38"
    ).annotations


def test_lower_case_record_gets_the_upper_case_label(tmp_path: Path) -> None:
    _cache(tmp_path, "GAC", "GAC")
    (finding,) = _analyze(tmp_path, "g", "gac", "1/1")
    assert finding.genotype_match == "GAC/GAC"
    assert finding.zygosity_code == "homozygous"


@pytest.mark.parametrize("gt", ["1/2", "2|1"])
def test_multi_alt_sequence_call_is_heterozygous(tmp_path: Path, gt: str) -> None:
    _cache(tmp_path, "GAC", "GT")
    (finding,) = _analyze(tmp_path, "G", "GAC,GT", gt)
    assert set(finding.genotype_match.split("/")) == {"GAC", "GT"}
    assert finding.zygosity_code == "heterozygous"


@pytest.mark.parametrize(
    ("a1", "a2", "label", "zygosity"),
    [("G", "A", "AG", "heterozygous"), ("A", "A", "AA", "homozygous")],
)
def test_chip_call_label_is_unchanged(
    tmp_path: Path, a1: str, a2: str, label: str, zygosity: str
) -> None:
    """A chip call has no REF from the input: two bases, joined and sorted, as before."""
    _cache(tmp_path, a1, a2)
    with SNPediaAnnotator(tmp_path) as source:
        (finding,) = source.annotate(Variant(RS, "1", 1000, a1, a2, build="GRCh38"))
    assert finding.genotype_match == label
    assert finding.zygosity_code == zygosity


@pytest.mark.parametrize("orientation", ["plus", "minus"])
def test_haploid_sequence_call_is_one_allele(tmp_path: Path, orientation: str) -> None:
    _cache(tmp_path, "GAC", "GAC", orientation)
    (finding,) = _analyze(tmp_path, "G", "GAC", "1")
    assert finding.genotype_match == "GAC"
    assert finding.zygosity_code == "hemizygous"

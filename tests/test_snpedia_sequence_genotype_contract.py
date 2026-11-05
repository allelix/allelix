"""O4: source page matching must not concatenate sequence alleles in reports."""

import contextlib
import json
import sqlite3

import pytest

from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.models import Variant, annotation_to_public_dict
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from allelix.reports.json_report import render_json

RS = "rs900000001"
CASES = [
    ("insertion-hom", "G", "GAC", "1/1", "GAC/GAC", "homozygous"),
    ("insertion-het", "G", "GAC", "0/1", "G/GAC", "heterozygous"),
    ("deletion-hom", "GAC", "G", "1/1", "G/G", "homozygous"),
    ("deletion-het", "GAC", "G", "0/1", "GAC/G", "heterozygous"),
    ("mnv-hom", "GAC", "ACT", "1/1", "ACT/ACT", "homozygous"),
    ("mnv-het", "GAC", "ACT", "0|1", "GAC/ACT", "heterozygous"),
    ("base-hom-control", "G", "A", "1/1", "AA", "homozygous"),
    ("base-het-control", "G", "A", "0|1", "AG", "heterozygous"),
]


def cache(path, a1, a2, orientation):
    if orientation == "minus":
        tr = str.maketrans("ACGT", "TGCA")
        a1, a2 = a1.translate(tr)[::-1], a2.translate(tr)[::-1]
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


@pytest.mark.parametrize("case", CASES, ids=[x[0] for x in CASES])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("orientation", ["plus", "minus"])
@pytest.mark.parametrize("path", ["single", "batch", "pipeline-json"])
def test_sequence_alleles_stay_separate(tmp_path, case, build, orientation, path):
    _, ref, alt, gt, expected_label, expected_zygosity = case
    a1, a2 = (ref if n == "0" else alt for n in gt.replace("|", "/").split("/"))
    cache(tmp_path, a1, a2, orientation)
    call = Variant(RS, "1", 1000, a1, a2, build=build, ref=ref, ref_from_input=True)
    with SNPediaAnnotator(tmp_path) as source:
        if path == "single":
            rows = source.annotate(call)
        elif path == "batch":
            rows = list(source.batch_annotate([call]))
        else:
            src = tmp_path / "input.vcf"
            src.write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                f"1\t1000\t{RS}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
            )
            result = run_analysis(src, VcfParser(), [source], build_override=build)
            rows = result.annotations
            out = tmp_path / "report.json"
            render_json(result, output_path=out, min_magnitude=0)
            public = json.loads(out.read_text())["annotations"]
            assert len(public) == 1, "the sequence page must actually match"
            assert public[0]["genotype_match"] == expected_label
            assert public[0]["zygosity"] == expected_zygosity
        assert len(rows) == 1, "the sequence page must actually match"
        # Source direct calls have not yet had pipeline carried_alleles patched.
        assert annotation_to_public_dict(rows[0])["genotype_match"] == expected_label
        assert rows[0].zygosity_code == expected_zygosity
        if path == "pipeline-json":
            assert sorted(rows[0].carried_alleles) == sorted((a1, a2))


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("padded", [False, True])
@pytest.mark.parametrize("haploid", [False, True])
def test_sequence_label_padding_and_ploidy(tmp_path, build, padded, haploid):
    ref, alt = ("CG", "CGAC") if padded else ("G", "GAC")
    cache(tmp_path, alt, alt, "plus")
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t1000\t{RS}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{'1' if haploid else '1/1'}\n"
    )
    result = run_analysis(src, VcfParser(), [SNPediaAnnotator(tmp_path)], build_override=build)
    assert len(result.annotations) == 1
    annotation = result.annotations[0]
    assert annotation.genotype_match == (alt if haploid else f"{alt}/{alt}")
    assert annotation.zygosity_code == ("hemizygous" if haploid else "homozygous")

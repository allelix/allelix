"""O1 writer boundary: primary, enrichment and REF INFO keep record provenance.

These hand-built annotations exercise every source-dispatch arm; the companion
pipeline test uses actual SNPedia output to prove the upstream eligibility leak.
"""

import pytest

from allelix.models import Annotation
from allelix.reports.vcf import render_vcf

SOURCES = ("clinvar", "snpedia", "pharmgkb", "gwas", "gnomad", "cadd", "alphamissense")


def finding(source, origin, *, rank=6, reference=False):
    enrichment = source in ("gnomad", "cadd", "alphamissense")
    return Annotation(
        source="clinvar" if enrichment else source,
        rsid="rs900000001",
        significance=f"row{origin}",
        category="clinical",
        magnitude=rank,
        description=f"record {origin}",
        attribution="Synthetic",
        genotype_match="AG",
        gene=f"GENE{origin}",
        trait=f"Trait{origin}",
        alt="G" if reference else "A",
        locus=("1", 1000),
        locus_ref="G",
        record_index=origin,
        reference_carries_alt=reference,
        allele_frequency=0.1 if source == "gnomad" else None,
        cadd_phred=20 if source == "cadd" else None,
        am_pathogenicity=0.9 if source == "alphamissense" else None,
        am_class="pathogenic" if source == "alphamissense" else None,
    )


def write(tmp_path, annotations, *, eligible, reference=False, gvcf=False):
    src = tmp_path / "input.vcf"
    alt = "A,<NON_REF>" if gvcf else "A"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + f"1\t1000\trs900000001\tG\t{alt}\t.\tPASS\t.\tGT\t0/1\n"
        * 2
    )
    key = ("1", 1000, "G", "G" if reference else "A")
    output = tmp_path / "output.vcf"
    render_vcf(
        input_path=src,
        output_path=output,
        build="GRCh38",
        annotators_used=[(s, "synthetic") for s in SOURCES],
        annotations_by_key={key: annotations},
        recovered_rsids={},
        record_rows={key: eligible},
        reference_rows={},
        record_recovered={},
    )
    return [
        line.split("\t")[7] for line in output.read_text().splitlines() if not line.startswith("#")
    ]


@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("origin", [0, 1])
@pytest.mark.parametrize("gvcf", [False, True])
def test_no_info_is_borrowed_from_another_records_finding(tmp_path, source, origin, gvcf):
    info = write(tmp_path, [finding(source, origin)], eligible={0, 1}, gvcf=gvcf)
    assert "ALLELIX_" in info[origin], "the originating record must still be stamped"
    assert "ALLELIX_" not in info[1 - origin], (
        "eligibility is not ownership of every allele finding"
    )


@pytest.mark.parametrize("source", SOURCES[:4])
@pytest.mark.parametrize("higher_rank", [0, 1])
def test_report_rank_is_chosen_within_each_record(tmp_path, source, higher_rank):
    annotations = [finding(source, i, rank=9 if i == higher_rank else 5) for i in range(2)]
    info = write(tmp_path, annotations, eligible={0, 1})
    token = "Trait" if source == "gwas" else "row"
    for i in range(2):
        assert f"{token}{i}" in info[i]
        assert f"{token}{1 - i}" not in info[i]


@pytest.mark.parametrize("origin", [0, 1])
@pytest.mark.parametrize("both_eligible", [False, True])
def test_reference_info_keeps_its_origin(tmp_path, origin, both_eligible):
    info = write(
        tmp_path,
        [finding("clinvar", origin, reference=True)],
        reference=True,
        eligible={0, 1} if both_eligible else {origin},
    )
    assert "ALLELIX_CLINVAR_REF=" in info[origin]
    assert "ALLELIX_CLINVAR_REF=" not in info[1 - origin]
    assert all("ALLELIX_CLINVAR=" not in item for item in info)

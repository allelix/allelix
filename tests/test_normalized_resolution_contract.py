"""Call normalization precedes resolution without rewriting original VCF records."""

import pytest
from test_annotation_parity_and_named_disclosures import _run
from validated_call_helpers import (
    BUILDS,
    RSID,
    SOURCES,
    analyze,
    call,
    homozygous_sources,
    setup,
    signatures,
)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["plain", "prefix", "suffix", "both"])
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("gt", ["0/1", "0|1", "1/1", "1", "./1"])
def test_resolution_and_interpretation_are_spelling_independent(tmp_path, build, shape, named, gt):
    data = tmp_path / "data"
    pos = setup(data, build)
    if gt in ("1/1", "1"):
        homozygous_sources(data)
    control = analyze(data, call(tmp_path / "control.vcf", pos, gt=gt), build)
    assert {a.source for a in control.annotations} >= {"clinvar", "gwas"}
    result = analyze(
        data,
        call(tmp_path / "input.vcf", pos, shape=shape, gt=gt, rsid=RSID if named else "."),
        build,
    )
    assert signatures(result) == signatures(control)
    assert all(a.haploid == (gt == "1") for a in result.annotations)


@pytest.mark.parametrize("shape", ["prefix", "suffix", "both"])
@pytest.mark.parametrize("build", BUILDS)
def test_annotated_output_preserves_original_spelling(tmp_path, shape, build):
    data = tmp_path / "data"
    pos = setup(data, build)
    source = call(tmp_path / "input.vcf", pos, shape=shape, rsid=".")
    out = tmp_path / "annotated.vcf"
    report, _ = _run(
        tmp_path,
        data,
        source,
        build,
        "json",
        extra=("--min-magnitude", "0", "--gwas-min-magnitude", "0", "--vcf-out", str(out)),
    )
    assert {a["source"] for a in report["annotations"]} == set(SOURCES)
    original = source.read_text().splitlines()[-1].split("\t")
    written = [
        line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")
    ]
    assert len(written) == 1
    assert written[0][:2] == original[:2]
    assert written[0][3:7] == original[3:7]
    assert written[0][8:] == original[8:]
    assert written[0][2] == RSID


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["prefix", "suffix"])
@pytest.mark.parametrize("gt", ["1/2", "2|1"])
def test_multiallelic_named_and_unnamed_padding_agree(tmp_path, build, shape, gt):
    from validated_call_helpers import identity, sql

    data = tmp_path / "data"
    pos = setup(data, build)
    identity(data, "AGT")
    sql(data, "pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='AT'")
    sql(data, "snpedia.sqlite", "UPDATE snpedia_genotypes SET allele1='A',allele2='T'")
    alt = "CA,CT" if shape == "prefix" else "AC,TC"
    control = analyze(data, call(tmp_path / "named.vcf", pos, shape=shape, gt=gt, alt=alt), build)
    assert {a.source for a in control.annotations} >= {"clinvar", "pharmgkb", "snpedia"}
    result = analyze(
        data, call(tmp_path / "unnamed.vcf", pos, shape=shape, gt=gt, alt=alt, rsid="."), build
    )
    assert signatures(result) == signatures(control)

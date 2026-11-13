"""A finding belongs to its input record, never to a locus-wide ploidy label."""

from collections import Counter

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import (
    BUILDS,
    RSID,
    analyze,
    call,
    homozygous_sources,
    setup,
    warning_count,
)

from allelix.reports import _pipeline


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("batch_size", [1, 1000])
@pytest.mark.parametrize("same_alt", [False, True])
@pytest.mark.parametrize("shape", ["plain", "suffix"])
@pytest.mark.parametrize("named", [False, True])
def test_diploid_record_retains_ploidy_unless_distinct_alts_require_abstention(
    tmp_path, monkeypatch, build, reverse, batch_size, same_alt, shape, named
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    ref, alt, other = ("G", "A", "T") if shape == "plain" else ("GC", "AC", "TC")
    diploid = (RSID if named else ".", "1", pos, ref, alt, "0/1", ".")
    haploid = (RSID if named else ".", "1", pos, ref, alt if same_alt else other, "1", ".")
    control = analyze(data, _vcf(tmp_path / "control.vcf", [diploid]), build)
    rows = [haploid, diploid]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    expected = Counter((a.source, a.genotype_match, a.zygosity_code) for a in control.annotations)
    found = Counter((a.source, a.genotype_match, a.zygosity_code) for a in result.annotations)
    assert expected, "The isolated diploid control must produce findings"
    if shape == "plain" and not same_alt:
        # The two different SNV ALTs share one REF, but mixed ploidy cannot
        # form the exact complementary diploid pair required in 2.7.0.
        assert not result.annotations
        assert warning_count(result, "split_site_calls") == 2
    else:
        assert not (expected - found)
        assert warning_count(result, "split_site_calls") == 0


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch_size", [1, 1000])
@pytest.mark.parametrize("reverse", [False, True])
def test_identical_allele_spelling_does_not_merge_haploid_and_diploid_origins(
    tmp_path, monkeypatch, build, batch_size, reverse
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    homozygous_sources(data)
    rows = [(RSID, "1", pos, "G", "A", gt, ".") for gt in ("1", "1/1")]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    for source in ("clinvar", "pharmgkb", "snpedia", "gwas"):
        annotations = [a for a in result.annotations if a.source == source]
        assert Counter(a.haploid for a in annotations) == {True: 1, False: 1}


@pytest.mark.parametrize("batch_size", [1, 1000])
@pytest.mark.parametrize("conflicting", [False, True])
def test_contradictory_records_are_disclosed_once_per_affected_record(
    tmp_path, monkeypatch, batch_size, conflicting
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    row = (RSID, "1", pos, "G", "A", "0/1", ".")
    second = (*row[:5], "1" if conflicting else "0/1", ".")
    result = analyze(data, _vcf(tmp_path / "input.vcf", [row, second]), "GRCh38")
    assert warning_count(result, "contradictory_calls") == (2 if conflicting else 0)


@pytest.mark.parametrize("gt,expected", [("1", True), ("1/1", False)])
def test_isolated_ploidy_control(tmp_path, gt, expected):
    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    homozygous_sources(data)
    result = analyze(data, call(tmp_path / "input.vcf", pos, gt=gt), "GRCh38")
    assert result.annotations
    assert all(a.haploid == expected for a in result.annotations)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("batch_size", [1, 1000])
def test_chip_calls_also_keep_record_local_ploidy(
    tmp_path, monkeypatch, build, reverse, batch_size
):
    from test_haploid_and_homozygous_warning import _array

    from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser

    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch_size)
    data = tmp_path / "data"
    pos = setup(data, build)
    rows = [(RSID, "1", pos, "T"), (RSID, "1", pos, "AG")]
    if reverse:
        rows.reverse()
    path = _array(tmp_path / "input.txt", rows)
    result = analyze(data, path, build, parser=TwentyThreeAndMeParser())
    diploid = [a for a in result.annotations if a.genotype_match == "AG"]
    assert diploid and all(a.zygosity_code == "heterozygous" for a in diploid)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
def test_annotations_retain_the_input_record_index(tmp_path, build, reverse):
    data = tmp_path / "data"
    pos = setup(data, build)
    homozygous_sources(data)
    genotypes = ["1", "1/1"]
    if reverse:
        genotypes.reverse()
    rows = [(RSID, "1", pos, "G", "A", gt, ".") for gt in genotypes]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build)
    assert len(result.annotations) == 8
    for annotation in result.annotations:
        expected = genotypes.index("1" if annotation.haploid else "1/1")
        assert getattr(annotation, "record_index", None) == expected

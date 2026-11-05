"""Each genotype source names its own sites for explicit VCF reference calls (ADR-0064)."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import RSID, analyze, reference_effect


@pytest.mark.parametrize("source", ["gwas", "pharmgkb", "snpedia"])
def test_one_source_alone_reads_an_explicit_reference_call(tmp_path, source):
    data = tmp_path / "data"
    pos = reference_effect(data, "GRCh38")
    path = _vcf(tmp_path / "input.vcf", [(".", "1", pos, "G", "A", "0/0", ".")])
    result = analyze(data, path, "GRCh38", (source,))
    assert {a.source for a in result.annotations} == {source}
    assert {a.rsid for a in result.annotations} == {RSID}


def test_other_annotators_name_no_reference_call_sites(tmp_path):
    from allelix.annotators.clinvar import ClinVarAnnotator

    data = tmp_path / "data"
    reference_effect(data, "GRCh38")
    with ClinVarAnnotator(data) as clinvar:
        assert clinvar.reference_call_rsids() == frozenset()


def test_one_rsid_at_two_placements_keeps_each_record(tmp_path):
    # GWAS Catalog and SNPedia abstain on an rsID with two placements.
    from test_annotation_parity_and_named_disclosures import _identity

    data = tmp_path / "data"
    pos = reference_effect(data, "GRCh38")
    (data / "rsid_identity.sqlite").unlink()
    _identity(data, [(RSID, "1", p, "G", "A", p, "G", "A", 0) for p in (pos, pos + 500)])
    rows = [(".", "1", pos, "G", "A", "0", "."), (".", "1", pos + 500, "G", "A", "0/0", ".")]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), "GRCh38", ("pharmgkb",))
    found = {(bool(a.haploid), a.record_index) for a in result.annotations}
    assert len({index for _haploid, index in found}) == 2
    assert {haploid for haploid, _index in found} == {True, False}
    assert len(found) == 2


_SOURCES = ("gwas", "pharmgkb", "snpedia")


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("order", ["reference_first", "call_first"])
@pytest.mark.parametrize(
    "other",
    [
        pytest.param(("G", "T", "1/1"), id="split_alt_call"),
        pytest.param(("G", "T", "0/1"), id="split_het"),
        pytest.param(("G", "T", "./."), id="split_no_call"),
    ],
)
def test_a_record_saying_more_than_reference_blocks_the_reference_call(
    tmp_path, build, order, other
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    ref, alt, gt = other
    rows = [(".", "1", pos, "G", "A", "0/0", "."), (".", "1", pos, ref, alt, gt, ".")]
    if order == "call_first":
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, _SOURCES)
    assert not [a for a in result.annotations if a.genotype_match in ("GG", "G/G", "G")]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("order", ["padded_first", "call_first"])
def test_a_padded_reference_record_over_a_called_site_gives_no_reference_call(
    tmp_path, build, order
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    rows = [(".", "1", pos - 1, "CG", "C", "0|0", "."), (".", "1", pos, "G", "A", "1|1", ".")]
    if order == "call_first":
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, _SOURCES)
    assert not [a for a in result.annotations if a.genotype_match in ("GG", "G/G", "G")]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_two_reference_records_at_one_site_still_give_the_reference_call(tmp_path, build):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    rows = [(".", "1", pos, "G", "A", "0/0", "."), (".", "1", pos, "G", "T", "0/0", ".")]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, _SOURCES)
    assert {a.source for a in result.annotations} == set(_SOURCES)

"""Direct and complement fits are evaluated together; known orientation survives."""

import pytest
from validated_call_helpers import (
    BUILDS,
    analyze,
    call,
    homozygous_sources,
    identity,
    setup,
    sql,
)

from allelix.utils.allele import forward_strand_reading


@pytest.mark.parametrize(
    "forward", ["AG", "CT", "AC", "GT", "AT", "CG", "ACG", "AGT", "ACT", "CGT", "ACGT"]
)
@pytest.mark.parametrize("allele", "ACGT")
def test_source_allele_requires_one_reading(forward, allele):
    complement = dict(zip("ACGT", "TGCA", strict=True))[allele]
    direct, reverse = allele in forward, complement in forward
    answer = forward_strand_reading((allele,), frozenset(forward))
    if direct != reverse:
        assert answer == ("direct" if direct else "complement")
    elif direct:
        assert answer not in ("direct", "complement")
    else:
        assert answer == "no_fit"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("alleles", ["AG", "AGT", "ACGT"])
@pytest.mark.parametrize("risk", ["A", "T"])
@pytest.mark.parametrize("gt", ["1/1", "1"])
def test_gwas_direct_and_reverse_paths_have_equal_ambiguity_rules(
    tmp_path, build, alleles, risk, gt
):
    data = tmp_path / "data"
    pos = setup(data, build)
    identity(data, alleles)
    sql(data, "gwas.sqlite", "UPDATE gwas_associations SET risk_allele=?", (risk,))
    result = analyze(data, call(tmp_path / "input.vcf", pos, gt=gt), build, ("gwas",))
    if "T" not in alleles:
        assert len(result.annotations) == 1
        assert result.annotations[0].alt == "A"
        assert not sum(result.strand_abstentions.get("gwas", {}).values())
    else:
        assert not result.annotations
        assert sum(result.strand_abstentions.get("gwas", {}).values()) == 1


@pytest.mark.parametrize("orientation", ["plus", "minus", ""])
@pytest.mark.parametrize("build", BUILDS)
def test_independently_oriented_snpedia_rows_are_not_discarded(tmp_path, orientation, build):
    data = tmp_path / "data"
    pos = setup(data, build)
    identity(data, "ACGT")
    homozygous_sources(data)
    allele = "T" if orientation == "minus" else "A"
    sql(
        data,
        "snpedia.sqlite",
        "UPDATE snpedia_genotypes SET orientation=?,allele1=?,allele2=?",
        (orientation, allele, allele),
    )
    result = analyze(data, call(tmp_path / "input.vcf", pos, gt="1/1"), build, ("snpedia",))
    assert len(result.annotations) == bool(orientation)
    if not orientation:
        assert sum(result.strand_abstentions.get("snpedia", {}).values()) == 1

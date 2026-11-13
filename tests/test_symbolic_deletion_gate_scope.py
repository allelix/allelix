"""A deletion span follows the selected allele and the actual quality gate."""

import pytest
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect

from allelix.parsers.vcf import VcfParser

PATHS = (("gwas", "snpedia", "pharmgkb"), ("clinvar",))


def _input(
    path,
    pos,
    *,
    alt="<DEL>",
    gt="0/0",
    filt="PASS",
    extra=(),
    reverse=False,
    reference=True,
    end=None,
    samples=None,
):
    fmt = ":".join(("GT", *(key for key, _ in extra)))
    value = ":".join((gt, *(value for _, value in extra)))
    names = ("S",) if samples is None else tuple(samples)
    values = (value,) if samples is None else tuple(samples.values())
    deletion = (
        f"1\t{pos - 10}\t.\tG\t{alt}\t.\t{filt}\tEND={pos + 10 if end is None else end}"
        f"\t{fmt}\t" + "\t".join(values)
    )
    ref = f"1\t{pos}\t{RSID}\tG\tA\t.\tPASS\t.\tGT\t" + "\t".join("0/0" for _ in names)
    rows = [deletion, ref] if reference else [deletion]
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(names)
        + "\n"
        + "\n".join(rows[::-1] if reverse else rows)
        + "\n"
    )
    return path


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize(
    "alt,gt",
    (("<DEL>,<DUP>", "0/2"), ("<DUP>,<DEL>", "0/1"), ("<DEL>,*", "0/2"), ("<DEL>,A", "0/2")),
)
def test_an_unselected_deletion_cannot_borrow_a_symbolic_no_call_reason(
    tmp_path, build, sources, reverse, alt, gt
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _input(tmp_path / "input.vcf", pos, alt=alt, gt=gt, reverse=reverse)
    result = analyze(data, source, build, sources)
    assert sorted(a.source for a in result.annotations) == sorted(sources), (
        "A passing call of the OTHER ALT does not carry the deletion; symbolic_allele "
        "and spanning_deletion are not quality-gate rejections"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
@pytest.mark.parametrize(
    "gt,filt,extra,options,blocked",
    (
        ("0/0", "PASS", (("GQ", "1"),), {}, True),
        ("0/0", "PASS", (("DP", "0"),), {}, True),
        ("0/0", "PASS", (("DP", "0"),), {"min_gq": 0, "keep_filtered": True}, True),
        ("0/0", "LowQual", (), {}, True),
        ("0/0", "PASS", (("FT", "LowQual"),), {}, True),
        ("0/0", "LowQual", (), {"keep_filtered": True}, False),
        ("0/0", "PASS", (("FT", "LowQual"),), {"keep_filtered": True}, False),
        ("0/0", "PASS", (("GQ", "1"),), {"min_gq": 0}, False),
        ("0/0", "PASS", (("GQ", "99"), ("DP", "10")), {}, False),
        ("./.", "PASS", (), {}, True),
        ("0/.", "PASS", (), {}, True),
        ("1", "PASS", (), {}, True),
    ),
    ids=(
        "low-gq",
        "zero-depth",
        "zero-depth-not-overridden",
        "site-filter",
        "sample-filter",
        "keep-site-filter",
        "keep-sample-filter",
        "relaxed-gq",
        "passing-reference",
        "unknown",
        "partial",
        "carried-haploid",
    ),
)
def test_deletion_quality_gate_and_unknown_alleles_apply_to_both_reference_paths(
    tmp_path, build, sources, gt, filt, extra, options, blocked
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _input(tmp_path / "input.vcf", pos, gt=gt, filt=filt, extra=extra)
    result = analyze(data, source, build, sources, parser=VcfParser(**options))
    assert sorted(a.source for a in result.annotations) == ([] if blocked else sorted(sources)), (
        "The symbolic END span follows the effective FILTER/FT/GQ/DP gate, not just GT"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
@pytest.mark.parametrize("extra", ((("GQ", "1"),), (("DP", "0"),)), ids=("low-gq", "zero-depth"))
def test_rejected_other_symbolic_alt_makes_the_deletion_possible(tmp_path, build, sources, extra):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _input(tmp_path / "input.vcf", pos, alt="<DEL>,<DUP>", gt="0/2", extra=extra)
    result = analyze(data, source, build, sources)
    assert not result.annotations, (
        "Quality rejection makes every ALT possible, even an unselected DEL"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
@pytest.mark.parametrize("sample,blocked", (("LOW", True), ("PASSING", False)))
def test_quality_of_another_sample_does_not_decide_the_deletion_span(
    tmp_path, build, sources, sample, blocked
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _input(
        tmp_path / "samples.vcf",
        pos,
        extra=(("GQ", "99"),),
        samples={"LOW": "0/0:1", "PASSING": "0/0:99"},
    )
    result = analyze(data, source, build, sources, parser=VcfParser(sample=sample))
    assert sorted(a.source for a in result.annotations) == ([] if blocked else sorted(sources))


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
@pytest.mark.parametrize("alt", ("<DEL>,A", "A,<DEL>"))
def test_an_uncarried_deletion_cannot_supply_reference_evidence_from_its_end(
    tmp_path, build, sources, alt
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _input(tmp_path / "only.vcf", pos, alt=alt, gt="0/0", reference=False)
    result = analyze(data, source, build, sources)
    assert not result.annotations, "The target has no explicit reference record or reference block"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", PATHS, ids=("reference-genotype", "clinvar-carried"))
def test_a_long_deletion_span_is_cleared_before_reusing_the_parser(tmp_path, build, sources):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    parser = VcfParser()
    long_span = _input(tmp_path / "long.vcf", pos, alt="<DEL:ME>", gt="0/1", end=10**12)
    assert not analyze(data, long_span, build, sources, parser=parser).annotations
    next_input = _input(tmp_path / "next.vcf", pos, gt="0/0")
    result = analyze(data, next_input, build, sources, parser=parser)
    assert sorted(a.source for a in result.annotations) == sorted(sources)

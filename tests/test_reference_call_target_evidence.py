"""Reference evidence describes the target allele, with trimming but no left-alignment."""

import pytest
from test_annotation_parity_and_named_disclosures import _identity, _vcf
from validated_call_helpers import BUILDS, POSITIONS, RSID, analyze, reference_effect, sql

from allelix.annotators.pharmgkb import PharmGKBAnnotator

TARGETS = {
    "snv": ("G", "A", "T"),
    "insertion": ("G", "GA", "GT"),
    "deletion": ("GC", "G", "C"),
    "mnv": ("GC", "AT", "TA"),
}


def _padded(pos, ref, alt):
    return pos - 1, "T" + ref + "C", "T" + alt + "C"


def _cache(data, build, ref, alt, *, padded=False):
    pos = reference_effect(data, build)
    sites = [(POSITIONS[b], ref, alt) for b in BUILDS]
    if padded:
        sites = [_padded(*site) for site in sites]
    (data / "rsid_identity.sqlite").unlink()
    _identity(data, [(RSID, "1", *sites[0], *sites[1], 0)])
    return pos


@pytest.fixture
def observed(monkeypatch):
    # Observe real source input, retaining its normal annotation behavior.
    # Multi-base genotypes may be unsupported by the source; eligibility
    # must still distinguish a matching record from unrelated evidence.
    seen = []
    original = PharmGKBAnnotator.batch_annotate

    def track(source, variants):
        variants = list(variants)
        seen.extend(variants)
        yield from original(source, variants)

    monkeypatch.setattr(PharmGKBAnnotator, "batch_annotate", track)
    return seen


def _reference_ids(seen):
    return [v.rsid for v in seen if v.allele1 == v.allele2 == v.ref and not v.is_no_call]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize(
    "evidence", ("matching", "other_alt", "anchor_snv", "no_alt", "multi_match", "multi_miss")
)
def test_only_evidence_for_the_target_reaches_the_source(
    tmp_path, observed, build, named, target, evidence
):
    data = tmp_path / "data"
    ref, alt, other = TARGETS[target]
    pos = _cache(data, build, ref, alt)
    record_ref, record_alt = ref, alt
    if evidence == "other_alt":
        record_alt = other
    elif evidence == "anchor_snv":
        record_ref, record_alt = ref[0], "T"
    elif evidence == "no_alt":
        record_alt = "."
    elif evidence == "multi_match":
        record_alt = f"{other},{alt}"
    elif evidence == "multi_miss":
        record_alt = "T,C" if target == "snv" else f"{other},{other}A"
    path = _vcf(
        tmp_path / "input.vcf",
        [(RSID if named else ".", "1", pos, record_ref, record_alt, "0/0", ".")],
    )
    analyze(data, path, build, ("pharmgkb",))
    allowed = target == "snv" or evidence in ("matching", "multi_match")
    assert _reference_ids(observed) == ([RSID] if allowed else [])
    if allowed:
        call = next(v for v in observed if v.rsid == RSID)
        assert (call.ref, call.allele1, call.allele2) == (record_ref,) * 3


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("cache_padded", (False, True))
@pytest.mark.parametrize("record_padded", (False, True))
def test_matching_target_and_record_are_compared_after_trimming(
    tmp_path, observed, build, named, target, cache_padded, record_padded
):
    data = tmp_path / "data"
    ref, alt, _other = TARGETS[target]
    pos = _cache(data, build, ref, alt, padded=cache_padded)
    if record_padded:
        pos, ref, alt = _padded(pos, ref, alt)
    path = _vcf(
        tmp_path / "input.vcf",
        [(RSID if named else ".", "1", pos, ref, alt, "0|0", ".")],
    )
    analyze(data, path, build, ("pharmgkb",))
    assert _reference_ids(observed) == [RSID]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", (False, True))
def test_repeat_equivalent_insertion_at_another_position_is_not_left_aligned(
    tmp_path, observed, build, named
):
    data = tmp_path / "data"
    pos = _cache(data, build, "G", "GG")
    path = _vcf(
        tmp_path / "input.vcf",
        [(RSID if named else ".", "1", pos - 1, "G", "GG", "0/0", ".")],
    )
    analyze(data, path, build, ("pharmgkb",))
    assert not _reference_ids(observed)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("shape", ("one", "other", "combined", "split", "carried_split"))
def test_multiallelic_reference_evidence_is_attributed_to_matching_targets(
    tmp_path, observed, build, reverse, shape
):
    data = tmp_path / "data"
    pos = _cache(data, build, "G", "GA")
    second = "rs900260902"
    sql(
        data,
        "pharmgkb.sqlite",
        "INSERT INTO pharmgkb_annotations "
        "(rsid,genotype,gene,drugs,phenotype,phenotype_category,annotation_text,"
        "level_of_evidence,score,pgkb_annotation_id,function_class,is_nonfinding) "
        "SELECT ?,genotype,gene,drugs,phenotype,phenotype_category,annotation_text,"
        "level_of_evidence,score,'synthetic-second',function_class,is_nonfinding "
        "FROM pharmgkb_annotations WHERE rsid=?",
        (second, RSID),
    )
    (data / "rsid_identity.sqlite").unlink()
    _identity(
        data,
        [
            (rsid, "1", POSITIONS["GRCh37"], "G", alt, POSITIONS["GRCh38"], "G", alt, 0)
            for rsid, alt in ((RSID, "GA"), (second, "GT"))
        ],
    )
    choices = {
        "one": [("GA", "0/0")],
        "other": [("GT", "0/0")],
        "combined": [("GA,GT", "0/0")],
        "split": [("GA", "0/0"), ("GT", "0/0")],
        "carried_split": [("GA", "0/0"), ("GT", "1/1")],
    }
    rows = [(".", "1", pos, "G", alt, gt, ".") for alt, gt in choices[shape]]
    if reverse:
        rows.reverse()
    analyze(data, _vcf(tmp_path / "input.vcf", rows), build, ("pharmgkb",))
    expected = {
        "one": [RSID],
        "other": [second],
        "combined": [RSID, second],
        "split": [RSID, second],
        "carried_split": [],
    }[shape]
    assert sorted(_reference_ids(observed)) == sorted(expected)

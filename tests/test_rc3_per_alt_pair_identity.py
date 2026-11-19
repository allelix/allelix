"""An exact pair resolves each ALT's own ID before choosing full or half calls."""

import pytest
from rc3_verification_helpers import (
    OTHER,
    POS,
    RS,
    analyze,
    caches,
    comparable,
    diploid_sources,
    extracted,
    findings,
    vcf,
    warning,
)

# Written IDs belong to one component only, including when the other component
# has a different allele placement. A non-rs label never supplies an rsID.
CASES = [
    pytest.param("AT", (".", "."), [(RS, "A")], RS, id="one-known"),
    pytest.param("AC", ("pos", "pos"), [(RS, "A")], RS, id="positional-ids"),
    pytest.param("AT", ("probe-a", "probe-t"), [(RS, "A")], RS, id="non-rs-labels"),
    pytest.param("AT", (RS, "."), [(RS, "A"), (OTHER, "T")], None, id="first-written"),
    pytest.param("AT", (".", OTHER), [(RS, "A"), (OTHER, "T")], None, id="second-written"),
    pytest.param(
        "AT", ("tag;RS880271201", "."), [(RS, "A"), (OTHER, "T")], None, id="canonical-written"
    ),
    pytest.param("AC", (".", "."), [(RS, "A"), (OTHER, "C")], None, id="two-resolved-control"),
    pytest.param("AT", (RS, OTHER), [(RS, "A"), (OTHER, "T")], None, id="two-written-control"),
    pytest.param("AC", (".", "."), [(RS, "A"), (RS, "C")], RS, id="same-rsid-control"),
    pytest.param("AT", (RS, "."), [(RS, "A")], RS, id="one-written-control"),
]


def _ids(ids, alts, build):
    return tuple(
        f"1:{POS[build]}:G:{a}" if x == "pos" else x for x, a in zip(ids, alts, strict=True)
    )


def _assert_counts(report, full):
    quality = report["input"]["vcf_quality"]
    assert report["input"]["total_variants"] == (3 if full else 4), (
        "Two unrelated flanks plus one full call or two assessed half calls"
    )
    assert quality["half_calls"] == 0
    assert quality["split_site_half_calls"] == (0 if full else 2)
    notices = warning(report, "split_site_calls")
    assert [w["count"] for w in notices] == ([] if full else [2])


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("auto", [False, True])
@pytest.mark.parametrize("alts,ids,placements,full_owner", CASES)
def test_per_alt_identity_findings_counts_and_annotated_round_trip(
    tmp_path, build, auto, alts, ids, placements, full_owner
):
    data = caches(tmp_path / "data", placements, placements, placements, placements)
    raw_ids = _ids(ids, alts, build)
    source = vcf(tmp_path / "input.vcf", build, alts, raw_ids, ("1/0", "0/1"))
    annotated = tmp_path / "annotated.vcf"
    mode = "auto" if auto else build
    report, _ = analyze(
        data, source, tmp_path / "report.json", mode, extra=("--vcf-out", annotated)
    )
    # Every shape has the unambiguous G>A risk row: it must not disappear.
    assert any(a["rsid"] == RS and a["alt"] == "A" for a in findings(report, "gwas")), (
        "The A association was present in 2.6.3; partial or component-local identity "
        "must retain it"
    )
    full = full_owner is not None
    _assert_counts(report, full)
    cv = findings(report, "clinvar")
    assert len(cv) == len(placements)
    if full:
        assert {a["rsid"] for a in cv} == {full_owner}
        assert {a["genotype_match"] for a in cv} == {"".join(sorted(alts))}
        assert {a["zygosity"] for a in cv} == {"heterozygous"}
    else:
        assert {(a["rsid"], a["alt"], a["genotype_match"], a["zygosity"]) for a in cv} == {
            (rid, alt, alt + "/.", "half_call") for rid, alt in placements
        }
        assert not findings(report, "snpedia") and not findings(report, "pharmgkb")
    # The writer must not spread one component's ID over both records.
    records = [
        line.split("\t") for line in annotated.read_text().splitlines() if not line.startswith("#")
    ]
    at_site = [row for row in records if row[1] == str(POS[build])]
    assert len(at_site) == 2
    expected_ids = (
        [full_owner, full_owner]
        if full
        else [dict((a, r) for r, a in placements)[a] for a in alts]
    )
    assert [row[2] for row in at_site] == expected_ids
    assert [row[4] for row in at_site] == list(alts)
    assert [row[-1].split(":")[0] for row in at_site] == ["1/0", "0/1"]
    for row, original in zip(at_site, raw_ids, strict=True):
        if row[2] != original:
            assert "ALLELIX_ORIGINAL_ID=" in row[7]
    reread, _ = analyze(data, annotated, tmp_path / "reread.json", mode)
    assert comparable(reread) == comparable(report), (
        "Findings, totals, quality and warnings must round-trip"
    )


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    "gts,phase",
    [
        (("1/0", "0/1"), False),
        (("1/.", "./1"), False),
        (("1|0", "0|1"), True),
        (("1|.", ".|1"), True),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_one_known_alt_keeps_the_full_pair_in_every_accepted_spelling(
    tmp_path, build, gts, phase, reverse
):
    pairs = [(RS, "A")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, pairs)
    alts, ids = "AT", (".", ".")
    if reverse:
        alts, gts = alts[::-1], gts[::-1]
    source = vcf(tmp_path / "input.vcf", build, alts, ids, gts, phase)
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    assert any(a["rsid"] == RS for a in findings(report, "gwas"))
    _assert_counts(report, True)
    assert {
        (a["rsid"], a["genotype_match"], a["zygosity"]) for a in findings(report, "clinvar")
    } == {(RS, "AT", "heterozygous")}


@pytest.mark.parametrize("build,indexed", [("GRCh37", False), ("GRCh38", False), ("GRCh38", True)])
@pytest.mark.parametrize("shape", ["one-known", "first-written", "second-written"])
def test_extract_uses_the_same_per_alt_decision(tmp_path, monkeypatch, build, indexed, shape):
    pairs = [(RS, "A")] if shape == "one-known" else [(RS, "A"), (OTHER, "T")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, pairs)
    ids = {"one-known": (".", "."), "first-written": (RS, "."), "second-written": (".", OTHER)}[
        shape
    ]
    source = vcf(tmp_path / "input.vcf", build, "AT", ids, ("1/0", "0/1"))
    wanted = [RS] if shape == "one-known" else [RS, OTHER]
    cells = extracted(data, source, build, wanted, monkeypatch, indexed=indexed)
    if shape == "one-known":
        assert cells[RS][3] in ("A/T", "T/A")
    else:
        assert cells[RS][3] in ("A/-", "-/A")
        assert cells[OTHER][3] in ("T/-", "-/T")


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("second_named", [False, True])
@pytest.mark.parametrize("second_known", [False, True])
def test_gnomad_only_identity_is_per_alt_and_remains_grch38_only(
    tmp_path, build, second_named, second_known
):
    second = [(OTHER, "T")] if second_known else []
    data = caches(tmp_path / "data", second, second, [(RS, "A"), *second], [(RS, None), *second])
    ids = (".", OTHER if second_named else ".")
    source = vcf(tmp_path / "input.vcf", build, "AT", ids, ("1/0", "0/1"))
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    two = build == "GRCh38" and (second_named or second_known)
    _assert_counts(report, not two)
    assert any(a["rsid"] == RS for a in findings(report, "gwas")) == (build == "GRCh38"), (
        "The risk-unspecified association exposes gnomAD-only rsID routing, "
        "without borrowing allele identity"
    )
    if second_known:
        assert any(a["rsid"] == OTHER for a in findings(report, "gwas"))


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_merged_multiallelic_record_and_unresolved_pair_are_unchanged(tmp_path, build):
    pairs = [(RS, "A"), (OTHER, "T")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, pairs)
    source = vcf(tmp_path / "input.vcf", build, ["A,T"], [RS], ["1/2"])
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    _assert_counts(report, True)
    assert {a["rsid"] for a in findings(report, "gwas")} == {RS}
    assert all(a["genotype_match"] == "AT" for a in findings(report, "clinvar"))
    # At an unrelated locus, no component resolves; it remains one full call.
    unresolved = vcf(tmp_path / "none.vcf", build, "AC", (".", "."), ("1/0", "0/1"))
    unresolved.write_text(
        unresolved.read_text().replace(f"\t{POS[build]}\t", f"\t{POS[build] + 10}\t")
    )
    none, _ = analyze(data, unresolved, tmp_path / "none.json", build)
    _assert_counts(none, True)
    assert not none["annotations"]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("two_owners", [False, True])
def test_diploid_sources_get_a_full_resolved_pair_but_never_its_component_hets(
    tmp_path, build, two_owners
):
    pairs = [(RS, "A"), (OTHER, "T")] if two_owners else [(RS, "A")]
    data = caches(tmp_path / "data", pairs, pairs, pairs, [])
    diploid_sources(data, RS)
    source = vcf(
        tmp_path / "input.vcf", build, "AT", (RS if two_owners else ".", "."), ("1/0", "0/1")
    )
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    if not two_owners:
        assert findings(report, "pharmgkb"), "The resolved full call reaches ClinPGx"
        # Eligibility is not a promise to match every source row: SNPedia's
        # palindromic-row guard also applies to the explicitly named control.
        named = vcf(tmp_path / "named.vcf", build, "AT", (RS, RS), ("1/0", "0/1"))
        control, _ = analyze(data, named, tmp_path / "named.json", build)
    for src in ("pharmgkb", "snpedia"):
        rows = findings(report, src)
        if two_owners:
            assert not rows, "Route S's half calls must not become AG or GT diploid calls"
        else:
            assert rows == findings(control, src)
            assert all(a["genotype_match"] == "AT" for a in rows)
    _assert_counts(report, not two_owners)

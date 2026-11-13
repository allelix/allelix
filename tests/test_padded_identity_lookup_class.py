"""Equivalent padded identity rows agree across coordinate and rsID consumers."""

import sqlite3
from contextlib import closing

import pytest
from test_annotation_parity_and_named_disclosures import _identity, _vcf
from test_padded_identity_source_parity import _input
from validated_call_helpers import (
    BUILDS,
    POSITIONS,
    RSID,
    analyze,
    reference_effect,
    setup,
    signatures,
)

from allelix.annotators.rsid_identity import RsidIdentityAnnotator

SHAPES = ((1, 0), (0, 1), (1, 1), (11, 0), (50, 1))
METHODS = (
    "forward37",
    "forward38",
    "alleles37",
    "alleles38",
    "placements",
    "iter",
    "loci",
    "reference_differing",
    "normal37",
    "normal38",
    "translation",
)


def _row(prefix, suffix):
    left, right = "C" * prefix, "T" * suffix
    return (
        RSID,
        "1",
        1000 - prefix,
        left + "A" + right,
        left + "G" + right,
        2000 - prefix,
        left + "G" + right,
        left + "A" + right,
        0,
    )


def _lookup(identity, method):
    if method.startswith("forward"):
        return identity.forward_sites({RSID}, "GRCh" + method[-2:])
    if method.startswith("alleles"):
        build = "GRCh" + method[-2:]
        pos = 1000 if build == "GRCh37" else 2000
        return identity.alleles_at_loci(build, {("1", pos)})
    if method == "placements":
        return identity.placements({RSID})
    if method == "iter":
        return list(identity.iter_placements([RSID]))
    if method == "loci":
        return identity.loci_of_rsids({RSID})
    if method == "reference_differing":
        return identity.reference_differing_sites()
    if method.startswith("normal"):
        build = "GRCh" + method[-2:]
        return identity.normal_alleles_on_other_build(
            {RSID: "A" if build == "GRCh37" else "G"}, build
        )
    assert method == "translation"
    return identity.translate_to_grch38({("1", 1000, "A", "G")})


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("method", METHODS)
def test_every_lookup_agrees_with_the_same_unpadded_identity_row(tmp_path, shape, method):
    plain, padded = tmp_path / "plain", tmp_path / "padded"
    plain.mkdir()
    padded.mkdir()
    _identity(plain, [_row(0, 0)])
    _identity(padded, [_row(*shape)])
    with RsidIdentityAnnotator(plain) as identity:
        expected = _lookup(identity, method)
        assert expected, "The equivalent unpadded row must decide this lookup"
    with RsidIdentityAnnotator(padded) as identity:
        assert _lookup(identity, method) == expected, (
            "Padding changes storage spelling, not the site or reference-strand decision"
        )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", (False, True))
@pytest.mark.parametrize("padding", (11, 50))
def test_source_findings_do_not_depend_on_ten_base_coordinate_window(
    tmp_path, build, named, padding
):
    data = tmp_path / "data"
    setup(data, build)
    source = _input(tmp_path / "input.vcf", build, "trimmed", named)
    before = analyze(data, source, build, ("gwas", "snpedia"))
    assert {a.source for a in before.annotations} == {"gwas", "snpedia"}
    (data / "rsid_identity.sqlite").unlink()
    prefix = "C" * padding
    _identity(
        data,
        [
            (
                RSID,
                "1",
                POSITIONS["GRCh37"] - padding,
                prefix + "G",
                prefix + "A",
                POSITIONS["GRCh38"] - padding,
                prefix + "G",
                prefix + "A",
                0,
            )
        ],
    )
    after = analyze(data, source, build, ("gwas", "snpedia"))
    assert signatures(after) == signatures(before), (
        "Recovery and findings must not stop at ten padding bases"
    )
    assert after.strand_abstentions == before.strand_abstentions


@pytest.mark.parametrize("ref,alt", (("CG", "CAT"), ("AC", "GT")))
def test_genuine_indel_and_mnv_placements_and_translation_stay_as_stored(tmp_path, ref, alt):
    _identity(tmp_path, [(RSID, "1", 1000, ref, alt, 2000, ref, alt, 0)])
    expected = ("1", 1000, ref, alt, 2000, ref, alt)
    with RsidIdentityAnnotator(tmp_path) as identity:
        assert identity.placements({RSID}) == {RSID: [expected]}
        assert list(identity.iter_placements([RSID])) == [(RSID, *expected)]
        assert identity.translate_to_grch38({("1", 1000, ref, alt)}) == {
            ("1", 1000, ref, alt): ("1", 2000, ref, alt)
        }


@pytest.mark.parametrize("palindromic", (False, True))
def test_translation_keeps_ambiguity_and_palindromic_guards(tmp_path, palindromic):
    row = list(_row(1, 0))
    row[-1] = int(palindromic)
    rows = [tuple(row)]
    if not palindromic:
        other = list(row)
        other[5] += 30
        rows.append(tuple(other))
    _identity(tmp_path, rows)
    with RsidIdentityAnnotator(tmp_path) as identity:
        assert identity.translate_to_grch38({("1", 1000, "A", "G")}) == {}


@pytest.mark.parametrize("build", BUILDS)
def test_padding_does_not_turn_a_normal_genotype_into_a_snpedia_finding(tmp_path, build):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    with closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.execute("UPDATE snpedia_genotypes SET allele1='A',allele2='A'")
    source = _vcf(tmp_path / "input.vcf", [(RSID, "1", pos, "G", "A", "1/1", ".")])
    before = analyze(data, source, build, ("snpedia",))
    assert not before.annotations, "The one-build normal-allele decision suppresses this row"
    with closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE rsid_identity SET pos37=pos37-1,pos38=pos38-1,"
            "ref37='C'||ref37,alt37='C'||alt37,ref38='C'||ref38,alt38='C'||alt38"
        )
    after = analyze(data, source, build, ("snpedia",))
    assert not after.annotations, (
        "Identity padding must not erase the normal allele and expose an AA row"
    )

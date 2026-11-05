"""O2: confident normal alleles are expressed on the input build's strand.

Synthetic real caches; original inverted multiallelic TT SNPedia seed retained.
No network, cache setup, shared data, or left-alignment is involved.
"""

import contextlib
import sqlite3

import pytest

from allelix.annotators import get_annotators
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.reference_fallback import ReferenceCarriedFallback
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.models import Variant
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from tests.test_annotation_parity_and_named_disclosures import _clinvar, _identity, _pharm, _stamp

RS = "rs900000001"
BUILDS = ("GRCh37", "GRCh38")


def identity(data, pairs):
    data.mkdir(exist_ok=True)
    _identity(
        data,
        [(RS, "1", 1000, r37, a37, 2000, r38, a38, flag) for r37, a37, r38, a38, flag in pairs],
    )


def gnomad(data, ref, alts):
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn, conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES ('1',2000,?,?,?,?)",
            [(ref, alt, RS, af) for alt, af in alts.items()],
        )
        _stamp(conn, "gnomad", len(alts), f"sv:{GNOMAD_SCHEMA_VERSION}")


def clinical(data, build, ref, alt):
    _clinvar(
        data,
        build,
        [
            (
                RS,
                "1",
                1000 if build == "GRCh37" else 2000,
                ref,
                alt,
                "TEST",
                "Pathogenic",
                "Pathogenic",
                1,
                1,
            )
        ],
    )


CASES = [
    # name, complete per-build identity, gnomAD frequencies, expected37, expected38
    ("inverted-alt", [("C", "A", "G", "T", 0)], {"T": 0.8}, "A", "T"),
    ("inverted-ref", [("C", "A", "G", "T", 0)], {"T": 0.1}, "C", "G"),
    (
        "outside-multiallelic",
        [("C", "A", "G", "T", 0), ("C", "T", "G", "A", 0)],
        {"T": 0.8, "A": 0.01},
        "A",
        "T",
    ),
    ("inverted-swapped", [("C", "A", "T", "G", 0)], {"G": 0.1}, "A", "T"),
    ("same-strand-swapped", [("C", "A", "A", "C", 0)], {"C": 0.1}, "A", "A"),
    ("tie", [("C", "A", "G", "T", 0)], {"T": 0.5}, None, None),
    ("unknown-af", [("C", "A", "G", "T", 0)], {"T": None}, None, None),
]


@pytest.mark.parametrize("case", CASES, ids=[c[0] for c in CASES])
@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("prefetch", [False, True])
def test_gnomad_normal_is_in_the_input_build(tmp_path, case, build, prefetch):
    _, pairs, frequencies, expected37, expected38 = case
    identity(tmp_path, pairs)
    gnomad(tmp_path, pairs[0][2], frequencies)
    pos = 1000 if build == "GRCh37" else 2000
    expected = expected37 if build == "GRCh37" else expected38
    with contextlib.ExitStack() as stack:
        cv = stack.enter_context(ClinVarAnnotator(tmp_path))
        ident = stack.enter_context(RsidIdentityAnnotator(tmp_path))
        gn = stack.enter_context(GnomadAnnotator(tmp_path))
        provider = ReferenceCarriedFallback(cv, ident, gn)
        if prefetch:
            provider.prefetch([Variant(RS, "1", pos, "A", "A", build=build)])
        assert provider.normal_allele_for(RS, build, ("1", pos)) == expected
        assert provider.reference_at(RS, build, ("1", pos)) == expected
        assert provider.normal_allele_for(RS, build, ("1", pos + 1)) is None
        assert provider.abstentions() == int(expected is None)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("inverted", [False, True])
@pytest.mark.parametrize("provider_kind", ["identity-only", "both-clinvar", "other-clinvar"])
@pytest.mark.parametrize("prefetch", [False, True])
def test_other_build_clinvar_and_identity_paths_translate(
    tmp_path, build, inverted, provider_kind, prefetch
):
    # Inversion with a real reference change, unlike C>A / G>T where REF is the same allele.
    pair = ("C", "A", "T", "G", 0) if inverted else ("C", "A", "A", "C", 0)
    identity(tmp_path, [pair])
    ref = pair[0] if build == "GRCh37" else pair[2]
    other = "GRCh38" if build == "GRCh37" else "GRCh37"
    other_ref = pair[2] if build == "GRCh37" else pair[0]
    other_alt = pair[3] if build == "GRCh37" else pair[1]
    expected = pair[1] if build == "GRCh37" else pair[3]
    if provider_kind != "other-clinvar":
        clinical(tmp_path, build, ref, ref)
    if provider_kind != "identity-only":
        clinical(tmp_path, other, other_ref, other_alt)
    with contextlib.ExitStack() as stack:
        ident = stack.enter_context(RsidIdentityAnnotator(tmp_path))
        cv = stack.enter_context(
            ClinVarAnnotator(tmp_path, other_build_refs=ident.normal_alleles_on_other_build)
        )
        gn = stack.enter_context(GnomadAnnotator(tmp_path))
        provider = ReferenceCarriedFallback(cv, ident, gn)
        pos = 1000 if build == "GRCh37" else 2000
        if prefetch:
            provider.prefetch([Variant(RS, "1", pos, ref, ref, build=build)])
        assert provider.normal_allele_for(RS, build, ("1", pos)) == expected
        assert provider.reference_at(RS, build, ("1", pos)) == expected


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("other_alt", ["A", "T"])
def test_same_variation_id_without_a_valid_allele_relation_abstains(tmp_path, build, other_alt):
    # Carried C: the matching other ALT must be C or G, never A or T.
    other = "GRCh38" if build == "GRCh37" else "GRCh37"
    clinical(tmp_path, build, "C", "C")
    clinical(tmp_path, other, "G", other_alt)
    with ClinVarAnnotator(tmp_path) as cv:
        assert cv.normal_allele_for(RS, build) is None


def sources(data):
    _pharm(data, RS, "TEST", ["TT", "AA"])
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        conn.executemany(
            "INSERT INTO snpedia_genotypes "
            "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
            "VALUES (?,?,?,6,'Bad','Synthetic','TEST','plus')",
            [(RS, a, a) for a in "TA"],
        )
        _stamp(conn, "snpedia", 2, f"pv:{_PARSER_VERSION}")


@pytest.mark.parametrize("source", ["snpedia", "pharmgkb"])
@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("clinvar_state", ["absent", "disabled-present"])
@pytest.mark.parametrize("path", ["single", "batch", "pipeline"])
def test_reference_fallback_does_not_hide_the_wrong_homozygote(
    tmp_path, source, build, clinvar_state, path
):
    identity(tmp_path, [("C", "A", "G", "T", 0), ("C", "T", "G", "A", 0)])
    gnomad(tmp_path, "G", {"T": 0.8, "A": 0.01})
    sources(tmp_path)
    if clinvar_state == "disabled-present":
        for b in BUILDS:
            clinical(tmp_path, b, "T", "A")
    off = {"clinvar"} if clinvar_state == "disabled-present" else set()
    pos, ref, alt = (1000, "C", "A,T") if build == "GRCh37" else (2000, "G", "T,A")
    with contextlib.ExitStack() as stack:
        registry = get_annotators(tmp_path, switched_off=off)
        for item in registry:
            stack.enter_context(item)
        annotator = next(a for a in registry if a.name == source)
        # T is normal only on GRCh38; A is normal only on GRCh37.
        for allele, expected in [("T", build == "GRCh37"), ("A", build == "GRCh38")]:
            call = Variant(RS, "1", pos, allele, allele, build=build, ref=ref, ref_from_input=True)
            if path == "single":
                found = annotator.annotate(call)
            elif path == "batch":
                found = list(annotator.batch_annotate([call]))
            else:
                src = tmp_path / "input.vcf"
                gt = str(alt.split(",").index(allele) + 1)
                src.write_text(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                    f"1\t{pos}\t{RS}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}/{gt}\n"
                )
                result = run_analysis(src, VcfParser(), [annotator], build_override=build)
                found = result.annotations
            assert len(found) == int(expected), (build, allele, source, path)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["named", "unnamed", "padded-named", "padded-unnamed", "gvcf"])
@pytest.mark.parametrize("clinvar_state", ["absent", "disabled-present"])
def test_outside_seed_survives_resolution_and_padding(tmp_path, build, shape, clinvar_state):
    identity(tmp_path, [("C", "A", "G", "T", 0), ("C", "T", "G", "A", 0)])
    gnomad(tmp_path, "G", {"T": 0.8, "A": 0.01})
    sources(tmp_path)
    if clinvar_state == "disabled-present":
        for b in BUILDS:
            clinical(tmp_path, b, "T", "A")
    off = {"clinvar"} if clinvar_state == "disabled-present" else set()
    pos, ref, alts = (1000, "C", ["A", "T"]) if build == "GRCh37" else (2000, "G", ["T", "A"])
    gt = alts.index("T") + 1
    if shape.startswith("padded"):
        pos -= 1
        ref = "G" + ref
        alts = ["G" + a for a in alts]
    if shape == "gvcf":
        alts.append("<NON_REF>")
    rsid = "." if "unnamed" in shape else RS
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t{pos}\t{rsid}\t{ref}\t{','.join(alts)}\t.\tPASS\t.\tGT\t{gt}/{gt}\n"
    )
    registry = get_annotators(tmp_path, switched_off=off)
    by_name = {a.name: a for a in registry}
    result = run_analysis(
        src,
        VcfParser(),
        [by_name["snpedia"]],
        identity=by_name["rsid_identity"],
        resource_owners=registry,
        build_override=build,
    )
    assert len(result.annotations) == int(build == "GRCh37")
    if result.annotations:
        assert result.annotations[0].rsid == RS


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("disabled", [False, True])
def test_no_frequency_evidence_never_suppresses_a_call(tmp_path, build, disabled):
    identity(tmp_path, [("C", "A", "G", "T", 0), ("C", "T", "G", "A", 0)])
    sources(tmp_path)
    if disabled:
        gnomad(tmp_path, "G", {"T": 0.8, "A": 0.01})
    with contextlib.ExitStack() as stack:
        registry = get_annotators(tmp_path, switched_off={"gnomad"} if disabled else set())
        for a in registry:
            stack.enter_context(a)
        source = next(a for a in registry if a.name == "snpedia")
        pos, ref = (1000, "C") if build == "GRCh37" else (2000, "G")
        call = Variant(RS, "1", pos, "T", "T", build=build, ref=ref, ref_from_input=True)
        assert len(source.annotate(call)) == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("ambiguity", ["palindromic", "disagree"])
def test_identity_only_normal_mapping_abstains_on_ambiguity(tmp_path, build, ambiguity):
    pairs = (
        [("T", "C", "C", "T", 1)]
        if ambiguity == "palindromic"
        else [("C", "A", "T", "G", 0), ("C", "G", "A", "C", 0)]
    )
    identity(tmp_path, pairs)
    carried = pairs[0][0] if build == "GRCh37" else pairs[0][2]
    with RsidIdentityAnnotator(tmp_path) as ident:
        assert ident.normal_alleles_on_other_build({RS: carried}, build) == {}


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("prefetch", [False, True])
def test_palindromic_gnomad_translation_is_not_a_confident_cross_build_normal(
    tmp_path, build, prefetch
):
    identity(tmp_path, [("T", "C", "C", "T", 1)])
    gnomad(tmp_path, "C", {"T": 0.01})
    with contextlib.ExitStack() as stack:
        cv = stack.enter_context(ClinVarAnnotator(tmp_path))
        ident = stack.enter_context(RsidIdentityAnnotator(tmp_path))
        gn = stack.enter_context(GnomadAnnotator(tmp_path))
        provider = ReferenceCarriedFallback(cv, ident, gn)
        pos = 1000 if build == "GRCh37" else 2000
        if prefetch:
            provider.prefetch([Variant(RS, "1", pos, "T", "T", build=build)])
        # gnomAD is already GRCh38; only the cross-build reading needs translation.
        expected = None if build == "GRCh37" else "C"
        assert provider.normal_allele_for(RS, build, ("1", pos)) == expected
        assert provider.reference_at(RS, build, ("1", pos)) == expected


def f5(data, both):
    from tests.test_annotation_parity_and_named_disclosures import _f5

    _f5(data, BUILDS if both else ("GRCh37",))
    # Real F5 shape: its identity rows are flagged; the same VariationID in
    # both ClinVar builds independently establishes the normal C allele.
    with contextlib.closing(sqlite3.connect(data / "rsid_identity.sqlite")) as conn, conn:
        conn.execute("UPDATE rsid_identity SET palindromic = 1")
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn, conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        conn.executemany(
            "INSERT INTO snpedia_genotypes "
            "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
            "VALUES ('rs6025',?,?,6,'Bad','Synthetic F5 genotype','F5','plus')",
            [(a, a) for a in "CT"],
        )
        _stamp(conn, "snpedia", 2, f"pv:{_PARSER_VERSION}")


@pytest.mark.parametrize("both", [False, True], ids=["one-build", "both-builds"])
@pytest.mark.parametrize("allele", ["T", "C"])
@pytest.mark.parametrize("path", ["single", "batch", "pipeline-json"])
def test_f5_one_build_abstention_preserves_carrier_and_shows_uncertain_normal_rows(
    tmp_path, both, allele, path
):
    data = tmp_path / "data"
    f5(data, both)
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data)
        for item in registry:
            stack.enter_context(item)
        by_name = {a.name: a for a in registry}
        call = Variant(
            "rs6025", "1", 169519049, allele, allele, build="GRCh37", ref="T", ref_from_input=True
        )
        selected = [by_name[n] for n in ("clinvar", "snpedia", "pharmgkb")]
        if path == "single":
            rows = [row for a in selected for row in a.annotate(call)]
        elif path == "batch":
            rows = [row for a in selected for row in a.batch_annotate([call])]
        else:
            import json

            from allelix.reports.json_report import render_json

            src = tmp_path / "f5.vcf"
            gt = "0/0" if allele == "T" else "1/1"
            src.write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                f"1\t169519049\trs6025\tT\tC\t.\tPASS\t.\tGT\t{gt}\n"
            )
            result = run_analysis(
                src,
                VcfParser(),
                selected,
                identity=by_name["rsid_identity"],
                build_override="GRCh37",
            )
            rows = result.annotations
            out = tmp_path / "report.json"
            render_json(
                result,
                output_path=out,
                min_magnitude=0,
                source_min_magnitudes={"snpedia": 0, "pharmgkb": 0},
            )
            public = json.loads(out.read_text())["annotations"]
            assert sorted(a["source"] for a in public) == sorted(a.source for a in rows)
        counts = {
            name: sum(a.source == name for a in rows)
            for name in ("clinvar", "snpedia", "pharmgkb")
        }
        assert counts["clinvar"] == int(allele == "T"), "ClinVar's F5 carrier assessment stays"
        # Without the other ClinVar build the flagged identity mapping abstains.
        # C/C therefore shows the cache's C/C genotype rows, not a ClinVar P/LP finding.
        expected = int(allele == "T" or not both)
        assert counts["snpedia"] == expected
        assert counts["pharmgkb"] == expected
        assert by_name["clinvar"].normal_allele_for("rs6025", "GRCh37") == ("C" if both else None)


@pytest.mark.parametrize("build", BUILDS)
def test_disabled_clinvar_cache_does_not_change_translated_fallback(tmp_path, build):
    import json

    from click.testing import CliRunner

    from allelix.cli import main

    data = tmp_path / "data"
    identity(data, [("C", "A", "G", "T", 0), ("C", "T", "G", "A", 0)])
    gnomad(data, "G", {"T": 0.8, "A": 0.01})
    sources(data)
    for b in BUILDS:
        clinical(data, b, "T", "A")
    (data / "config.toml").write_text("[sources]\nclinvar = false\n")
    pos, ref, alt, gt = (
        (1000, "C", "A,T", "2/2") if build == "GRCh37" else (2000, "G", "T,A", "1/1")
    )
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t{pos}\t{RS}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
    )
    out = tmp_path / "output.json"
    args = [
        "analyze",
        str(src),
        "--build",
        build.lower(),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(out),
        "--min-magnitude",
        "0",
        "--snpedia-min-magnitude",
        "0",
    ]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    rows = json.loads(out.read_text())["annotations"]
    assert not any(a["source"] == "clinvar" for a in rows)
    for source in ("snpedia", "pharmgkb"):
        assert sum(a["source"] == source for a in rows) == int(build == "GRCh37")

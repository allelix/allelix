# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar variants whose allele a build's reference carries (F5 rs6025).

On GRCh37 the reference base at rs6025 is the Factor V Leiden allele, so
ClinVar writes the real variant (c.1601G>A) as T>T there and the
normal-allele record (c.1601=) as T>C. On GRCh38 they are C>T and C>C.
The same genotype must give the same findings on both builds, on the
array and VCF paths, for ClinVar, ClinPGx and SNPedia.

The ClinVar caches are built through the TSV loader from synthetic rows,
so the loader's record classification is part of what is tested.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.manager import (
    _SS_REQUIRED_COLUMNS,
    _VS_REQUIRED_COLUMNS,
    load_clinvar_tsv,
)
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import PHARMGKB_SCHEMA
from allelix.models import Annotation, AnnotationInvariantError, Variant
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

if TYPE_CHECKING:
    from pathlib import Path

GRCH37 = "GRCh37"
GRCH38 = "GRCh38"
BUILDS = (GRCH37, GRCH38)
F5 = "rs6025"
POS = {GRCH37: 169_519_049, GRCH38: 169_549_811}
REF = {GRCH37: "T", GRCH38: "C"}

# (VariationID, HGVS name, {build: (ref, alt)}, [per-SCV significance])
_F5_RECORDS = (
    (
        "642",
        "NM_000130.4(F5):c.1601G>A (p.Arg534Gln)",
        {GRCH37: ("T", "T"), GRCH38: ("C", "T")},
        ["drug response"],
    ),
    (
        "226007",
        "NM_000130.5(F5):c.1601=",
        {GRCH37: ("T", "C"), GRCH38: ("C", "C")},
        ["Pathogenic", "Benign", "Benign", "Pathogenic", "Benign"],
    ),
)
# A substitution REF == ALT on both builds (like LRRK2 rs2256408): the
# reference carries the variant allele on both, and no normal allele is known.
_BOTH_BUILDS_RECORD = (
    "39134",
    "NM_198578.3(LRRK2):c.149A>G (p.His50Arg)",
    {GRCH37: ("G", "G"), GRCH38: ("G", "G")},
    ["Benign"],
)
_BOTH_POS = {GRCH37: 40_619_082, GRCH38: 40_225_280}
# A normal-allele record with a risk-factor classification (shaped like
# CASP12 rs497116): written G>G on GRCh37, whose reference is the normal
# allele, and A>G on GRCh38, whose reference carries the other allele.
NORMAL_RSID = "rs900231001"
NORMAL_POS = {GRCH37: 104_761_000, GRCH38: 104_890_000}
_NORMAL_RECORD = (
    "900231",
    "NM_000000.1(CASP12):c.100=",
    {GRCH37: ("G", "G"), GRCH38: ("A", "G")},
    ["risk factor", "Pathogenic"],
)

# What each genotype must give, on every build and path (default filters).
EXPECTED = {
    "CC": [],
    "CT": [("clinvar_drug_response", "T", "Heterozygous")],
    "TT": [("clinvar_drug_response", "T", "Homozygous")],
}


def _vs_line(variation_id: str, name: str, build: str, ref: str, alt: str, rs: str, pos: int):
    cols = [""] * 43
    cols[0] = variation_id
    cols[1] = "single nucleotide variant"
    cols[2] = name
    cols[4] = name.split("(")[1].split(")")[0]
    cols[9] = rs.removeprefix("rs")
    cols[16] = build
    cols[18] = "1" if rs == F5 else "12"
    cols[24] = "reviewed by expert panel"
    cols[30] = variation_id
    cols[31] = str(pos)
    cols[32] = ref
    cols[33] = alt
    return "\t".join(cols)


def _ss_line(variation_id: str, significance: str, scv: int) -> str:
    cols = [""] * 16
    cols[0] = variation_id
    cols[1] = significance
    cols[5] = "C0015499:Factor V deficiency"
    cols[6] = "criteria provided, single submitter"
    cols[10] = f"SCV{scv:09d}.1"
    cols[15] = "yes"
    return "\t".join(cols)


def _header(columns: dict[int, str], width: int) -> str:
    """A header row carrying the names the loader checks, placeholders elsewhere."""
    return "#" + "\t".join(columns.get(i, f"col{i}") for i in range(width)) + "\n"


def _build_clinvar(
    data: Path,
    builds: tuple[str, ...],
    *,
    both_builds_record: bool = False,
    normal_allele_record: bool = False,
):
    records = [(*r, F5, POS) for r in _F5_RECORDS]
    if both_builds_record:
        records.append((*_BOTH_BUILDS_RECORD, "rs2256408", _BOTH_POS))
    if normal_allele_record:
        records.append((*_NORMAL_RECORD, NORMAL_RSID, NORMAL_POS))
    vs = data / "vs.tsv"
    ss = data / "ss.tsv"
    vs.write_text(
        _header(_VS_REQUIRED_COLUMNS, 43)
        + "\n".join(
            _vs_line(vid, name, b, *alleles[b], rs, pos[b])
            for vid, name, alleles, _sigs, rs, pos in records
            for b in BUILDS
        )
        + "\n"
    )
    scv = iter(range(1, 1000))
    ss.write_text(
        _header(_SS_REQUIRED_COLUMNS, 16)
        + "\n".join(
            _ss_line(vid, sig, next(scv))
            # submission_summary is sorted by VariationID, as ClinVar ships it.
            for vid, _n, _a, sigs, _r, _p in sorted(records, key=lambda r: int(r[0]))
            for sig in sigs
        )
        + "\n"
    )
    for build in builds:
        db = data / clinvar_db_filename(build)
        load_clinvar_tsv(vs, ss, db, build, source_url="test://f5", remote_signal="md5:f5")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute(
                "UPDATE database_versions SET name = ? WHERE name = 'clinvar'",
                (clinvar_record_name(build),),
            )
            conn.commit()


def _build_identity(data: Path) -> None:
    write_identity_db(
        [
            (F5, "1", POS[GRCH37], "T", "C", POS[GRCH38], "C", "C", 0),
            (F5, "1", POS[GRCH37], "T", "T", POS[GRCH38], "C", "T", 0),
        ],
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://f5",
        downloaded_at="2026-10-13T00:00:00Z",
    )


def test_reference_only_factor_v_leiden_cli(f5_data: Path, tmp_path: Path) -> None:
    import json

    from click.testing import CliRunner

    from allelix.cli import main

    path = tmp_path / "f5.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{POS[GRCH37]}\t{F5}\tT\tC\t.\tPASS\t.\tGT:GQ\t0/0:99\n"
    )
    out = tmp_path / "f5.json"
    result = CliRunner().invoke(
        main,
        ["analyze", str(path), "--no-update", "--data-dir", str(f5_data), "--output", str(out)],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out.read_text())
    assert any(a["rsid"] == F5 and a["genotype_match"] == "TT" for a in payload["annotations"])
    assert "This file has no non-reference calls." in result.output


def test_factor_v_leiden_reference_call_is_stamped_in_vcf_out(
    f5_data: Path, tmp_path: Path
) -> None:
    """GRCh37 0/0 at rs6025 carries the Leiden finding in --vcf-out, as in JSON."""
    from click.testing import CliRunner

    from allelix.cli import main

    pysam = pytest.importorskip("pysam")
    path = tmp_path / "f5.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{POS[GRCH37]}\t{F5}\tT\tC\t.\tPASS\t.\tGT:GQ\t0/0:99\n"
    )
    out = tmp_path / "f5.annotated.vcf"
    result = CliRunner().invoke(
        main,
        ["analyze", str(path), "--no-update", "--data-dir", str(f5_data), "--vcf-out", str(out)],
    )
    assert result.exit_code == 0, result.output
    assert "1 stamped row(s)" in result.output
    with pysam.VariantFile(str(out)) as stream:
        record = next(stream)
        assert (record.ref, record.alts) == ("T", ("C",))
        assert record.samples["SAMPLE"]["GT"] == (0, 0)
        (value,) = record.info["ALLELIX_CLINVAR"]
        assert value.split("|")[1] == F5


def _build_pharmgkb(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.executescript(PHARMGKB_SCHEMA)
        conn.executemany(
            "INSERT INTO pharmgkb_annotations (rsid, genotype, gene, drugs, phenotype, "
            "phenotype_category, annotation_text, level_of_evidence, score, "
            "pgkb_annotation_id, allele_function, function_class, is_nonfinding) "
            "VALUES (?, ?, 'F5', 'hormonal contraceptives', 'Thrombosis', 'Toxicity', ?, "
            "'2B', 1.0, ?, '', '', 0)",
            [
                (F5, geno, f"Patients with the {F5} {geno} genotype", f"PA{n}")
                for n, geno in enumerate(("CC", "CT", "TT"))
            ],
        )
        conn.commit()


def _build_snpedia(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.execute(
            "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
            "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
            "orientation TEXT, scraped_at TEXT)"
        )
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, ?, ?, 'F5', 'plus', '2026-05-20')",
            [
                (F5, "C", "C", 0.0, "Good", "normal risk of thrombosis"),
                (F5, "C", "T", 4.1, "Bad", "higher risk of thrombosis"),
                (F5, "T", "T", 7.0, "Bad", "much higher risk of thrombosis"),
            ],
        )
        conn.commit()


@pytest.fixture
def f5_data(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    _build_clinvar(data, BUILDS, both_builds_record=True)
    _build_identity(data)
    _build_pharmgkb(data)
    _build_snpedia(data)
    return data


def _findings(annotations: list[Annotation], source: str = "clinvar") -> list[tuple]:
    return sorted((a.significance, a.alt, a.zygosity) for a in annotations if a.source == source)


# ── The loader keeps exactly the rows the rule needs ─────────────


def test_cache_rows(f5_data: Path) -> None:
    rows = {}
    for build in BUILDS:
        with contextlib.closing(sqlite3.connect(f5_data / clinvar_db_filename(build))) as conn:
            rows[build] = sorted(
                conn.execute(
                    "SELECT allele_id, ref, alt, clinical_significance FROM clinvar_variants "
                    "WHERE rsid = ?",
                    (F5,),
                )
            )
    # GRCh37: the real variant T>T is kept; the normal-allele record keeps
    # its Benign rows and loses its Pathogenic ones.
    assert rows[GRCH37] == [
        (642, "T", "T", "drug response"),
        (226007, "T", "C", "Benign"),
        (226007, "T", "C", "Benign"),
        (226007, "T", "C", "Benign"),
    ]
    # GRCh38: the normal-allele record is written REF == ALT; since iv:9 it
    # is kept like on GRCh37, with its Benign rows and without the P/LP ones.
    assert rows[GRCH38] == [
        (642, "C", "T", "drug response"),
        (226007, "C", "C", "Benign"),
        (226007, "C", "C", "Benign"),
        (226007, "C", "C", "Benign"),
    ]


# ── ClinVar: the same genotype gives the same findings on both builds ──


def _chip(path: Path, build: str, genotype: str, rsid: str = F5) -> Path:
    path.write_text(
        "# MyHappyGenes [TEMPUS]\n"
        "# Sample ID\tPROBE000003\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        f"{rsid}\t1\t{POS[build]}\t{genotype[0]}\t{genotype[1]}\n"
    )
    return path


@pytest.mark.parametrize("with_identity", [False, True])
@pytest.mark.parametrize("genotype", ["CC", "CT", "TT"])
@pytest.mark.parametrize("build", BUILDS)
def test_array_call(
    f5_data: Path,
    tmp_path: Path,
    cm_stack: contextlib.ExitStack,
    build: str,
    genotype: str,
    with_identity: bool,
) -> None:
    identity = cm_stack.enter_context(RsidIdentityAnnotator(f5_data)) if with_identity else None
    result = run_analysis(
        _chip(tmp_path / "chip.txt", build, genotype),
        MyHappyGenesParser(),
        [cm_stack.enter_context(ClinVarAnnotator(f5_data))],
        build_override=build,
        identity=identity,
    )
    assert _findings(result.annotations) == EXPECTED[genotype]


_GT = {
    GRCH37: {"CC": ("C", "1/1"), "CT": ("C", "0/1"), "TT": ("C", "0/0")},
    GRCH38: {"CC": ("T", "0/0"), "CT": ("T", "0/1"), "TT": ("T", "1/1")},
}


@pytest.mark.parametrize("genotype", ["CC", "CT", "TT"])
@pytest.mark.parametrize("build", BUILDS)
def test_vcf_record(
    f5_data: Path, tmp_path: Path, cm_stack: contextlib.ExitStack, build: str, genotype: str
) -> None:
    alt, gt = _GT[build][genotype]
    vcf = tmp_path / "f5.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"1\t{POS[build]}\t{F5}\t{REF[build]}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n"
    )
    result = run_analysis(
        vcf,
        VcfParser(),
        [cm_stack.enter_context(ClinVarAnnotator(f5_data))],
        build_override=build,
    )
    # TT is the GRCh37 reference here, so its record is a 0/0 call: the
    # pipeline watches the reference-carried site and assesses the call as TT.
    assert _findings(result.annotations) == EXPECTED[genotype]


def test_normal_allele_record_never_reported_as_pathogenic(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    """Benign rows are hidden by default; with them included, a C carrier on
    GRCh37 still sees no Pathogenic row from the normal-allele record."""
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data, include_benign=True))
    for build in BUILDS:
        v = Variant(F5, "1", POS[build], "C", "C", build=build, ref=REF[build])
        assert not [a for a in clinvar.annotate(v) if "pathogenic" in a.significance]


# ── ClinPGx and SNPedia: homozygous reference means the normal allele ──


def test_reference_for_returns_the_normal_allele(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    assert clinvar.reference_for(F5, GRCH37) == "C"
    assert clinvar.reference_for(F5, GRCH38) == "C"


# ClinPGx and SNPedia: what each genotype must give, on both builds. The alt
# is the Leiden allele (T) whichever build's reference carries it.
EXPECTED_GENOTYPE_ROWS = {
    "CC": {"pharmgkb": [], "snpedia": []},
    "CT": {
        "pharmgkb": [("pharmgkb_loe_2b", "T", "Heterozygous")],
        "snpedia": [("snpedia_bad", "T", "Heterozygous")],
    },
    "TT": {
        "pharmgkb": [("pharmgkb_loe_2b", "T", "Homozygous")],
        "snpedia": [("snpedia_bad", "T", "Homozygous")],
    },
}


@pytest.mark.parametrize("genotype", ["CC", "CT", "TT"])
@pytest.mark.parametrize("build", BUILDS)
def test_clinpgx_and_snpedia_parity(
    f5_data: Path, cm_stack: contextlib.ExitStack, build: str, genotype: str
) -> None:
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    pharmgkb = cm_stack.enter_context(
        PharmGKBAnnotator(
            f5_data,
            clinvar_ref_provider=clinvar.reference_at,
            normal_allele_provider=clinvar.normal_allele_for,
        )
    )
    snpedia = cm_stack.enter_context(
        SNPediaAnnotator(
            f5_data,
            clinvar_ref_provider=clinvar.reference_at,
            normal_allele_provider=clinvar.normal_allele_for,
        )
    )
    v = Variant(F5, "1", POS[build], genotype[0], genotype[1], build=build, ref=REF[build])
    for source, annotator in (("pharmgkb", pharmgkb), ("snpedia", snpedia)):
        single = annotator.annotate(v)
        assert _findings(single, source) == EXPECTED_GENOTYPE_ROWS[genotype][source]
        assert _findings(list(annotator.batch_annotate([v])), source) == _findings(single, source)
        for a in single:
            # The build's reference base is kept as locus_ref (the enrichment
            # key); the exemption is set exactly when alt is that base.
            assert a.locus_ref == REF[build]
            assert a.reference_carries_alt == (a.alt == REF[build])


# ── Fallback: no normal allele known, nothing is suppressed ──────


def test_one_build_install_suppresses_nothing(
    tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    data = tmp_path / "one"
    data.mkdir()
    _build_clinvar(data, (GRCH37,))
    _build_pharmgkb(data)
    clinvar = cm_stack.enter_context(ClinVarAnnotator(data, builds=(GRCH37,)))
    assert clinvar.reference_for(F5, GRCH37) is None
    pharmgkb = cm_stack.enter_context(
        PharmGKBAnnotator(
            data,
            clinvar_ref_provider=clinvar.reference_at,
            normal_allele_provider=clinvar.normal_allele_for,
        )
    )
    for genotype in ("CC", "TT"):
        v = Variant(F5, "1", POS[GRCH37], genotype[0], genotype[1], build=GRCH37, ref="T")
        assert [a.genotype_match for a in pharmgkb.annotate(v)] == [genotype]
    # ClinVar itself does not need the other build.
    v = Variant(F5, "1", POS[GRCH37], "T", "T", build=GRCH37, ref="T")
    assert _findings(clinvar.annotate(v)) == EXPECTED["TT"]


@pytest.mark.parametrize("build", BUILDS)
def test_ref_equal_alt_on_both_builds(
    f5_data: Path, cm_stack: contextlib.ExitStack, build: str
) -> None:
    """No normal allele is known: no REF for the non-finding check, and the
    ClinVar row is reported to carriers of its allele like any other."""
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data, include_benign=True))
    assert clinvar.reference_for("rs2256408", build) is None
    v = Variant("rs2256408", "12", _BOTH_POS[build], "G", "G", build=build, ref="G")
    assert _findings(clinvar.annotate(v)) == [("clinvar_benign", "G", "Homozygous")]
    v = Variant("rs2256408", "12", _BOTH_POS[build], "A", "A", build=build)
    assert clinvar.annotate(v) == []


# ── The invariant exception is narrow ───────────────────────────


def test_invariant_still_rejects_alt_equal_to_reference() -> None:
    kwargs = {
        "source": "clinvar",
        "rsid": F5,
        "significance": "clinvar_drug_response",
        "category": "clinical",
        "magnitude": 5.0,
        "description": "",
        "attribution": "ClinVar",
        "alt": "T",
        "locus_ref": "T",
        "carried_alleles": ("T", "T"),
        "genotype_match": "TT",
    }
    with pytest.raises(AnnotationInvariantError):
        Annotation(**kwargs)
    Annotation(**kwargs, reference_carries_alt=True)


def test_exception_marks_only_ref_equal_alt_rows(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data, include_benign=True))
    seen = []
    for build in BUILDS:
        for genotype in ("CC", "CT", "TT"):
            v = Variant(F5, "1", POS[build], genotype[0], genotype[1], build=build, ref=REF[build])
            seen.extend(clinvar.annotate(v))
    assert any(a.reference_carries_alt for a in seen)
    assert any(not a.reference_carries_alt for a in seen)
    for a in seen:
        assert a.reference_carries_alt == (a.alt == a.locus_ref)


def test_old_cache_keeps_the_locus_ref(f5_data: Path, cm_stack: contextlib.ExitStack) -> None:
    """Until `db update` rebuilds an iv:4 cache, it still holds normal-allele
    rows written REF == ALT (226007 C>C on GRCh38), so REF == ALT does not
    mean the reference carries a variant. ``reference_for`` then gives the
    locus REF, as 2.4.1 did, and never borrows from the other build."""
    with contextlib.closing(sqlite3.connect(f5_data / clinvar_db_filename(GRCH38))) as conn:
        conn.execute(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, '1', ?, 'C', 'C', 'Pathogenic', 'F5-related disorders', 'F5', '', 226007)",
            (F5, POS[GRCH38]),
        )
        conn.commit()
    for build in BUILDS:
        with contextlib.closing(sqlite3.connect(f5_data / clinvar_db_filename(build))) as conn:
            conn.execute("UPDATE database_versions SET local_version_tag = 'iv:4'")
            conn.commit()
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    assert not clinvar.is_ready()
    assert clinvar.reference_for(F5, GRCH37) == "T"
    assert clinvar.reference_for(F5, GRCH38) == "C"


def test_reference_carried_row_gets_the_grch38_enrichment_key(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    """The GRCh37 key (T, T) of a reference-carried row translates to the real
    GRCh38 variant (C>T), so gnomAD, AlphaMissense and CADD are looked up
    for it as on GRCh38. A row whose alt is the reference base for any other
    reason still gets no key. (The real rs6025 identity rows are flagged
    palindromic, all four bases, and abstain under ADR-0040; this fixture's
    rows are not.)"""
    from allelix.reports._pipeline import _enrichment_keys

    identity = cm_stack.enter_context(RsidIdentityAnnotator(f5_data))
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    v = Variant(F5, "1", POS[GRCH37], "T", "T", build=GRCH37, ref="T")
    (row,) = clinvar.annotate(v)
    assert row.reference_carries_alt
    assert _enrichment_keys([row], GRCH37, identity) == [(row, ("1", POS[GRCH38], "C", "T"))]
    plain = Annotation(
        source="gwas",
        rsid=F5,
        significance="gwas_association",
        category="trait",
        magnitude=1.0,
        description="",
        attribution="GWAS Catalog",
        genotype_match="TT",
        alt="T",
        locus_ref="T",
        locus=("1", POS[GRCH37]),
    )
    assert _enrichment_keys([plain], GRCH37, identity) == []


def test_differing_clinvar_ref_without_reference_carried_row_keeps_variant_ref(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    """Only a reference-carried rsID takes its normal allele from the other
    build. Where ClinVar's REF merely differs from ``Variant.ref`` (no REF ==
    ALT row), ClinPGx derives alt against ``Variant.ref`` as before."""
    other = "rs900000123"
    with contextlib.closing(sqlite3.connect(f5_data / clinvar_db_filename(GRCH37))) as conn:
        conn.execute(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, '2', 2000000, 'A', 'G', 'Benign', 'x', 'GENEX', '', 999)",
            (other,),
        )
        conn.commit()
    with contextlib.closing(sqlite3.connect(f5_data / "pharmgkb.sqlite")) as conn:
        conn.execute(
            "INSERT INTO pharmgkb_annotations (rsid, genotype, gene, drugs, phenotype, "
            "phenotype_category, annotation_text, level_of_evidence, score, "
            "pgkb_annotation_id, allele_function, function_class, is_nonfinding) "
            "VALUES (?, 'AG', 'GENEX', 'd', 'p', 'Toxicity', 't', '3', 1.0, 'PA9', '', '', 0)",
            (other,),
        )
        conn.commit()
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    assert clinvar.reference_for(other, GRCH37) == "A"
    assert clinvar.normal_allele_for(other, GRCH37) is None
    assert clinvar.normal_allele_for(F5, GRCH37) == "C"
    pharmgkb = cm_stack.enter_context(
        PharmGKBAnnotator(
            f5_data,
            clinvar_ref_provider=clinvar.reference_at,
            normal_allele_provider=clinvar.normal_allele_for,
        )
    )
    # Variant.ref G disagrees with ClinVar's A; alt stays the non-G allele.
    v = Variant(other, "2", 2_000_000, "A", "G", build=GRCH37, ref="G")
    (row,) = pharmgkb.annotate(v)
    assert (row.alt, row.locus_ref, row.reference_carries_alt) == ("A", "G", False)


# ── One ClinVar build: the identity table gives the normal allele ──


def _one_build_data(tmp_path: Path, build: str, *, identity: bool = True) -> Path:
    data = tmp_path / f"one-{build}"
    data.mkdir()
    _build_clinvar(data, (build,))
    if identity:
        _build_identity(data)
    _build_pharmgkb(data)
    _build_snpedia(data)
    return data


def _one_build_annotators(data: Path, build: str, stack: contextlib.ExitStack):
    from allelix.annotators import get_annotators

    annotators = {a.name: a for a in get_annotators(data, clinvar_builds=(build,))}
    for name in ("clinvar", "pharmgkb", "snpedia", "rsid_identity"):
        stack.enter_context(annotators[name])
    return annotators


@pytest.mark.parametrize("genotype", ["CC", "CT", "TT"])
def test_one_build_install_uses_the_identity_table(
    tmp_path: Path, cm_stack: contextlib.ExitStack, genotype: str
) -> None:
    """GRCh37 ClinVar only: CC at rs6025 is the normal genotype, as with both builds."""
    data = _one_build_data(tmp_path, GRCH37)
    annotators = _one_build_annotators(data, GRCH37, cm_stack)
    clinvar = annotators["clinvar"]
    assert clinvar.reference_for(F5, GRCH37) == "C"
    assert clinvar.normal_allele_for(F5, GRCH37) == "C"
    v = Variant(F5, "1", POS[GRCH37], genotype[0], genotype[1], build=GRCH37, ref="T")
    for source in ("pharmgkb", "snpedia"):
        rows = annotators[source].annotate(v)
        assert _findings(rows, source) == EXPECTED_GENOTYPE_ROWS[genotype][source]
        assert _findings(list(annotators[source].batch_annotate([v])), source) == _findings(
            rows, source
        )


def test_one_build_install_on_the_array_path(tmp_path: Path, cm_stack: contextlib.ExitStack):
    """The pipeline on a GRCh37 chip with the GRCh37 ClinVar cache only."""
    data = _one_build_data(tmp_path, GRCH37)
    annotators = _one_build_annotators(data, GRCH37, cm_stack)
    for genotype in ("CC", "CT", "TT"):
        result = run_analysis(
            _chip(tmp_path / f"chip-{genotype}.txt", GRCH37, genotype),
            MyHappyGenesParser(),
            [annotators["clinvar"], annotators["pharmgkb"], annotators["snpedia"]],
            build_override=GRCH37,
            identity=annotators["rsid_identity"],
        )
        for source in ("pharmgkb", "snpedia"):
            assert (
                _findings(result.annotations, source) == EXPECTED_GENOTYPE_ROWS[genotype][source]
            ), (genotype, source)


def test_one_build_install_without_the_identity_table_is_unchanged(
    tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    """No identity cache: nothing is known, nothing is suppressed (as before)."""
    data = _one_build_data(tmp_path, GRCH37, identity=False)
    annotators = _one_build_annotators(data, GRCH37, cm_stack)
    assert annotators["clinvar"].reference_for(F5, GRCH37) is None


def test_both_builds_do_not_ask_the_identity_table(
    f5_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    def refuse(carried: dict[str, str], build: str) -> dict[str, str]:
        raise AssertionError("the other build's ClinVar cache answers")

    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data, other_build_refs=refuse))
    assert clinvar.normal_allele_for(F5, GRCH37) == "C"


def _identity(tmp_path: Path, rows: list[tuple]) -> RsidIdentityAnnotator:
    data = tmp_path / "identity"
    data.mkdir()
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://f5",
        downloaded_at="2026-10-13T00:00:00Z",
    )
    return RsidIdentityAnnotator(data)


def test_identity_normal_alleles_abstain_unless_unambiguous(
    tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    rows = [
        # rs6025 as dbSNP lists it: every base, flagged palindromic. The
        # flag concerns translating alleles, not reading a REF.
        (F5, "1", POS[GRCH37], "T", "C", POS[GRCH38], "C", "C", 1),
        (F5, "1", POS[GRCH37], "T", "A", POS[GRCH38], "C", "A", 1),
        # Same REF on both builds: no normal allele.
        ("rs900953001", "1", 100, "G", "A", 200, "G", "A", 0),
        # GRCh37 REF is not the carried allele.
        ("rs900953002", "1", 300, "A", "G", 400, "C", "G", 0),
        # Two GRCh38 REFs for one rsID.
        ("rs900953003", "1", 500, "T", "C", 600, "C", "C", 0),
        ("rs900953003", "1", 500, "T", "C", 700, "G", "C", 0),
        # A multi-base GRCh38 REF.
        ("rs900953004", "1", 800, "T", "C", 900, "CA", "C", 0),
    ]
    identity = cm_stack.enter_context(_identity(tmp_path, rows))
    carried = {
        F5: "T",
        "rs900953001": "G",
        "rs900953002": "T",
        "rs900953003": "T",
        "rs900953004": "T",
        "rs900953099": "T",
    }
    assert identity.normal_alleles_on_other_build(carried, GRCH37) == {F5: "C"}
    # From GRCh38 the other way round: GRCh38's reference carries C at rs6025
    # only in this made-up case, and the normal allele is GRCh37's T.
    assert identity.normal_alleles_on_other_build({F5: "C"}, GRCH38) == {F5: "T"}
    assert identity.normal_alleles_on_other_build({}, GRCH37) == {}


# ── A normal-allele record reads the same on both builds ────────


@pytest.fixture
def normal_data(tmp_path: Path) -> Path:
    data = tmp_path / "normal-data"
    data.mkdir()
    _build_clinvar(data, BUILDS, normal_allele_record=True)
    return data


def test_normal_allele_record_is_kept_on_both_builds(normal_data: Path) -> None:
    """G>G on GRCh37 is kept (it was dropped), flagged; P/LP is dropped on both."""
    rows = {}
    for build in BUILDS:
        with contextlib.closing(sqlite3.connect(normal_data / clinvar_db_filename(build))) as conn:
            rows[build] = conn.execute(
                "SELECT ref, alt, clinical_significance, normal_allele FROM clinvar_variants "
                "WHERE rsid = ?",
                (NORMAL_RSID,),
            ).fetchall()
    assert rows[GRCH37] == [("G", "G", "risk factor", 1)]
    assert rows[GRCH38] == [("A", "G", "risk factor", 1)]


def _normal_chip(path: Path, build: str, genotype: str) -> Path:
    path.write_text(
        "# MyHappyGenes [TEMPUS]\n"
        "# Sample ID\tPROBE000004\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        f"{NORMAL_RSID}\t12\t{NORMAL_POS[build]}\t{genotype[0]}\t{genotype[1]}\n"
    )
    return path


_NORMAL_EXPECTED = {
    "GG": [("clinvar_risk_factor", "G", "Homozygous")],
    "AG": [("clinvar_risk_factor", "G", "Heterozygous")],
    "AA": [],
}


@pytest.mark.parametrize("genotype", ["GG", "AG", "AA"])
@pytest.mark.parametrize("build", BUILDS)
def test_normal_allele_record_array_call(
    normal_data: Path, tmp_path: Path, cm_stack: contextlib.ExitStack, build: str, genotype: str
) -> None:
    result = run_analysis(
        _normal_chip(tmp_path / "chip.txt", build, genotype),
        MyHappyGenesParser(),
        [cm_stack.enter_context(ClinVarAnnotator(normal_data))],
        build_override=build,
    )
    assert _findings(result.annotations) == _NORMAL_EXPECTED[genotype]


@pytest.mark.parametrize("genotype", ["GG", "AG", "AA"])
@pytest.mark.parametrize("build", BUILDS)
def test_normal_allele_record_vcf_call(
    normal_data: Path, tmp_path: Path, cm_stack: contextlib.ExitStack, build: str, genotype: str
) -> None:
    """GG on GRCh37 is a 0/0 call: the watched site gives the same row as GRCh38's 1/1."""
    ref, alt = {GRCH37: ("G", "A"), GRCH38: ("A", "G")}[build]
    gt = {"GG": {GRCH37: "0/0", GRCH38: "1/1"}, "AG": {GRCH37: "0/1", GRCH38: "0/1"},
          "AA": {GRCH37: "1/1", GRCH38: "0/0"}}[genotype][build]  # fmt: skip
    vcf = tmp_path / "n.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=12>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"12\t{NORMAL_POS[build]}\t{NORMAL_RSID}\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n"
    )
    result = run_analysis(
        vcf,
        VcfParser(),
        [cm_stack.enter_context(ClinVarAnnotator(normal_data))],
        build_override=build,
    )
    assert _findings(result.annotations) == _NORMAL_EXPECTED[genotype]


def test_normal_allele_record_gives_no_ref_where_ref_equals_alt(
    normal_data: Path, cm_stack: contextlib.ExitStack
) -> None:
    """A normal-allele REF == ALT row is not a reference-carried variant, and
    gives no REF for the ClinPGx and SNPedia non-finding check (as before
    iv:9, when it was not loaded). Written A>G on GRCh38, it keeps its REF."""
    clinvar = cm_stack.enter_context(ClinVarAnnotator(normal_data))
    assert clinvar.reference_for(NORMAL_RSID, GRCH37) is None
    assert clinvar.reference_for(NORMAL_RSID, GRCH38) == "A"
    assert clinvar.normal_allele_for(NORMAL_RSID, GRCH37) is None
    sites = clinvar.reference_carried_sites(GRCH37)
    assert sites[("12", NORMAL_POS[GRCH37])] == ("G", NORMAL_RSID, True)
    assert sites[("1", POS[GRCH37])] == ("T", F5, False)


# ── ClinPGx and SNPedia at a normal-allele record's site ─────────
#
# GRCh37's reference carries the normal allele (G), which the record
# classifies. A homozygous-reference GG call carries it, so it is not a
# non-finding there: ClinPGx and SNPedia show their GG rows, as before iv:9.
# On GRCh38 the record is written A>G, and AA (homozygous reference there)
# stays a non-finding, as before.


def _build_normal_sources(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.executescript(PHARMGKB_SCHEMA)
        conn.executemany(
            "INSERT INTO pharmgkb_annotations (rsid, genotype, gene, drugs, phenotype, "
            "phenotype_category, annotation_text, level_of_evidence, score, "
            "pgkb_annotation_id, allele_function, function_class, is_nonfinding) "
            "VALUES (?, ?, 'CASP12', 'placeholder drug', 'Sepsis', 'Efficacy', ?, "
            "'3', 1.0, ?, '', '', 0)",
            [
                (NORMAL_RSID, geno, f"Patients with the {geno} genotype", f"PA9{n}")
                for n, geno in enumerate(("AA", "AG", "GG"))
            ],
        )
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.execute(
            "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
            "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
            "orientation TEXT, scraped_at TEXT)"
        )
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, ?, ?, 'CASP12', 'plus', "
            "'2026-05-20')",
            [
                (NORMAL_RSID, "A", "A", 2.0, "Good", "placeholder summary one"),
                (NORMAL_RSID, "A", "G", 2.5, "Bad", "placeholder summary two"),
                (NORMAL_RSID, "G", "G", 3.0, "Bad", "placeholder summary three"),
            ],
        )
        conn.commit()


@pytest.fixture
def normal_sources(normal_data: Path) -> Path:
    _build_normal_sources(normal_data)
    return normal_data


# Which genotypes ClinPGx and SNPedia report, per build (default filters).
_NORMAL_SHOWN = {
    GRCH37: {"GG": True, "AG": True, "AA": True},
    GRCH38: {"GG": True, "AG": True, "AA": False},
}


@pytest.mark.parametrize("genotype", ["GG", "AG", "AA"])
@pytest.mark.parametrize("build", BUILDS)
def test_normal_allele_site_clinpgx_and_snpedia(
    normal_sources: Path, cm_stack: contextlib.ExitStack, build: str, genotype: str
) -> None:
    clinvar = cm_stack.enter_context(ClinVarAnnotator(normal_sources))
    kwargs = {
        "clinvar_ref_provider": clinvar.reference_at,
        "normal_allele_provider": clinvar.normal_allele_for,
    }
    pharmgkb = cm_stack.enter_context(PharmGKBAnnotator(normal_sources, **kwargs))
    snpedia = cm_stack.enter_context(SNPediaAnnotator(normal_sources, **kwargs))
    ref = {GRCH37: "G", GRCH38: "A"}[build]
    v = Variant(
        NORMAL_RSID, "12", NORMAL_POS[build], genotype[0], genotype[1], build=build, ref=ref
    )
    for annotator in (pharmgkb, snpedia):
        single = annotator.annotate(v)
        batch = list(annotator.batch_annotate([v]))
        assert bool(single) == _NORMAL_SHOWN[build][genotype], (annotator.name, single)
        assert [a.genotype_match for a in batch] == [a.genotype_match for a in single]


def test_normal_allele_site_risk_factor_and_snpedia_agree(
    normal_sources: Path, tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    """The case found on real files: ClinVar's risk-factor row for a GRCh37
    homozygous-reference call is shown, and so is SNPedia's row for it."""
    clinvar = cm_stack.enter_context(ClinVarAnnotator(normal_sources))
    snpedia = cm_stack.enter_context(
        SNPediaAnnotator(
            normal_sources,
            clinvar_ref_provider=clinvar.reference_at,
            normal_allele_provider=clinvar.normal_allele_for,
        )
    )
    result = run_analysis(
        _normal_chip(tmp_path / "chip.txt", GRCH37, "GG"),
        MyHappyGenesParser(),
        [clinvar, snpedia],
        build_override=GRCH37,
    )
    assert _findings(result.annotations) == [("clinvar_risk_factor", "G", "Homozygous")]
    assert [a.genotype_match for a in result.annotations if a.source == "snpedia"] == ["GG"]


@pytest.mark.parametrize("build", BUILDS)
def test_outdated_cache_keeps_the_locus_ref_of_ref_equal_alt_rows(
    f5_data: Path, cm_stack: contextlib.ExitStack, build: str
) -> None:
    """Before `db update`, a cache older than iv:9 still serves ClinPGx and
    SNPedia. Its REF == ALT rows are reference-carried variants (it held no
    normal-allele ones), and each keeps its locus REF as before iv:9."""
    with contextlib.closing(sqlite3.connect(f5_data / clinvar_db_filename(build))) as conn:
        conn.execute(
            "UPDATE database_versions SET local_version_tag = 'iv:5' WHERE name = ?",
            (clinvar_record_name(build),),
        )
        conn.commit()
    clinvar = cm_stack.enter_context(ClinVarAnnotator(f5_data))
    assert not clinvar._cache_is_current(build)
    # Written G>G on both builds, with no other single-base row.
    assert clinvar.reference_for("rs2256408", build) == "G"

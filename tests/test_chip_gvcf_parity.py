# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Chip / gVCF annotation parity harness (release 2.4 gate).

The same genotype at the same position must produce the same annotation
set whether it arrives as a consumer-array row (rsID-keyed path) or as an
rsID-less GRCh38 HaplotypeCaller gVCF record (position-resolved path):
alt for alt, magnitude for magnitude, frequency for frequency.

Every case below is one synthetic site. The harness writes a MyHappyGenes
chip file (GRCh38 coordinates, rsID column populated) and a gVCF (``ID=.``,
``<NON_REF>`` reference blocks, ``chr``-prefixed contigs) carrying the
identical genotype, runs both through ``run_analysis`` against the same
synthetic ClinVar (and optionally gnomAD) caches, and compares the
complete annotation multisets. Each case is also checked against a
first-principles oracle — the ClinVar rows whose alt the user carries —
so that two paths agreeing on an empty set cannot pass vacuously.

Three annotator configurations are exercised because they take different
routes through rsID-less resolution: ClinVar alone (the position
resolver is the only way an rsID-less record gets an rsID), ClinVar
plus gnomAD (gnomAD's second-pass resolver and ``Variant.ref``
population for array data both fire), and ClinVar plus the cross-build
rsID identity table (its locus resolver and ``Variant.ref`` source fire).

Sites, rsIDs and genotypes are synthetic; they reproduce the shapes of
the field-reported discrepancies (multi-allelic sites, several ClinVar
submissions per rsID) without carrying any genotype from a real file.

Known defects on the current code are listed in ``KNOWN_DEFECTS`` and
run as strict xfails: the harness stays green today, and the change that
fixes a defect must delete its entry (a strict xfail that starts passing
fails the suite).
"""

from __future__ import annotations

import contextlib
import sqlite3
from collections import Counter
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import (
    _BENIGN_CLNSIGS,
    ClinVarAnnotator,
    _magnitude,
    _normalize_clnsig,
    clinvar_db_filename,
    clinvar_record_name,
)
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.rsid_identity_loader import (
    RSID_IDENTITY_DB_FILENAME,
    pair_alleles,
    write_identity_db,
)
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.parsers.base import GenotypeParser

BUILD = "GRCh38"
_SINGLE = "criteria_provided,_single_submitter"
_MULTI = "criteria_provided,_multiple_submitters,_no_conflicts"


@dataclass(frozen=True)
class ClinVarRow:
    """One per-SCV ClinVar cache row at a parity site."""

    alt: str
    significance: str
    condition: str
    review_status: str = _SINGLE


@dataclass(frozen=True)
class ParitySite:
    """A synthetic site and the genotype both fixtures carry there.

    ``clinvar`` rows share the site's rsID, chromosome, position and REF,
    as per-SCV rows from the variant_summary / submission_summary loader
    do. ``gnomad`` maps alt to allele frequency; empty means the site is
    absent from gnomAD.
    """

    case: str
    rsid: str
    chrom: str
    pos: int
    ref: str
    genotype: tuple[str, str]
    gene: str
    clinvar: tuple[ClinVarRow, ...]
    gnomad: tuple[tuple[str, float], ...] = ()


SITES: tuple[ParitySite, ...] = (
    ParitySite(
        case="snv_het_single_alt",
        rsid="rs900178001",
        chrom="1",
        pos=1_780_010,
        ref="A",
        genotype=("A", "G"),
        gene="PARA",
        clinvar=(ClinVarRow("G", "Pathogenic", "Parity_condition_A"),),
        gnomad=(("G", 0.012),),
    ),
    ParitySite(
        case="snv_hom_alt_single_alt",
        rsid="rs900178002",
        chrom="1",
        pos=1_780_020,
        ref="C",
        genotype=("T", "T"),
        gene="PARB",
        clinvar=(ClinVarRow("T", "Likely_pathogenic", "Parity_condition_B", _MULTI),),
        gnomad=(("T", 0.21),),
    ),
    ParitySite(
        # Multi-allelic in ClinVar; the user carries exactly one of the alts.
        case="multiallelic_het_one_alt_carried",
        rsid="rs900178003",
        chrom="2",
        pos=2_780_030,
        ref="C",
        genotype=("C", "T"),
        gene="PARC",
        clinvar=(
            ClinVarRow("A", "Pathogenic", "Parity_condition_C_A"),
            ClinVarRow("G", "Pathogenic", "Parity_condition_C_G"),
            ClinVarRow("T", "Pathogenic", "Parity_condition_C_T"),
        ),
        gnomad=(("A", 0.0001), ("G", 0.0002), ("T", 0.03)),
    ),
    ParitySite(
        # Homozygous reference at a multi-allelic pathogenic site: the chip
        # reports REF/REF, the gVCF covers the position with a confident
        # <NON_REF> block. Neither may produce an annotation.
        case="multiallelic_hom_ref",
        rsid="rs900178004",
        chrom="2",
        pos=2_780_040,
        ref="C",
        genotype=("C", "C"),
        gene="PARD",
        clinvar=(
            ClinVarRow("A", "Pathogenic", "Parity_condition_D"),
            ClinVarRow("G", "Pathogenic", "Parity_condition_D"),
            ClinVarRow("T", "Pathogenic", "Parity_condition_D"),
        ),
        gnomad=(("A", 0.0003), ("T", 0.02)),
    ),
    ParitySite(
        # Two different alts carried (gVCF GT 1/2).
        case="multiallelic_het_two_alts",
        rsid="rs900178005",
        chrom="10",
        pos=10_780_050,
        ref="A",
        genotype=("G", "T"),
        gene="PARE",
        clinvar=(
            ClinVarRow("G", "Pathogenic", "Parity_condition_E_G"),
            ClinVarRow("T", "Uncertain_significance", "Parity_condition_E_T"),
        ),
        gnomad=(("G", 0.004), ("T", 0.001)),
    ),
    ParitySite(
        # Several submissions for the same alt plus a second alt at the
        # site: the shape shared by the sites the gVCF path dropped in the
        # field comparison.
        case="multi_scv_shared_alt",
        rsid="rs900178006",
        chrom="10",
        pos=10_780_060,
        ref="T",
        genotype=("T", "G"),
        gene="PARF",
        clinvar=(
            ClinVarRow("G", "Pathogenic", "Parity_condition_F1"),
            ClinVarRow("G", "Likely_pathogenic", "Parity_condition_F2"),
            ClinVarRow("C", "Uncertain_significance", "Parity_condition_F3"),
        ),
        gnomad=(("G", 0.35), ("C", 0.0005)),
    ),
    ParitySite(
        # Benign rows are suppressed by default on both paths; the
        # pathogenic row for the same alt must survive on both.
        case="benign_suppression",
        rsid="rs900178007",
        chrom="10",
        pos=10_780_070,
        ref="G",
        genotype=("G", "A"),
        gene="PARG",
        clinvar=(
            ClinVarRow("A", "Benign", "Parity_condition_G1"),
            ClinVarRow("A", "Pathogenic", "Parity_condition_G2"),
        ),
        gnomad=(("A", 0.06),),
    ),
    ParitySite(
        # ClinVar reference-allele record (HGVS "=", REF == ALT in the VCF
        # columns) at a multi-allelic pathogenic site; the user is
        # homozygous reference. Neither path may report a finding.
        case="reference_allele_record_hom_ref",
        rsid="rs900178008",
        chrom="10",
        pos=10_780_080,
        ref="C",
        genotype=("C", "C"),
        gene="PARH",
        clinvar=(
            ClinVarRow("T", "Pathogenic", "Parity_condition_H"),
            ClinVarRow("C", "Pathogenic", "Parity_condition_H"),
        ),
        gnomad=(("T", 0.0004),),
    ),
    ParitySite(
        # Reference-allele record at a site where the user is het for a
        # real alt: only the real alt is reported.
        case="reference_allele_record_het",
        rsid="rs900178009",
        chrom="10",
        pos=10_780_090,
        ref="A",
        genotype=("A", "C"),
        gene="PARI",
        clinvar=(
            ClinVarRow("C", "Pathogenic", "Parity_condition_I1"),
            ClinVarRow("A", "Pathogenic", "Parity_condition_I2"),
        ),
        gnomad=(("C", 0.4),),
    ),
)

# Configurations: ClinVar alone, plus the gnomAD cache, or plus the cross-build
# rsID identity table (ADR-0040), which adds its own locus resolver and
# ``Variant.ref`` source.
CONFIGS: tuple[str, ...] = ("clinvar", "clinvar+gnomad", "clinvar+identity")

# (case, config) -> reason. Strict xfail: the fix must remove the entry.
KNOWN_DEFECTS: dict[tuple[str, str], str] = {}


# ── Fixture writers ─────────────────────────────────────────────


def _write_clinvar_cache(data_dir: Path, sites: tuple[ParitySite, ...]) -> None:
    rows = []
    allele_id = 178_000
    for site in sites:
        for row in site.clinvar:
            allele_id += 1
            rows.append(
                (
                    site.rsid,
                    site.chrom,
                    site.pos,
                    site.ref,
                    row.alt,
                    row.significance,
                    row.condition,
                    site.gene,
                    row.review_status,
                    allele_id,
                )
            )
    with contextlib.closing(sqlite3.connect(data_dir / clinvar_db_filename(BUILD))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://parity', '20261010', '2026-10-10', ?, NULL, ?)",
            (clinvar_record_name(BUILD), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _write_gnomad_cache(data_dir: Path, sites: tuple[ParitySite, ...]) -> None:
    rows = [
        (site.chrom, site.pos, site.ref, alt, site.rsid, af)
        for site in sites
        for alt, af in site.gnomad
    ]
    with contextlib.closing(sqlite3.connect(data_dir / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) "
            "VALUES ('gnomad', 'test://parity', '4.1', '2026-10-10', ?, 'sv:1')",
            (len(rows),),
        )
        conn.commit()


def _write_identity_cache(data_dir: Path, sites: tuple[ParitySite, ...]) -> None:
    """Identity rows for every non-reference allele at each site (same strand)."""
    rows = []
    for site in sites:
        alts = {row.alt for row in site.clinvar} | set(site.genotype)
        alts.discard(site.ref)
        palindromic = int(bool(pair_alleles(site.ref, alts, site.ref, alts)[1]))  # type: ignore[index]
        rows += [
            (
                site.rsid,
                site.chrom,
                site.pos - 1000,
                site.ref,
                alt,
                site.pos,
                site.ref,
                alt,
                palindromic,
            )
            for alt in sorted(alts)
        ]
    write_identity_db(
        rows,
        data_dir / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP parity",
        source_url="test://parity",
        downloaded_at="2026-10-10T00:00:00Z",
    )


def _write_chip(path: Path, sites: tuple[ParitySite, ...]) -> None:
    lines = [
        "# MyHappyGenes [TEMPUS]",
        "# Synthetic parity fixture. Positions are GRCh38.",
        "# Sample ID\tPARITY000001",
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward",
    ]
    lines += [f"{s.rsid}\t{s.chrom}\t{s.pos}\t{s.genotype[0]}\t{s.genotype[1]}" for s in sites]
    path.write_text("\n".join(lines) + "\n")


_GVCF_HEADER = """\
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP in block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr10,length=133797422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPARITY000001
"""


def _gvcf_record(site: ParitySite) -> str:
    """Render the site as a HaplotypeCaller gVCF line with ``ID=.``."""
    chrom = f"chr{site.chrom}"
    a1, a2 = site.genotype
    if a1 == site.ref and a2 == site.ref:
        # Confident single-base reference block covering the site.
        return (
            f"{chrom}\t{site.pos}\t.\t{site.ref}\t<NON_REF>\t.\t.\tEND={site.pos}"
            "\tGT:DP:GQ:MIN_DP:PL\t0/0:40:99:40:0,120,1800"
        )
    alts: list[str] = []
    for allele in (a1, a2):
        if allele != site.ref and allele not in alts:
            alts.append(allele)
    alleles = [site.ref, *alts]
    gt = "/".join(sorted(str(alleles.index(a)) for a in (a1, a2)))
    ad = ",".join(["20" if a in (a1, a2) else "0" for a in alleles] + ["0"])
    return (
        f"{chrom}\t{site.pos}\t.\t{site.ref}\t{','.join([*alts, '<NON_REF>'])}\t512.64\t.\t."
        f"\tGT:AD:DP:GQ:PL\t{gt}:{ad}:40:99:0,90,1200"
    )


def _write_gvcf(path: Path, sites: tuple[ParitySite, ...]) -> None:
    """Write records in coordinate order with reference blocks in between."""
    body: list[str] = []
    covered_to: dict[str, int] = {}
    for site in sorted(sites, key=lambda s: (int(s.chrom), s.pos)):
        chrom = f"chr{site.chrom}"
        block_start = covered_to.get(chrom, site.pos - 500) + 1
        if block_start < site.pos:
            body.append(
                f"{chrom}\t{block_start}\t.\tN\t<NON_REF>\t.\t.\tEND={site.pos - 1}"
                "\tGT:DP:GQ:MIN_DP:PL\t0/0:35:90:30:0,90,1350"
            )
        body.append(_gvcf_record(site))
        covered_to[chrom] = site.pos
    path.write_text(_GVCF_HEADER + "\n".join(body) + "\n")


# ── Harness ─────────────────────────────────────────────────────


AnnotationKey = tuple[object, ...]


def _canonical(a: Annotation) -> AnnotationKey:
    """Every user-visible Annotation field, in a hashable comparison key."""
    return (
        a.source,
        a.rsid,
        a.alt,
        a.significance,
        a.magnitude,
        a.condition,
        a.gene,
        a.review_status,
        a.genotype_match,
        a.zygosity,
        tuple(a.references),
        a.allele_frequency,
        a.am_pathogenicity,
        a.am_class,
        a.cadd_phred,
        a.description,
        a.trait,
        a.p_value,
        a.phecode,
    )


def _expected(site: ParitySite) -> Counter[tuple[str, str, str, float]]:
    """Oracle: ClinVar rows whose non-reference alt the user carries, benign suppressed."""
    carried = set(site.genotype)
    return Counter(
        (
            site.rsid,
            row.alt,
            f"clinvar_{_normalize_clnsig(row.significance)}",
            _magnitude(row.significance),
        )
        for row in site.clinvar
        if row.alt != site.ref
        and row.alt in carried
        and _normalize_clnsig(row.significance) not in _BENIGN_CLNSIGS
    )


@dataclass
class ParityRun:
    chip: list[Annotation]
    gvcf: list[Annotation]


def _run(
    data_dir: Path,
    file_path: Path,
    parser: GenotypeParser,
    config: str,
    cm_stack: contextlib.ExitStack,
) -> list[Annotation]:
    assert parser.can_parse(file_path), f"{parser.name} rejected the parity fixture"
    clinvar = cm_stack.enter_context(ClinVarAnnotator(data_dir, builds=(BUILD,)))
    gnomad = (
        cm_stack.enter_context(GnomadAnnotator(data_dir)) if config == "clinvar+gnomad" else None
    )
    identity = (
        cm_stack.enter_context(RsidIdentityAnnotator(data_dir))
        if config == "clinvar+identity"
        else None
    )
    result = run_analysis(
        file_path, parser, [clinvar], build_override=BUILD, gnomad=gnomad, identity=identity
    )
    return result.annotations


def run_parity(
    tmp_path: Path,
    sites: tuple[ParitySite, ...],
    config: str,
    cm_stack: contextlib.ExitStack,
) -> ParityRun:
    """Build caches and the chip / gVCF pair for ``sites``; annotate both."""
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    _write_clinvar_cache(data_dir, sites)
    if config == "clinvar+gnomad":
        _write_gnomad_cache(data_dir, sites)
    if config == "clinvar+identity":
        _write_identity_cache(data_dir, sites)
    chip_path = tmp_path / "parity_chip.txt"
    gvcf_path = tmp_path / "parity.g.vcf"
    _write_chip(chip_path, sites)
    _write_gvcf(gvcf_path, sites)
    return ParityRun(
        chip=_run(data_dir, chip_path, MyHappyGenesParser(), config, cm_stack),
        gvcf=_run(data_dir, gvcf_path, VcfParser(), config, cm_stack),
    )


def assert_parity(run: ParityRun) -> None:
    """Fail with a readable symmetric diff when the two paths disagree."""
    chip = Counter(_canonical(a) for a in run.chip)
    gvcf = Counter(_canonical(a) for a in run.gvcf)
    if chip == gvcf:
        return
    only_chip = "\n  ".join(repr(k) for k in sorted((chip - gvcf).elements(), key=repr))
    only_gvcf = "\n  ".join(repr(k) for k in sorted((gvcf - chip).elements(), key=repr))
    pytest.fail(
        "chip / gVCF annotation sets differ\n"
        f"only on chip path:\n  {only_chip or '(none)'}\n"
        f"only on gVCF path:\n  {only_gvcf or '(none)'}"
    )


def assert_matches_oracle(annotations: list[Annotation], sites: tuple[ParitySite, ...]) -> None:
    expected: Counter[tuple[str, str, str, float]] = Counter()
    for site in sites:
        expected += _expected(site)
    actual = Counter((a.rsid, a.alt, a.significance, a.magnitude) for a in annotations)
    assert actual == expected, (
        f"expected {sorted(expected.elements())}, got {sorted(actual.elements())}"
    )


# ── Tests ───────────────────────────────────────────────────────


def _cases() -> list[object]:
    params: list[object] = []
    for site in SITES:
        for config in CONFIGS:
            reason = KNOWN_DEFECTS.get((site.case, config))
            marks = [pytest.mark.xfail(reason=reason, strict=True)] if reason else []
            params.append(pytest.param(site, config, id=f"{site.case}-{config}", marks=marks))
    return params


@pytest.mark.parametrize(("site", "config"), _cases())
def test_site_parity(
    site: ParitySite, config: str, tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    """One site per fixture pair: both paths match each other and the oracle."""
    run = run_parity(tmp_path, (site,), config, cm_stack)
    assert_matches_oracle(run.chip, (site,))
    assert_parity(run)
    assert_matches_oracle(run.gvcf, (site,))


@pytest.mark.parametrize("config", CONFIGS)
def test_panel_parity(config: str, tmp_path: Path, cm_stack: contextlib.ExitStack) -> None:
    """All sites without a known defect in one file pair (no cross-site interference)."""
    sites = tuple(s for s in SITES if (s.case, config) not in KNOWN_DEFECTS)
    run = run_parity(tmp_path, sites, config, cm_stack)
    assert run.chip, "parity panel produced no annotations; the harness is vacuous"
    assert_matches_oracle(run.chip, sites)
    assert_parity(run)


def test_gvcf_fixture_shape(tmp_path: Path) -> None:
    """The gVCF writer emits rsID-less records with <NON_REF> on every line."""
    path = tmp_path / "shape.g.vcf"
    _write_gvcf(path, SITES)
    records = [ln.split("\t") for ln in path.read_text().splitlines() if not ln.startswith("#")]
    assert records
    assert all(cols[2] == "." for cols in records)
    assert all(cols[4].endswith("<NON_REF>") for cols in records)
    assert all(cols[0].startswith("chr") for cols in records)
    by_pos = {int(cols[1]): cols for cols in records}
    two_alts = next(s for s in SITES if s.case == "multiallelic_het_two_alts")
    assert by_pos[two_alts.pos][4] == "G,T,<NON_REF>"
    assert by_pos[two_alts.pos][9].startswith("1/2:")

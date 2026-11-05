"""Small synthetic caches and call shapes for validated-call contracts."""

import sqlite3
from collections import Counter
from contextlib import closing

from test_annotation_parity_and_named_disclosures import _identity, _mthfr, _vcf

from allelix.annotators import get_annotators
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

BUILDS = ("GRCh37", "GRCh38")
SOURCES = ("clinvar", "pharmgkb", "snpedia", "gwas")
RSID = "rs1801133"
POSITIONS = {"GRCh37": 11856378, "GRCh38": 11796321}


def sql(data, name, statement, args=()):
    with closing(sqlite3.connect(data / name)) as conn:
        conn.execute(statement, args)
        conn.commit()


def setup(data, build):
    return _mthfr(data, build)


def call(path, pos, *, rsid=RSID, chrom="1", shape="plain", gt="0/1", alt=None):
    ref, default_alt = {
        "plain": ("G", "A"),
        "prefix": ("CG", "CA"),
        "suffix": ("GC", "AC"),
        "both": ("CGC", "CAC"),
    }[shape]
    return _vcf(
        path,
        [(rsid, chrom, pos - (shape in ("prefix", "both")), ref, alt or default_alt, gt, ".")],
    )


def analyze(data, path, build, sources=SOURCES, *, parser=None, **kwargs):
    registry = get_annotators(data, clinvar_builds=(build,))
    by = {a.name: a for a in registry}
    return run_analysis(
        path,
        parser or VcfParser(),
        [by[n] for n in sources],
        identity=by["rsid_identity"],
        build_override=build,
        resource_owners=registry,
        **kwargs,
    )


def signatures(result):
    # Source identity and interpretation, excluding spelling/ploidy on purpose.
    return Counter(
        (a.source, a.rsid, a.significance, a.condition, a.magnitude, a.description)
        for a in result.annotations
    )


def warning_count(result, code):
    return sum(w.count or 0 for w in result.report_warnings if w.code == code)


def identity(data, alleles="AG"):
    (data / "rsid_identity.sqlite").unlink()
    _identity(
        data,
        [
            (
                RSID,
                "1",
                POSITIONS["GRCh37"],
                "G",
                a,
                POSITIONS["GRCh38"],
                "G",
                a,
                int(alleles == "ACGT"),
            )
            for a in alleles
            if a != "G"
        ],
    )


def homozygous_sources(data):
    sql(data, "pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='AA'")
    sql(data, "snpedia.sqlite", "UPDATE snpedia_genotypes SET allele1='A',allele2='A'")


def reference_effect(data, build):
    pos = setup(data, build)
    sql(data, f"clinvar.{build}.sqlite", "UPDATE clinvar_variants SET alt='G'")
    sql(data, "pharmgkb.sqlite", "UPDATE pharmgkb_annotations SET genotype='GG'")
    sql(data, "snpedia.sqlite", "UPDATE snpedia_genotypes SET allele1='G',allele2='G'")
    sql(data, "gwas.sqlite", "UPDATE gwas_associations SET risk_allele='G'")
    (data / "rsid_identity.sqlite").unlink()
    refs = ("G", "A", "A", "G") if build == "GRCh37" else ("A", "G", "G", "A")
    _identity(
        data, [(RSID, "1", POSITIONS["GRCh37"], *refs[:2], POSITIONS["GRCh38"], *refs[2:], 0)]
    )
    return pos

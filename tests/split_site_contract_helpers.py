"""Synthetic source caches and raw records for the split-site contract."""

import json
import re
import sqlite3
from contextlib import closing

from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity, _pharm, _stamp

from allelix.cli import main
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION
from allelix.databases.schema import GWAS_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA

BUILDS = ("GRCh37", "GRCh38")
RSID = "rs999000555"
OTHER = "rs999000556"
POS = {"GRCh37": 200400, "GRCh38": 200500}
LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}
SOURCES = ("clinvar", "gwas", "pharmgkb", "snpedia")


def data_dir(path, *, identity=True, rsid=RSID, pathogenic=True, alts="AC"):
    path.mkdir()
    for build in BUILDS:
        _clinvar(
            path,
            build,
            [
                (rsid, "1", POS[build], "G", alt, "SYNGENE", sig, sig, 1, 1)
                for alt in "ATC"
                for sig in ["Pathogenic" if pathogenic else "Benign"]
            ],
        )
    _identity(path, [(rsid, "1", POS["GRCh37"], "G", a, POS["GRCh38"], "G", a, 0) for a in alts])
    genotypes = sorted(
        {
            "".join(sorted(pair))
            for pair in [
                ("G", "G"),
                (alts[0], alts[0]),
                (alts[1], alts[1]),
                ("G", alts[0]),
                ("G", alts[1]),
                tuple(alts),
            ]
        }
    )
    _pharm(path, rsid, "SYNGENE", genotypes)
    with closing(sqlite3.connect(path / "snpedia.sqlite")) as conn, conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        for gt in genotypes:
            conn.execute(
                "INSERT INTO snpedia_genotypes "
                "(rsid,allele1,allele2,magnitude,repute,summary,gene,orientation) "
                "VALUES (?,?,?,6,'Bad',?,'SYNGENE','plus')",
                (rsid, *gt, "Synthetic " + gt),
            )
        _stamp(conn, "snpedia", len(genotypes), f"pv:{_PARSER_VERSION}")
    with closing(sqlite3.connect(path / "gwas.sqlite")) as conn, conn:
        conn.executescript(GWAS_SCHEMA)
        for a in alts:
            conn.execute(
                "INSERT INTO gwas_associations "
                "(rsid,risk_allele,trait,p_value,gene,study_accession,trait_category) "
                "VALUES (?,?,?,1e-30,'SYNGENE',?,'disease')",
                (rsid, a, "Synthetic risk " + a, "GCST90027100" + a),
            )
        _stamp(conn, "gwas", len(alts), f"cv:{_CATEGORIZER_VERSION}")
    # G>A,C gives GWAS one unambiguous A risk row. G>A,T is used
    # separately for AT-person guards; both risk bases there are strand-ambiguous.
    (path / "config.toml").write_text(
        "[sources]\nclinvar=true\ngwas=true\npharmgkb=true\nsnpedia=true\n"
        "gnomad=false\nalphamissense=false\ncadd=false\n"
        f"rsid_identity={str(identity).lower()}\n"
    )
    return path


def record(alt, gt, *, rsid=RSID, filt="PASS", gq=99, dp=30, ps=".", ft="PASS", pos=None):
    return dict(alt=alt, gt=gt, rsid=rsid, filt=filt, gq=gq, dp=dp, ps=ps, ft=ft, pos=pos)


def write_vcf(path, build, rows):
    text = (
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={LENGTH[build]}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    for r in rows:
        text += (
            f"1\t{r['pos'] or POS[build]}\t{r['rsid']}\tG\t{r['alt']}\t50\t{r['filt']}\t.\t"
            f"GT:GQ:DP:PS:FT\t{r['gt']}:{r['gq']}:{r['dp']}:{r['ps']}:{r['ft']}\n"
        )
    path.write_text(text)
    return path


def run(data, source, output, build, *, auto=False, extra=()):
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(output),
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--snpedia-min-magnitude",
        "0",
        "--include-benign",
        "--gwas-all",
        *extra,
    ]
    if not auto:
        args += ["--build", build.lower()]
    result = CliRunner(env={"COLUMNS": "240"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, repr(result.exception))
    return json.loads(output.read_text()), result.output


def by_source(report):
    return {s: [a for a in report["annotations"] if a["source"] == s] for s in SOURCES}


def signature(report):
    """All public finding fields, preserving multiplicity; input provenance is separate."""
    return sorted(json.dumps(a, sort_keys=True) for a in report["annotations"])


def assert_assessed_halves(text, count):
    """Allow either the existing half-call note or a truthful split-allele note."""
    flat = " ".join(text.split())
    pattern = (
        rf"\b{count}\s+(?:(?:split[- ]site\s+)?half call(?:\(s\)|s)?|"
        rf"allele(?:\(s\)|s)?\s+(?:of|from)\s+split sites)"
    )
    match = re.search(pattern, flat)
    assert match, f"The assessed count {count} must be disclosed: {flat}"
    note = flat[match.start() : match.start() + 650]
    assert "ClinVar" in note and "GWAS" in note and "assess" in note
    assert "ClinPGx" in note and "SNPedia" in note

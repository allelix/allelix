"""Shared synthetic RC3 verification fixtures; no network or private data."""

import contextlib
import json
import sqlite3

import pysam
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main, utility
from allelix.databases._versions import (
    CLINVAR_INTERPRETER_VERSION,
    GNOMAD_SCHEMA_VERSION,
    PHARMGKB_INTERPRETER_VERSION,
)
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION
from allelix.databases.rsid_identity_loader import write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA, GWAS_SCHEMA, PHARMGKB_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA

RS = "rs880271201"
OTHER = "rs880271202"
FOREIGN = "rs880271299"
POS = {"GRCh37": 777101, "GRCh38": 777201}
LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}


def stamp(c, name, n, tag):
    c.execute(
        (
            "INSERT INTO database_versions(name,source_url,version,downlo"
            "aded_at,record_count,local_version_tag) VALUES (?,'test://of"
            "fline','synthetic','2026-01-01',?,?)"
        ),
        (name, n, tag),
    )


def caches(d, cv, ident, gn, risks):
    d.mkdir()
    for build, pos in POS.items():
        with contextlib.closing(sqlite3.connect(d / clinvar_db_filename(build))) as c, c:
            c.executescript(CLINVAR_SCHEMA)
            for i, (rid, alt) in enumerate(cv):
                c.execute(
                    (
                        "INSERT INTO clinvar_variants(rsid,chromosome,position,ref,al"
                        "t,gene,clinical_significance,aggregate,condition,review_stat"
                        "us,allele_id) VALUES (?,'1',?,'G',?,'VERIFYGENE','Pathogenic"
                        "','Pathogenic','Synthetic verification','criteria provided, "
                        "single submitter',?)"
                    ),
                    (rid, pos, alt, 88101 + i),
                )
            stamp(c, clinvar_record_name(build), len(cv), f"iv:{CLINVAR_INTERPRETER_VERSION}")
    write_identity_db(
        [(rid, "1", POS["GRCh37"], "G", alt, POS["GRCh38"], "G", alt, 0) for rid, alt in ident],
        d / "rsid_identity.sqlite",
        version="synthetic",
        source_url="test://offline",
        downloaded_at="2026-01-01",
    )
    with contextlib.closing(sqlite3.connect(d / "gnomad.sqlite")) as c, c:
        c.executescript(GNOMAD_SCHEMA)
        for rid, alt in gn:
            c.execute(
                (
                    "INSERT INTO gnomad_frequencies(chrom,pos,ref,alt,rsid,af) VA"
                    "LUES ('1',?,'G',?,?,0.01)"
                ),
                (POS["GRCh38"], alt, rid),
            )
        stamp(c, "gnomad", len(gn), f"sv:{GNOMAD_SCHEMA_VERSION}")
    with contextlib.closing(sqlite3.connect(d / "gwas.sqlite")) as c, c:
        c.executescript(GWAS_SCHEMA)
        for rid, alt in risks:
            c.execute(
                (
                    "INSERT INTO gwas_associations(rsid,risk_allele,trait,p_value"
                    ",gene,study_accession,trait_category) VALUES (?,?,?,1e-30,'V"
                    "ERIFYGENE',?,'disease')"
                ),
                (rid, alt, "Verification risk " + str(alt), "GCST880271" + str(alt)),
            )
        stamp(c, "gwas", len(risks), f"cv:{_CATEGORIZER_VERSION}")
    cfg = {
        "clinvar": bool(cv),
        "rsid_identity": bool(ident),
        "gnomad": bool(gn),
        "gwas": bool(risks),
        "pharmgkb": False,
        "snpedia": False,
        "cadd": False,
        "alphamissense": False,
    }
    (d / "config.toml").write_text(
        "[sources]\n" + "".join(f"{k}={str(v).lower()}\n" for k, v in cfg.items())
    )
    return d


def vcf(p, build, alts, ids, gts, phase=False):
    text = (
        "##fileformat=VCFv4.2\n"
        + f"##contig=<ID=1,length={LENGTH[build]}>\n"
        + (
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
            '\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality"'
            '>\n##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase s'
            'et">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
        )
    )
    # Unrelated flanks keep contradiction-only inputs analyzable on the old release.
    rows = [
        (100001, ".", "T", "0/1"),
        *[(POS[build], r, a, g) for r, a, g in zip(ids, alts, gts, strict=True)],
        (999991, ".", "T", "0/1"),
    ]
    for pos, rid, alt, gt in rows:
        text += (
            f"1\t{pos}\t{rid}\tG\t{alt}\t99\tPASS\t.\tGT:GQ:PS\t{gt}:99:{7 if phase else '.'}\n"
        )
    p.write_text(text)
    return p


def diploid_sources(d, rid):
    genotypes = ["AA", "AG", "GG", "CG", "CC", "AC", "AT", "TT", "GT", "CT"]
    with contextlib.closing(sqlite3.connect(d / "pharmgkb.sqlite")) as c, c:
        c.executescript(PHARMGKB_SCHEMA)
        for i, gt in enumerate(genotypes):
            c.execute(
                (
                    "INSERT INTO pharmgkb_annotations(rsid,genotype,gene,drugs,ph"
                    "enotype,phenotype_category,annotation_text,level_of_evidence"
                    ",score,pgkb_annotation_id,function_class,is_nonfinding) VALU"
                    "ES (?,?,'VERIFYGENE','Synthetic drug','Synthetic response','"
                    "Efficacy','Verification','1A',10,?,'decreased',0)"
                ),
                (rid, gt, str(i + 89101)),
            )
        stamp(c, "pharmgkb", len(genotypes), f"iv:{PHARMGKB_INTERPRETER_VERSION}")
    with contextlib.closing(sqlite3.connect(d / "snpedia.sqlite")) as c, c:
        c.executescript(_STRUCTURED_SCHEMA)
        for gt in genotypes:
            c.execute(
                (
                    "INSERT INTO snpedia_genotypes(rsid,allele1,allele2,magnitude"
                    ",repute,summary,gene,orientation) VALUES (?,?,?,6,'Bad',?,'V"
                    "ERIFYGENE','plus')"
                ),
                (rid, gt[0], gt[1], "Verification " + gt),
            )
        stamp(c, "snpedia", len(genotypes), f"pv:{_PARSER_VERSION}")
    p = d / "config.toml"
    p.write_text(
        p.read_text()
        .replace("pharmgkb=false", "pharmgkb=true")
        .replace("snpedia=false", "snpedia=true")
    )


def write_vcf(p, build, rows):
    vcf(p, build, "", (), ())
    text = p.read_text()
    headers = text.split("#CHROM")[0]
    headers += (
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n#'
        '#FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filt'
        'er">\n##FILTER=<ID=LowQual,Description="Rejected">\n'
    )
    text = headers + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    full = [
        dict(pos=100001, alt="T", gt="0/1", id="."),
        *rows,
        dict(pos=999991, alt="T", gt="0/1", id="."),
    ]
    for r in full:
        text += (
            f"1\t{r.get('pos', POS[build])}\t{r.get('id', RS)}\tG\t{r['alt']}\t99\t"
            f"{r.get('filter', 'PASS')}\t.\tGT:GQ:DP:PS:FT\t{r['gt']}:"
            f"{r.get('gq', 99)}:{r.get('dp', 30)}:{r.get('ps', '.')}:{r.get('ft', 'PASS')}\n"
        )
    p.write_text(text)
    return p


def analyze(data, source, output, build, *, extra=()):
    result = CliRunner(env={"COLUMNS": "500", "TERM": "dumb"}).invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            build,
            "--gwas-all",
            "--gwas-min-magnitude",
            "0",
            "--min-magnitude",
            "0",
            "--include-benign",
            "--output",
            str(output),
            *map(str, extra),
        ],
    )
    assert result.exit_code == 0, (result.output, repr(result.exception))
    return json.loads(output.read_text()), " ".join(result.output.split())


def extracted(data, source, build, rsids, monkeypatch, *, indexed=False):
    if indexed:
        zipped = source.with_suffix(".vcf.gz")
        pysam.tabix_compress(str(source), str(zipped), force=True)
        pysam.tabix_index(str(zipped), preset="vcf", force=True)
        source = zipped
    trace = []
    original = utility._try_tabix_extract

    def traced(*args, **kwargs):
        result = original(*args, **kwargs)
        trace.append(result is not None)
        return result

    with monkeypatch.context() as m:
        m.setattr(utility, "_try_tabix_extract", traced)
        result = CliRunner(env={"COLUMNS": "500", "TERM": "dumb"}).invoke(
            main,
            [
                "extract",
                str(source),
                "--data-dir",
                str(data),
                "--build",
                build,
                "--snps",
                ",".join(rsids),
            ],
        )
    assert result.exit_code == 0, (result.output, repr(result.exception))
    if indexed:
        assert any(trace), "The indexed control must actually take the tabix path"
    cells = {}
    for line in result.output.splitlines():
        if "│" in line:
            row = [cell.strip() for cell in line.split("│")[1:-1]]
            if row and row[0] in rsids:
                cells[row[0]] = row
    assert set(cells) == set(rsids), result.output
    return cells


def warning(report, code):
    return [item for item in report.get("report_warnings", []) if item["code"] == code]


def findings(report, source):
    return [item for item in report["annotations"] if item["source"] == source]


def comparable(report):
    return (
        sorted(json.dumps(a, sort_keys=True) for a in report["annotations"]),
        report["input"]["total_variants"],
        report["input"]["vcf_quality"],
        report.get("report_warnings", []),
    )

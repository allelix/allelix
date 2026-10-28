# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Synthetic annotation parity and named disclosures without changing clinical scores."""

from __future__ import annotations

import contextlib
import html
import json
import re
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators import get_annotators
from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import (
    CLINVAR_INTERPRETER_VERSION,
    GNOMAD_SCHEMA_VERSION,
    PHARMGKB_INTERPRETER_VERSION,
)
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA, GWAS_SCHEMA, PHARMGKB_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION, _STRUCTURED_SCHEMA
from allelix.models import Variant

BUILDS = ("GRCh37", "GRCh38")
SURFACES = ("json", "html", "terminal")


def _stamp(conn, name, count, tag):
    conn.execute(
        "INSERT INTO database_versions (name,source_url,version,"
        "downloaded_at,record_count,local_version_tag) "
        "VALUES (?, 'test://synthetic', 'synthetic', '2001-01-01', ?, ?)",
        (name, count, tag),
    )


def _clinvar(data, build, rows):
    """Rows: rsid, chrom, pos, ref, alt, gene, classification, aggregate, same, total."""
    data.mkdir(exist_ok=True)
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        for i, (rsid, chrom, pos, ref, alt, gene, sig, aggregate, same, total) in enumerate(rows):
            conn.execute(
                "INSERT INTO clinvar_variants (rsid,chromosome,position,"
                "ref,alt,gene,clinical_significance,"
                "condition,review_status,allele_id,aggregate,submitters_same,submitters_total) "
                "VALUES (?,?,?,?,?,?,?,'Synthetic condition',"
                "'criteria provided, single submitter',?,?,?,?)",
                (rsid, chrom, pos, ref, alt, gene, sig, 900252800 + i, aggregate, same, total),
            )
        _stamp(conn, clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}")
        conn.commit()


def _identity(data, rows):
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://identity",
        downloaded_at="2001-01-01",
    )


def _pharm(data, rsid, gene, genotypes, *, level="1A"):
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.executescript(PHARMGKB_SCHEMA)
        for i, genotype in enumerate(genotypes):
            conn.execute(
                "INSERT INTO pharmgkb_annotations (rsid,genotype,gene,"
                "drugs,phenotype,phenotype_category,"
                "annotation_text,level_of_evidence,score,"
                "pgkb_annotation_id,function_class,is_nonfinding) "
                "VALUES (?,?,?,'Synthetic drug','Synthetic response',"
                "'Efficacy','Synthetic finding',?,10,?,"
                "'decreased',0)",
                (rsid, genotype, gene, level, str(900252800 + i)),
            )
        _stamp(conn, "pharmgkb", len(genotypes), f"iv:{PHARMGKB_INTERPRETER_VERSION}")
        conn.commit()


def _vcf(path, rows, *, block=False):
    header = "##fileformat=VCFv4.2\n"
    if block:
        header += '##ALT=<ID=NON_REF,Description="Any non-reference allele">\n'
    path.write_text(
        header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(
            f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{gt}\n"
            for rsid, chrom, pos, ref, alt, gt, info in rows
        )
    )
    return path


def _run(tmp_path, data, source, build, surface, *, command="analyze", extra=(), label="report"):
    args = [command, str(source), "--build", build.lower(), "--data-dir", str(data), "--no-update"]
    output = tmp_path / f"{label}.{surface}"
    if surface != "terminal":
        args += ["--output", str(output)]
    result = CliRunner(env={"COLUMNS": "500"}).invoke(main, [*args, *extra])
    assert result.exit_code == 0, (result.output, result.exception)
    if surface == "json":
        return json.loads(output.read_text()), result.output
    if surface == "html":
        raw = output.read_text()
        # A static guide or a hidden JSON payload is not a visible named disclosure.
        visible = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", raw, flags=re.S | re.I)
        return html.unescape(re.sub(r"<[^>]+>", " ", visible)), result.output
    return result.output, result.output


def _gnomad_site(data, rsid, chrom, pos, ref, alt, af):
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            (chrom, pos, ref, alt, rsid, af),
        )
        _stamp(conn, "gnomad", 1, f"sv:{GNOMAD_SCHEMA_VERSION}")
        conn.commit()


def _f5(data, installed, *, with_gnomad=False):
    data.mkdir()
    _identity(data, [("rs6025", "1", 169519049, "T", "C", 169549811, "C", "T", 0)])
    _pharm(data, "rs6025", "F5", ("CC", "CT", "TT"), level="2B")
    if with_gnomad:
        _gnomad_site(data, "rs6025", "1", 169549811, "C", "T", 0.01)
    for build in installed:
        _clinvar(
            data,
            build,
            [
                (
                    "rs6025",
                    "1",
                    169519049 if build == "GRCh37" else 169549811,
                    "T" if build == "GRCh37" else "C",
                    "T",
                    "F5",
                    "Pathogenic",
                    "Pathogenic",
                    1,
                    1,
                )
            ],
        )


@pytest.mark.parametrize(
    "installed,with_gnomad",
    [((), False), ((), True), (("GRCh37",), False), (("GRCh38",), False), (BUILDS, False)],
    ids=["none-no-gnomad", "none-with-gnomad", "only37", "only38", "both"],
)
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("genotype", ["CC", "CT"])
def test_f5_noncarrier_filter_does_not_depend_on_installed_clinvar_builds(
    tmp_path, installed, with_gnomad, batch, genotype
):
    data = tmp_path / "data"
    _f5(data, installed, with_gnomad=with_gnomad)
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, clinvar_builds=installed or BUILDS)
        for source in registry:
            stack.enter_context(source)
        pharm = next(a for a in registry if a.name == "pharmgkb")
        assert pharm.is_ready()
        variant = Variant(
            "rs6025", "1", 169519049, *genotype, ref="T", build="GRCh37", ref_from_input=True
        )
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    normal_known = bool(installed) or with_gnomad
    assert len(rows) == (0 if genotype == "CC" and normal_known else 1), rows
    if rows:
        assert rows[0].significance == "pharmgkb_loe_2b"
        if normal_known:
            assert rows[0].alt == "T"


@pytest.mark.parametrize(
    "installed,with_gnomad",
    [((), False), ((), True), (("GRCh37",), False), (("GRCh38",), False), (BUILDS, False)],
    ids=["none-no-gnomad", "none-with-gnomad", "only37", "only38", "both"],
)
@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics"])
def test_cli_f5_noncarrier_is_suppressed_only_with_normal_allele_evidence(
    tmp_path, installed, with_gnomad, command
):
    data = tmp_path / "data"
    _f5(data, installed, with_gnomad=with_gnomad)
    source = _vcf(tmp_path / "input.vcf", [("rs6025", "1", 169519049, "T", "C", "1/1", ".")])
    report, _ = _run(
        tmp_path, data, source, "GRCh37", "json", command=command, extra=("--min-magnitude", "0")
    )
    rows = [row for row in report["annotations"] if row["source"] == "pharmgkb"]
    assert len(rows) == (0 if installed or with_gnomad else 1), rows


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("genotype", ["CC", "CT", "TT"])
def test_grch38_reference_is_not_guessed_normal_when_gnomad_says_it_is_minor(
    tmp_path, build, batch, genotype
):
    data = tmp_path / "data"
    data.mkdir()
    rsid = "rs900252841"
    _identity(data, [(rsid, "1", 200, "T", "C", 1000, "C", "T", 0)])
    _pharm(data, rsid, "SYNTHETIC", ("CC", "CT", "TT"), level="2B")
    _gnomad_site(data, rsid, "1", 1000, "C", "T", 0.9)
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, clinvar_builds=BUILDS)
        for source in registry:
            stack.enter_context(source)
        pharm = next(a for a in registry if a.name == "pharmgkb")
        variant = Variant(
            rsid,
            "1",
            200 if build == "GRCh37" else 1000,
            *genotype,
            ref="T" if build == "GRCh37" else "C",
            build=build,
            ref_from_input=True,
        )
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    assert len(rows) == (0 if genotype == "TT" else 1), rows
    if rows:
        assert rows[0].alt == "C"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("genotype", ["AA", "AG", "GG"])
def test_equal_reference_sites_keep_existing_clinpgx_rows_without_clinvar(
    tmp_path, build, batch, genotype
):
    data = tmp_path / "data"
    data.mkdir()
    rsid = "rs900252842"
    _identity(data, [(rsid, "1", 200, "A", "G", 1000, "A", "G", 0)])
    _pharm(data, rsid, "SYNTHETIC", ("AA", "AG", "GG"))
    _gnomad_site(data, rsid, "1", 1000, "A", "G", 0.9)
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, clinvar_builds=BUILDS)
        for source in registry:
            stack.enter_context(source)
        pharm = next(a for a in registry if a.name == "pharmgkb")
        variant = Variant(
            rsid,
            "1",
            200 if build == "GRCh37" else 1000,
            *genotype,
            ref="A",
            build=build,
            ref_from_input=True,
        )
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    assert len(rows) == 1, rows


def _mthfr(data, build):
    data.mkdir()
    p37, p38 = 11856378, 11796321
    pos = p37 if build == "GRCh37" else p38
    _identity(data, [("rs1801133", "1", p37, "G", "A", p38, "G", "A", 0)])
    _clinvar(
        data, build, [("rs1801133", "1", pos, "G", "A", "MTHFR", "Pathogenic", "Pathogenic", 1, 1)]
    )
    _pharm(data, "rs1801133", "MTHFR", ("AG",))
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.executescript(_STRUCTURED_SCHEMA)
        conn.execute(
            "INSERT INTO snpedia_genotypes (rsid,allele1,allele2,magnitude,"
            "repute,summary,gene,orientation) "
            "VALUES ('rs1801133','A','G',6,'Bad','Synthetic interpretation','MTHFR','plus')"
        )
        _stamp(conn, "snpedia", 1, f"pv:{_PARSER_VERSION}")
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.executescript(GWAS_SCHEMA)
        conn.execute(
            "INSERT INTO gwas_associations (rsid,risk_allele,trait,p_value,"
            "gene,study_accession,trait_category) "
            "VALUES ('rs1801133','A','Synthetic condition',1e-30,'MTHFR',"
            "'GCST900252800','disease')"
        )
        _stamp(conn, "gwas", 1, f"cv:{_CATEGORIZER_VERSION}")
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,"
            "af) VALUES ('1',?,'G','A','rs1801133',1e-5)",
            (p38,),
        )
        _stamp(conn, "gnomad", 1, f"sv:{GNOMAD_SCHEMA_VERSION}")
        conn.commit()
    return pos


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["suffix", "prefix"])
@pytest.mark.parametrize("source_name", ["clinvar", "pharmgkb", "snpedia", "gwas"])
def test_padded_vcf_alleles_get_the_same_source_findings(tmp_path, build, shape, source_name):
    data = tmp_path / "data"
    pos = _mthfr(data, build)
    normal = _vcf(tmp_path / "normal.vcf", [("rs1801133", "1", pos, "G", "A", "0/1", ".")])
    padded = _vcf(
        tmp_path / "padded.vcf",
        [
            (
                "rs1801133",
                "1",
                pos if shape == "suffix" else pos - 1,
                "GC" if shape == "suffix" else "CG",
                "AC" if shape == "suffix" else "CA",
                "0/1",
                ".",
            )
        ],
    )
    extra = ("--min-magnitude", "0", "--gwas-min-magnitude", "0")
    expected, _ = _run(tmp_path, data, normal, build, "json", extra=extra, label="normal")
    actual, _ = _run(tmp_path, data, padded, build, "json", extra=extra, label="padded")
    before = [r for r in expected["annotations"] if r["source"] == source_name]
    after = [r for r in actual["annotations"] if r["source"] == source_name]
    assert len(before) == 1, expected["annotations"]
    fields = (
        "source",
        "rsid",
        "alt",
        "condition",
        "description",
        "magnitude",
        "allele_frequency",
        "plausibility",
    )
    assert [tuple(r.get(k) for k in fields) for r in after] == [
        tuple(r.get(k) for k in fields) for r in before
    ]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["1/2", "2/1", "1|2", "2|1", "1/1"])
def test_a_second_alt_equal_to_the_trimmed_alt_is_not_a_second_copy(tmp_path, build, gt):
    data = tmp_path / "data"
    pos = _mthfr(data, build)
    source = _vcf(tmp_path / "input.vcf", [("rs1801133", "1", pos, "GC", "AC,A", gt, ".")])
    report, _ = _run(tmp_path, data, source, build, "json", extra=("--min-magnitude", "0"))
    rows = [r for r in report["annotations"] if r["source"] == "clinvar" and r["alt"] == "A"]
    assert len(rows) == 1, rows
    assert rows[0]["allele_frequency"] == pytest.approx(1e-5)
    assert rows[0]["zygosity"] == ("homozygous" if gt == "1/1" else "heterozygous")
    assert rows[0]["plausibility"] == ("implausible" if gt == "1/1" else "plausible"), rows[0]


@pytest.mark.parametrize("surface", SURFACES)
@pytest.mark.parametrize("present", [False, True], ids=["absent", "covered"])
def test_missing_reference_carried_high_value_site_is_named(tmp_path, surface, present):
    data = tmp_path / "data"
    data.mkdir()
    rsid, gene, pos = "rs776746", "CYP3A5", 99270539
    _identity(data, [(rsid, "7", pos, "C", "T", 99672916, "T", "C", 0)])
    _clinvar(data, "GRCh37", [(rsid, "7", pos, "C", "C", gene, "Pathogenic", "Pathogenic", 1, 1)])
    rows = [("rs900252801", "1", 1000, "A", "G", "0/1", ".")]
    if present:
        rows.append((rsid, "7", pos, "C", "T", "0/1", "."))
    source = _vcf(tmp_path / "input.vcf", rows)
    report, _ = _run(tmp_path, data, source, "GRCh37", surface)
    if surface == "json":
        assert report["input"]["vcf_quality"]["reference_carried_not_in_file"] == (
            0 if present else 1
        )
        text = json.dumps({k: v for k, v in report.items() if k != "annotations"})
    else:
        text = report
    if not present:
        assert rsid in text and gene in text, text
        assert re.search(r"not (?:in|covered)|missing|no record|absent|could not", text, re.I), (
            text
        )
    elif surface == "json":
        warnings = json.dumps(
            report.get("report_warnings", []) + report.get("high_value_no_calls", [])
        )
        assert rsid not in warnings


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("surface", SURFACES)
@pytest.mark.parametrize("floor", [None, "0"], ids=["default-floor", "unfiltered"])
def test_conflicting_capped_variant_is_named_with_split_without_rescoring(
    tmp_path, build, surface, floor
):
    data = tmp_path / "data"
    data.mkdir()
    rsid, gene = "rs75391579", "GALT"
    pos = 1000 if build == "GRCh37" else 2000
    _clinvar(
        data,
        build,
        [
            (
                rsid,
                "9",
                pos,
                "A",
                "G",
                gene,
                "Pathogenic",
                "Conflicting classifications of pathogenicity",
                40,
                41,
            )
        ],
    )
    source = _vcf(tmp_path / "input.vcf", [(rsid, "9", pos, "A", "G", "1/1", ".")])
    report, _ = _run(
        tmp_path,
        data,
        source,
        build,
        surface,
        extra=() if floor is None else ("--min-magnitude", floor),
    )
    if surface == "json":
        assert report["clinvar"]["pathogenic_rows_scored_down"] == 1
        if floor is None:
            assert report["annotations"] == []
        else:
            assert len(report["annotations"]) == 1
            row = report["annotations"][0]
            assert row["magnitude"] == 4.0 and row["aggregate_capped"] is True
        text = json.dumps(report)
    else:
        text = report
    assert rsid in text and gene in text, text
    assert "Conflicting" in text or "conflicting" in text
    assert re.search(r"\b40\b", text) and re.search(r"\b41\b", text), text


# Public variant labels and synthetic calls, using the previously probed chip envelopes.
INDELS = (
    ("rs80359550", "BRCA2", "6174delT", "13", 32914438, 32914437, "AT", "A", "II"),
    ("rs80357906", "BRCA1", "5382insC", "17", 41209083, 41209079, "A", "AC", "DD"),
    ("rs113993960", "CFTR", "F508del", "7", 117199646, 117199644, "ACTT", "A", "II"),
)


def _chip(path, calls, kind):
    if kind == "23andme":
        path.write_text(
            "# This data file generated by 23andMe\n# rsid\tchromosome\tposition\tgenotype\n"
            + "".join(f"{rsid}\t{chrom}\t{pos}\t{gt}\n" for rsid, chrom, pos, gt in calls)
        )
    else:
        path.write_text(
            "#AncestryDNA raw data download\nrsid\tchromosome\tposition\tallele1\tallele2\n"
            + "".join(
                f"{rsid}\t{chrom}\t{pos}\t{gt[0]}\t{gt[1]}\n" for rsid, chrom, pos, gt in calls
            )
        )
    return path


def _named_unchecked(text, rsid, gene, label, call):
    assert rsid in text and gene in text and label in text, text
    # Alleles may be formatted with the input delimiter or a display slash, but not reinterpreted.
    assert re.search(rf"\b{call[0]}\s*[/,|\t ]?\s*{call[1]}\b", text), text
    assert re.search(
        r"could not|cannot|unable|not (?:assess|check|interpret)|uninterpret", text, re.I
    ), text


def _high_value_section(tmp_path, surface, report):
    if surface == "json":
        return json.dumps({k: v for k, v in report.items() if k.startswith("high_value")})
    if surface == "html":
        raw = (tmp_path / "report.html").read_text()
        sections = []
        for match in re.finditer(r"<(div|section|aside)\b[^>]*>(.*?)</\1>", raw, re.S | re.I):
            body = match.group(2)
            if re.search(r"<(?:strong|h[1-6])\b[^>]*>[^<]*high[- ]value", body, re.I):
                sections.append(html.unescape(re.sub(r"<[^>]+>", " ", body)))
        return " ".join(sections)
    # Console disclosures follow the high-value heading; unrelated earlier warnings do not count.
    match = re.search(r"high[- ]value[^\n]*\n", report, re.I)
    return report[match.start() :] if match else ""


@pytest.mark.parametrize("kind", ["23andme", "ancestrydna"])
@pytest.mark.parametrize("coordinate", ["chip", "identity"])
@pytest.mark.parametrize("carrier", [True, False], ids=["carrier", "reference-control"])
@pytest.mark.parametrize("surface", SURFACES)
def test_uninterpretable_indel_calls_are_named_in_warnings_and_high_value_section(
    tmp_path, kind, coordinate, carrier, surface
):
    data = tmp_path / "data"
    data.mkdir()
    cache_rows, calls = [], []
    for rsid, gene, _label, chrom, chip_pos, identity_pos, ref, alt, control in INDELS:
        cache_rows.append(
            (rsid, chrom, identity_pos, ref, alt, gene, "Pathogenic", "Pathogenic", 1, 1)
        )
        calls.append(
            (
                rsid,
                chrom,
                chip_pos if coordinate == "chip" else identity_pos,
                "DI" if carrier else control,
            )
        )
    _clinvar(data, "GRCh37", cache_rows)
    source = _chip(tmp_path / "input.txt", calls, kind)
    report, console = _run(tmp_path, data, source, "GRCh37", surface)
    if surface == "json":
        assert report["annotations"] == []
        text = json.dumps(report.get("report_warnings", []))
        high_value = _high_value_section(tmp_path, surface, report)
    else:
        text = report
        high_value = _high_value_section(tmp_path, surface, report)
        assert re.search(r"high.value", text, re.I), text
    for rsid, gene, label, _chrom, _chip_pos, _identity_pos, _ref, _alt, control in INDELS:
        call = "DI" if carrier else control
        _named_unchecked(text, rsid, gene, label, call)
        _named_unchecked(high_value, rsid, gene, label, call)
        _named_unchecked(console, rsid, gene, label, call)


@pytest.mark.parametrize("kind", ["23andme", "ancestrydna"])
@pytest.mark.parametrize("surface", SURFACES)
@pytest.mark.parametrize("selection", ["high-value-only", "plp-only", "unlisted-non-plp"])
def test_indel_naming_uses_the_union_of_high_value_and_plp_sites(
    tmp_path, kind, surface, selection
):
    data = tmp_path / "data"
    data.mkdir()
    if selection == "high-value-only":
        # No ClinVar installation: the built-in high-value list is sufficient to name this call.
        rsid, gene, label = "rs80359550", "BRCA2", "6174delT"
        _pharm(data, "rs900252899", "GENE", ("AG",))
    else:
        rsid, gene, label = "rs900252898", "GENE", "rs900252898"
        aggregate = "Pathogenic" if selection == "plp-only" else "Uncertain significance"
        _clinvar(data, "GRCh37", [(rsid, "2", 1000, "AT", "A", gene, aggregate, aggregate, 1, 1)])
    source = _chip(tmp_path / "input.txt", [(rsid, "2", 1000, "DI")], kind)
    report, console = _run(tmp_path, data, source, "GRCh37", surface)
    if surface == "json":
        assert report["annotations"] == []
        text = json.dumps(report.get("report_warnings", []))
        high_value = _high_value_section(tmp_path, surface, report)
    else:
        text = report
        high_value = _high_value_section(tmp_path, surface, report)
    if selection == "unlisted-non-plp":
        # The aggregate D/I count may remain, but an unrelated site is not a high-value alert.
        assert rsid not in text
    else:
        _named_unchecked(text, rsid, gene, label, "DI")
        _named_unchecked(high_value, rsid, gene, label, "DI")
        _named_unchecked(console, rsid, gene, label, "DI")

"""A risk base belongs to SNPS, before lowest-p selection; synthetic TSV only."""

import contextlib
import csv
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _identity

from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.cli import main
from allelix.databases import gwas_loader as loader
from allelix.databases.gwas_loader import (
    gwas_source_hash,
    iter_gwas_records,
    load_gwas_tsv,
    retained_source_stamp,
)

RS = "rs900541001"
OTHER = "rs900541002"
TRAIT = "Synthetic vascular disease"
FIELDS = (
    "SNPS",
    "STRONGEST SNP-RISK ALLELE",
    "DISEASE/TRAIT",
    "P-VALUE",
    "OR or BETA",
    "95% CI (TEXT)",
    "MAPPED_GENE",
    "STUDY ACCESSION",
    "PUBMEDID",
)


def _row(risk, *, snp=RS, trait=TRAIT, p="1e-20", study="GCST900541001"):
    return dict(
        zip(
            FIELDS,
            (snp, risk, trait, p, "1.8", "[1.5-2.1]", "SYNGENE", study, "900541001"),
            strict=True,
        )
    )


def _tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.mark.parametrize(
    "prefix",
    (
        OTHER,
        "RS900541002",
        "not-an-rsid",
        "",
        RS + ";" + OTHER,
        RS + "-A",
        "rs900541001suffix",
        "chr1:1000",
    ),
)
@pytest.mark.parametrize("base", ("A", "C", "G", "T"))
@pytest.mark.parametrize("surface", ("iterator", "sqlite"))
def test_foreign_or_unattributed_risk_base_is_not_an_association(tmp_path, prefix, base, surface):
    source = _tsv(tmp_path / "source.tsv", [_row(prefix + "-" + base)])
    if surface == "iterator":
        assert list(iter_gwas_records(source)) == [], "Never attach another SNP's base to SNPS"
    else:
        target = tmp_path / "gwas.sqlite"
        count = load_gwas_tsv(source, target)
        with contextlib.closing(sqlite3.connect(target)) as conn:
            rows = conn.execute("SELECT rsid, risk_allele FROM gwas_associations").fetchall()
        assert rows == [] and count == 0, "Do not turn a rejected row into a riskless finding"


@pytest.mark.parametrize("prefix", (RS, RS.upper(), " " + RS + " ", " " + RS.upper() + " "))
@pytest.mark.parametrize("base", ("A", "c", " G ", "t"))
def test_own_rsid_is_compared_trimmed_and_case_insensitively(tmp_path, prefix, base):
    source = _tsv(tmp_path / "source.tsv", [_row(prefix + "-" + base, snp=" " + RS + " ")])
    rows = list(iter_gwas_records(source))
    assert len(rows) == 1
    assert (rows[0]["rsid"], rows[0]["risk_allele"]) == (RS, base.strip().upper())


@pytest.mark.parametrize("bad", (OTHER + "-G", "other-G", RS + "-A-G"))
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("valid_risk", (RS + "-A", RS + "-?"))
def test_rejected_lowest_p_row_cannot_hide_the_valid_fallback(tmp_path, bad, reverse, valid_risk):
    rows = [_row(bad, p="1e-50", study="BAD"), _row(valid_risk, p="1e-20", study="GOOD")]
    source = _tsv(tmp_path / "source.tsv", rows[::-1] if reverse else rows)
    records = list(iter_gwas_records(source))
    assert len(records) == 1
    assert records[0]["study_accession"] == "GOOD", "Filter attribution before deduplicating"
    assert records[0]["risk_allele"] == ("A" if valid_risk.endswith("-A") else None)
    assert records[0]["p_value"] == 1e-20


@pytest.mark.parametrize(
    "risk", ("", "?", RS, "A", RS + "-?", OTHER + "-?", OTHER + "-AT", "other-", "other-N")
)
def test_missing_or_non_base_risk_keeps_the_existing_riskless_contract(tmp_path, risk):
    source = _tsv(tmp_path / "source.tsv", [_row(risk)])
    records = list(iter_gwas_records(source))
    assert len(records) == 1 and records[0]["risk_allele"] is None


@pytest.mark.parametrize("snp,trait", ((RS + ";" + OTHER, TRAIT), (RS, ""), ("not-rsid", TRAIT)))
def test_existing_haplotype_identifier_and_no_trait_filters_stay(tmp_path, snp, trait):
    source = _tsv(tmp_path / "source.tsv", [_row(RS + "-A", snp=snp, trait=trait)])
    assert list(iter_gwas_records(source)) == []


def test_removed_rows_do_not_count_toward_the_publication_floor(tmp_path):
    target = tmp_path / "gwas.sqlite"
    old = _tsv(tmp_path / "old.tsv", [_row(RS + "-A")])
    load_gwas_tsv(old, target)
    before = target.read_bytes()
    source = _tsv(
        tmp_path / "new.tsv", [_row(RS + "-A"), _row(OTHER + "-G", trait="Other disease")]
    )
    with pytest.raises(OSError, match="only 1"):
        load_gwas_tsv(source, target, min_rows=2)
    assert target.read_bytes() == before


@pytest.mark.parametrize("retained", (False, True))
def test_a_released_cv3_cache_is_outdated_and_reingests_its_owned_table(
    tmp_path, monkeypatch, retained
):
    from allelix.annotators import gwas as module

    source = _tsv(
        tmp_path / "gwas_catalog_associations.tsv",
        [_row(OTHER + "-G"), _row(RS + "-A", p="1e-10", study="GOOD")],
    )
    target = tmp_path / "gwas.sqlite"
    digest = gwas_source_hash(source)
    load_gwas_tsv(source, target, remote_signal="etag:synthetic-release", content_hash=digest)
    with contextlib.closing(sqlite3.connect(target)) as conn, conn:
        conn.execute("UPDATE database_versions SET local_version_tag='cv:3'")
        conn.execute("UPDATE gwas_associations SET risk_allele='G', study_accession='BAD'")
    if not retained:
        source.unlink()
    before = target.read_bytes()
    monkeypatch.setattr(module, "GWAS_MIN_ROWS", 0)

    def no_network(*args, **kwargs):
        pytest.fail("The retained, hash-matching source needs no download")

    monkeypatch.setattr(module, "download", no_network)
    with GWASCatalogAnnotator(tmp_path) as annotator:
        assert not annotator.is_ready(), "The version must invalidate cv:3's old parse"
        assert target.read_bytes() == before, "Readiness is read-only"
        if not retained:
            assert not annotator.prepare_cache()
            assert not annotator.is_ready()
            assert target.read_bytes() == before
            return
        assert annotator.prepare_cache()
        assert annotator.is_ready()
    with contextlib.closing(sqlite3.connect(target)) as conn:
        rows = conn.execute(
            "SELECT risk_allele, study_accession FROM gwas_associations"
        ).fetchall()
        tag = conn.execute("SELECT local_version_tag FROM database_versions").fetchone()[0]
    assert rows == [("A", "GOOD")]
    assert tag != "cv:3" and loader._CATEGORIZER_VERSION > 3
    assert retained_source_stamp(target) == ("etag:synthetic-release", digest)


@pytest.mark.parametrize(
    "build,length,pos", (("grch37", 249250621, 1000), ("grch38", 248956422, 2000))
)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("bad", (OTHER + "-A", "other-A"))
def test_rejected_rows_never_reach_a_gwas_finding_on_either_build(
    tmp_path, build, length, pos, explicit, bad
):
    data = tmp_path / "data"
    data.mkdir()
    source = _tsv(
        tmp_path / "source.tsv",
        [
            _row(bad, trait="Rejected synthetic disease"),
            _row(RS + "-A", trait="Kept synthetic disease"),
        ],
    )
    load_gwas_tsv(source, data / "gwas.sqlite")
    _identity(data, [(RS, "1", 1000, "G", "A", 2000, "G", "A", 0)])
    vcf = tmp_path / "call.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n1\t{pos}\t{RS}\tG\tA\t.\tPASS\t.\tGT\t1/1\n"
    )
    output = tmp_path / "result.json"
    args = [
        "analyze",
        str(vcf),
        "--data-dir",
        str(data),
        "--no-update",
        "--gwas-all",
        "--gwas-min-magnitude",
        "0",
        "--min-magnitude",
        "0",
        "--output",
        str(output),
    ]
    if explicit:
        args += ["--build", build]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    report = json.loads(output.read_text())
    assert report["input"]["build"].lower() == build
    rows = [r for r in report["annotations"] if r["source"] == "gwas"]
    assert len(rows) == 1 and rows[0]["trait"] == "Kept synthetic disease"

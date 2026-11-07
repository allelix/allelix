# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The ClinPGx/SNPedia hom-ref filter takes ClinVar's REF from the row at the call's locus.

ADR-0023 (amendment 2026-10-16) and ADR-0039 (2.6.0 amendment): ClinVar's REF
suppresses a ClinPGx or SNPedia row only where ClinVar places the rsID at the
call's own position. `reference_at` checked that a row is there, but returned
the REF of whichever one-base row of the rsID it read first, an insertion
anchor included. The B-01 tests and their fixtures are reviewer B's (outside
review of 2.6.0 RC4, finding B-01), unchanged. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import (
    ClinVarAnnotator,
    clinvar_db_filename,
    clinvar_record_name,
)
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION, PHARMGKB_INTERPRETER_VERSION
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION, GWAS_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GWAS_SCHEMA, PHARMGKB_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION

if TYPE_CHECKING:
    from pathlib import Path

# ── synthetic caches ─────────────────────────────────────────────────────

_GRCH37 = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}  # fmt: skip
_GRCH38 = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
    "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
    "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}  # fmt: skip


def _cv(rsid, chrom, pos, ref, alt, *, sig="Pathogenic", agg="Pathogenic", allele_id=1):
    return (
        rsid, chrom, pos, ref, alt, sig, "Synthetic_disease", "SYNGENE",
        "criteria_provided,_single_submitter", allele_id, agg, 1, 1, 0, 0,
    )  # fmt: skip


def _write_clinvar(data: Path, build: str, rows) -> None:
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total, outside_aggregate, normal_allele) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) VALUES (?, 'test://cv', "
            "'20260101', '2026-06-08', ?, 'md5:x', ?)",
            (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _write_identity(data: Path, rows) -> None:
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2026-10-16T00:00:00Z",
    )


def _write_pharm(data: Path, rows) -> None:
    """rows: (rsid, genotype, annotation id)."""
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.executescript(PHARMGKB_SCHEMA)
        conn.executemany(
            "INSERT INTO pharmgkb_annotations (rsid, genotype, gene, drugs, phenotype, "
            "phenotype_category, annotation_text, level_of_evidence, score, "
            "pgkb_annotation_id, function_class, is_nonfinding) "
            "VALUES (?, ?, 'SYNGENE', 'Synthetic drug', 'Synthetic response', 'Efficacy', "
            "'Synthetic finding', '1A', 10, ?, 'decreased', 0)",
            rows,
        )
        conn.execute(
            "INSERT INTO pharmgkb_allele_function (rsid, allele, function_class) "
            "VALUES ('rs900990000', 'A', 'normal')"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) VALUES ('pharmgkb', "
            "'test://pharmgkb', 'test', '2026-10-10T00:00:00Z', ?, 'pgkb:etag:x|cpic:2026', ?)",
            (max(len(rows), 1), f"iv:{PHARMGKB_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _write_snpedia(data: Path, rows) -> None:
    """rows: (rsid, allele1, allele2, magnitude, summary, orientation)."""
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.executescript(
            "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
            "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
            "orientation TEXT, scraped_at TEXT);"
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
            "version TEXT, downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
            "remote_signal TEXT, local_version_tag TEXT);"
        )
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, 'Bad', ?, 'SYNGENE', ?, "
            "'2026-10-10T00:00:00')",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions VALUES ('snpedia', 'test://snpedia', 'test', "
            "'2026-10-10T00:00:00Z', ?, '', ?)",
            (len(rows), f"pv:{_PARSER_VERSION}"),
        )
        conn.commit()


def _write_empty_gwas(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / GWAS_DB_FILENAME)) as conn:
        conn.executescript(GWAS_SCHEMA)
        conn.execute(
            "INSERT INTO gwas_associations (rsid, risk_allele, trait, p_value, trait_category) "
            "VALUES ('rs900990001', 'A', 'filler', 1e-9, 'disease')"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) VALUES ('gwas', 'test://gwas', "
            "'test', '2026-10-16T00:00:00Z', 1, '', ?)",
            (f"cv:{_CATEGORIZER_VERSION}",),
        )
        conn.commit()


def _vcf(path: Path, records, *, build: str) -> Path:
    lengths = _GRCH38 if build == "GRCh38" else _GRCH37
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={c},length={n}>" for c, n in lengths.items()]
    lines += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    lines += ['##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">']
    lines += ["\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                         "FORMAT", "S1"])]  # fmt: skip
    lines += ["\t".join(str(x) for x in r) for r in records]
    path.write_text("\n".join(lines) + "\n")
    return path


def _analyze(vcf: Path, data: Path, *extra: str) -> dict:
    out = vcf.with_suffix(".json")
    result = CliRunner().invoke(
        main,
        ["analyze", str(vcf), "--data-dir", str(data), "--no-update", "--output", str(out),
         "--min-magnitude", "0", *extra],
        catch_exceptions=False,
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    return json.loads(out.read_text())


def _vcf_out(vcf: Path, data: Path) -> list[list[str]]:
    out = vcf.with_suffix(".annotated.vcf")
    result = CliRunner().invoke(
        main,
        ["analyze", str(vcf), "--data-dir", str(data), "--no-update", "--vcf-out", str(out)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    return [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]


# ── B-01 ─────────────────────────────────────────────────────────────────
# ADR-0023 amendment 2026-10-16 / ADR-0039 amendment 2.6.0: the ClinPGx and
# SNPedia non-finding filter uses ClinVar's REF "only where ClinVar places the
# rsID at the call's own (chromosome, position)" (CHANGELOG 2.5: "Before, an
# rsID placed elsewhere could hide a finding"). `reference_at` checks that a
# row exists at the locus, but returns `reference_for(rsid)`, the REF of the
# first single-base row of the rsID anywhere on the build. When ClinVar holds
# the rsID at two loci with different REFs, a call that is homozygous for the
# ALT at one locus equals the other locus's REF and is suppressed.

B01 = "rs900910001"
_B01_OTHER = _cv(B01, "1", 5000, "T", "C", sig="Uncertain_significance",
                 agg="Uncertain significance", allele_id=91)  # fmt: skip
_B01_HERE = _cv(B01, "1", 1000, "C", "T", sig="Uncertain_significance",
                agg="Uncertain significance", allele_id=92)  # fmt: skip


@pytest.fixture
def b01_data(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    # ClinVar lists the rsID at two GRCh38 loci; the other locus's row is read first.
    _write_clinvar(data, "GRCh38", [_B01_OTHER, _B01_HERE])
    _write_clinvar(data, "GRCh37", [_cv("rs900990002", "2", 10, "A", "G")])
    # The identity table places the rsID at the call's locus (no ADR-0063 conflict).
    _write_identity(data, [(B01, "1", 900, "C", "T", 1000, "C", "T", 0)])
    _write_pharm(data, [(B01, "TT", "pgkb-b01-tt")])
    _write_snpedia(data, [(B01, "T", "T", 4.0, "synthetic TT page", "plus")])
    _write_empty_gwas(data)
    return data


def test_b01_reference_at_gives_the_ref_at_the_calls_own_locus(b01_data):
    """ClinVar's REF at 1:1000 is C; the T is the REF of the same rsID at 1:5000."""
    with ClinVarAnnotator(b01_data) as clinvar:
        assert clinvar.reference_at(B01, "GRCh38", ("1", 1000)) == "C"
        assert clinvar.reference_at(B01, "GRCh38", ("1", 5000)) == "T"


def test_b01_homozygous_alt_call_keeps_its_clinpgx_and_snpedia_findings(b01_data, tmp_path):
    """C>T 1/1 at 1:1000 is homozygous for the ALT, not homozygous reference."""
    vcf = _vcf(
        tmp_path / "b01.vcf",
        [("1", 1000, B01, "C", "T", 50, "PASS", ".", "GT:GQ", "1/1:99")],
        build="GRCh38",
    )
    report = _analyze(vcf, b01_data)
    found = {(a["source"], a["rsid"], a["genotype_match"]) for a in report["annotations"]}
    assert ("clinvar", B01, "TT") in found  # control: ClinVar sees the same call
    assert ("pharmgkb", B01, "TT") in found
    assert ("snpedia", B01, "TT") in found


B01_INS = "rs900910002"


def test_b01_insertion_anchor_ref_does_not_hide_the_snv_finding(tmp_path):
    """One rsID with an insertion anchored at 1:999 (REF C, one base) and an SNV A>C at 1:1000.

    dbSNP rsIDs often cover an SNV and an indel at one site; ClinVar writes the
    insertion VCF-style, anchored one base earlier with a one-base REF. A call
    homozygous for the SNV's ALT (C/C) equals the anchor base and is hidden.
    """
    data = tmp_path / "data"
    data.mkdir()
    _write_clinvar(
        data,
        "GRCh38",
        [
            _cv(B01_INS, "1", 999, "C", "CT", sig="Uncertain_significance",
                agg="Uncertain significance", allele_id=94),
            _cv(B01_INS, "1", 1000, "A", "C", sig="Uncertain_significance",
                agg="Uncertain significance", allele_id=95),
        ],
    )  # fmt: skip
    _write_clinvar(data, "GRCh37", [_cv("rs900990002", "2", 10, "A", "G")])
    _write_identity(data, [(B01_INS, "1", 900, "A", "C", 1000, "A", "C", 0)])
    _write_pharm(data, [(B01_INS, "CC", "pgkb-b01-cc")])
    _write_snpedia(data, [(B01_INS, "C", "C", 4.0, "synthetic CC page", "plus")])
    _write_empty_gwas(data)
    with ClinVarAnnotator(data) as clinvar:
        assert clinvar.reference_at(B01_INS, "GRCh38", ("1", 1000)) == "A"
    vcf = _vcf(
        tmp_path / "b01i.vcf",
        [("1", 1000, B01_INS, "A", "C", 50, "PASS", ".", "GT:GQ", "1/1:99")],
        build="GRCh38",
    )
    report = _analyze(vcf, data)
    found = {(a["source"], a["genotype_match"]) for a in report["annotations"]}
    assert ("clinvar", "CC") in found  # control
    assert ("pharmgkb", "CC") in found
    assert ("snpedia", "CC") in found


def test_b01_control_single_locus_is_not_suppressed(tmp_path):
    """The same call with ClinVar holding only the call's own locus is reported (passes on RC4)."""
    data = tmp_path / "data"
    data.mkdir()
    _write_clinvar(data, "GRCh38", [_B01_HERE])
    _write_clinvar(data, "GRCh37", [_cv("rs900990002", "2", 10, "A", "G")])
    _write_identity(data, [(B01, "1", 900, "C", "T", 1000, "C", "T", 0)])
    _write_pharm(data, [(B01, "TT", "pgkb-b01-tt")])
    _write_snpedia(data, [(B01, "T", "T", 4.0, "synthetic TT page", "plus")])
    _write_empty_gwas(data)
    vcf = _vcf(
        tmp_path / "b01c.vcf",
        [("1", 1000, B01, "C", "T", 50, "PASS", ".", "GT:GQ", "1/1:99")],
        build="GRCh38",
    )
    report = _analyze(vcf, data)
    found = {(a["source"], a["genotype_match"]) for a in report["annotations"]}
    assert {("pharmgkb", "TT"), ("snpedia", "TT")} <= found


def test_a_true_hom_ref_call_is_still_suppressed(b01_data, tmp_path):
    """C/C at 1:1000, where ClinVar's REF is C, is homozygous reference, so its
    ClinPGx and SNPedia rows are hidden. Before, the REF read was the other
    locus's T, and this hom-ref call was reported as a finding."""
    import contextlib as _contextlib
    import sqlite3 as _sqlite3

    with _contextlib.closing(_sqlite3.connect(b01_data / "pharmgkb.sqlite")) as conn:
        conn.execute("UPDATE pharmgkb_annotations SET genotype = 'CC' WHERE rsid = ?", (B01,))
        conn.commit()
    with _contextlib.closing(_sqlite3.connect(b01_data / "snpedia.sqlite")) as conn:
        conn.execute(
            "UPDATE snpedia_genotypes SET allele1 = 'C', allele2 = 'C' WHERE rsid = ?", (B01,)
        )
        conn.commit()
    vcf = _vcf(
        tmp_path / "b01r.vcf",
        [("1", 1000, B01, "C", "T", 50, "PASS", ".", "GT:GQ", "0/0:99")],
        build="GRCh38",
    )
    report = _analyze(vcf, b01_data)
    sources = {a["source"] for a in report["annotations"] if a["rsid"] == B01}
    assert "pharmgkb" not in sources
    assert "snpedia" not in sources

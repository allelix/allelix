# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""--vcf-out stamps ALLELIX_CLINVAR_REF only on the record whose call produced the finding.

At a site where the reference carries ClinVar's variant allele, a hom-ref
record (a 0/0 duplicate, or a 0/0 deletion whose REF covers the site) beside a
call carrying REF gave no call of its own, so it borrows nothing
(docs/formats.md, Annotated VCF; ADR-0036, ADR-0059, ADR-0062).

The ``test_f1_*`` tests and the ``test_r1_lone_*`` control are reviewer A's
(outside review of 2.6.0 RC4, finding F1), with the fixture builders they use
from A's ``rc4_fixtures.py`` inlined, unchanged but for line wrapping. The
``.``-ALT case is reviewer B's B-02. Synthetic data only.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.rsid_identity_loader import write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"


def clinvar_db(data: Path, build: str, rows: list[dict]) -> None:
    """rows: dicts with rsid, chromosome, position, ref, alt, sig and optional extras."""
    data.mkdir(parents=True, exist_ok=True)
    path = data / clinvar_db_filename(build)
    path.unlink(missing_ok=True)
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        for i, r in enumerate(rows):
            conn.execute(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
                "submitters_same, submitters_total, outside_aggregate, normal_allele) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    r["rsid"],
                    r["chromosome"],
                    r["position"],
                    r["ref"],
                    r["alt"],
                    r.get("sig", "Pathogenic"),
                    r.get("condition", "Synthetic condition"),
                    r.get("gene", "TEST"),
                    r.get("review", "criteria_provided,_single_submitter"),
                    r.get("allele_id", 500000 + i),
                    r.get("aggregate", r.get("sig", "Pathogenic")),
                    1,
                    1,
                    0,
                    int(r.get("normal_allele", 0)),
                ),
            )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) VALUES (?, 'test://synthetic', "
            "'synthetic', '2026-11-04', ?, '', ?)",
            (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def identity_db(data: Path, rows: list[tuple]) -> None:
    """rows: (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic)."""
    data.mkdir(parents=True, exist_ok=True)
    write_identity_db(
        rows,
        data / "rsid_identity.sqlite",
        version="synthetic",
        source_url="test://synthetic",
        downloaded_at="2026-11-04",
    )


def write_vcf(path: Path, rows: list[str], meta: str = "") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("##fileformat=VCFv4.2\n" + meta + HEADER + "".join(r + "\n" for r in rows))
    return path


def info(cols: list[str], field: str) -> str | None:
    for token in cols[7].split(";"):
        if token.startswith(field + "="):
            return token.split("=", 1)[1]
    return None


def analyze_cli(
    data: Path, vcf: Path, out: Path, *, build: str | None = None, extra: tuple[str, ...] = ()
) -> dict:
    """Run ``allelix analyze`` end to end: JSON report and annotated VCF."""
    out.mkdir(parents=True, exist_ok=True)
    report, annotated = out / "report.json", out / "annotated.vcf"
    for p in (report, annotated):
        p.unlink(missing_ok=True)
    args = [
        "analyze",
        str(vcf),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(report),
        "--vcf-out",
        str(annotated),
        "--min-magnitude",
        "0",
        "--snpedia-min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--include-benign",
        *extra,
    ]
    if build:
        args += ["--build", build]
    res = CliRunner(env={"COLUMNS": "200"}).invoke(main, args, catch_exceptions=True)
    payload = json.loads(report.read_text()) if report.exists() else None
    rows = (
        [ln.split("\t") for ln in annotated.read_text().splitlines() if not ln.startswith("#")]
        if annotated.exists()
        else None
    )
    return {
        "exit_code": res.exit_code,
        "output": res.output,
        "exception": repr(res.exception),
        "json": payload,
        "vcf": rows,
    }


GRCH37_CHR1 = "##contig=<ID=1,length=249250621>\n"
F5 = "rs900006025"


# ── F1 (P1): an explicit hom-ref record borrows another call's ALLELIX_CLINVAR_REF ──


def _f5_like_data(data: Path) -> None:
    """GRCh37's reference carries ClinVar's variant allele T at 1:5000 (as at F5 rs6025)."""
    clinvar_db(
        data,
        "GRCh37",
        [dict(rsid=F5, chromosome="1", position=5000, ref="T", alt="T", allele_id=7)],
    )
    clinvar_db(
        data,
        "GRCh38",
        [dict(rsid=F5, chromosome="1", position=6000, ref="C", alt="T", allele_id=7)],
    )
    identity_db(data, [(F5, "1", 5000, "T", "C", 6000, "C", "T", 0)])


def _ref_slot(cols: list[str]) -> str | None:
    return info(cols, "ALLELIX_CLINVAR_REF")


HOM_REF = f"1\t5000\t{F5}\tT\tC\t.\tPASS\t.\tGT\t0/0"


@pytest.mark.parametrize("het_gt", ["0/1", "0/."])
@pytest.mark.parametrize("hom_ref_first", [True, False])
def test_f1_hom_ref_duplicate_does_not_borrow_the_het_calls_ref_finding(
    tmp_path, het_gt, hom_ref_first
):
    """The analysis reads the site from the het call (one CT/half-call finding);
    the 0/0 gets nothing."""
    _f5_like_data(tmp_path / "data")
    het = f"1\t5000\t{F5}\tT\tC\t.\tPASS\t.\tGT\t{het_gt}"
    rows = [HOM_REF, het] if hom_ref_first else [het, HOM_REF]
    vcf = write_vcf(tmp_path / "in.vcf", rows, GRCH37_CHR1)
    r = analyze_cli(tmp_path / "data", vcf, tmp_path / "out")
    assert r["exit_code"] == 0, r["output"]
    findings = [a for a in r["json"]["annotations"] if a["source"] == "clinvar"]
    # The analysis produced one finding, from the call carrying REF; none is homozygous.
    assert len(findings) == 1 and findings[0]["zygosity"] != "homozygous", findings
    by_gt = {cols[9]: cols for cols in r["vcf"]}
    assert _ref_slot(by_gt[het_gt]) is not None
    # formats.md: a hom-ref duplicate that gave no call "borrows nothing",
    # ALLELIX_CLINVAR_REF included.
    assert _ref_slot(by_gt["0/0"]) is None, by_gt["0/0"]


def test_f1_overlapping_hom_ref_deletion_record_does_not_borrow_the_snv_finding(tmp_path):
    """A 0/0 deletion-site record whose REF covers the site is not the record the
    het SNV finding came from."""
    _f5_like_data(tmp_path / "data")
    rows = ["1\t4999\t.\tAT\tA\t.\tPASS\t.\tGT\t0/0", f"1\t5000\t{F5}\tT\tC\t.\tPASS\t.\tGT\t0/1"]
    vcf = write_vcf(tmp_path / "in.vcf", rows, GRCH37_CHR1)
    r = analyze_cli(tmp_path / "data", vcf, tmp_path / "out")
    assert r["exit_code"] == 0, r["output"]
    assert [(a["source"], a["genotype_match"]) for a in r["json"]["annotations"]] == [
        ("clinvar", "CT")
    ]
    deletion, snv = r["vcf"]
    assert _ref_slot(snv) is not None
    assert _ref_slot(deletion) is None, deletion


def test_r1_lone_hom_ref_record_keeps_its_reference_finding(tmp_path):
    """Control for F1 (passes on RC4): a 0/0 record the analysis used is stamped (ADR-0059)."""
    _f5_like_data(tmp_path / "data")
    vcf = write_vcf(tmp_path / "in.vcf", [HOM_REF], GRCH37_CHR1)
    r = analyze_cli(tmp_path / "data", vcf, tmp_path / "out")
    assert [(a["genotype_match"], a["zygosity"]) for a in r["json"]["annotations"]] == [
        ("TT", "homozygous")
    ]
    assert _ref_slot(r["vcf"][0]) == f"clinvar_pathogenic|{F5}|TEST"


def test_b02_hom_ref_duplicate_with_dot_alt_borrows_nothing(tmp_path):
    """Reviewer B's B-02 case: the 0/0 duplicate writes ALT as ".", not C."""
    _f5_like_data(tmp_path / "data")
    rows = [
        f"1\t5000\t{F5}\tT\tC\t.\tPASS\t.\tGT\t0/1",
        f"1\t5000\t{F5}\tT\t.\t.\tPASS\t.\tGT\t0/0",
    ]
    vcf = write_vcf(tmp_path / "in.vcf", rows, GRCH37_CHR1)
    r = analyze_cli(tmp_path / "data", vcf, tmp_path / "out")
    assert r["exit_code"] == 0, r["output"]
    het, hom_ref = r["vcf"]
    assert _ref_slot(het) is not None
    assert _ref_slot(hom_ref) is None, hom_ref


def test_identical_hom_ref_duplicates_both_keep_the_reference_finding(tmp_path):
    """Control (documented rule): identical 0/0 rows the reference call was
    rebuilt from are both stamped."""
    _f5_like_data(tmp_path / "data")
    vcf = write_vcf(tmp_path / "in.vcf", [HOM_REF, HOM_REF], GRCH37_CHR1)
    r = analyze_cli(tmp_path / "data", vcf, tmp_path / "out")
    assert r["exit_code"] == 0, r["output"]
    assert [(a["genotype_match"], a["zygosity"]) for a in r["json"]["annotations"]] == [
        ("TT", "homozygous")
    ]
    assert [_ref_slot(cols) for cols in r["vcf"]] == [f"clinvar_pathogenic|{F5}|TEST"] * 2

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
#!/usr/bin/env python3
"""Generate per-build synthetic ClinVar VCFs for tests.

ADR-0021 + ADR-0015: every record carries authoritative GRCh37 and
GRCh38 coordinates. The generator emits one VCF per build with
build-correct positions and (for strand-inverted variants like NIPA1
rs104894490) build-correct REF/ALT. Mock-as-spec is enforced by
`tests/test_mock_data_invariants.py::test_mock_clinvar_positions_match_declared_build`.

Usage:
    python tests/generate_clinvar_fixture.py                 # both builds
    python tests/generate_clinvar_fixture.py --build grch37  # one build
    python tests/generate_clinvar_fixture.py --build grch38

Default output paths: tests/fixtures/mock_clinvar_grch37.vcf and
tests/fixtures/mock_clinvar_grch38.vcf.
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Each record carries:
#   chrom, GRCh37 pos, GRCh38 pos,
#   GRCh37 (ref, alt), GRCh38 (ref, alt) — usually identical but DIFFERENT for
#       variants whose strand inverted between assemblies (rs104894490 NIPA1),
#   rsid (digits only), CLNSIG, condition, gene, review status, allele_id
#
# Positions are sourced from NCBI dbSNP (cross-checked against the live
# ClinVar GRCh37 VCF). Conditions use underscores per real ClinVar VCF
# format. If a position appears in `allelix.utils.build_detect.KNOWN_SNP_POSITIONS`,
# the values here MUST agree — the invariant test fails the build otherwise.
RECORDS = [
    # MTHFR C677T — mock MHG fixture is G/A (heterozygous carrier).
    {
        "chrom": "1",
        "pos_grch37": 11856378,
        "pos_grch38": 11796321,
        "ref_alt_grch37": ("G", "A"),
        "ref_alt_grch38": ("G", "A"),
        "rsid": "1801133",
        "clnsig": "Pathogenic",
        "condition": "MTHFR_deficiency",
        "gene": "MTHFR",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100001,
    },
    # MTHFR A1298C — mock MHG fixture is T/T (no carrier; ALT=G unmatched).
    {
        "chrom": "1",
        "pos_grch37": 11854476,
        "pos_grch38": 11794419,
        "ref_alt_grch37": ("T", "G"),
        "ref_alt_grch38": ("T", "G"),
        "rsid": "1801131",
        "clnsig": "Likely_pathogenic",
        "condition": "Hyperhomocysteinemia",
        "gene": "MTHFR",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100002,
    },
    # MTRR A66G — mock MHG fixture is G/G (homozygous ALT carrier).
    {
        "chrom": "5",
        "pos_grch37": 7870860,
        "pos_grch38": 7870973,
        "ref_alt_grch37": ("A", "G"),
        "ref_alt_grch38": ("A", "G"),
        "rsid": "1801394",
        "clnsig": "Likely_benign",
        "condition": "Folate_metabolism_disorder",
        "gene": "MTRR",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100003,
    },
    # COMT V158M — mock MHG fixture is A/A (homozygous ALT carrier).
    {
        "chrom": "22",
        "pos_grch37": 19951271,
        "pos_grch38": 19963748,
        "ref_alt_grch37": ("G", "A"),
        "ref_alt_grch38": ("G", "A"),
        "rsid": "4680",
        "clnsig": "Drug_response",
        "condition": "Methylphenidate_response",
        "gene": "COMT",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100004,
    },
    # CYP2C9*2 — mock MHG fixture is C/T (het carrier).
    {
        "chrom": "10",
        "pos_grch37": 96702047,
        "pos_grch38": 94942290,
        "ref_alt_grch37": ("C", "T"),
        "ref_alt_grch38": ("C", "T"),
        "rsid": "1799853",
        "clnsig": "Drug_response",
        "condition": "Warfarin_response",
        "gene": "CYP2C9",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100005,
    },
    # SLCO1B1 — mock MHG fixture is T/C (het carrier).
    {
        "chrom": "12",
        "pos_grch37": 21331549,
        "pos_grch38": 21178615,
        "ref_alt_grch37": ("T", "C"),
        "ref_alt_grch38": ("T", "C"),
        "rsid": "4149056",
        "clnsig": "Drug_response",
        "condition": "Statin-induced_myopathy",
        "gene": "SLCO1B1",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100006,
    },
    # BRCA1 — mock MHG fixture is G/A (het carrier — pathogenic).
    # Position corrected against NCBI dbSNP: GRCh37 41209080, GRCh38 43057063.
    # The pre-Round 23 fixture used 41245466, which doesn't match either build.
    {
        "chrom": "17",
        "pos_grch37": 41209080,
        "pos_grch38": 43057063,
        "ref_alt_grch37": ("G", "A"),
        "ref_alt_grch38": ("G", "A"),
        "rsid": "80357906",
        "clnsig": "Pathogenic",
        "condition": "Hereditary_breast_and_ovarian_cancer_syndrome",
        "gene": "BRCA1",
        "review": "criteria_provided,_multiple_submitters,_no_conflicts",
        "allele_id": 100007,
    },
    # TP53 — mock MHG fixture is G/G (homozygous reference, must NOT trigger).
    {
        "chrom": "17",
        "pos_grch37": 7577538,
        "pos_grch38": 7674222,
        "ref_alt_grch37": ("G", "T"),
        "ref_alt_grch38": ("G", "T"),
        "rsid": "121918506",
        "clnsig": "Pathogenic",
        "condition": "Li-Fraumeni_syndrome",
        "gene": "TP53",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100008,
    },
    # rsID not in the MHG fixture — used for negative-match assertions.
    # Build-invariant synthetic; same coordinates in both files.
    {
        "chrom": "1",
        "pos_grch37": 100,
        "pos_grch38": 100,
        "ref_alt_grch37": ("A", "T"),
        "ref_alt_grch38": ("A", "T"),
        "rsid": "999999999",
        "clnsig": "Benign",
        "condition": "Synthetic_test_only",
        "gene": "TESTGENE",
        "review": "no_assertion_provided",
        "allele_id": 100009,
    },
    # CFTR ΔF508 indel — mock MHG fixture is no-call (-/-, array can't call indels).
    # Same REF/ALT in both builds; only position changes.
    {
        "chrom": "7",
        "pos_grch37": 117199644,
        "pos_grch38": 117559590,
        "ref_alt_grch37": ("CTT", "C"),
        "ref_alt_grch38": ("CTT", "C"),
        "rsid": "113993960",
        "clnsig": "Pathogenic",
        "condition": "Cystic_fibrosis",
        "gene": "CFTR",
        "review": "criteria_provided,_multiple_submitters,_no_conflicts",
        "allele_id": 100010,
    },
    # ── NIPA1 rs104894490: the strand-inverted regression case (ADR-0021) ──
    # Between GRCh37 and GRCh38 the reference strand for chr15 region inverted.
    # GRCh37 ClinVar VCF: chr15:23060816 REF=C ALT=G (Pathogenic).
    # GRCh38 ClinVar VCF: chr15:22812251 REF=G ALT=A (Pathogenic).
    # A user with G/G on GRCh38 (homozygous reference) was historically flagged
    # as a pathogenic carrier because the analyzer compared against the
    # GRCh37 entry where G IS the ALT. The carrier check at clinvar.py:179
    # was correct; the cross-build comparison wasn't. This row exists so the
    # regression test pins the dispatch contract per build.
    {
        "chrom": "15",
        "pos_grch37": 23060816,
        "pos_grch38": 22812251,
        "ref_alt_grch37": ("C", "G"),
        "ref_alt_grch38": ("G", "A"),
        "rsid": "104894490",
        "clnsig": "Pathogenic",
        "condition": "Hereditary_spastic_paraplegia_6",
        "gene": "NIPA1",
        "review": "criteria_provided,_single_submitter",
        "allele_id": 100011,
    },
]

# Multi-allelic rows test the parser's per-ALT split. Each ALT in a comma
# list pairs with its corresponding pipe-separated CLNSIG/CLNDN/ALLELEID.
MULTI_ALLELIC_RECORDS = [
    # rs1065852 (CYP2D6) — mock MHG fixture is G/A (carries the A allele).
    # Two ALTs: A → Drug_response (should match), C → Benign (should not match).
    {
        "chrom": "22",
        "pos_grch37": 42526694,
        "pos_grch38": 42130692,
        "ref": "G",
        "alts": "A,C",
        "rsid": "1065852",
        "clnsigs": "Drug_response|Benign",
        "clndns": "Codeine_response|Synthetic_benign_pair",
        "gene": "CYP2D6",
        "review": "criteria_provided,_single_submitter",
        "allele_ids": "100020|100021",
    },
]

HEADER_TEMPLATE = """##fileformat=VCFv4.1
##fileDate=20260101
##source=ClinVarSyntheticForTests
##reference={reference}
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ClinVar allele id">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease name">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene symbol:gene id">
##INFO=<ID=RS,Number=1,Type=String,Description="dbSNP ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def write_vcf(output: Path, build: str) -> None:
    """Emit a build-specific synthetic ClinVar VCF."""
    pos_key = "pos_grch37" if build == "GRCh37" else "pos_grch38"
    ref_alt_key = "ref_alt_grch37" if build == "GRCh37" else "ref_alt_grch38"

    lines = [HEADER_TEMPLATE.format(reference=build)]
    for rec in RECORDS:
        ref, alt = rec[ref_alt_key]
        info = (
            f"ALLELEID={rec['allele_id']};CLNDN={rec['condition']};"
            f"CLNSIG={rec['clnsig']};CLNREVSTAT={rec['review']};"
            f"GENEINFO={rec['gene']}:{rec['allele_id']};RS={rec['rsid']}"
        )
        lines.append(
            f"{rec['chrom']}\t{rec[pos_key]}\t{rec['allele_id']}\t{ref}\t{alt}\t.\t.\t{info}\n"
        )
    for r in MULTI_ALLELIC_RECORDS:
        info = (
            f"ALLELEID={r['allele_ids']};CLNDN={r['clndns']};CLNSIG={r['clnsigs']};"
            f"CLNREVSTAT={r['review']};GENEINFO={r['gene']}:{r['rsid']};RS={r['rsid']}"
        )
        lines.append(
            f"{r['chrom']}\t{r[pos_key]}\t{r['rsid']}\t{r['ref']}\t{r['alts']}\t.\t.\t{info}\n"
        )
    output.write_text("".join(lines), encoding="utf-8")
    total = len(RECORDS) + len(MULTI_ALLELIC_RECORDS)
    print(f"Wrote {total} {build} ClinVar VCF lines to {output}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--build",
        choices=("grch37", "grch38", "both"),
        default="both",
        help="Which build's VCF to emit. Default: both.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "fixtures",
        help="Output directory for fixture files.",
    )
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    builds = (
        ("GRCh37", "GRCh38")
        if args.build == "both"
        else ("GRCh37" if args.build == "grch37" else "GRCh38",)
    )
    for build in builds:
        suffix = build.lower()  # "grch37" / "grch38"
        write_vcf(args.output_dir / f"mock_clinvar_{suffix}.vcf", build)


if __name__ == "__main__":
    main()

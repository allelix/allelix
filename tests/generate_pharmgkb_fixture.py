# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
#!/usr/bin/env python3
"""Generate a tiny synthetic ClinPGx clinical-annotations dump for tests.

Produces two TSVs that mirror real ClinPGx structure:
  - clinical_annotations.tsv: one row per clinical annotation
  - clinical_ann_alleles.tsv: per-genotype rows linked by Clinical Annotation ID

The synthetic genotypes line up with the MyHappyGenes mock fixture's known
SNPs, so end-to-end `analyze` tests on the MHG fixture will fire on
known-carrier rows.

Usage:
    python tests/generate_pharmgkb_fixture.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Each annotation: id, rsid (or non-rsid for skip-tests), gene, drugs,
# phenotype, phenotype_category, level_of_evidence, score
ANNOTATIONS = [
    (
        "PA-001",
        "rs1801133",
        "MTHFR",
        "methotrexate",
        "Increased toxicity risk",
        "Toxicity",
        "2A",
        "0.85",
    ),
    ("PA-002", "rs4680", "COMT", "tramadol", "Reduced response", "Efficacy", "3", "0.50"),
    ("PA-003", "rs1799853", "CYP2C9", "warfarin", "Lower dose required", "Dosage", "1A", "0.95"),
    (
        "PA-004",
        "rs4149056",
        "SLCO1B1",
        "simvastatin",
        "Increased myopathy risk",
        "Toxicity",
        "1A",
        "0.95",
    ),
    # Star allele — loader must skip
    ("PA-005", "CYP2D6*1/*2", "CYP2D6", "codeine", "Normal metabolism", "Efficacy", "2A", "0.80"),
    # Multi-rsid composite — loader must skip
    ("PA-006", "rs1065852, rs16947", "CYP2D6", "tramadol", "Composite", "Efficacy", "3", "0.70"),
    # Single rsid + non-SNV per-allele rows: exercises the inner
    # `_normalize_genotype is None: continue` skip in iter_pharmgkb_records.
    # The annotation row passes _is_single_rsid, but every allele row below
    # is rejected by _normalize_genotype, so no records are yielded.
    (
        "PA-007",
        "rs99999999",
        "TESTGENE",
        "testdrug",
        "Synthetic — for inner-skip coverage",
        "Efficacy",
        "4",
        "0.10",
    ),
    # ADR-0020 synthetic non-finding: the per-allele function lookup
    # provided to the loader (see conftest.mock_cpic_lookup) classifies
    # rs900000010 G as Normal function and A as Decreased. The GG row
    # below must store is_nonfinding=1 in the cache; the AG row must
    # store is_nonfinding=0.
    (
        "PA-008",
        "rs900000010",
        "CFTR",
        "ivacaftor",
        "Synthetic non-finding — for ADR-0020 coverage",
        "Efficacy",
        "1A",
        "0.95",
    ),
]

# Per-genotype annotations: clinical_annotation_id, genotype, annotation_text,
# allele_function. Genotypes are aligned with MHG mock fixture diploid calls
# so the end-to-end analyze tests trigger expected matches.
# v0.6.1 fixture correction (ADR-0017): real ClinPGx populates `Allele
# Function` ONLY on haplotype rows (*1, *2, *15:02). SNV rows have it empty
# on every row. Earlier fixture revisions inverted this and concealed the
# v0.6.0 regression — see CHANGELOG. The fixture now models real shape:
#   - SNV genotype rows (AA, AG, etc.): Allele Function = ""
#   - Haplotype rows (*1/*2 etc.): Allele Function populated (but those
#     rows are dropped by `_normalize_genotype`, so they never reach the
#     classifier — they exist here only to model ClinPGx's actual dump).
ALLELES = [
    # rs1801133: MHG fixture has G/A → normalized "AG"
    ("PA-001", "AG", "Heterozygous carriers may have elevated toxicity risk", ""),
    ("PA-001", "AA", "Homozygous variant carriers have high toxicity risk", ""),
    ("PA-001", "GG", "Normal toxicity profile", ""),
    # rs4680: MHG fixture has A/A
    ("PA-002", "AA", "Reduced opioid response possible", ""),
    ("PA-002", "GG", "Normal opioid response", ""),
    # rs1799853: MHG fixture has C/T → normalized "CT"
    ("PA-003", "CT", "Consider lower warfarin starting dose", ""),
    ("PA-003", "TT", "Significantly lower starting dose recommended", ""),
    # rs4149056: MHG fixture has T/C → normalized "CT"
    ("PA-004", "CT", "Increased myopathy risk on simvastatin", ""),
    ("PA-004", "CC", "Normal statin tolerance", ""),
    # PA-005 (star allele) — real ClinPGx populates Allele Function for
    # haplotype rows. Loader rejects via _normalize_genotype.
    ("PA-005", "*1/*2", "Normal metabolizer", "Normal function"),
    # PA-006 (multi-rsid composite) — rejected at annotation level.
    ("PA-006", "AG", "Composite haplotype", ""),
    # PA-007 — non-SNV per-allele rows: rejected by _normalize_genotype.
    ("PA-007", "*1/*1", "Wild type haplotype", "Normal function"),
    ("PA-007", "ins", "Insertion variant", ""),
    # PA-008: ADR-0020 non-finding via structured lookup. The mock CPIC
    # lookup classifies rs900000010 G as Normal, A as Decreased. GG is
    # therefore a reference-homozygote → non-finding. AG is a carrier →
    # finding. The annotator's filter is a pure join — no description text
    # is consulted to make the decision.
    ("PA-008", "GG", "Reference homozygote row (non-finding via CPIC lookup)", ""),
    ("PA-008", "AG", "Heterozygous carrier row (finding via CPIC lookup)", ""),
    ("PA-008", "AA", "Homozygous variant row (finding via CPIC lookup)", ""),
]

CLINICAL_ANN_HEADER = (
    "Clinical Annotation ID\tVariant/Haplotypes\tGene\tLevel of Evidence\t"
    "Score\tPhenotype Category\tDrug(s)\tPhenotype(s)\n"
)
CLINICAL_ANN_ALLELES_HEADER = (
    "Clinical Annotation ID\tGenotype/Allele\tAnnotation Text\tAllele Function\n"
)


def write_tsvs(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ann_path = out_dir / "clinical_annotations.tsv"
    alleles_path = out_dir / "clinical_ann_alleles.tsv"

    with ann_path.open("w", encoding="utf-8") as fh:
        fh.write(CLINICAL_ANN_HEADER)
        for ann_id, variant, gene, drugs, phenotype, phen_cat, loe, score in ANNOTATIONS:
            fh.write(
                f"{ann_id}\t{variant}\t{gene}\t{loe}\t{score}\t{phen_cat}\t{drugs}\t{phenotype}\n"
            )

    with alleles_path.open("w", encoding="utf-8") as fh:
        fh.write(CLINICAL_ANN_ALLELES_HEADER)
        for ann_id, geno, text, func in ALLELES:
            fh.write(f"{ann_id}\t{geno}\t{text}\t{func}\n")

    print(f"Wrote {len(ANNOTATIONS)} annotations + {len(ALLELES)} allele rows to {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "fixtures" / "mock_pharmgkb",
    )
    args = ap.parse_args()
    write_tsvs(args.output_dir)


if __name__ == "__main__":
    main()

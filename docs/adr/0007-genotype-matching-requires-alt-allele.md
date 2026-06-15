# ADR-0007: Genotype matching requires the user to carry the ALT allele

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Reference databases like ClinVar publish (rsid, REF, ALT, significance) tuples — "at this position, when the ALT allele is present, here's what's known." A naive annotator that triggers on any rsID match would flag everyone as a carrier of every annotated variant in the database, including users who carry only the reference (normal) allele.

This is the difference between "you have a position where someone else's pathogenic variant lives" (true for everyone, meaningless) and "you carry the pathogenic allele at this position" (true for some users, actionable).

A variant being in ClinVar as pathogenic means nothing if the person has the reference (normal) allele.

## Decision

An annotation triggers if and only if at least one of the user's two alleles equals the database entry's ALT allele:

```python
if variant.allele1 == clinvar_row.alt or variant.allele2 == clinvar_row.alt:
    yield Annotation(...)
```

Specifically:

- **Homozygous reference** (e.g., user G/G when REF=G ALT=A) → no annotation.
- **Heterozygous carrier** (user G/A when REF=G ALT=A) → annotation.
- **Homozygous variant** (user A/A when REF=G ALT=A) → annotation.
- **No-call** (user `-/-`) → no annotation. We don't speculate from missing data.
- **rsID not in database** → no annotation, regardless of genotype.

Multiple ClinVar rows can share an rsID (multi-allelic sites). Each row is evaluated independently against the user's two alleles.

Strand orientation is assumed consistent (forward) for v0.2.0. Strand-flipping for ambiguous SNPs is deferred; see ADR-0010.

## Consequences

- Allelix only flags variants the user actually carries — the report length is realistic, not pathological.
- Tests pin both directions: known carriers must trigger, known reference-homozygotes must not. `tests/annotators/test_clinvar.py::TestGenotypeMatching` is the regression suite.
- Future annotators (PharmGKB, GWAS, SNPedia) inherit the same rule. The Annotator ABC docstring states it; per-annotator implementation enforces it.
- The `Annotation.genotype_match` field records which ALT allele triggered the annotation, so reports can show the user *why* a row appeared.
- Edge case: A/T and G/C SNPs are strand-ambiguous (the complement is the same letter pair). We accept the small risk of false matches and document the strand-flipping work in ADR-0010.

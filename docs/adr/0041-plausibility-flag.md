# ADR-0041: Plausibility flag from zygosity and population frequency

**Status:** Accepted
**Date:** 2026-10-10
**Builds on:** ADR-0039 (every annotation names its allele and carries the
user's alleles in the same orientation), ADR-0040 (gnomAD frequency per
allele by position, on GRCh37 as well as GRCh38).

## Context

Genotyping arrays sometimes miscall very rare variants. The typical
pattern is a homozygous call for an allele that population data say almost
nobody carries even once. Reported without comment, such a call reads as a
homozygous pathogenic finding. Allelix already has both inputs needed to
spot it: the call's copy count of the allele, and that allele's gnomAD
frequency.

## Decision

Each annotation gets a `plausibility` value, computed after gnomAD
enrichment from the copies of the annotation's allele that the user
carries and that allele's gnomAD frequency (`allele_frequency`):

| Copies | Frequency | Value |
|---|---|---|
| 2 | below `homozygous_max_af` (default 0.0001) | `implausible` |
| 1 | at or below `heterozygous_max_af` (default 0) | `suspicious` |
| 1 or 2 | otherwise | `plausible` |
| any | no gnomAD value | empty (not assessed) |

- **The allele is the one the frequency describes.** For ClinVar, ClinPGx
  and SNPedia that is the annotation's `alt`. For a GWAS row without a
  reported risk allele, it is the single non-reference allele the user
  carries, the same choice enrichment makes. The copy count comes from the
  user's alleles in the annotation's orientation, so strand-flipped
  matches count correctly. A separate risk-allele field is not needed.
- **A missing frequency is not zero.** The gnomAD cache holds exome sites
  that carry an rsID. Most non-coding variants are missing from it for
  that reason alone, so absence is never read as "no carriers". The
  heterozygous rule fires only when gnomAD reports the allele with a
  frequency of 0 (or at or below a raised threshold).
- **Chromosomes.** The homozygous rule applies to autosomes only: on X a
  male's hemizygous call looks homozygous, and Allelix does not infer sex.
  The heterozygous rule applies to autosomes and X. Y and MT are not
  assessed.
- **A flag, never a filter.** Flagged rows stay in every report. The HTML
  report shows a warning icon next to the genotype, with an explanation on
  hover and in the detail panel. The terminal table marks the genotype and
  prints a legend. JSON carries the value on every annotation.
- **Tunable.** Both thresholds live in `[plausibility]` in `config.toml`
  and can be set with `allelix config set plausibility.<key> <value>`.
  Values outside [0, 1] are rejected.

The new JSON field bumps the report schema to 7 (ADR-0033). The annotated
VCF output is unchanged (ADR-0036).

## Consequences

- A very rare homozygous call on an autosome is visibly marked for review,
  for both GRCh37 chips and GRCh38 input, wherever gnomAD has the allele.
- Variants without a gnomAD value (palindromic GRCh37 sites, non-exome
  sites, sources without an allele) are never flagged. The flag can miss
  artifacts; it is not designed to clear calls.
- A homozygous rare X-linked call is never flagged, even in a female.
  Inferring sex from the X calls would close that gap and is a separate
  decision.

*Amended by ADR-0050:* "not assessed" is the value `not_assessed` with a reason (`plausibility_reason`), every annotation gets a value, a call the file declares haploid counts as one copy, and the report warns on three or more rare homozygous ClinVar Pathogenic / Likely pathogenic calls.

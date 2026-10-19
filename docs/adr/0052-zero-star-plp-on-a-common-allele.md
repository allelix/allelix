# ADR-0052: A zero-star ClinVar P/LP submission on a common allele leaves the default report

**Status:** Accepted
**Date:** 2026-10-17
**Amends:** ADR-0045 (ClinVar rows scored against ClinVar's own
classification). **Builds on:** ADR-0040 (gnomAD keys and the identity
table), ADR-0046 (ClinVar rows without an rsID).

## Context

ADR-0045 caps a submission by ClinVar's aggregate classification. It does
not help when the aggregate itself is Pathogenic or Likely pathogenic but
rests on a single submission without assertion criteria, and the allele is
common.

CDKN2B-AS1 rs1063192 is the case that showed it. ClinVar's aggregate is
"Likely pathogenic; protective", built from two zero-star submissions:
- "Likely pathogenic" for a breast tumour;
- "protective" for coronary disease.

The Likely pathogenic allele has a gnomAD frequency of about 0.65, and it
is the major allele in most populations. The HG002 benchmark sample is
homozygous for it, and Allelix 2.4.2 and 2.5 both show it as a homozygous
Likely pathogenic finding at magnitude 7.0 in the default report. An allele
that common cannot cause a rare disease. ACMG/AMP's BA1 criterion treats a
population frequency of 5% or more as stand-alone evidence of benign
impact.

## Decision

A ClinVar row is scored 4.0, below the default display floor of 5.0, when
all of these hold:
1. its classification is Pathogenic or Likely pathogenic (low penetrance
   and the combined Pathogenic/Likely pathogenic included);
2. its own submission has zero stars: review status "no assertion criteria
   provided" or "no classification provided";
3. gnomAD's frequency for ClinVar's allele is 0.05 or more (the cache
   holds the gnomAD v4.1 exome overall AF);
4. its magnitude is above 4.0, so a row already capped by ADR-0045 is
   left as it is.

### The frequency is the ClinVar allele's, not the call's

The frequency enrichment gives a call (ADR-0040) does not serve here.
- **On GRCh37 input it abstains at palindromic sites,** to protect the
  orientation of the user's call.
- **rs1063192 is such a site:** dbSNP lists G>A, G>C and G>T there. So
  HG002 GRCh37 and every GRCh37 array get no frequency for it.

This rule describes ClinVar's variant, not the call, so it looks up
ClinVar's allele directly:
- **The locus:** the rsID's GRCh38 locus from the identity table.
- **The key:** ClinVar's own forward-strand REF and ALT at that locus, in
  the gnomAD cache.
- **No palindromic abstention:** no call is translated.

**On any doubt the row is unchanged.** A doubt is any of these:
- the row has no rsID;
- the table places the rsID at more than one locus on either build;
- either ClinVar build has no readable cache;
- either build lacks a record with the same REF and ALT at the table's
  locus;
- the table has no row with that REF and ALT on both builds;
- the row is not at the table's locus on the analysis build;
- the build's reference carries the variant allele (REF equals ALT);
- gnomAD has no row for the key.

Every row this lookup reaches carries the frequency in `clinvar_allele_af`,
whether or not it is scored down.

### What is unchanged

- A row without that frequency keeps its magnitude. This covers no gnomAD,
  no identity table, one ClinVar build only, and every doubt above.
- A submission with one star or more keeps its magnitude. HFE C282Y
  (rs1800562, gnomAD AF 0.059) and HbS (rs334) are such variants.
- Rows of other kinds (drug response, protective, risk factor) and other
  sources keep their magnitude.

The row says why, and the run says how many:
- JSON: `common_allele_capped: true` and `clinvar_allele_af` on the row,
  and `clinvar.common_allele_rows_scored_down` for the run.
- Terminal: `(1 of 2; zero-star; common in gnomAD)` after the significance.
- HTML: a "Scored" line in the detail panel, a neutral repute, and a line in
  the ClinVar notice.
- Console, HTML and JSON each give the count: "N ClinVar Pathogenic or
  Likely pathogenic rows not shown: zero-star classification on a common
  allele". The rows stay in the report at `--min-magnitude 4` or lower.

The rule runs after gnomAD enrichment, in the pipeline, because the
frequency is not known when the ClinVar annotator scores a row. The ClinVar
cache does not change, so the interpreter version does not change.

## Consequences

- HG002 no longer opens with a Likely pathogenic finding on an allele most
  people carry.
- The rule depends on gnomAD, the identity table and both ClinVar builds.
  Without any of them, the report is as before. The
  `multiple_homozygous_plp` warning (ADR-0050) already ignores common
  alleles, so it is unaffected.
- A zero-star P/LP submission on a common allele that is a real low-risk
  finding is still in the report below the default floor, with its reason.

# ADR-0045: ClinVar rows are scored against ClinVar's own classification of the variant

**Status:** Accepted; amended by ADR-0052
**Date:** 2026-10-13
**Amends:** ADR-0008 (ClinVar magnitude from the classification). A row's
magnitude is still its submission's classification, now capped by
ClinVar's aggregate classification of the variant.

## Context

Since 2.2 the ClinVar cache holds one row per submission (SCV) and each row
is scored alone. ClinVar also publishes an aggregate classification per
variant (variant_summary `ClinicalSignificance`), built from the
submissions in the strongest review tier present. When those disagree it
says "Conflicting classifications of pathogenicity".

Scoring submissions alone let a single Pathogenic submission headline a
variant that ClinVar itself calls Conflicting. In the 2026-10-10 data,
28,678 GRCh37 variants with an rsID have a contributing P/LP submission
under a Conflicting aggregate. On three real 23andMe files, 10 or 11 of the
13 P/LP variants in the 2.4.1 default report were of this kind (TYR
rs1042602, BBS2 rs4784677 "1 of 10 submitters", SLC3A1 rs698761), so a
healthy person's report opened with about ten Pathogenic 9.0 rows.

The loader keeps only submissions ClinVar marks as contributing to the
aggregate (`ContributesToAggregateClassification = yes`). submission_summary
holds current submissions only; "no" means the submission was outranked by
a higher review tier. That is right when the higher tier makes a disease
call: an expert panel's VUS or Benign should win over single submitters.
It is wrong when the higher tier classifies something else. For Factor V
Leiden (rs6025, VariationID 642) the only contributing submission is an
expert panel's "drug response" for hormonal contraceptives; its 25 P/LP
submissions (thrombophilia, factor V deficiency) were dropped,
so a carrier saw only "drug response". 46 GRCh37 variants with an rsID lose
every P/LP submission this way, among them CFTR, DPYD *2A, RYR1 and CACNA1S
variants. The 2.4.x count was 53, because 2.4.x still loaded submissions
ClinVar has flagged; 2.5 leaves flagged submissions out first.

## Decision

### The cache records ClinVar's view

Each row carries the variant's aggregate classification, verbatim, and the
submitter split: how many of the submissions counted for the aggregate
gave this row's classification, and how many were counted. Flagged
submissions and placeholder values ("not specified", "-") are not counted
and not loaded.

### Capping

A row whose own classification is a disease classification (P, LP and
their low-penetrance forms, VUS and its subdivisions, LB, B) is scored at
its own magnitude, capped by the aggregate's germline term (the text
before the first ";", so "Conflicting classifications of pathogenicity;
risk factor" counts as Conflicting):

| Aggregate | Cap |
|---|---|
| Conflicting classifications of pathogenicity | 4.0 |
| Likely benign | 2.0 |
| Benign/Likely benign | 1.5 |
| Benign | 1.0 |

Every other aggregate leaves the row's magnitude unchanged. Rows of
another kind (drug response, risk factor, protective, Affects, the risk
allele terms) are never capped: under "Benign; risk factor" the germline
term judges the disease calls, not the risk factor call. No threshold
of Allelix's own applies: HFE rs1799945 (H63D, 43 of 45 submitters
Pathogenic or Likely pathogenic) is capped because ClinVar calls it
Conflicting. The split counts each classification on its own, so its
Pathogenic row shows 40 of 45.

### Outranked submissions where the aggregate is not a disease classification

When the aggregate's germline term is not a disease classification (drug
response, risk factor, protective, association, other, Affects, the risk
allele terms, not provided, or none), the outranked P/LP submissions are
kept, judged the way ClinVar judges an aggregate:

1. Take the outranked submissions that are not flagged and give a disease
   classification (P, LP and their low-penetrance forms, VUS and its
   subdivisions, LB, B).
2. Keep only the strongest review tier among them: practice guideline, then
   expert panel, then criteria provided, then no assertion criteria.
3. If that tier holds no P/LP submission, nothing is kept. If every
   submission in it is P/LP, its P/LP submissions are scored as submitted.
   If it also holds VUS, LB or B, its P/LP submissions are kept but capped
   at 4.0, as Conflicting.

The rows are marked as outside the aggregate, and their split counts that
tier. Factor V Leiden's criteria-provided tier is all P/LP (22 Pathogenic,
1 Pathogenic low penetrance), so those rows are shown as submitted beside
the drug response row.
MTHFR rs1801133 (3 of 10 in its top tier) is capped. VKORC1 rs9923231,
whose only P/LP is a no-assertion submission below criteria-provided Benign
ones, gets no row. When the aggregate is a disease classification,
outranked submissions stay out, as before.

### Display and disclosure

- JSON (schema 9): every annotation carries `clinvar_aggregate`,
  `submitters_same`, `submitters_total`, `aggregate_capped` and
  `outside_aggregate`. A top-level `clinvar` block counts the P/LP rows
  scored down and the rows hidden for an unscored classification.
- Terminal: a scored-down or outranked row's significance reads
  `pathogenic (1 of 10; ClinVar: conflicting)` or
  `pathogenic (22 of 23; outranked; ClinVar: drug response)`.
- HTML: the detail panel shows ClinVar's classification and the split, and
  a capped row's repute is neutral.
- Console, HTML and JSON each give the number of P/LP rows scored down,
  since those rows leave the default report.

### Vocabulary

The per-submission terms beside the five-tier scale get magnitudes:
Pathogenic, low penetrance 7.5; Likely pathogenic, low penetrance 6.5;
Established risk allele and confers sensitivity 6.0; Likely risk allele,
Affects and protective 5.0; VUS-high 4.5; Uncertain risk allele and VUS-mid
4.0; VUS-low 3.5. A classification with no magnitude is no longer scored
5.0: the row is hidden and the run counts it. A test fails when the
vocabulary snapshot holds a term with no decision, and the live drift
sentinel fails on a new term or on a Conflicting or Benign aggregate the
cap does not know.

## Consequences

- Most Pathogenic rows in a typical array report leave the default view;
  they remain at `--min-magnitude 4` or lower, with the split, and the
  report says how many.
- A Factor V Leiden carrier sees the thrombophilia submissions again.
- Allelix follows ClinVar's aggregate even where a large majority of
  submitters say Pathogenic; the split makes that visible.
- The cache gains columns, so the ClinVar interpreter version is bumped
  (6 to 7). A cache built before cannot be read with this rule and is
  rebuilt by the next `db update`.

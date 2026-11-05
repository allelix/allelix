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
- **Amended 2026-10-25 (2.5.2): the scored-down variants are named.** Each
  one is listed with its gene, rsID, submitter split and ClinVar's own
  classification ("GALT rs75391579: 40 of 41 submitters Pathogenic;
  ClinVar's classification of the variant: Conflicting classifications of
  pathogenicity"). The list is on the console, in the HTML ClinVar notice
  and in JSON `clinvar.scored_down_variants` (schema 11). It is taken from
  the unfiltered findings, so it is the same at the default floor, where
  the rows themselves are hidden, and at `--min-magnitude 0`. This is
  reporting only: the cap and the magnitudes are unchanged.

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

## Amendment 2.6.0 (2026-10-30): outranked P/LP under "Likely risk allele" or "other"

Supervisor ruling, 2026-10-30.

**The gap.** An outranked P/LP submission, from a lower review tier whose own submissions all agree, was stored as agreed and scored in full (9.0 Pathogenic, 7.0 Likely pathogenic) when ClinVar's aggregate was "Likely risk allele" or "other". Neither term is in the cap table. The top tier did not call the variant P/LP, so a lower tier's P/LP is a conflict with ClinVar's view, as under Conflicting.

**The rule.** Under those two aggregate terms, compound forms included (for example "Likely risk allele; drug response"), an outranked P/LP submission takes the Conflicting cap (4.0).
- The row stays, with its significance, aggregate, outside-aggregate marker and split.
- Like other scored-down rows, it is named in `clinvar.scored_down_variants`.
- Contributing submissions, non-disease submissions and every other aggregate are unchanged.

**No rebuild.** The rule is applied when the cache is read. The stored cache does not change, so the interpreter version and existing caches stay as they are.

**Not decided here.** A display tier for lopsided Conflicting splits was left to the owner; the next amendment records the decision (no tier).

## Amendment 2.6.0 (2026-10-30): no tier by submitter count

Owner decision.

**The question.** Some Conflicting variants have a lopsided split, where nearly every submitter says Pathogenic or Likely pathogenic. The question was whether Allelix should score those above the Conflicting cap, or show them in a tier of their own, by counting the submitters.

**Decision.** No. ClinVar's aggregate classification stays the rule. A variant ClinVar calls Conflicting is capped at 4.0 whatever its split, and the split is shown beside it.

**Why.**
- There is no accepted standard for weighing submitter votes. Neither ClinVar nor ClinGen publishes one, so any threshold would be Allelix's own.
- The dissent can be the newer expert view. One recent submission with stronger evidence can be right against many older ones.
- A vote count ignores who submitted and when. Review status, the evidence each submitter used and the submission date all count for more than the number of votes.
- Allelix reports what ClinVar says; it does not classify variants (Regulatory Posture in the README). Overriding ClinVar's aggregate with a vote count would be a classification of Allelix's own.

**The facts.** In the GRCh38 ClinVar data used for 2.6:
- 2,161 Conflicting variants have at least 5 submitters with 80% or more of them P/LP.
- 258 have at least 10 submitters with 90% or more.
- Examples: GALT rs75391579 (40 of 41 submitters Pathogenic) and BTD rs13078881 (51 of 52).

These are capped like every other Conflicting variant. The report names each one with its split (`clinvar.scored_down_variants`, 2.5.2 amendment above), so a reader can see how lopsided it is and look it up in ClinVar.

**Revisit when** ClinVar or ClinGen publishes a weighting standard for submissions, or an expert panel resolves a given variant. A resolved variant needs no change in Allelix: ClinVar's aggregate then follows the panel, and the next `db update` picks it up.

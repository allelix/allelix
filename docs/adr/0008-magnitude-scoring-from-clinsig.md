# ADR-0008: Magnitude scoring is Allelix-derived from source classifications

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

`Annotation.magnitude` is a 0-10 importance score, originally inspired by SNPedia (where contributors assign magnitudes directly). For SNPedia, magnitude is source-attributed — it came from SNPedia.

For ClinVar, there is no native magnitude. ClinVar publishes categorical labels (Pathogenic, Likely_pathogenic, Uncertain_significance, Benign, etc.). To sort and filter annotations across sources, Allelix needs a numeric score. That means *Allelix* is mapping ClinVar's categories to numbers.

ADR-0003 says Allelix never asserts. This ADR is the carve-out: synthesizing a magnitude from a source category is an Allelix decision, and it must be visible as such.

## Decision

ClinVar's `CLNSIG` is mapped to a magnitude via a fixed lookup table in `allelix.annotators.clinvar._CLNSIG_MAGNITUDE`. Default for unrecognized values is 5.0 (mid-scale, neutral).

The mapping (subject to refinement; pin via tests):

| CLNSIG | Magnitude |
|---|---|
| Pathogenic | 9.0 |
| Pathogenic/Likely_pathogenic | 8.5 |
| Likely_pathogenic | 7.0 |
| Drug_response | 6.5 |
| Risk_factor | 6.0 |
| Uncertain_significance | 4.0 |
| Conflicting_interpretations_of_pathogenicity | 4.0 |
| Likely_benign | 2.0 |
| Benign/Likely_benign | 1.5 |
| Benign | 1.0 |

The score is presented in reports alongside the source attribution, so a row reads: "ClinVar / clinvar_pathogenic / magnitude 9.0." The user sees that ClinVar said "Pathogenic" and that Allelix's mapping turned it into 9.0.

Documentation in the README and `Annotator.attribution` field make clear that magnitude in non-SNPedia rows is a derived sort key, not an external claim.

## Benign suppression (amended 2026-05-19)

ClinVar rows with CLNSIG in `{Benign, Likely_benign, Benign/Likely_benign}` are
suppressed by default in `ClinVarAnnotator.annotate()`. The `--include-benign`
CLI flag restores them.

**Rationale.** On real data, benign/likely_benign rows are the majority of
ClinVar matches. They carry magnitudes 1.0–2.0 and provide no actionable
information in a report. Emitting them inflates annotation counts and buries
clinically relevant findings. The suppression is at the annotator level (not the
renderer) so downstream count snapshots, JSON output, and HTML reports all
reflect only non-benign annotations by default.

**`--include-benign` semantics.** Passing `--include-benign` on `analyze`,
`methylation`, or `pharmacogenomics` sets `ClinVarAnnotator(include_benign=True)`.
Benign rows then emit at their natural magnitude and are subject to the same
`--min-magnitude` filter as everything else. This is the "full dump" mode for
researchers who want to see everything ClinVar says about their variants.

## Consequences

- A single magnitude scale lets reports rank annotations across heterogeneous sources.
- Tests pin the table values (e.g., Pathogenic → 9.0) so a refactor that breaks the mapping is visible at PR time.
- SNPedia magnitudes are passed through as-is from the source wiki. They are editorial scores assigned by wiki contributors, not subject to Allelix's magnitude tier methodology. Values like 8.8 are valid SNPedia scores that do not correspond to Allelix's tier boundaries. The `source` field on each annotation distinguishes SNPedia-originated magnitudes from Allelix-derived ones.
- When PharmGKB and GWAS Catalog land, this ADR is the precedent: each source-to-magnitude mapping is documented, table-driven, and visible to the user.
- **SNPedia/ClinVar overlap:** SNPedia/ClinVar rsID overlap on real-data input: ~11% (18/157 in a reference run). Overlap is complementary, not redundant — ClinVar assigns a generic classification (risk factor / drug response / association); SNPedia provides the explicit phenotype direction and effect size at the user's actual genotype. Examples: rs602662 ClinVar "confers sensitivity" → SNPedia "Higher vitamin B12 levels"; rs9264942 ClinVar "risk factor" → SNPedia "90% reduction in HIV viral load". The two sources are kept side-by-side in default output for this reason.
- If the ClinVar table needs revision (e.g., regulators publish guidance on weighting), update this ADR and the table together — don't mutate one without the other.
- Benign suppression reduces default ClinVar output to actionable findings only; `--include-benign` is the escape hatch.

# ADR-0024: GWAS Catalog magnitude scoring from p-value + effect size

- **Date:** 2026-05-19
- **Status:** Accepted
- **Precedent:** ADR-0008 (each source-to-magnitude mapping documented here)

## Table of contents

1. [Context](#context) — why GWAS needs its own scoring
2. [Decision](#decision) — p-value tiers, effect size modifier, unknown risk allele handling
3. [Per-source magnitude floor](#per-source-magnitude-floor-amended-2026-05-19) — `--gwas-min-magnitude` default 9.0
4. [Trait-category filtering](#trait-category-filtering-amended-2026-05-19) — excluded measurement/behavioral categories
5. [Structural noise detection](#structural-noise-detection--step-15-amended-2026-05-20) — UKB body fields, NMR ratios, uncharacterized analytes
6. [Must-include rsID allowlist](#must-include-rsid-allowlist-amended-2026-05-20) — clinically significant associations that bypass the magnitude floor
7. [Consequences](#consequences)
8. [Cache invalidation on categorizer change](#cache-invalidation-on-categorizer-change-amended-2026-05-20) — `|cv:N` stamp + auto-reingest
9. [MTAG and PheCode rollup](#mtag-and-phecode-rollup-amended-2026-05-21) — duplicate collapsing

## Context

The GWAS Catalog publishes trait–SNP associations with p-values and optional
effect sizes (odds ratios or beta coefficients). Like ClinVar's CLNSIG→magnitude
mapping (ADR-0008), Allelix derives a numeric magnitude from these structured
fields to enable cross-source sorting and filtering.

Two design questions arise:

1. How to map continuous p-values to the discrete 0–10 magnitude scale.
2. What to do when the GWAS Catalog row doesn't specify a risk allele.

## Decision

### P-value magnitude tiers

The mapping uses genome-wide significance conventions from statistical genetics:

| P-value threshold | Base magnitude | Rationale |
|---|---|---|
| p < 5×10⁻¹⁰⁰ | 8.0 | Hyper-significant (huge meta-analyses) |
| p < 5×10⁻²⁰ | 7.0 | Strong GWAS signal, well-replicated |
| p < 5×10⁻⁸ | 6.0 | Genome-wide significant (standard threshold) |
| p < 5×10⁻⁶ | 4.0 | Suggestive (sub-threshold) |
| p < 5×10⁻⁴ | 3.0 | Nominal |
| p ≥ 5×10⁻⁴ or None | 2.0 | Weak or unparseable |

The 5×10⁻⁸ threshold is the Bonferroni-corrected boundary for ~1M independent
tests, universally accepted as "genome-wide significant" since Hoggart et al.
2008. Thresholds above and below are half-decade steps matching how the field
informally grades evidence strength.

### Effect size modifier

When an odds ratio (OR) or beta is available, large effect sizes bump magnitude:

| OR range | Modifier |
|---|---|
| OR ≥ 3.0 or OR ≤ 0.33 | +1.0 (large effect) |
| OR ≥ 2.0 or OR ≤ 0.5 | +0.5 (moderate effect) |
| Otherwise | +0.0 |

Final magnitude is capped at 9.0. The protective direction (OR < 1) is treated
symmetrically — a protective variant with OR = 0.2 is as noteworthy as a risk
variant with OR = 5.0.

### Comparison with other sources

| Source | Max magnitude | Scoring basis |
|---|---|---|
| ClinVar | 9.0 | Categorical CLNSIG (ADR-0008) |
| PharmGKB | 8.0 | Level of Evidence (LoE 1A → 8.0) |
| GWAS Catalog | 9.0 | Continuous p-value + OR (this ADR) |

ClinVar Pathogenic and GWAS hyper-significant + large-effect associations both
reach 9.0. This is intentional: a GWAS hit with p < 10⁻¹⁰⁰ and OR > 3 is
genuinely high-impact evidence.

### Unknown risk allele handling

GWAS Catalog rows where STRONGEST SNP-RISK ALLELE is `rs123-?` or absent don't
specify which allele is the risk allele. Without this information, the carrier
rule (ADR-0007) cannot be applied — we don't know whether the user carries the
risk allele or the protective allele.

This is the same false-positive pattern that ClinVar v0.4.x exhibited:
triggering on rsID match alone regardless of genotype.

**Decision:** Emit the annotation on rsID match alone, but cap magnitude at 3.0
regardless of p-value or effect size. This ensures unknown-risk-allele
associations:

- Still appear in comprehensive reports (--min-magnitude 0).
- Don't appear in standard reports using typical thresholds (--min-magnitude 5).
- Are clearly marked in the description: "(risk allele not specified in study)".

The cap value 3.0 is below the genome-wide significance tier (6.0) and below
typical user-facing filtering thresholds (5.0). The constant
`_UNKNOWN_RISK_ALLELE_MAG_CAP` is pinned by tests.

**Alternative rejected:** Penalty-based approach (subtract 1.0 from base
magnitude). Rejected because a genome-wide significant hit (base 6.0) would
still score 5.0, passing standard filters even when we can't verify the user
carries the risk allele.

### Strand normalization

The GWAS Catalog FAQ states that risk alleles are reported on the forward
reference strand. The MHG parser also reports forward strand. The carrier rule
therefore uses direct base comparison without strand flipping.

This assumption is not validated against pre-2014 GWAS Catalog entries, which
may have reported risk alleles on the gene-coding strand. When R-1 (strand-aware
genotype comparison) ships, the GWAS annotator should wire in the same
strand-flip fallback used by ClinVar (see `allelix/utils/allele.py`).

## Per-source magnitude floor (amended 2026-05-19)

On real data, the GWAS Catalog produces ~88,000 trait–SNP associations at
genome-wide significance or better. Even with a default `--min-magnitude 5.0`,
thousands of GWAS annotations pass the threshold (base magnitude 6.0 for the
standard genome-wide tier). This overwhelms the report with low-evidence
trait associations and buries ClinVar/PharmGKB findings.

**Decision.** A per-source magnitude floor of 9.0 is applied to GWAS Catalog
annotations via `--gwas-min-magnitude` (default 9.0). The floor is implemented
in `AnalysisResult.filter()` as a `source_min_magnitudes` dict that renderers
pass through. GWAS annotations must meet `max(min_magnitude, gwas_min_magnitude)`
to appear in the report.

The initial floor of 7.0 (v0.7.x) was insufficient: on real data, 30,000+
GWAS rows have mag 7 (genome-wide significant common-trait territory: Height,
BMI, hematological counts). The 9.0 floor restricts output to hyper-significant
associations (p < 5x10^-100) or strong signals with large effect sizes
(p < 5x10^-20 + OR >= 3.0). Verified against the full GWAS Catalog (795k
records): mock MHG produces 331 raw GWAS annotations, 7 pass at floor 9.0.

**Focused reports (methylation, pharmacogenomics)** exclude GWAS entirely by
default. Methylation biology is interpreted from ClinVar + PharmGKB, not GWAS
trait associations. The `--include-gwas` flag opts in when completeness is
desired.

**`--gwas-min-magnitude 0` semantics.** Setting the floor to 0 disables the
per-source filter and subjects GWAS annotations to only the global
`--min-magnitude` threshold. Combined with `--min-magnitude 0`, this produces
the full GWAS dump.

## Trait-category filtering (amended 2026-05-19)

The per-source magnitude floor (9.0) reduced GWAS output from thousands to
hundreds, but common-trait polygenic noise still dominated: body height (159
rows at mag 9), cholesterol panel (~100 rows), hematological counts, etc.
These are statistically robust GWAS associations with large meta-analyses,
but they are not clinically actionable for individual report consumers.

**Decision.** Each GWAS Catalog row is classified into a trait category at
load time using the `MAPPED_TRAIT` text (EFO ontology label) and
`MAPPED_TRAIT_URI` (ontology URI). The classification is stored in the
`trait_category` column of `gwas_associations`.

### Category taxonomy

| Category | Examples | Default |
|---|---|---|
| `disease` | MONDO_ URI entries, generic disease/disorder/syndrome | Included |
| `cancer` | carcinoma, neoplasm, melanoma, lymphoma | Included |
| `drug_response` | warfarin dose, response to, drug-induced | Included |
| `immune` | rheumatoid arthritis, lupus, psoriasis, asthma | Included |
| `cardiovascular` | coronary artery disease, atrial fibrillation, stroke | Included |
| `metabolic` | type 2 diabetes, obesity, insulin resistance | Included |
| `neurological` | Alzheimer's, Parkinson's, schizophrenia, epilepsy | Included |
| `other` | unclassified traits | Included |
| `body_measurement` | body height, BMI, waist circumference, lean mass | **Excluded** |
| `lipid_measurement` | cholesterol, triglyceride, VLDL, phospholipids | **Excluded** |
| `hematological_measurement` | platelet count, red cell, hemoglobin | **Excluded** |
| `other_measurement` | blood pressure, QT interval, hair color, FEV/FVC | **Excluded** |
| `behavioral` | educational attainment, smoking, alcohol consumption | **Excluded** |

### Classification priority

Disease categories are checked before measurement categories so that
multi-trait entries (e.g., "breast cancer, body height") are classified as
disease rather than measurement. The priority chain:

1. `MONDO_` URI prefix -> `disease` (structured ontology, highest confidence)
2. `OBA_` URI prefix -> `other_measurement` (biological attribute ontology)
3. Cancer keywords -> `cancer`
4. Drug response keywords -> `drug_response`
5. Immune keywords -> `immune`
6. Cardiovascular keywords -> `cardiovascular`
7. Metabolic keywords -> `metabolic`
8. Neurological keywords -> `neurological`
9. Disease catch-all keywords -> `disease`
10. Body measurement keywords -> `body_measurement`
11. Lipid keywords -> `lipid_measurement`
12. Hematological keywords -> `hematological_measurement`
13. Other measurement keywords -> `other_measurement`
14. Behavioral keywords -> `behavioral`
15. "measurement" suffix catch-all -> `other_measurement`
16. Default -> `other`

### Field source (amended 2026-05-20)

Keyword matching runs against the concatenation of `MAPPED_TRAIT` (EFO-mapped
ontology label) and `DISEASE/TRAIT` (raw study label). GWAS Catalog's EFO
mapping is inconsistent — UKB data-field rows frequently have empty
`MAPPED_TRAIT` but informative `DISEASE/TRAIT` (e.g., "Impedance of whole body
(UKB data field 23106)"). Either field populates the substring keyword space,
with no preference between them.

### Implementation

- `classify_gwas_trait()` in `gwas_loader.py` performs keyword matching
  against the lowercased concatenation of `MAPPED_TRAIT` and `DISEASE/TRAIT`.
  Keyword lists are module-level tuples (`_CANCER_KW`, `_LIPID_KW`, etc.).
- `_EXCLUDED_TRAIT_CATEGORIES` frozenset in `gwas.py` defines the excluded
  set. The annotator's `filter_traits` flag controls whether filtering is
  applied (default: True via `get_annotators()`).
- `--gwas-all` CLI flag sets `filter_traits=False`, disabling category
  filtering. Combined with `--gwas-min-magnitude 0`, this produces the
  complete unfiltered GWAS dump.

### Verified against real data

Full GWAS Catalog (930k raw rows, 796k after dedup):
- 77.9% classified as excluded categories
- At mag >= 9.0: 18,675 total, 605 pass trait filter
- Remaining "other" at mag-9: 88 rows (health traits, rare conditions)
- Combined with ClinVar (~139) and PharmGKB (~21): ~765 total annotations
  in default-invocation report

## Structural noise detection — step 1.5 (amended 2026-05-20)

The keyword-based classifier missed three structural noise patterns that
leaked 11 UKB body-composition and metabolite-ratio rows into the "other"
category at mag-9:

1. **UKB bioimpedance traits.** "whole body water mass", "impedance of arm",
   "impedance of trunk", etc. — UK Biobank data-field traits that are body
   composition measurements but don't contain the existing `_BODY_MEASUREMENT_KW`
   keywords. Fixed by adding 7 keywords to `_BODY_MEASUREMENT_KW`.

2. **NMR metabolite ratios.** "cholesterol-to-phospholipid ratio in small VLDL",
   "triglyceride-to-phospholipid ratio in large HDL", etc. — Nightingale NMR
   metabolomics traits. These contain lipid substrings but the "-to-…ratio"
   structure is a more reliable signal. Fixed by `_is_metabolite_ratio()`:
   returns True when both `-to-` and ` ratio` appear in the trait.

3. **Uncharacterized analytes.** "X-12345 level" — Metabolon platform traits
   with no assigned identity. Fixed by `_is_uncharacterized_analyte()`:
   returns True when the trait starts with `x-` and contains ` level`.

Both helper functions are checked after URI prefix detection and before the
keyword chain, so MONDO_ diseases are never misrouted. The helpers' structural
patterns (`-to-` + ` ratio`, `x-` prefix + ` level`) are narrow enough that
false positives on disease traits are not a concern.

**Tests:** `TestStructuralNoiseDetection` in `test_gwas.py` with parametrized
tests for all three noise patterns plus non-misrouting assertions for diseases.

## Must-include rsID allowlist (amended 2026-05-20)

The per-source magnitude floor (`--gwas-min-magnitude`, default 9.0) correctly
suppresses the long tail of common-trait GWAS noise, but it also hides a small
set of clinically significant GWAS associations whose p-values fall in the
genome-wide-significant tier (mag 6.0–7.0) rather than hyper-significant (mag
8.0–9.0). These are well-replicated associations with clear clinical utility
that should appear in reports regardless of the GWAS magnitude floor.

**Decision.** A `_MUST_INCLUDE_RSIDS` frozenset in `gwas.py` lists rsIDs
whose GWAS annotations bypass the per-source magnitude floor
(`source_min_magnitudes`). The global `--min-magnitude` floor still applies.
Trait-category filtering still applies. The carrier rule still applies.

### Initial allowlist

| rsID | Gene | Condition |
|---|---|---|
| rs10737680 | CFH | Age-related macular degeneration |
| rs11209026 | IL23R | Inflammatory bowel disease |
| rs9271366 | HLA-DRB1 | Multiple sclerosis |

### Implementation

- `Annotation` dataclass gains `is_must_include: bool = False`.
- `GWASCatalogAnnotator.annotate()` sets `is_must_include=True` when
  `variant.rsid in _MUST_INCLUDE_RSIDS`.
- `AnalysisResult.filter()` skips the `source_min_magnitudes` elevation
  for annotations where `is_must_include` is True. The base `min_magnitude`
  floor is still applied.

### What it does NOT do

- Does not override the global `--min-magnitude` threshold.
- Does not override trait-category filtering (`--gwas-all` is separate).
- Does not override the carrier rule — the user must carry the risk allele.
- Does not affect non-GWAS sources.

**Tests:** `TestMustInclude` in `test_gwas.py` with constant check, carrier
flag test, source floor bypass test, global min-magnitude test, and trait
filter test.

## Consequences

- Tests pin each p-value tier boundary and the OR modifier via
  `TestMagnitudeScoring` in `tests/annotators/test_gwas.py`.
- `_UNKNOWN_RISK_ALLELE_MAG_CAP` is pinned at 3.0 by
  `test_unknown_risk_allele_cap_value_is_3_0`.
- Unknown-risk-allele rows are always identifiable by their description text
  and by their magnitude ≤ 3.0.
- Future strand normalization (R-1) is a latent gap, not a blocking issue.
  Documented here so R-1 implementation can reference this ADR.
- If the magnitude tiers need revision, update this ADR and the
  `_magnitude()` function together.
- The per-source floor keeps GWAS output tractable without hiding data —
  `--gwas-min-magnitude 0` is the escape hatch for full dumps.
- Trait-category filtering reduces default GWAS output to disease-relevant
  associations. `--gwas-all` is the escape hatch for full trait coverage.
- `TestTraitClassifier` and `TestTraitFiltering` in `test_gwas.py` pin
  the classification logic and the excluded category set.
- Adding keywords to the classifier is a maintenance task: update the
  keyword tuple, add a test, amend this ADR.

## Cache invalidation on categorizer change (amended 2026-05-20)

The `trait_category` column is populated at ingest time. Code changes to
`classify_gwas_trait()` are invisible to `schema_is_current()` unless
explicitly signaled — the column set doesn't change, so the PRAGMA check
passes on stale caches.

**Decision.** `_CATEGORIZER_VERSION` is an integer constant in
`gwas_loader.py`. `load_gwas_tsv()` appends `|cv:N` to the
`database_versions.remote_signal` value on ingest. `schema_is_current()`
verifies the cached signal contains `|cv:{current}`. Mismatch returns
False, causing `db update` to rebuild the cache.

Increment `_CATEGORIZER_VERSION` whenever `classify_gwas_trait()` content
semantics change (new keywords, new helpers, field source changes).
Schema-only changes (new columns) are still caught by the existing
`_REQUIRED_GWAS_COLUMNS` check.

**Tests.** `TestCategorizerVersion` in `test_gwas_loader.py` pins the
accept/reject behavior for matching, missing, and stale `cv:` markers.

### Auto-reingest from cached TSV (amended 2026-05-22)

The raw GWAS Catalog TSV (~200 MB) is retained on disk after the first
successful ingest. When `_CATEGORIZER_VERSION` bumps and the cached TSV is
still present, `is_ready()` auto-reingests from the local file without
re-downloading from EBI.

This eliminates the "code fix shipped but cache not rebuilt" failure mode:
users who have already downloaded the GWAS data get the updated
classification automatically on first run, without network access.

**Implementation.** `GWASCatalogAnnotator.setup()` deletes only the
downloaded ZIP archive; the extracted TSV at
`data_dir/gwas_catalog_associations.tsv` is kept. `is_ready()` checks
`schema_is_current()` first; on failure, checks for the TSV and calls
`load_gwas_tsv()` if found. The cached `remote_signal` is preserved
through the reingest so freshness comparison still works.

If the TSV has been manually deleted, `is_ready()` returns False and the
user sees the standard "run `db update`" message — same as before.

**Tests.** `TestAutoReingest` in `test_gwas.py` with two tests:
auto-reingest when categorizer bumped with TSV present, and False when
TSV is missing.

## MTAG and PheCode rollup (amended 2026-05-21)

GWAS Catalog publishes the same finding under multiple labels: MTAG
re-analyses appear as "(MTAG)" suffix variants of plain traits, and PheCode
hierarchical codes (411 / 411.4 / 411.8) appear as sibling rows for the same
underlying disease family. Allelix collapses both classes after the
magnitude/category filter, before rendering.

**MTAG dedup.** When (rsid, base-trait) has both plain and (MTAG) versions, drop
the (MTAG) row. MTAG is multi-trait analysis of the same underlying GWAS
data; the plain single-trait row is the canonical finding. A solo (MTAG) row
with no plain twin is kept.

**PheCode hierarchical rollup.** When the same rsid has multiple PheCode rows in
the same parent group (numeric prefix before the dot), keep the row with
the strongest p-value and drop the rest.

**Must-include exemption.** Rows flagged `is_must_include=True` are never dropped
by rollup, only by global magnitude.

Rollup is conservative — distinct concepts (different base traits, different
PheCode parent numbers) are never collapsed. Substring detection only, no
regex (ADR-0016).

**Implementation:** `rollup_gwas_duplicates()` in `allelix/reports/_pipeline.py`.
Called by all three renderers (terminal, JSON, HTML) immediately after
`AnalysisResult.filter()`.

**Tests.** `test_rollup.py` in `tests/reports/` covers MTAG twin collapse,
solo MTAG keep, PheCode parent/child collapse (strongest p wins), distinct
parents kept, must-include exemption, non-GWAS passthrough, and the
rs10455872 8→5 real-data case.

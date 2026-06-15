# ADR-0016: Data Classification Principle

- **Date:** 2026-05-13
- **Status:** Accepted (non-negotiable architectural directive)

## Context

Earlier prose-regex approaches to PharmGKB non-finding detection missed phrasings with each release; this directive ends that arms race.

The fundamental error: classifying database records by regexing English prose written by curators for human readers. Prose is not a stable interface. The patterns drift with editorial style. The same regex that catches one wave of phrasings misses the next.

PharmGKB itself stores the answer in a structured field. The `clinical_ann_alleles.tsv` `Allele Function` column has values "Normal function" / "Decreased function" / "No function" / "Increased function" / empty. That is the classification input. The annotation text is for the user to read.

The same principle applies to every annotator. Structured fields are authoritative; prose is for display.

## Decision

**Classification uses structured database fields. Regex on prose is forbidden in production code.**

Concretely:

1. **Each annotator's ingestion layer extracts structured classification fields.** PharmGKB: `Allele Function`. ClinVar: `CLNSIG`, `CLNREVSTAT`. GWAS Catalog: `p_value`, `effect_size`. SNPedia: magnitude/genotype structured entries.

2. **Filtering and scoring decisions reference structured columns.** The annotator's SQL WHERE clause and the report's filter dataclass operate on enumerated values, not text matches.

3. **If the source has no structured field for a needed decision, that decision is not automated.** Document the gap; surface the rows with a flag the user can interpret. Do not approximate structure with regex.

4. **Regex on description text is permitted only in the test suite as a safety net.** A safety-net test loads the canonical mock fixture, runs the regex over each row's text, and asserts the structured classifier reached the same conclusion. A safety-net failure means the structured classifier missed a case the regex caught — fix the ingestion to extract the structured signal, NOT add another production regex.

### How this applies to the existing PharmGKB filters

- **Non-finding suppression (ADR-0013):** the `_NONFINDING_PATTERNS` regex is replaced by reading PharmGKB's `Allele Function` column. A row whose matched genotype has `Allele Function = "Normal function"` is the non-finding; classify at load time, filter at query time. The old regex moves to `tests/test_pharmgkb_safety_net.py` and runs as a smoke check.

- **Somatic suppression (ADR-0014):** PharmGKB exposes no structured germline/somatic flag at the row level. The regex-based suppression is removed entirely. In practice the common case is handled by non-finding suppression (somatic-context reference rows are typically `Allele Function = "Normal function"`). Rare residual cases surface to the user; ADR-0014 is amended to document the gap and the residual.

### Schema migration

A new `function_class TEXT` column is added to `pharmgkb_annotations` storing the normalized classification (`normal` / `decreased` / `no_function` / `increased` / `unknown`). v0.5.x caches lack this column, so `schema_is_current()` returns False and the next `db update` automatically refreshes into the v0.6.0 schema with structured classification.

## Consequences

- The PharmGKB filter is stable against editorial style changes. PharmGKB can rewrite annotation text without our filter degrading.
- The category of "regex-pattern-arms-race" bugs is closed for the entire codebase. Every future annotator added to the project is bound by the same principle.
- Some decisions become unautomatable. Somatic-vs-germline is the visible example. The trade is correct: surfacing slightly more rows that the user evaluates is better than silently classifying based on regex hits that drift.
- The safety-net regex serves a real role — it's the canary that warns us when the structured ingestion missed a field. It is not, and cannot be, a fallback when the structured field is empty.
- Audit of the existing codebase identified exactly two violations (`_NONFINDING_PATTERNS`, `_SOMATIC_PATTERN`), both in `pharmgkb_loader.py`. Both refactored in v0.6.0.
- The principle is a top-level project policy. PRs that violate it are rejected.

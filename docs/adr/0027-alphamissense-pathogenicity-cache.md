# ADR-0027: AlphaMissense Pathogenicity Cache

**Status:** Accepted
**Date:** 2026-06-09

## Context

Allelix v1.4.0 adds missense variant pathogenicity predictions from
DeepMind's AlphaMissense (Cheng et al., Science 2023). AlphaMissense scores
71M missense variants on a 0–1 scale using protein structure predictions.
This complements gnomAD (ADR-0026), which answers "how common is this
variant?" — AlphaMissense answers "how likely is this missense change to
damage the protein?"

AlphaMissense source data is coordinate-keyed (chrom/pos/ref/alt), not
rsID-keyed. Consumer genotype arrays report rsIDs. A build-time join
against the gnomAD cache provides the coordinate-to-rsID mapping.

## Decision

Ship the AlphaMissense cache as a pre-built SQLite asset distributed via
HuggingFace (`dial481/allelix-alphamissense`). Schema uses composite
primary key `(chrom, pos, ref, alt)` — same pattern as gnomAD (ADR-0026)
— to preserve multi-allelic variants. An index on `rsid` supports the
annotator's lookup path.

The build script (`scripts/build_alphamissense_cache.py`) streams the
source TSV from Zenodo over HTTPS by default (never saves the ~3.6 GB
source file to disk) or reads a pre-downloaded local file. The gnomAD
cache must be built first to provide the rsID mapping; `--no-gnomad`
builds without rsIDs (not useful for allelix but unblocks testing).

AlphaMissense is an enrichment-only annotator — it stamps
`am_pathogenicity` and `am_class` on existing annotations via
`bulk_lookup()`, same as gnomAD stamps `allele_frequency`. It does not
participate in the per-variant annotation loop.

### PharmGKB caveat

AlphaMissense measures protein structure impact. PharmGKB annotations
describe drug-response phenotypes, which depend on pharmacokinetic
pathways, not just protein folding. An AM score on a PharmGKB row is
technically accurate but clinically misleading. All three renderers
mark PharmGKB AM scores as neutral:

- HTML: grey `am-score` class with "protein structure impact only" tooltip
- Terminal: dimmed score with `*` marker and footnote
- JSON: `am_caveat: "protein structure impact only"` field

### MAX aggregation caveat

Both gnomAD and AlphaMissense use `MAX()` aggregation with `GROUP BY
rsid` for their `bulk_lookup()` queries. At multi-allelic sites (same
rsid, different alt alleles), this returns the highest score regardless
of which alt the user actually carries. The proper fix requires threading
the alt allele through the enrichment lookup path, which touches the
Annotation model, pipeline keying, and both enrichment annotators. This
is deferred to a post-v1.4.0 patch.

## Consequences

**Positive:**

- Adds computational pathogenicity context to every missense variant in
  the report, especially valuable for the thousands of variants ClinVar
  hasn't reviewed
- Same distribution pattern as gnomAD — familiar to users, reuses
  `install_prebuilt_cache` infrastructure
- Composite PK preserves multi-allelic sites (no silent data loss)
- Streaming build mode avoids 3.6 GB disk requirement for source data
- CC BY 4.0 license is compatible with commercial use

**Negative:**

- Large download (~8 GB on disk after decompression)
- Build-time gnomAD dependency adds ordering constraint
- MAX aggregation over-reports at multi-allelic sites (documented, deferred)

**Future:**

- CADD variant deleteriousness scores (v1.5.0) — same enrichment pattern,
  non-commercial license requiring config system gating
- Alt-aware enrichment lookups — proper fix for MAX aggregation caveat

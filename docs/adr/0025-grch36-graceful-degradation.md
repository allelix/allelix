# ADR-0025: GRCh36 Graceful Degradation (No Full Build Support)

**Status:** Accepted
**Date:** 2026-06-07

## Context

Real genotype files from ~2010–2012 (early FTDNA, 23andMe) use GRCh36 (hg18)
coordinates. Allelix detects these files correctly via position-based build
detection (P-B, shipped v1.1.0) and warns the user that ClinVar annotations
will be skipped (v1.0.0+).

The original roadmap (R-12 sub-item 2) proposed shipping a GRCh36 ClinVar
cache — either by downloading NCBI's GRCh36-mapped ClinVar VCF or by building
a liftover pipeline (GRCh36 → GRCh37 coordinates). Both approaches have
significant costs:

- **NCBI no longer publishes ClinVar VCFs on GRCh36 coordinates.** The last
  GRCh36 ClinVar VCF was retired years ago.
- **Liftover** requires a new coordinate-mapping layer, chain file management,
  and testing against edge cases (indels near rearrangement breakpoints,
  split/merge regions). This is a substantial architectural addition for a
  build that no active genotyping service exports.
- **User population is tiny.** GRCh36 files come from services that stopped
  using build 36 over a decade ago. Users with these files can re-export on a
  newer build or run a one-time liftover themselves.

## Decision

Allelix will **not** ship a GRCh36 ClinVar cache or a built-in liftover
pipeline. Instead, it implements graceful degradation:

1. **Detect GRCh36 accurately.** The 3-way probe table (11 SNPs × 3 builds)
   identifies GRCh36 files even without header metadata. Shipped in v1.1.0.

2. **Run rsID-based pipelines normally.** PharmGKB, GWAS Catalog, SNPedia, and
   gnomAD all use rsID lookups — build-independent. These fire on GRCh36 files
   exactly as they do on GRCh37/38.

3. **Skip ClinVar position matching.** ClinVar queries by chromosome+position
   against per-build caches. With no GRCh36 cache, `annotate()` returns `[]`
   for GRCh36 variants. This is safe — no false matches, no silent
   misannotation.

4. **Emit an explicit warning.** The CLI tells the user what happened and why,
   and points to documentation with liftover instructions.

5. **Provide liftover documentation.** A docs page with copy-paste one-liners
   for UCSC liftOver and CrossMap, so users can convert their files to GRCh38
   before re-running allelix for full ClinVar coverage.

## Consequences

- R-12 sub-item 2 (GRCh36 ClinVar cache) is removed from the roadmap entirely.
  It is not deferred — it is cancelled.
- No liftover dependency is added to allelix's install footprint.
- Users with GRCh36 files get partial but accurate results, with a clear path
  to full results via external liftover.
- The CLI warning must link to the liftover docs page (v1.3.1).

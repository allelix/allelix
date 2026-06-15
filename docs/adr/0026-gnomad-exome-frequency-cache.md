# ADR-0026: gnomAD Exome Frequency Cache

**Status:** Accepted
**Date:** 2026-06-07

## Context

Allelix v1.3.0 adds population allele frequencies from gnomAD v4.1 exomes
(730K individuals, protein-coding regions). This data enriches existing
annotations with allele frequency context.

The cache must store per-allele variant data (not just per-rsID) to support
future joins with coordinate-keyed databases like AlphaMissense and CADD.

## Decision

Ship the full gnomAD v4.1 exome cache as a pre-built SQLite asset distributed
via HuggingFace (`dial481/allelix-gnomad`). Schema uses composite primary key
`(chrom, pos, ref, alt)` to preserve multi-allelic variants. An index on `rsid`
supports the current annotator lookup path.

The build script (`scripts/build_gnomad_cache.py`) supports HTTP streaming and
local file modes for reproducibility.

## Consequences

**Positive:**

- Full exome coverage — all ~16M variants with rsIDs
- Per-allele storage preserves multi-allelic sites (no silent data loss)
- Coordinate columns and composite PK enable future AlphaMissense/CADD joins
  without rebuilding
- HuggingFace distribution has no file size limits and supports ETag-based
  freshness checking

**Negative:**

- Larger download than a filtered subset
- HuggingFace is an external dependency

**Future:**

- AlphaMissense pathogenicity scores (coordinate-keyed) — shipped in v1.4.0
  (ADR-0027)
- CADD scores (coordinate-keyed) — planned for v1.5.0
- gnomAD whole-genome frequencies (non-coding regions) — distant future, no
  timeline

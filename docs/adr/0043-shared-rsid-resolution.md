# ADR-0043: rsID-less resolution as a shared step

**Status:** Accepted
**Date:** 2026-10-10
**Builds on:** ADR-0021 and ADR-0037 (build detection), ADR-0039
(`resolve_rsid_at_locus`), ADR-0040 (identity table), ADR-0042 (high-value
positions on VCF input).

## Context

Variant callers leave the VCF ID column empty. `analyze` recovered rsIDs
by position inside the pipeline's batch flush. Every other command that
compared `Variant.rsid` against an rsID-keyed set never reached that code:

- `stats` (the high-value list);
- `export plink` (gnomAD ref/alt coding);
- build detection (the seed SNP table).

`extract` carried its own copy of a gnomAD coordinate match. Each new
entry point was exposed to the same miss.

## Decision

### 1. One resolution step

`allelix.annotators._resolution` holds the whole chain. `resolve_batch`
works on a list and `resolve_rsidless` streams an iterator in batches.
The chain runs in this order:

1. `Variant.ref` from the identity table, for rsID-carrying calls only.
2. ClinVar position resolution.
3. Identity-table position resolution.
4. gnomAD position resolution (GRCh38).
5. The gnomAD REF fallback, checked against the variant's own position.

The order matters. REF population comes first so ClinVar's REF check can
confirm array calls. A vendor probe without an rsID never gets a REF from
position alone (the ADR-0040 guard), so no resolver touches it.

Users of the step:

- **The pipeline's batch flush** calls `resolve_batch`.
- **`stats`, `extract` and `export plink`** go through `_resolved_variants`
  in the CLI helpers. It runs the same build detection as `analyze` first,
  because resolution needs the build. Sources whose caches are missing are
  skipped.
- **`stats`** also watches the high-value positions as `analyze` does
  (ADR-0042). An rsID-less `./.` record at a high-value position, which no
  allele-based resolver can identify, is still reported. Resolved calls
  count as called cluster partners.
- **`extract`'s tabix path** fetches rows at gnomAD's coordinates for each
  wanted rsID. It then passes them through `resolve_batch` instead of its
  own (ref, alt) comparison.

### 2. Build detection votes by position

Build detection cannot call the resolver. Resolution needs the build being
detected, and stamping rsIDs under a guessed build would leave wrong rsIDs
behind whenever the guess lost.

The seed table already holds every seed's coordinate on all three builds,
and no coordinate repeats (pinned by a test). So a record without an `rs`
ID at a seed's coordinate on build B votes for B, through
`known_snp_seed`. A record that names a different rsID at that coordinate
is a different variant and does not vote. The ADR-0037 thresholds are
unchanged: three concordant seeds and an 80% ratio.

A gVCF reference block covering a seed is not counted. Blocks tile the
genome, so they cover every build's seed coordinates alike, and their REF
gives only the block's first base.

## Consequences

- `stats` reports high-value no-calls and incomplete clusters on rsID-less
  VCFs. `export plink` gives gnomAD-resolvable rsID-less calls a real
  A1/A2 instead of monomorphic coding. `extract` and `analyze` use the same
  rule.
- An rsID-less GRCh37 file is detected from its positions. That includes a
  UCSC-style `chr`-prefixed hg19 file, which the `chr`-prefix heuristic
  alone would call GRCh38.
- A future command that needs rsIDs calls `_resolved_variants` (or
  `resolve_rsidless`) instead of adding another copy.
- `stats` now reads the data directory and gains `--data-dir`.

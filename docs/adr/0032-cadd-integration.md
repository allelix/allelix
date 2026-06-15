# ADR-0032: CADD v1.7 Integration

**Status:** Accepted
**Date:** 2026-06-10

## Context

CADD (Combined Annotation Dependent Depletion) ranks how deleterious any
single-nucleotide variant is, using over 100 annotation tracks. A
PHRED-scaled score where 10 = top 10% most deleterious, 20 = top 1%,
30 = top 0.1%. Unlike AlphaMissense (protein-structure only) or gnomAD
(frequency only), CADD covers coding, non-coding, and regulatory
variants — making it a complementary enrichment signal.

CADD is the first non-commercial data source in allelix. The CADD
license (University of Washington / Hudson River Biotechnology) permits
free non-commercial use; commercial use requires a separate license.
This breaks the assumption that all sources can be downloaded silently
via `db update`.

The full CADD dataset is ~8.8 billion SNVs (~81 GB compressed) plus a
separate indel prescored file (~1.2 GB, ~320M rows). The allelix
enrichment pipeline resolves rsIDs to genomic coordinates via gnomAD,
so only positions present in allelix's existing databases are ever
queried. The pre-built cache filters to this position universe.

## Decision

### Enrichment-only annotator

CADD follows the gnomAD / AlphaMissense enrichment pattern:
`annotate()` returns `[]`. The pipeline calls `bulk_lookup()` after all
annotators have run and stamps each annotation's `cadd_phred` field.
Coordinate resolution and strand normalization go through gnomAD's
`bulk_resolve_coordinates()` and `resolve_strand()`.

### Disabled by default

`sources.cadd = false` in the default config. Users opt in via:
- `allelix config set sources.cadd true` (persistent), or
- `allelix db update --cadd` (one-shot download)

The first download shows a non-commercial license confirmation prompt.

### Two modes

**Cache mode** (default): pre-built SQLite from HuggingFace filtered
to positions present in gnomAD, AlphaMissense, and ClinVar (GRCh38).
~5 GB on disk, ~120M keys (117M SNV + 3M indel). The build script uses
int64 packing for SNV keys and a separate tuple set for indel keys.
Fast lookups.

**Full mode** (`options.cadd_full = true`): queries the complete CADD
tabix file locally via pysam. Covers every scored position. Requires
GRCh38 — the pipeline skips enrichment with a warning if the detected
build is not GRCh38. pysam is an optional dependency
(`pip install allelix[cadd]`).

### When full mode matters

Cache mode covers the large majority of variants present in gnomAD,
AlphaMissense, and ClinVar — nearly every position allelix can annotate
from its other databases. For genotyping chip data (23andMe, AncestryDNA,
MyHappyGenes, etc.), cache and full mode produce effectively identical
results because chip probes overwhelmingly target known, cataloged
variants. Full mode adds coverage for novel or private
variants that appear only in whole-genome or whole-exome sequencing data
and are not in any pre-computed database. If your input is a genotyping
chip file, cache mode is all you need.

### Strand normalization

Array genotype data may report alleles on the minus strand. CADD scores
are reference-forward. `resolve_strand()` maps array alleles to
reference orientation using gnomAD ref/alt as ground truth. Palindromic
SNPs (A/T, C/G ref/alt pairs) cannot be resolved and return None — the
variant gets no CADD score rather than a wrong one.

### Indel coverage

The CADD indel prescored file (`gnomad.genomes.r4.0.indel.tsv.gz`,
~1.2GB) covers gnomAD-observed indels. The build script processes it
in a second pass after the SNV file. In full mode, indel lookups route
to the indel tabix file (optional — graceful degradation if absent).

### License descriptor

```python
license = LicenseDescriptor(
    spdx="LicenseRef-CADD",
    license_url="https://cadd.gs.washington.edu/license",
    commercial_ok=False,
    licensable=True,
    purchase_url="https://els2.comotion.uw.edu/product/cadd-scores",
)
```

`commercial_ok=False` marks CADD as non-commercial. `licensable=True`
with a `purchase_url` means a commercial license can be purchased.
The three-state permission model resolves each source to one of:

- **ALLOW** — source is usable (non-commercial mode, or commercial with
  license held).
- **BLOCK_PURCHASABLE** — blocked in commercial mode, but a license can
  be purchased. `config show` displays the purchase URL.
- **BLOCK_FINAL** — blocked in commercial mode with no license available
  (e.g., SNPedia).

Users assert a held license via `allelix config set license.cadd true`.

## Related ADRs

- [ADR-0028](0028-local-version-tag-convention.md): local version tag convention — CADD cache uses `sv:` tag for schema versioning
- [ADR-0029](0029-download-integrity-verification.md): download integrity verification — CADD download uses SHA-256 hash verification
- [ADR-0030](0030-data-source-update-cadence.md): data source update cadence — CADD is a code-driven source with pinned URL
- [ADR-0031](0031-centralized-license-descriptors.md): centralized license descriptors — CADD's `commercial_ok=False` is enforced by the `__init_subclass__` mechanism

## Consequences

- CADD is the first opt-in source. The pattern (disabled-by-default +
  confirmation prompt + `commercial_ok=False`) is reusable for future
  non-commercial sources.
- The `cadd_phred` field on Annotation is always None when CADD is not
  enabled. Report renderers conditionally show the column only when at
  least one annotation carries a score.
- The build script depends on existing gnomAD, AlphaMissense, and
  ClinVar (GRCh38) caches for coordinate filtering. gnomAD is required;
  AlphaMissense and ClinVar are optional. Building CADD without gnomAD
  is not supported.
- JSON schema version stays at "3" — `cadd_phred` is a new optional
  field that `asdict()` serializes automatically.

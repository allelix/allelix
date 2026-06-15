# ADR-0010: Strand-flip helpers ship in v0.4.0; liftover is deferred

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Two related coordinate-system problems block clean cross-database matching:

1. **Strand orientation.** Different sources publish variants on different strands. The same biological variant `rs123 G→A` on the forward strand reads as `C→T` on the reverse. For most allele pairs a strand flip is reversible by complementing each base. But `A↔T` and `C↔G` SNPs are palindromic — their complement equals their alternative, so strand cannot be inferred from sequence alone.

2. **Reference build conversion.** Same SNP, different coordinates in GRCh37 (hg19) vs GRCh38 (hg38). Most consumer DNA tests are GRCh37; ClinVar and PharmGKB now publish both. WGS providers and academic studies increasingly use GRCh38 only. Cross-build matching needs liftover (UCSC chain-file conversion).

Both problems are real, but they have very different costs.

## Decision

**Ship in v0.4.0:** strand-flip helpers in `allelix.utils.allele`:

- `complement(allele)` — single base or multi-base indel (reverse-complement).
- `flip_genotype(allele1, allele2)` — diploid call flipped to the other strand.
- `is_strand_ambiguous(ref, alt)` — true for A/T and C/G pairs.

These are stdlib-only, ~20 lines, and useful immediately for parser/annotator authors who need to normalize cross-strand calls. None of the current annotators *invoke* them yet — the existing fixtures align by construction. They're available for the next parser (23andMe / VCF) which is more likely to need them.

**Deferred:** liftover. Deferred until a real GRCh38-only data source forces the work; tracked as a known limitation in CHANGELOG. The reasons:

- Liftover requires UCSC `.chain.gz` files (~50MB hg19↔hg38) that we don't want to bundle (privacy/license-wise neutral, but bandwidth and version-tracking are non-trivial).
- Real implementations (`pyliftover`, `pysam`'s liftover) are heavy dependencies pulled in for a feature most users won't need until they swap consumer DNA test for WGS, or until a GRCh38-only database lands.
- Until a real-world need surfaces, liftover would be untested-against-real-data infrastructure that ages badly.

The first user with a GRCh38 input file or a GRCh38-only database is the forcing function. A future ADR (0011 or later) will pick the implementation when that arrives.

## Consequences

- v0.4.0 ships strand helpers; annotators can opt-in to using them when needed (no global change to matching logic in this release).
- ambiguous-SNP false matches remain possible but rare and small in count.
- ADR-0021 subsequently added dual-build support (GRCh37 + GRCh38) with auto-detection. Liftover remains deferred.
- A future ADR will document the liftover decision (data dir, chain-file source, UI surface) once a real case forces it.
- The v0.4.0 utility tests (`tests/utils/test_allele.py`) lock in the small contract: complement is reversible, ambiguity detection covers exactly A/T and C/G, indels are reverse-complemented base-by-base.

# ADR-0009: PharmGKB matches the user's exact normalized diploid call

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

ADR-0007 established the carrier rule for ClinVar: an annotation triggers if the user carries at least one ALT allele at a position. That rule fits ClinVar's REF/ALT VCF schema cleanly.

PharmGKB does not publish in a REF/ALT model. It publishes per-genotype clinical predictions: for each (variant, drug) pair, separate annotations describe what's known for each diploid genotype (e.g., "rs1801133 + AG → consider lower methotrexate dose; rs1801133 + GG → normal toxicity profile"). The relevant question is not "does the user carry the ALT?" but "what does PharmGKB say about the *exact* genotype the user has?"

Naively reusing the carrier rule would either fire too many annotations (matching any allele in any PharmGKB row) or none (PharmGKB rows don't have a single "ALT" to compare against).

## Decision

For PharmGKB, an annotation triggers if and only if the user's normalized 2-letter SNV genotype equals the annotation's normalized genotype:

```python
user_geno = normalize(variant.allele1 + variant.allele2)
# row triggers iff user_geno == row.genotype (also normalized at load time)
```

Normalization rules (`_normalize_genotype`):

- Strip `:`, `;`, `/` separators and uppercase.
- Reject anything that isn't exactly two A/C/G/T characters (no indels, star alleles, complex variants in v0.3.0).
- Sort alphabetically so `AG` and `GA` collide at the same key.

Both sides — the loader writing rows into SQLite and the annotator querying — apply the same normalization. SQLite's index on `rsid` makes the lookup O(log n).

What this excludes from v0.3.0:

- **Star alleles** (`CYP2D6*1/*2`): require haplotype reconstruction from multiple SNPs. Deferred.
- **Multi-rsid composite annotations** (`rs1801133, rs1801131`): require multi-SNP genotype reasoning. Deferred.
- **Indel genotypes** (`CTT/C`): rare in PharmGKB clinical annotations and ambiguous in their per-genotype model.

Skipped rows are silently dropped at parse time — the annotator never sees them.

## Consequences

- The carrier-rule contract from ADR-0007 still holds at the annotator-interface level: only annotations for variants the user actually carries are emitted. The implementation differs (exact match vs ALT-presence) because the source models differ.
- PharmGKB never fires for a no-call (`-/-`) — `_normalize_genotype` returns None.
- PharmGKB never fires for indels in v0.3.0 — the same rule. A future ADR will lift this restriction when haplotype support lands.
- Star alleles being skipped means CYP2D6, CYP2C19, and similar genes that PharmGKB primarily annotates by haplotype are underserved in v0.3.0. Documented in the README's Supported Databases note.
- The synthetic fixture (`tests/fixtures/mock_pharmgkb/`) deliberately includes a star-allele row (PA-005) and a multi-rsid row (PA-006) so tests pin both skip behaviors.

# ADR-0040: Cross-build rsID identity table

**Status:** Accepted
**Date:** 2026-10-10
**Builds on:** ADR-0035 (allele model), ADR-0039 (locus-keyed selection).
**Replaces:** the 2.3.1 rule that gnomAD, AlphaMissense and CADD are skipped
for every GRCh37 analysis. ADR-0010's deferral of general liftover still
stands: this table translates only the rsIDs Allelix can report.

## Context

gnomAD v4.1, AlphaMissense and CADD v1.7 are GRCh38 datasets. Release 2.3.1
closed a wrong-locus hazard by skipping all three for any analysis that was
not on GRCh38: an rsID or coordinate collision could attach a valid GRCh38
value to the wrong allele on another assembly. That was the right emergency
call, but every consumer chip format (23andMe, AncestryDNA, FTDNA, LivingDNA,
MyHeritage) is GRCh37. As a result, those inputs lost population frequency,
missense pathogenicity and deleteriousness scores, and `Variant.ref` was
never populated for them, so strand-aware matching had no reference allele.
rsID-less GRCh37 VCFs could also recover rsIDs only from ClinVar.

The hazard is real only where the reference strand was inverted between
assemblies (the GRCh37 forward allele is the complement of the GRCh38 one)
or where the strand relation cannot be decided. That set can be enumerated.

## Decision

### 1. One table, built from dbSNP

`rsid_identity.sqlite` holds one row per allele of every rsID the
annotators can report (ClinVar, ClinPGx, GWAS Catalog, SNPedia, plus the
high-value and build-detection seed lists):

```
rsid | chrom | pos37 | ref37 | alt37 | pos38 | ref38 | alt38 | palindromic
```

`alt37` and `alt38` are the same allele written on each build's forward
strand. Either can equal its build's REF where the reference allele changed
between assemblies. `scripts/build_rsid_identity_cache.py` streams the
GRCh37 and GRCh38 dbSNP VCFs once each and keeps only primary-assembly
placements. An rsID is included only when it has exactly one placement per
build on the same chromosome and its allele sets are equal on the same strand
or, for SNVs, after complementing. `palindromic` marks allele sets equal to
their own complement (A/T, C/G, or all four bases), where the strand
relation is undecidable. The build streams each VCF through a Bloom filter
of the rsID universe and pairs placements inside SQLite, so it runs in
tens of MB of memory. The cache
ships like the other prebuilt caches: gzipped on HuggingFace and verified by
SHA-256. dbSNP is in the public domain.

### 2. Three uses, all keyed by locus on the input's own build

- **`Variant.ref`** for rsID-carrying array data on GRCh37 and GRCh38.
  The REF at the variant's locus is used when it is unique for the call's
  shape. Vendor probes without an rsID get none (section 3). On GRCh38,
  gnomAD stays as a fallback, now accepted only when its coordinate is the
  variant's own locus.
- **rsID recovery** for rsID-less VCF records on both builds, using the shared
  `resolve_rsid_at_locus` rule (ADR-0039).
- **Enrichment translation.** Every annotation carries its locus, the REF at
  that locus and its allele. GRCh38 input is looked up at
  `(chrom, pos, ref, alt)` directly. GRCh37 input is translated to its GRCh38
  key only when exactly one non-palindromic row matches and the allele is not
  the GRCh38 reference. gnomAD, AlphaMissense and CADD are all queried by that
  key. No enrichment is keyed by rsID, and the separate GRCh38-only gates and
  rsID-resolution side channel are gone. An annotation without its own allele
  (a GWAS row with no reported risk allele) uses the single non-reference
  allele the user carries, and gets nothing when there are two.

### 3. Guards on locus matching (refines ADR-0039)

Locus-keyed ClinVar selection needs evidence that the row at a position is
the user's variant:

- A call with **neither a reference allele nor an rsID** (array vendor
  probes such as 23andMe `i5000123`) is never matched by locus, and the
  ClinVar, identity and gnomAD position resolvers never stamp an rsID on
  it. The identity table does not supply a REF for such a call either.
  The REF at a position says nothing about the strand a vendor reported
  the probe's alleles on. With an inferred REF, the strand-aware match
  would read the complement of a homozygous-reference call as a
  homozygous pathogenic genotype. A call without an rsID that has a REF
  is therefore always a VCF record whose REF came from the file.
- When the input's **rsID is placed by ClinVar elsewhere** on this build and
  not at the input's locus, the rows at the locus belong to a different
  variant and are not reported under the input's rsID. An rsID that ClinVar
  does not hold at all (a retired or merged ID) still matches by locus. The
  rsID is never a positive selection key; it is used only to reject a
  contradicting locus.

Amended in 2.4.2 (see ADR-0039, Amendment): the REF comes only from the
rows of the call's own rsID, and a call with no REF matches by locus only
when ClinVar holds its rsID.

## Consequences

- GRCh37 chip reports regain gnomAD, AlphaMissense and CADD values wherever
  the identity table proves the allele, and get none where it cannot. Without
  the table installed, GRCh37 behaves exactly as in 2.3.1.
- Strand-aware matching (ADR-0035) gets its reference allele on GRCh37
  arrays. Plausibility flagging can now be built on GRCh37 input.
- rsID-less GRCh37 VCFs recover non-ClinVar rsIDs, which matters for panel
  coverage and high-value no-call detection.
- The dbSNP b157 build holds about 4.9M rows for about 3.0M rsIDs
  (about 440 MB on disk, 150 MB gzipped). About a sixth of those rsIDs
  are palindromic and abstain from enrichment. Most are multi-allelic SNVs
  whose dbSNP allele set includes all four bases, not only A/T and C/G
  sites. Fewer than 4,000 rsIDs pair on the complemented strand. The
  abstention is deliberate. Using REF agreement between builds to decide
  orientation at those sites would be a separate decision.
- The table covers only annotator rsIDs. A full dbSNP resolver for arbitrary
  positions remains optional future work.
- Rebuilding the table (for a new dbSNP release or a larger annotator rsID
  set) means re-running the build script and re-pinning the cache URL and
  SHA-256.

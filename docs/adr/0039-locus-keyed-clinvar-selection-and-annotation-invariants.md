# ADR-0039: Locus-keyed ClinVar selection and Annotation allele-identity invariants

**Status:** Accepted
**Date:** 2026-10-10
**Completes:** ADR-0035 (the allele model: `Variant.ref` plus per-Annotation
`alt`). ADR-0007's carrier rule and ADR-0021's per-build dispatch are unchanged.

## Context

ADR-0035 gave `Variant` a reference allele and every `Annotation` an `alt`,
but two paths still identified a ClinVar variant by rsID alone:

- The array (chip) path queried ClinVar `WHERE rsid = ?` and accepted any row
  whose `alt` the user carried. ClinVar publishes **reference-allele records**
  (HGVS `=`, e.g. `c.1601=`) whose VCF columns carry REF == ALT. A
  homozygous-reference call matched such a row and was reported as a
  pathogenic finding with `alt` equal to the reference base. The rsID query
  also matched rows at another position when the input's build was wrong.
- The rsID-less VCF path resolved an rsID from the rows at a position, then
  annotated by that rsID. The resolver counted duplicate per-SCV rows as
  separate candidates, and it needed both user alleles inside one row's
  `{ref, alt}`. Either case made it abstain, and the site disappeared from
  the report without any signal. With gnomAD installed, its second-pass
  resolver sometimes recovered the rsID and hid the drop.

The same genotype could therefore produce different findings depending on
whether it arrived as a chip row or as a gVCF record.

## Decision

### 1. ClinVar rows are selected by `(chrom, pos, ref, alt)`

Every ClinVar path (`annotate`, `batch_annotate`, `bulk_resolve_rsids`) fetches
the cache rows at the variant's locus (chromosome normalized to ClinVar's bare
form, build-specific cache per ADR-0021). A row annotates only when:

- REF != ALT (reference-allele records never annotate);
- its REF equals `Variant.ref` when the input supplies one, directly or as the
  complement for reverse-stranded input;
- the user carries its ALT (ADR-0007, strand-aware per ADR-0035), with
  ADR-0011 indel-anchor protection.

The rsID is never a selection key. An annotation takes the input's rsID when
it has one; otherwise it takes the matching row's rsID. ClinVar therefore
annotates rsID-less calls without depending on resolution.

### 2. One position-resolver rule

`resolve_rsid_at_locus` is shared by the ClinVar and gnomAD resolvers. A call
resolves when exactly one REF's consistent rows explain every allele the user
carries and those rows name exactly one rsID. Duplicate rows are one
candidate. A `1/2` call resolves when both alts belong to the same rsID.

### 3. Construction invariants on `Annotation`

`Annotation.__post_init__` raises `AnnotationInvariantError` when:

1. `alt` is set and the user carries zero copies of it. This applies to every
   source.
2. `alt` equals the reference allele, for ClinVar, ClinPGx and SNPedia, where
   `alt` names a non-reference allele.

GWAS is exempt from (2). For GWAS rows `alt` is the study's effect (risk)
allele, which is often the reference allele.

The context comes from two internal fields, `locus_ref` and `carried_alleles`
(the user's alleles in `alt`'s orientation). They are excluded from public
serialization, equality and repr, so the output contract and schema version
are unchanged. Annotators filter before constructing; a violation is a
programming error.

## Consequences

- Reference-allele false findings and the rsID-less silent drops are gone,
  and the chip / gVCF parity harness (`tests/test_chip_gvcf_parity.py`) pins
  both paths to the same annotation set.
- An input whose positions disagree with its effective build no longer
  receives ClinVar rows matched by rsID. The cross-build false positive that
  ADR-0021 described for NIPA1 cannot recur through ClinVar even under a wrong
  `--build`.
- ClinVar's REF oracle for ClinPGx / SNPedia non-finding suppression
  (`reference_for`, ADR-0023) remains rsID-keyed. It reads a reference base,
  not a finding; per-build positional identity for array data arrives with
  the cross-build rsID identity table.

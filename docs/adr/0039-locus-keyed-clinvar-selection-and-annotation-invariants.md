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

## Amendment (2.4.2): array calls with an rsID nobody holds

ADR-0040 let an array call whose rsID ClinVar does not hold still match
ClinVar by locus, and the identity table filled its REF from whatever rsID's
rows sat at that position. The withdrawn probe ID rs62508729 reads TT at
the position of rs62507269, a pathogenic PAH variant, so array files that
carry it (for example 23andMe v5 exports) got a false homozygous
pathogenic phenylketonuria finding. From 2.4.2:

- The identity table fills `Variant.ref` only from the rows of the call's
  own rsID. An array call whose rsID the table does not hold gets no REF
  from another rsID at the same position.
- A call with no REF is matched by locus only when ClinVar holds its rsID
  on this build. VCF records keep their REF from the file and are
  unchanged.
- Merged rsIDs (for example rs10735810, merged into rs2228570) also lose
  their position match. They were not matched before 2.4.0 either. A
  merge-aware rule (dbSNP merge history in the identity table) would
  restore them; it is not part of this release.

## Amendment (2.4.2): REF == ALT, and records for the normal allele

This ADR's Context treats every REF == ALT row as a reference-allele record
and names `c.1601=` as the example. That is wrong for Factor V Leiden
(rs6025) on GRCh37, whose reference base is the Leiden allele. ClinVar
writes the Leiden variant (c.1601G>A) there as T>T, and the normal-allele
record (`c.1601=`) as T>C. Dropping REF == ALT rows dropped the real
variant, and the normal-allele record's Pathogenic submissions matched
every C carrier. On GRCh38 (reference C) the same records are C>T and C>C,
and were handled correctly. From 2.4.2 (ClinVar cache `iv:5`):

- The loader reads ClinVar's HGVS name. A record whose nucleotide
  description ends in `=` describes the normal allele. Its REF == ALT rows
  are dropped (no earlier release reported them). Where its REF differs
  from its ALT it keeps every submission except Pathogenic and Likely
  pathogenic ones (every value of that family, including "low
  penetrance"): carrying the normal allele is not a pathogenic finding,
  but some other classifications of a normal allele are real claims. The
  CASP12 normal allele (rs497116, `c.373=`) is classified a risk factor,
  and is reported on GRCh38, where ClinVar writes it A>G. On GRCh37 the
  same record is written G>G and is dropped, so the two builds differ
  there; a rule for such records on both builds is not part of this
  change. REF == ALT rows that are not single-nucleotide variants (indels
  ClinVar cannot place on that build) are dropped too.
- Every REF == ALT row left in the cache therefore means that this build's
  reference carries the variant allele. It annotates a user who carries
  that allele, counting copies of it, like any other row. This is the one
  exception to invariant 2 (`alt` is never the reference base): the
  Annotation is marked `reference_carries_alt` and only such rows are.
- `reference_for` (the non-finding check of ClinPGx and SNPedia) returns,
  at an rsID with a REF == ALT row, the REF of the same ClinVar variant on
  the other build when that row is a single-base REF != ALT (C for
  rs6025). So CC, not TT, is the non-finding genotype on both builds.
  Without such a row (only one build's cache installed, or REF == ALT on
  both builds) it returns nothing and nothing is suppressed, as on GRCh37
  before this release. An `iv:4` cache keeps the locus REF until
  `db update` rebuilds it.
- ClinPGx and SNPedia derive their `alt` against the same normal allele,
  so CT and TT at rs6025 carry `alt` T on both builds. `locus_ref` stays
  the build's reference base, so enrichment uses the same key as the
  ClinVar row; where `alt` is that base the Annotation is marked
  `reference_carries_alt` as well.
- VCF input on GRCh37 has no variant record where a person is homozygous
  for the allele the reference carries (Factor V Leiden TT). **Amended
  2026-10-16 (2.5):** the pipeline watches ClinVar's reference-carried
  sites (one-base REF == ALT rows with an rsID, on the effective build)
  in the VCF stream. A `0/0` call or a covering gVCF reference block there
  that passes the GQ gate becomes a REF/REF call, and every source
  assesses it (`input.vcf_quality.reference_carried_calls`). A low-quality
  one stays a no-call. A plain VCF with no record there cannot tell a
  homozygous carrier from no coverage: those sites are counted
  (`reference_carried_not_in_file`) and the report says so. A gVCF
  position no record covers is already a no-call. Array input was never
  affected.
- **Amended 2026-10-16 (2.5), cache `iv:9`:** a normal-allele record (HGVS
  `c.N=`) written REF == ALT on a build whose reference is the normal
  allele is kept, no longer dropped. CASP12 rs497116 is G>G on GRCh37 and
  A>G on GRCh38; its risk-factor classification was shown on GRCh38 only.
  Such rows carry `normal_allele = 1`. They keep the locus REF in
  `reference_for`: they name the reference's own allele, not a variant it
  carries. They annotate carriers of their allele like any row, and their
  P/LP classifications are dropped at load, as on the other build. The
  site watch above includes them; the plain-VCF count does not, since
  missing a normal genotype there is not a missed variant.

## Amendment 2.6.0 (2026-11-02): the other build's REF is written on this build's strand

`reference_for` took the other build's REF of the same ClinVar variant as
written. It now takes it on this build's strand: that row's ALT, the same
variant allele, must be the carried allele (same strand) or its complement
(inverted, the REF is complemented); any other ALT, or a result equal to
the carried allele, gives nothing. With one build installed, the identity
table's other-build REF is translated the same way and abstains where the
strand relation cannot be told (ADR-0023, amendment of this date, which
also states the one-build rs6025 change).

## Amendment 2.6.0 (2026-11-03): the REF oracle is checked at the call's locus

ClinPGx and SNPedia no longer read the rsID-keyed `reference_for`
directly. Their non-finding filter reads
`ClinVarAnnotator.reference_at(rsid, build, locus)`, which gives
`reference_for`'s answer only where ClinVar places the rsID at the call's
own position on that build. An rsID ClinVar places elsewhere, or a call
with no position, suppresses nothing. The cross-build rsID identity table
(ADR-0040) is in place.

## Amendment 2.6.1 (2026-11-05): REF and normal allele from the call's own locus

The 2.6.0 amendment checked that ClinVar places the rsID at the call's
position, and then took `reference_for`'s rsID-wide answer, and
`normal_allele_for` ignored the locus. Where ClinVar holds one rsID at more
than one position, another placement could decide a call: an insertion's
one-base anchor or another row's REF hid a homozygous-ALT ClinPGx or SNPedia
finding (or showed a homozygous-reference one), and a reference-carried
site's normal allele elsewhere replaced an ordinary call's REF, so ClinPGx
and SNPedia derived an empty ALT. Two reference-carried positions with
different normal alleles pooled into no answer.

Both lookups now read only the rows at the call's own locus
(`reference_at`; `normal_allele_for` with a locus). The locus's single-base
REF and ALT rows count; an insertion anchor does not. Where the reference
carries the variant allele there (a REF == ALT row), the normal allele is
that locus's own VariationIDs on the other build, written on this build's
strand as before (`_normal_at`). A locus without ClinVar's row of the rsID
gives no REF and no normal allele. A caller that has no locus keeps the
rsID's single answer, as before. An outdated cache (before iv:9) keeps its
locus REF until `db update`, the 2.4.1 behaviour.

## Amendment 2.7.0 (2026-11-09): every placement of a VariationID

ClinVar's variant_summary can place one VariationID more than once on a
build: a pseudo-autosomal variant is listed on X and on Y. The loader's
temporary index was keyed on the VariationID alone and kept the first row,
so the other chromosome's placement had no row and a call there no match.
It is now keyed on (VariationID, chromosome, position): every distinct
placement is stored, each with its own copy of the variant's submissions
and the same per-submission submitter split; an identical repeated row is
kept once. Selection stays locus-keyed as above, so each placement is
matched, and its REF and normal allele are answered, at its own locus. The
cross-build normal-allele lookup reads the same VariationID's rows on the
other build, and abstains where its placements disagree. The cache's
record count counts the stored placement rows. CLINVAR_INTERPRETER_VERSION
is 10, so an iv:9 cache is rebuilt.

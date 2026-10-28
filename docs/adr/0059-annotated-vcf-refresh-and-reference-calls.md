# ADR-0059: Annotated VCF refresh and explicit reference calls

**Status:** Accepted
**Date:** 2026-10-24

## Context

`--vcf-out` appended its INFO to whatever the input carried. On an input
Allelix had already annotated, a value ClinVar had since withdrawn stayed
in the file, beside or instead of the current one. Lowercase REF/ALT missed
the pipeline's uppercase keys and got no values. And a ClinVar finding for
an allele the reference itself carries (ADR-0039) reached the JSON and HTML
reports but not the VCF, because the writer only looked at ALT alleles.

## Decision

Before writing, remove every `ALLELIX_*` INFO token from each data row, the
`##INFO=<ID=ALLELIX_...` declarations and the `##ALLELIX_...` header lines,
then write the current result. This applies to every row, including filtered
rows and reference blocks, which still get no new values. Other INFO tokens,
header lines and sample columns are written byte for byte.

Look up REF/ALT in upper case; write the columns as the input has them.

On an explicit reference call (`0/0`, `0|0` or haploid `0`), write the
pipeline's reference-carried ClinVar finding for that exact locus, REF and
ploidy in `ALLELIX_CLINVAR`. The field stays `Number=A`: the value repeats
in each non-symbolic ALT slot. Its description gains one sentence: "On a
homozygous-reference or haploid reference record, the finding is for the
allele the reference itself carries." Only ClinVar gets this; gnomAD,
AlphaMissense and CADD values for the reference allele are not written,
since they would read as properties of the ALT.

The parser records which explicit records gave a usable reference call at a
watched position: CHROM, POS, REF, ALT, FILTER, FORMAT and the sample columns
as written, with the positions each covers. The writer stamps exactly those
records, at the REF base covering the finding's position. It does not decide
eligibility itself, so the analysis's sample, quality and mitochondrial rules
apply unchanged: an excluded `chrM` record is never a source, and a `0/0`
record with a multi-base REF that covers the site is.

Only a run that writes `--vcf-out` keeps these records: one entry per
explicit reference record at a watched position, the size of that record's
columns. An analysis without `--vcf-out` keeps only the bounded base and
ploidy counts per watched position.

gVCF reference blocks are not stamped (one line covers many positions), and
a record with ALT `.` has no slot. Findings from either stay in the JSON and
HTML reports, and the docs say so.

`VCF_SCHEMA_VERSION` goes from 0.4.0 to 0.5.0. No INFO key, Number, Type or
value layout changes, but `ALLELIX_CLINVAR` on a reference record now
describes the allele the reference carries rather than the ALT, which is a
change in meaning for anyone reading the file.

## Amendment 2026-10-25 (2.5.2): a field of its own, on every record that carries REF

Two faults in the decision above:
- A het or half call (`0/1`, `0/.`) at a site where the reference carries a
  ClinVar allele had that finding in the JSON but not in the VCF. Only
  `0/0` records were stamped.
- The finding on a `0/0` record was repeated in every ALT slot of a
  `Number=A` field, where it read as a finding for the ALT.

The decision now:
- **A field of its own.** VCF schema 0.6.0 adds `ALLELIX_CLINVAR_REF`:
  `Number=1`, `Type=String`, value `significance|rsID|gene`, with the
  escaping and the report-rank choice of `ALLELIX_CLINVAR`. It holds the
  highest-ranked ClinVar finding for the allele the reference itself
  carries, on that record. `ALLELIX_CLINVAR` stays `Number=A` and carries
  only findings for the record's carried ALT alleles. A record can have
  both.
- **Every record whose call carries REF.** That is a het, a half call, a
  hom-ref or haploid reference call, and a single-position record with ALT
  `.` and no `END`. With `Number=1` the last needs no ALT slot, and its site
  is exact: an all-sites `0/0` at F5 on GRCh37 is a homozygous Leiden carrier,
  and a VCF reader must see it. A multi-base REF is stamped where it covers
  the site, as before.
- **Blocks excluded.** gVCF `END` reference blocks stay unstamped: a block
  spans many positions and has no single site.
- **Only the record the call was read from.** This applies to every field:
  ALT findings, enrichment and recovered rsIDs as well. The parser numbers
  data records, and each call keeps its number. The writer stamps a key's
  findings only on the records whose own call carries that allele. A
  duplicate at the same CHROM, POS, REF and ALT whose call is a no-call,
  low quality, hom-ref or another ALT borrows nothing. A recovered rsID goes
  on the record whose own call it was recovered for.
- **Re-annotation** removes stale values of both fields.
- **Half calls with ALT `.`.** A single-position record with ALT `.` and GT
  `0/.` is now read as a half call of REF, like the same GT on a record with
  an ALT. Before, it gave no call, and its finding was lost from the JSON
  too.

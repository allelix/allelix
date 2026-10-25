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

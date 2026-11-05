# ADR-0061: An rsID-less padded call is resolved in its trimmed form

**Status:** Accepted
**Date:** 2026-10-29

## Context

An rsID-less VCF record gets its rsID from position resolution: ClinVar,
the rsID identity table and gnomAD, in `annotators/_resolution.py`. Most
of their rows are in VCF minimal form, so a record written with padding
bases does not equal them. Examples are `GC>AC`, or `CG>CA` one base
earlier, for the site `G>A`. Through 2.5.2 such a record kept no rsID. Not
every row is minimal, though: the identity table keeps dbSNP's spelling,
and some rows are stored padded.

The 2.5.2 trimming for every source (ADR-0047's 2.5.2 amendment) runs
later, when the annotators read the call. That helped only a record that
carried its rsID. The sources that match by rsID (ClinPGx, SNPedia, the
GWAS Catalog) found nothing for the unnamed one; ClinVar, which matches by
coordinates and trims itself, still did. The same call gave four sources
named and one source unnamed.

## Decision

- **Trim when the written form does not resolve.** `resolve_batch` first
  looks up each rsID-less call as written, so a cache row stored in the
  padded spelling (dbSNP keeps some) still matches. A call left unnamed
  that has an input REF is then looked up in its trimmed form. Each
  carried non-reference allele is trimmed against REF by
  `trim_allele_pair`: the shared suffix, then the shared prefix, keeping
  an indel's anchor base, with no left-alignment. The trimmed form is used
  when every carried allele lands on one position and REF, and that
  differs from the record's.
- **Every resolver, every command.** The batch path and the streamed path
  (`resolve_rsidless`) both use it, so `analyze`, `stats`, `extract`
  (sequential and indexed) and `export plink` resolve the same call the
  same way.
- **The record keeps its spelling.** The trimmed form exists only for the
  lookup. The record keeps its own position, REF and alleles and takes
  only the recovered rsID (`original_rsid` keeps the old ID, for debugging
  only). Reports show the genotype as written, and `--vcf-out` writes the
  record's own columns, with the recovered rsID in ID and the old one in
  `ALLELIX_ORIGINAL_ID`, as for any recovered call.
- **Extract finds it at its trimmed site too.** An unnamed padded record
  whose call carries no ALT (a no-call, a REF-only half call, `0/0`, or a
  rejected call) cannot be named by resolution. `extract` counts it at
  every position its record stands at: its POS, and the position each of
  its ALTs trims to. That is also how `extract` decides whether a call is
  at a placement for its choice (ADR-0054), on the sequential and indexed
  paths alike. The call shown keeps its record's own POS, REF and alleles.
- **No left-alignment.** A repeat indel written at another anchor (`TT>T`
  one base on) is not moved; that needs the reference sequence.

## Consequences

An unnamed padded call gets the rows of its named or trimmed form, and so
do half, phased, haploid and two-ALT calls. A call mixing two frames, a
call whose REF did not come from the input (arrays), and a reference-only
call are not renamed. `extract` still finds a reference-only or rejected
padded call at its trimmed site. Annotating explicit reference calls is
not part of this decision.

`compare` reads both files as written and does not resolve rsIDs, so it is
not changed here. Matching a padded record against a trimmed one across two
files is left to a later change to `compare`.

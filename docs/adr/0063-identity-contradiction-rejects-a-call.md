# ADR-0063: A call whose locus contradicts its rsID's identity placement reaches no source

**Status:** Accepted
**Date:** 2026-10-29

## Context

A named call (one whose ID is an rsID) can sit at a position other than
the one its rsID has on the analysis build. A stale rsID, coordinates
converted from another build, or a faulty merge can cause this. Through
2.5.2, the sources treated it differently:
- ClinPGx matched by rsID and genotype alone (ADR-0009). ADR-0023 says a
  call outside ClinVar's placement is not suppressed.
- ClinVar, SNPedia and the GWAS Catalog rejected the call or abstained,
  each by its own rule.

So one call could get a pharmacogenomic finding and nothing else, with no
word in the report. The rsID identity table (ADR-0040) knows where each
rsID sits on GRCh37 and GRCh38.

## Decision

- **Contradiction is not absence.** When the identity table places the
  rsID on the analysis build and none of those placements is the call's
  locus, the call contradicts its identity.
  - Absence is never a contradiction: no identity cache, a build the table
    does not cover (GRCh36), an rsID the table does not hold, or an rsID
    with no placement on this build.
  - Any one of several placements is enough.
- **Only base calls against base placements** (supervisor ruling,
  2026-10-29). The check compares a single-base call (one base per allele,
  as written or trimmed) only with the rsID's single-base placements on the
  build (as stored or trimmed). An indel call, a chip D/I call and an
  identity indel row are identity unavailable, never a contradiction:
  their positions legitimately differ by representation and in repeats,
  and ADR-0061 forbids left-alignment. An rsID with only indel rows
  therefore leaves every call unchecked; an indel row beside a base row
  does not excuse a base call that stands at no base placement.
- **The call's locus** is its position as written or its trimmed position
  (a padded record, ADR-0047). A padded call is judged by its single-base
  form (`CG>CA` is `G>A`) whatever its genotype, heterozygous or
  homozygous, phased or not. A padded call that names no ALT (a
  reference-only `0/0` with a multi-base REF) has no single-base form and is
  not checked.
- **One check, before every source.** The pipeline checks each batch's
  named calls before annotation. A contradicting call is given to no
  source: ClinVar, ClinPGx, SNPedia and the GWAS Catalog alike, even when
  a source's own cache agrees with the call's position. Enrichment
  (gnomAD, AlphaMissense, CADD) works on findings, so the call gets none.
  No-calls are not checked, since no source reads them.
- **Counted once per call.** The report warning `identity_locus_conflict`
  counts each such call once, not once per source. It appears in JSON
  `report_warnings`, the terminal and the HTML warnings.

## Real-data counts

On 2.6 dev with V1, V2 and V7, against the same tree without V3:

| Input | Conflicts, any call | Conflicts, base only | Findings lost |
|---|---|---|---|
| 23andMe v5 (960,613 calls) | 504 (463 D/I, 458 of them within 10 bp) | 40 | 0 |
| 23andMe (960,614 calls) | 514 | 40 | 0 |
| AncestryDNA (701,478 calls) | 2 | 2 | 1 (GWAS) |
| HG002 GRCh37 benchmark | 17 (15 indels in repeats) | 2 | 0 |
| HG002 GRCh38 benchmark | 0 | 0 | 0 |
| 14 more array files from six vendor formats | — | 1 to 103 per file | 1 GWAS finding in each of 7 files, 0 in the others |
| HG002 child gVCF, HG002 GRCh38 chr22 | — | 0 among calls the file names (scan) | not run |

Every row except the last comes from a full analysis, so it counts calls named
by position resolution too. The last row is a scan of the calls as the file names
them, before resolution; calls it would name are not covered. HG00187 has no
detected build and is not checked.

The 8 findings lost across all files come from one rsID: a single-base
array call at a position that none of the rsID's single-base placements
gives. GRCh36 files and files without a detected build are not checked.

A first version that also checked indels rejected 6 findings from two
23andMe D/I calls at rsIDs with only indel identity rows. It also hid most
D/I calls at P/LP rsIDs from the "D/I calls not checked" disclosure (376 to
23 on one file). With base-only checking, the disclosure stays at 376.

## Consequences

ADR-0009's matching and ADR-0023's "not suppressed outside ClinVar's
placement" now apply only to calls their rsID does not contradict. The
call still counts in the file's totals; only its annotations are
withheld. A call without an rsID that position resolution names
(ADR-0061) passes through the same check after resolution, and can
contradict too: ClinVar and gnomAD are resolution providers separate from
the identity table, so a ClinVar row outside the rsID's identity placement
names a call that the identity table then contradicts.

## Amendment 2.6.0 (2026-11-02): a record takes only its own call's findings

A rejected call produced no finding, and its record gets no Allelix finding
field in the annotated VCF either. Before, the writer could stamp it with a
finding of another record at the same CHROM, POS, REF and ALT, because it
was still listed among the records carrying that allele (ADR-0059,
amendment 2026-11-02). Its recovered rsID, if any, is handled as before.

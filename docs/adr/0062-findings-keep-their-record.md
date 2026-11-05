# ADR-0062: A finding keeps its own record's ploidy and record

**Status:** Accepted
**Date:** 2026-10-29

## Context

ADR-0050 shows a call the file declares haploid as hemizygous. Through
2.5.2 the pipeline marked findings haploid after annotating a whole batch
(`_mark_haploid`). It matched each finding to a haploid call by chromosome,
position and REF, by locus, or by rsID; none of these names the record or
its ALT. So a haploid `G>T` record and a separate diploid `G>A` record at
the same position made the diploid call's findings hemizygous, while their
genotype stayed `AG`. A chip file listing one position twice, once as one
letter and once as two, did the same. Nothing in the report said that two
records disagreed.

## Decision

- **Findings are traced to their call.** A batch is split, in file order,
  into runs in which no two calls share a locus (as written, or trimmed,
  ADR-0047) or an rsID. A batch without such duplicates stays one run. Each
  run is annotated on its own, so the existing match (record key, locus,
  rsID) finds exactly one call per finding.
- **Each finding takes that call's ploidy.** It is hemizygous only when its
  own call is haploid. `Annotation.record_index` names the VCF record the
  call was read from (`Variant.record_index`, zero-based in file order;
  None for other inputs). It is internal: not in JSON, HTML or `--diff`.
  A padded call's finding that carries only its trimmed locus (a GWAS row
  without a risk allele) finds it by that locus. A reference call rebuilt
  at a ClinVar reference-carried site from an explicit `0/0` or haploid
  `0` record keeps the number of the first record whose base the site
  accepts (a record with another REF there gives its number to no call);
  one rebuilt from a gVCF block has none.
- **Contradictions are counted.** When the calls at one position disagree
  on ploidy (a haploid call and a diploid call), every call at that
  position is counted. The report says so in `report_warnings` with the
  code `contradictory_calls` and the count of calls, and in the terminal
  and HTML warnings. Identical duplicates are not counted, and no-calls
  take no part. Every call's position and ploidy is kept (one integer per
  call; a position outside the compact range, such as a negative POS,
  is kept as a plain integer) and compared when the stream ends, so file order and batching do
  not change the count. A padded call is at its trimmed site (ADR-0061).

## Consequences

The genotype, zygosity and plausibility of a finding come from one record
however the file is batched, ordered or named, for VCF and chip input. Only
calls at the same locus or rsID are annotated in separate runs, so the
order of findings in a batch changes only where such duplicates occur. Two
diploid records that disagree (for example `0/1` and `1/1` for the same
ALT) are not counted here; each keeps its own genotype.

## Amendment 2.6.0 (2026-10-29): explicit reference records (ADR-0064)

`contradictory_calls` also counts explicit reference records (`0/0`, `0|0`,
haploid `0`) that pass the quality gate, but only at a target site of an
enabled source (ADR-0064: a ClinPGx, SNPedia or GWAS Catalog placement,
a single-base record naming one of their rsIDs, or a ClinVar
reference-carried site): there a haploid `0` beside a diploid `0/0` or
`1/1` is a contradiction. Non-reference records count everywhere, as
above. A reference record elsewhere is not counted and makes no other
record contradictory. The count is of records, not positions. A reference record is at the site its calls stand at: a padded
record whose ALTs all trim to one base at one position is there (`CG>CA`
at P-1 is at P, as `G>A` at P is), and an indel record, ALT `.` or a
single-base REF is at its written position. Reference blocks are not counted. A reference
call rebuilt at a ClinVar reference-carried site from such a record is
counted once, with the record; one rebuilt from a block is counted as
before.

## Amendment 2.6.0 (2026-11-02): a record takes only its own call's findings

`Annotation.record_index` is now also how the annotated VCF picks a
record's findings: an annotation is filed under its own record's
`(chrom, pos, ref, alt)`, and a record is stamped only with the annotations whose
`record_index` is its own (ADR-0059, amendment 2026-11-02). It stays
internal: not in JSON, HTML or `--diff`.

## Amendment 2.6.0 (2026-11-03): positions that are not counted

A call at POS 0 or below is not counted in `contradictory_calls`: its
position is unknown (an unplaced rsID-list row, ADR-0069) and is shared
with no other call. So a negative POS is not kept as a plain integer. The
plain-integer fallback is only for a position past the 64-bit range of
the compact array.

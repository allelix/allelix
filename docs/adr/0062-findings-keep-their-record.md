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

## Amendment 2.7.0 (2026-11-09): one site's split records are one genotype

A two-ALT call split over two records (`G>A 1/0` and `G>T 0/1`) gave each
record's own genotype (`AG` and `GT`) for an `A/T` person, so SNPedia and
ClinPGx rows for `AG` and `GT` were reported and the `AT` row missed. The
owner ruled (D4) that only the exact complementary pair is combined, and
anything else abstains, counted.

The analysis pipeline now reads the consecutive single-base records at one
CHROM and POS together (`VcfParser.combine_split_sites`, set only by the
analysis; `compare`, `extract`, `stats` and `export plink` are unchanged):
records that all name one ALT (repeats, or disagreeing calls of one allele,
counted as above) and records that disagree on REF give their calls as
before; among records of one REF that name different ALTs, exactly two
records with the same ID and REF, one different ALT each, both called,
with the two ALTs on different haplotypes (`1/0` and `0/1` unphased, or
`1|0` and `0|1` in one phase set: equal `PS`, or none on either, the
default phase set; SUP's ruling within D4), give one call with the two
ALTs, the merged record's call. Each allele keeps its own record
(`Variant.allele_records`): a finding is traced to, and stamped in
`--vcf-out` on, the record of its own allele, and the call's
`record_index` is the first of the two. Any other mix gives no call, its
records counted in `report_warnings` as `split_site_calls`: both ALTs on
one haplotype, differing or unrelatable phase sets, different IDs, half
calls, no-calls and records the quality gate rejects (which still count
as part of the site, though they give no call). A site is the records
that follow one another at one position as written: any record at
another position ends it. Indel, padded and symbolic records and
reference records are not grouped.

## Amendment 2.7.1 (2026-11-13): split records give their alleles unless they contradict

The 2.7.0 rule withheld every mix that was not the exact pair, so a carried
Pathogenic allele on a passing record went unreported where 2.6.3 reported
it, and the 1000 Genomes files, which write one positional ID per allele
(`22:11924846:G:C`), never paired at all: their raw IDs differed. The owner
amended D4 (2.7.1 decision D1; the outside review's Addendum 1, section 2).

Among records of one REF that name different ALTs:

- The exact pair is two records, both called, one different single-base
  ALT each, on different haplotypes: `1/0` and `0/1`, or `1/.` and `./1`
  (as `vt decompose` writes `1/2`), unphased, or the same phased in one
  phase set. The records' IDs are compared as canonical rsIDs: only two
  different rsIDs block the pair; a missing, positional or other ID never
  does. The call takes the one rsID, or none, and is then resolved by
  position. It is one call for every source, as before.
- Records that contradict each other give no call to any source: more
  than two ALT copies across the passing records (a homozygote with a
  heterozygote, a third ALT), a haploid call beside another call, or an
  ALT on the same haplotype of one phase set twice. The site is named in
  `high_value_unchecked` (reason `split_site`) at a high-value rsID, or
  where ClinVar holds a Pathogenic or Likely pathogenic row at its
  position for an ALT its records carry, whatever rsIDs they name (SUP's
  ruling); it is named under its records' rsIDs or, with none, under that
  row's rsID. A P/LP row for an ALT the records do not carry does not
  count.
- Any other mix agrees but is not one genotype (`0/1` and `0/1`
  unphased, phase sets that differ, phased with unphased, a half call
  beside a heterozygote, two different rsIDs, a passing call beside a
  rejected or `./.` record). Each passing record gives each ALT it
  carries as a half call `ALT/.` (ADR-0047), given only to the sources
  that read one allele (ClinVar and the GWAS Catalog, where they run);
  ClinPGx and SNPedia, which need the genotype, get nothing. A rejected
  record gives nothing and cancels nothing.

`split_site_calls` still counts every record of every site that is not the
exact pair. Half calls are counted after the decision: JSON schema 15
(decision D5) counts the split-site alleles given as `ALT/.` in
`input.vcf_quality.split_site_half_calls`, and `half_calls` keeps counting
standalone half calls only. In `--vcf-out`, both records of a merged pair
carry its rsID and their own `ALLELIX_ORIGINAL_ID` where their ID columns
did not both name it, so the annotated file pairs again when it is read.

## Amendment 2.7.1 RC2 (2026-11-14): two rsIDs found by position are no pair

An exact pair whose records name no rsID is one call resolved by position.
Where that lookup finds the two ALTs under two different rsIDs, the call
had no rsID at all (resolution needs one rsID to explain every carried
allele), so the sources that look a call up by rsID never saw it; 2.6.3
gave the GWAS Catalog finding (the outside review's A-01). SUP ruled, for
2.7.1, that two rsIDs found by position are point 1's two different rsIDs:
the pair is replaced by one half call per record, `ALT/.` under that ALT's
own rsID, the same result as the records with those rsIDs written. ClinVar
and the GWAS Catalog assess each allele; ClinPGx and SNPedia get nothing.
Both records count in `split_site_calls`, both alleles in
`split_site_half_calls`, and in `--vcf-out` each record carries its own
rsID, so the annotated file reads back the same. The cost: ClinVar shows
two half calls at such a site instead of the person's two-ALT genotype.
Keeping that genotype while binding each ALT to its own rsID, for split
pairs and merged `1/2` records alike, is left for 2.8.0.

## Amendment 2.7.1 RC3 (2026-11-14): the pair's rsID is decided per ALT

The RC2 amendment split a pair only when both of its ALTs were found by
position, under two different rsIDs. Two siblings still lost findings 2.6.3
gave (the outside reviews' A-01b, B-01 and B-02): a pair of which only one
ALT is known (resolution needs one rsID to explain every carried allele, so
the call kept no rsID and the GWAS Catalog never saw it), and a pair of
which one record names its rsID while the other ALT belongs to another
rsID found by position (the written rsID was taken for the call and the
other ALT was never looked up). SUP ruled, for 2.7.1, one rule for every
accepted exact pair: each record keeps the canonical rsID it names, and a
record naming none has its ALT looked up alone, by the usual ClinVar,
identity-table and gnomAD routes and their build limits. Two different
rsIDs: one half call per record, `ALT/.` under its own rsID, as the RC2
amendment says. Exactly one: the full two-ALT call takes it, as when one
record names it (point 1). None: the pair is unchanged. Each record's own
written ID is kept for this, not the merged call's. Sites whose records
only agree or contradict each other are not pairs and are not affected.
The annotated VCF reads back the same in each case.

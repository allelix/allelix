# ADR-0064: An explicit VCF reference call reaches the sources that read a reference genotype

**Status:** Accepted
**Date:** 2026-10-29

## Context

A chip file writes every genotype, a homozygous-reference one included.
ClinPGx, SNPedia and the GWAS Catalog read that genotype like any other:
for example, a GWAS association whose effect allele is the reference base,
or a ClinPGx row for `GG`. A VCF writes the same call as `0/0`, `0|0` or a
haploid `0`. Through 2.5.2 the VCF reader yielded no call for it. Only
ClinVar's reference-carried sites rebuilt one (ADR-0039, ADR-0059). So the
same person got a row from a chip file and nothing from a VCF. A half call
`0/.` was annotated while the more complete `0/0` was not.

## Decision

- **Which sites.** The sites are the placements, on the analysis build, of
  every rsID that ClinPGx, SNPedia or the GWAS Catalog holds
  (`Annotator.reference_call_rsids`). They come from the rsID identity
  table, each REF/ALT pair trimmed (ADR-0047, no left-alignment): a
  single-base change at its base, any other variant (an insertion, a
  deletion, an MNV) at its first base. A single-base record that names
  one of those rsIDs is a site too, at its own position, with or without
  the identity table. They are watched for a VCF input only.
- **Which records.** An explicit record whose selected call is homozygous
  reference or haploid reference, and that passes the quality gate,
  counts. The gate is FILTER, the sample's FT (ADR-0060), GQ and DP 0;
  `--keep-filtered` applies as usual.
  - ALT may be a base or `.`.
  - A padded record counts at each watched position its REF covers, with
    its REF base there.
  - gVCF reference blocks (`END`) are not expanded.
  - A no-call, a low-quality call or a filtered call gives nothing.
  - A position that another record changes gives nothing. A record
    changes the positions of each ALT its genotype carries, and of every
    ALT when an allele is unknown (`./.`, `0/.`). A record the quality
    gate rejects is a no-call there, as for every other check: a
    filtered, low-GQ or DP 0 record changes the positions of every ALT,
    whatever its genotype (`0/0`, `0/1` or `1/1`), since any of them is
    possible (`--keep-filtered` and `--min-gq` apply as usual; DP 0
    always rejects). Shared padding is no change (ADR-0047, ADR-0061): an
    equal-length ALT changes the bases that differ (`CGT>CGA` changes
    only the third), a deletion the bases it removes, and a symbolic or
    `*` ALT all of REF. An insertion changes its anchor base: the base
    survives, but the allele there is not the reference allele, so a
    `G>GT` `0/0` split from a carried `G>GA` gives no `G/G`. A split
    multi-allelic site (`0/0` for one ALT, `1/1` for another) is the
    common case. The record order does not matter.
  - So does a site the file does not mention: absence is never read as
    reference.
- **Evidence for the target** (supervisor ruling, 2026-10-29). An
  explicit `0/0` is evidence only for the variant its record describes.
  - A single-base change (`G>A`) counts from any record at the position
    with that REF base, as before.
  - Any other target counts only from a record with one ALT whose trimmed
    REF/ALT pair is the target's, padded on either side or not; no
    left-alignment, so a repeat-equivalent insertion one base away is not
    evidence. An SNV-only record at an insertion's anchor, an unrelated
    insertion or deletion, ALT `.`, or another ALT of a split record says
    nothing about it.
  - A matching record gives the full REF/REF call at the target's trimmed
    position (`GC/GC` for `GC>G`), haploid when the record is; the sources
    read it as they read any call. No call is made where another record
    changes any position of the target's trimmed REF.
  - Indel targets are watched like any other; none is excluded wholesale.
- **Which rsID.** A single-base record naming an rsID stands for that
  rsID alone, at its own position, as the same chip call does, unless the
  table places that rsID as another kind of variant. Identity then decides
  as for every call (ADR-0063): a placement elsewhere on the build is a
  counted `identity_locus_conflict`; no cache, an rsID the table does not
  hold or does not place leaves the call standing. A padded record naming
  an rsID stands for it only where the table places it and the record is
  evidence for it. An unnamed record stands for every rsID placed at its
  position that it is evidence for. Each rsID and position is called once, from the
  first record, and not again where ClinVar's reference-carried path
  already made the call.
- **Which sources.** These calls go to ClinPGx, SNPedia and the GWAS
  Catalog, as REF/REF calls (haploid when the record is). ClinVar keeps its
  own reference-carried path. The findings are in the JSON, HTML and
  terminal reports. `--vcf-out` does not stamp them: the ADR-0059 rule for
  reference records is unchanged.
- **Shared checks.** The calls pass the identity check (ADR-0063) like
  every call. They are annotated in runs that share no locus or rsID, and
  each finding keeps its own record's ploidy and number (ADR-0062).
- **Ploidy** (supervisor ruling, 2026-10-30; ADR-0062 amendment). A
  reference record is kept only at a target site, for annotation and for
  the ploidy count alike. `contradictory_calls` is the number of records
  in scope at positions where records in scope disagree on ploidy. In
  scope are:
  - every non-reference record (one carrying an ALT, or a half call),
    anywhere, as V7 counts them;
  - an explicit reference record (`0/0`, `0|0`, haploid `0`) that passes
    the quality gate, only at a target site of an enabled source: a
    position a ClinPGx, SNPedia or GWAS Catalog placement watches, a
    single-base record naming one of their rsIDs, or a ClinVar
    reference-carried site.

  A reference record anywhere else is not counted and makes no other
  record contradictory: it changes no finding. `0` and `0/0`, or `0` and
  `1/1`, at a target site are two contradictory records; three records
  of mixed ploidy are three. A record's site is where a call of it
  stands, as for every call (ADR-0061): a padded record whose ALTs all
  trim to one base at one position is there (`CG>CA` at P-1 is at P);
  an indel record, ALT `.` or a single-base REF is at its written
  position. Reference blocks do not count as reference records. A reference call
  ClinVar's path rebuilds from such a record is counted once, with the
  record; one rebuilt from a block is counted as before (ADR-0062).
- **Memory.** The targets are held as sorted position and rsID-number
  arrays while the file streams. A target's alleles are read from the
  identity table only where a reference record was kept, and nothing is
  kept for a reference record at any other site.

## Consequences

The same call gives the same rows from a chip file and from a VCF. An
all-sites or multi-sample VCF gains rows wherever a reference genotype has
a finding: on one 1000 Genomes chr22 sample, the confident `0/0` records
gave about 900 GWAS and 20 ClinPGx rows, the same rows the chip form of the
same calls gets. The default display floors still apply. gVCF reference blocks
stay a separate question.

## Amendment 2.7.0 (2026-11-09): a symbolic deletion's END span

A symbolic ALT changed only the bases of its REF, its padding base, so a
carried `<DEL>` with `END` beyond a `0/0` SNV left the SNV's reference call
standing, where a literal deletion spanning it removed it. A symbolic
deletion (`<DEL>` or `<DEL:...>`) the selected sample's genotype carries
now also changes POS + 1 through its INFO `END`, and so does one it may
carry: an allele is unknown, or the quality gate rejects the call (FILTER
or FT unless `--keep-filtered`, GQ below `--min-gq`, DP 0), whatever its
genotype, as above. The decision reads the selected sample's own genotype
and gate (`_deletion_span_for`): a passing call of another ALT of the
record, symbolic or `*`, says nothing about the deletion's bases. Those spans are
kept as intervals, not one entry per base, so a long deletion costs no
more than a short one. The same span counts for ClinVar's
reference-carried calls (ADR-0039), through the parser's covering-record
watch: there, past its REF, a deletion covers a position only for a
sample that carries it or may carry it, by the same decision. Other symbolic ALTs (`<DUP>`, `<INV>`, `<CNV>`,
`<NON_REF>`), a deletion the genotype confidently does not carry, and an
`END` that is missing, unreadable or not after POS change nothing more than
before.

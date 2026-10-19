# ADR-0047: VCF calls that were dropped silently are matched or counted

**Status:** Accepted
**Date:** 2026-10-14
**Builds on:** ADR-0007 (carrier rule), ADR-0011 (indel-anchor protection),
ADR-0039 (locus-keyed ClinVar selection), ADR-0041 (plausibility flag),
ADR-0046 (REF read from the input).

## Context

Four kinds of VCF call could lose a ClinVar finding with no warning:

- **Untrimmed alleles.** A caller or a merge can write an allele with bases
  it shares with REF: a deletion `GTT>GT,G`, an SNV `CA>TA`. ClinVar writes
  the trimmed form (`GT>G`, `C>T`). The match needs REF and ALT to be
  equal, so these never matched. A compound call (`ATG>ATGTG,A` 2/1) needs
  each ALT trimmed on its own: one trim of REF and both alleles together
  leaves `ATGTG` as it is.
- **Half calls.** GT `./1` was a no-call: the called allele was lost.
- **Homozygous deletions to one base** (`AC>A` 1/1). The ADR-0011 rule, a
  ClinVar indel never matches a call with no multi-base allele, is meant for
  array calls, which have no REF. It also fired here. A heterozygous call
  kept the multi-base REF allele and matched.
- **Spanning deletions** (`*`). A call with `*` on a haplotype was counted
  only when the `*` was on one haplotype, and its position was not kept.

Other rows the parser dropped (malformed, no genotype, a genotype index
outside the ALT list, a genotype naming a symbolic ALT such as `<NON_REF>`
next to a real one, a genotype naming an ALT on a record whose ALT is `.`)
were not counted at all.

## Decision

- **Trimming, in the ClinVar annotator, per carried allele.** For a call
  whose REF came from the input (`Variant.ref_from_input`), each carried ALT
  is trimmed against REF: the bases they share at the end, then at the
  start, while both keep one base, so an indel keeps its anchor. A dropped
  leading base moves the position on by one. The row equal to the trimmed
  `(position, REF, ALT)` matches, looked up at the trimmed position.
  - The annotation carries ClinVar's position and alleles, so gnomAD,
    AlphaMissense and CADD are looked up by ClinVar's key.
  - The genotype is shown as the file wrote it (`ATGTG/A`).
  - `--vcf-out` stamps the record and ALT as written.
  - Calls with a REF from a reference table (arrays) are never trimmed.
  - The run counts the calls that matched only after trimming.
- **No left-alignment.** It needs the reference sequence. An indel the
  file places at a different position from ClinVar's normalized one (in a
  repeat) still does not match.
- **The anchor rule applies to arrays only.** A call whose REF came
  from the input skips the ADR-0011 multi-base check: REF equality already
  pins the variant. Array calls keep it, also when a REF was filled in from
  the identity table.
- **Half calls.** A VCF call with one allele called (`./1`, `1/.`, also
  `./0`) is yielded with the called allele.
  - ClinVar and the GWAS Catalog, which need only a carried allele, assess
    it. ClinPGx and SNPedia, which need both alleles, do not.
  - The genotype is shown as `T/.` and the zygosity as "Half call", never
    "Homozygous" or "Heterozygous".
  - The plausibility flag (ADR-0041) is not assessed: the copy count is
    unknown.
  - The position-based rsID resolvers (ClinVar, identity table, gnomAD)
    include it.
  - A half call below the quality floor is a no-call, like any other call.
  - An array call with one allele missing stays a no-call.
- **Dropped calls are counted, with their reason**, in `input.vcf_quality`
  (JSON) and in the input-quality lines on the console and in HTML:
  - `spanning_deletion_calls` counts a `*` on one or both haplotypes, and
    `spanning_deletion_sites` lists each one as `chrom:pos REF>ALT GT`;
  - `half_calls`, `malformed_rows`, `no_genotype_rows`,
    `unusable_genotype_calls` (a genotype index outside the ALT list, or an
    ALT named on a record whose ALT is `.`), `symbolic_allele_calls` (a
    genotype naming a symbolic ALT such as `<NON_REF>` on a record that
    also has a real ALT), and `trimmed_matches` (informational);
  - with the existing filtered, low-quality and symbolic counts.
  - A call written as a no-call (`./.`) is reported as a no-call, as
    before. It is not a dropped call and is not in these counts.
  - A gVCF reference block (`<NON_REF>` or `<*>` as its only ALT) whose
    genotype names that ALT is not counted. Checking every block's genotype
    made parsing a 19.2-million-block gVCF measurably slower (about 10% to
    30% more CPU time over repeated runs on a 2-vCPU machine), and no such
    call was found in the two real gVCFs checked.
- No ClinVar cache change and no interpreter-version bump. The new JSON
  fields belong to schema 9, the one schema bump of 2.5.

### Amendment 2026-10-16 (2.5)

- A multi-allelic or untrimmed record can match rows of rsIDs other than
  the one in its ID column. When the record's ID is itself one of ClinVar's
  rsIDs among those rows, each matched row now carries its own rsID; before,
  every row was labelled with the record's ID. A record whose ID ClinVar
  does not hold there (a merged or alias rsID) keeps its own ID, as before.
- A call from an indel record keeps a slash in reports. A homozygous
  one-base deletion `AC>A` 1/1 is `A/A`, not `AA`, which read like an SNV
  genotype. This applies to ClinVar, ClinPGx and GWAS alike
  (`models.diploid_label`).
- `--vcf-out` does not stamp a match at a reference-carried row (REF ==
  ALT, ADR-0039). The matched allele is the record's REF, and the stamp is
  per ALT allele. The report and JSON carry the row.

## Consequences

- VCF users get ClinVar findings for these calls. Evidence for the HG002
  GRCh37 benchmark VCF is in the CHANGELOG.
- Array users see no change.
- The report says how many calls were not assessed and why.
- Indels written at a different position from ClinVar's remain unmatched
  until Allelix can left-align against a reference sequence.
- A call with one haplotype deleted (`*`) still has no Variant shape; that
  stays with the single-copy representation planned for 3.0.

# ADR-0051: Count the array calls ClinVar cannot match, and panel coverage in four states

**Status:** Accepted
**Date:** 2026-10-16
**Supersedes:** the premise of ADR-0011 that arrays do not call indels (its
decision stands). **Builds on:** ADR-0039 (locus-keyed ClinVar selection),
ADR-0045 (ClinVar's aggregate), ADR-0050 (`is_pathogenic_aggregate`).

## Context

- **Two kinds of array call are never matched against ClinVar, and the
  report did not say so.**
  - A vendor probe without an rsID (23andMe `i5000123`) is not
    locus-matched (ADR-0039): nothing confirms that the probe tests the
    variant ClinVar records at that position.
  - An array writes an insertion or deletion as `D` or `I`. Those are not
    bases, so they never equal a ClinVar ALT.
  - A user could read "no pathogenic findings" as covering those sites. In
    one public 23andMe v5 file (openSNP), about 6,000 i-probes sit on loci
    with ClinVar Pathogenic / Likely pathogenic rows, and 970 calls are D/I.
- **ADR-0011's premise is false.** It says arrays do not call indels;
  23andMe does, as D/I codes. Its decision (a ClinVar indel row never
  matches a single-base array call) is still right: a single base cannot
  confirm an indel, and a D/I code cannot say which indel.
- **Array input had no parse counts.** VCF input has `input.vcf_quality`;
  an array file had nowhere to report what its parse saw.
- **Panel coverage folded no-calls into "found".** `panel_coverage` was
  `found` / `missing` / `no_findings`, so a no-call at a requested rsID
  read as "genotyped, no findings".

## Decision

- **Annotator run counts.** `Annotator.run_counts()` returns what the
  source counted during the run (empty by default). ClinVar returns
  `unrecognised_rows` (ADR-0045), `trimmed_matches` (ADR-0047),
  `i_probes_at_plp_loci` and `di_calls_at_plp_rsids`. The pipeline reads
  these in place of the annotator's attributes.
- **Probe classes.** P/LP is ClinVar's aggregate classification
  (`is_pathogenic_aggregate`, ADR-0050). Counted by call shape, whatever
  the vendor:
  - `i_probes_at_plp_loci`: called probes whose ID starts with `i`, at a
    locus where ClinVar has a row with a P/LP aggregate. A no-call is not
    counted.
  - `di_calls_at_plp_rsids`: calls with a `D` or `I` allele whose rsID has
    a ClinVar row with a P/LP aggregate, wherever ClinVar places it.
  - JSON `input.probe_classes`; null for VCF input and when ClinVar did
    not run on the input's build.
- **Array input stats.** `ArrayInputStats` in `parsers/base.py`, the array
  counterpart of `VcfQualityStats`, in JSON `input.array_quality` (null
  for VCF input):
  - `di_calls`: every call with a `D` or `I` allele, counted by the
    pipeline;
  - `unreadable_genotypes`, `not_genotyped`, `wrong_column_rows`: for
    parsers that can tell (a genotype field no rule reads, a reference
    call a converter filled in, a row with the wrong column count). A
    parser keeps them through `GenotypeParser.input_stats()`; the others
    leave them at 0.
- **Disclosure.** The console and the HTML input-quality notice say how
  many i-probes sit at P/LP loci and how many D/I calls there are, with
  the number at P/LP rsIDs.
- **Panel coverage.** Each requested rsID is in exactly one state:
  - `annotated`: at least one shown annotation (after the display filter).
    It wins over a no-call, so a half call with a finding is annotated.
  - `typed_no_annotation`: in the input and called, nothing shown.
  - `no_call`: in the input with either allele missing.
  - `not_on_chip`: not in the input. For VCF input it is shown as "not in
    file": the position may be reference, or have no variant call.

  JSON `panel_coverage` is `{requested, counts, variants}`. `counts` has
  the four states and sums to `requested`; each `variants` entry is
  `{rsid, gene, state, genotype, annotation_count}`, sorted by state, then
  rsID. A VCF half call's genotype is shown as the annotations show it
  (`T/.`). HTML lists every rsID with its state; the terminal report
  prints the four counts and flags the no-call and not-on-chip rsIDs.

## Consequences

- Schema 9 gains `input.array_quality` and `input.probe_classes`, and
  `panel_coverage` changes shape (`found`, `missing` and `no_findings` are
  gone). The CHANGELOG lists each.
- Nothing is matched that was not matched before; the counts only say what
  was not assessed.
- Matching i-probes or D/I calls would need a vendor probe map, which
  Allelix does not have; that is a separate decision.

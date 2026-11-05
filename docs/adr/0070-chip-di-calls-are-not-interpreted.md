# ADR-0070: A chip D/I call is interpreted by no source

**Status:** Accepted
**Date:** 2026-11-01
**Builds on:** ADR-0051 (the array calls ClinVar cannot match), ADR-0024 (GWAS magnitude scoring), ADR-0069 (rsID and genotype lists)

## Context

An array writes an insertion or deletion call as `D` or `I`. A `D` or `I`
code says that some deletion or insertion was called. It does not say
which one. Allelix has no per-probe map that would tell it. ADR-0051
counts these calls. Its 2.5.2 amendment names them at high-value and P/LP
rsIDs as not checked, "neither a finding nor a normal result". The
owner's decision is that a chip D/I call is not interpreted.

ClinPGx already kept to that: it reads only a genotype of two of A, C, G
and T. Three sources still read the call's letters:

- **ClinVar** matches by locus, REF and ALT, so a pure D/I call (`DD`,
  `DI`, `II`) never matched. A mixed call, one base and one code (`AD`,
  `DA`, `AI`, `IA`), carried its base as an ALT and got that ALT's ClinVar
  row.
- **SNPedia** compares a call to its page genotypes by equality. The public
  cache holds 1,076 genotype rows written in `D` and `I`, at 1,017 IDs (924
  rs, 93 i). So a `DI` call matched a page's `(D;I)` row, as written,
  whatever the probe measured.
- **The GWAS Catalog**: a row without a risk allele fires on the rsID
  alone, for any call that is not a no-call (ADR-0024), so a D/I call got
  it. The cache holds 225,951 such rows at 162,238 rsIDs.

An rsID and genotype list's D/I row is unplaced (`di_call`, ADR-0069) but
still reached SNPedia and the GWAS Catalog.

## Decision

**A call with a `D` or `I` allele reaches no source.** `Variant.is_di_call`
is the one predicate; it is true for a mixed call. ClinVar's row check,
SNPedia's candidate check and the GWAS Catalog's single and batch paths
skip such a call before matching it to a finding, from an array file or
an rsID list. ClinVar counts it at a P/LP rsID first, so the count and
the naming below are unchanged.

The call stays where ADR-0051 puts it:

- counted in `array_quality.di_calls`;
- named in `high_value_unchecked` and the `di_calls_not_checked` warning at
  a high-value or P/LP rsID, in JSON, HTML and the terminal.

A sequence call is not a D/I call. A VCF's alleles are bases, so an indel
call such as `C/CT` still matches a SNPedia page written in those
sequences, as before.

## Consequences

- A D/I call no longer shows ClinVar (a mixed call), SNPedia or GWAS
  findings. The
  `di_calls_not_checked` message, "cannot check against its sources", is
  now true of every source.
- A D/I call is no longer counted in SNPedia's or GWAS's
  `strand_abstentions`, since it never reaches their strand check. On a
  file with D/I calls those counts fall.
- A GWAS row with a risk allele was never matched by a D/I call: the loader
  keeps only A, C, G and T risk alleles. Only rows without one change.
- Interpreting D/I calls would still need a vendor probe map. That is a
  separate decision.

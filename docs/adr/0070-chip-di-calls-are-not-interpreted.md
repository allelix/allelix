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

## Amendment 2.6.1 (2026-11-05): 23andMe's own i-probe pages

Owner's decision: SNPedia reads its page for a 23andMe `i` ID
again, for a D/I call in a file the 23andMe parser read. The `i` ID names
23andMe's own probe, and SNPedia's page for it is written for that probe,
so the page says which indel the `D` or `I` names: the question this ADR
is about. In 2.6.0 a carrier call at 23andMe's CFTR probe (`i4000292`,
`DI`) lost its SNPedia finding. The pipeline turns this on per run
(`_read_vendor_probe_pages`, `SNPediaAnnotator.read_vendor_probe_pages`)
only when the parser is the 23andMe parser itself, not a subclass: Genes
for Good writes 23andMe-style files, but its `i` IDs are not 23andMe's
probes.

Everywhere else this ADR stands: a D/I call at an rsID, a D/I call from
any other parser, and every other source (ClinVar, the GWAS Catalog,
ClinPGx). The call is still counted in `array_quality.di_calls`. At a
locus where ClinVar classifies a variant as Pathogenic or Likely
pathogenic, or near one (an indel whose REF span covers the call or ends
up to five bases before it), the call is still named, because ClinVar
does not match a probe to that variant (ADR-0051).

## Amendment 2.6.2 (2026-11-07): what is named where, and the warning's count

The Decision's "named in `high_value_unchecked` and the
`di_calls_not_checked` warning ... in JSON, HTML and the terminal" holds
for JSON. Since 2.6.2 the terminal and HTML name a D/I call only at a
high-value rsID and count the rest, one line per kind, with the full list
behind `--show-unchecked` and in a collapsed HTML section (ADR-0051).
The warning's count is the number of D/I calls listed in JSON
`high_value_unchecked` (reason `di_call`), so it grew in
2.6.1 with the vendor-probe and positional naming; the count of D/I
calls at P/LP rsIDs (`di_calls_at_plp_rsids`) did not change (ADR-0051).

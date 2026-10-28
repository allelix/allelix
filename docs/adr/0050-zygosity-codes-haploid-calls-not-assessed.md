# ADR-0050: Zygosity codes, file-declared haploid calls, a named "not assessed", and the multi-homozygote warning

**Status:** Accepted
**Date:** 2026-10-15
**Amends:** ADR-0041 (plausibility). **Builds on:** ADR-0047 (half calls),
ADR-0049 (`report_warnings`).

## Context

- **Haploid calls read as two copies.** A male X call written as one
  letter (23andMe), an MT or Y call in FamFinder's one-allele shape, and a
  VCF one-index GT (`1`) were doubled at parse time and shown as
  Homozygous, so the reports said "two copies" where the file said one.
- **"Not assessed" had no value of its own.** Plausibility was an empty
  string when it was not assessed, whatever the reason: no gnomAD
  frequency, gnomAD not run, a chromosome without a rule, a half call, or
  an annotation with no carried allele. HTML and the terminal showed
  nothing, so "not assessed" and "fine" looked the same.
- **No warning for many rare homozygous P/LP calls.** A file with several
  homozygous calls of rare Pathogenic alleles is usually a genotyping or
  file problem, and nothing said so.
- **Display strings in JSON.** JSON carried display strings ("Homozygous",
  "No Call", "Half call"); 2.5 is the schema 9 window to make them codes.

## Decision

- **Zygosity codes.** JSON `zygosity` is one of `homozygous`,
  `heterozygous`, `hemizygous`, `half_call` or `no_call`. HTML and the
  terminal show `Homozygous`, `Heterozygous`, `Hemizygous`, `Half call`
  and `No call`. The missing-allele glyph is `.` (`T/.`).
- **File-declared haploid only.** `Variant.haploid` is set where the file
  says one copy:
  - a one-letter array call (`split_genotype` / `is_haploid_call`, used by
    23andMe (and Genes for Good, which reuses its reader), deCODEme,
    LivingDNA, MyHeritage, FTDNA and FTDNA Illumina);
  - FamFinder's one-allele MT/Y line;
  - a VCF GT with one index.

  A no-call is never haploid. Matching still sees the called allele twice,
  so no finding changes. After each batch, the pipeline marks the
  annotations of a haploid call (`Annotation.haploid`, internal) and shows
  the one allele (`T`), with zygosity `hemizygous`. Sex is not inferred
  and there is no `--sex` option. A male X call written as two letters
  (AncestryDNA, FTDNA) stays Homozygous.
- **Plausibility "not assessed" is a value with a reason.** `plausibility`
  is `plausible`, `suspicious`, `implausible` or `not_assessed`. Every
  annotation gets one. `plausibility_reason`, for `not_assessed` only, is:
  - `gnomad_not_run`: gnomAD did not run;
  - `half_call`: one allele is missing, so the copy count is unknown;
  - `no_allele`: the annotation names no carried allele to look up;
  - `chromosome`: no rule for the copy count on that chromosome (Y, MT, a
    homozygous X call);
  - `no_frequency`: gnomAD has no frequency for the allele.

  A hemizygous call is one copy, so ADR-0041's heterozygous rule applies
  to it on X. `not_assessed` ranks below `plausible` when several
  annotations share a call. HTML shows a muted mark with the reason, and a
  detail-panel row. The terminal shows `○`, with a legend line.
- **The multi-homozygote warning.** `run_analysis` counts, before any
  display filter (so `--min-magnitude` and `--category` cannot hide it),
  distinct autosomal variants (1 to 22) that are:
  - called homozygous (a hemizygous call is not);
  - from a ClinVar row whose aggregate is Pathogenic or Likely pathogenic.
    Every `/`-separated part of the aggregate's germline term (before
    `;`) is Pathogenic, Likely pathogenic, or one of those with "low
    penetrance". "Pathogenic/Likely pathogenic; risk factor" counts;
    "Pathogenic/Likely risk allele" and Conflicting do not;
  - rare in gnomAD (frequency below `homozygous_max_af`) or without a
    gnomAD frequency.

  At `MULTI_HOMOZYGOUS_PLP_THRESHOLD` (3) or more, `report_warnings` gets
  `multiple_homozygous_plp` with the count, and the report is still
  written.

## Consequences

- The JSON values change within schema 9: `zygosity` codes,
  `plausibility` `not_assessed` in place of `""`, and the new
  `plausibility_reason`. The CHANGELOG lists each.
- Reading a two-letter male X call as hemizygous needs sex inference or a
  declared sex: a separate decision (ADR-0041), and not in 2.5.
- `stats` het/hom counts and the PLINK export are unchanged for haploid
  calls.

## Amendment 2026-10-25 (2.5.2): only assessed calls count as rare

`multiple_homozygous_plp` counts a homozygous P/LP call when gnomAD assessed
it as rare: either a frequency below the threshold, or an allele looked up in
gnomAD and absent (`no_frequency`). Before, any missing frequency counted.
A call gnomAD was not asked about no longer counts: gnomAD not run, no
identity table, or a GRCh37 translation that abstained (ADR-0040's 2.5.2
amendment). On GRCh37 input such calls are common, and nothing says they are
rare. The plausibility reasons are unchanged.

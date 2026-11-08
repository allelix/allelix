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

### Amendment 2026-10-29 (2.6.0): one record per finding

A finding is marked haploid only by the call it came from, never by
another record at the same position; calls of different ploidy at one
position are counted in the report warning `contradictory_calls`
(ADR-0062).

## Amendment 2.6.0 (2026-11-02): SNPedia sequence genotypes

The zygosity codes are unchanged. `zygosity_code` reads a slashed label by
its two alleles and an unslashed one by its bases, which is right only for
two one-base alleles. Every source now writes a sequence allele with a slash
(ADR-0047, amendment 2026-11-02), so a SNPedia finding on a homozygous
multi-base call (`GAC/GAC`) is `homozygous`; before, its run-together label
made it `heterozygous`.

## Amendment 2.6.0 (2026-11-03): a one-letter call in an rsID list

`Variant.haploid` is also set for a one-letter call in an rsID and
genotype list (`allelix/parsers/rsidlist.py`, ADR-0069): `A` is a haploid
call of A. A no-call is still never haploid.

## Amendment 2.6.2 (2026-11-07): plausibility is shown only when it warns

The HTML report and the terminal show plausibility only when it warns. A
call that is `implausible` or `suspicious` keeps its mark (HTML `⚠` with
its hover text and its detail-panel row; terminal `⚠` or `?`). A call that
was `not_assessed` gets no mark and no detail-panel row, and the terminal's
`○` legend line is gone. The muted `○` was barely visible and explained
itself only on hover, which a phone cannot do; the panel row appeared on
some findings and not others, which read as random.

The warning is explained once, in plain words, in the HTML "Reading This
Report" section, in the terminal's legend line (shown when a call is
flagged) and in the README: the call looks doubtful because the variant is
rarer in the population (gnomAD) than Allelix's limit for that call. The
wording names no input type, since VCF and chip calls are assessed alike,
and it gives the default limits as defaults: a configured
`heterozygous_max_af` above 0 flags alleles gnomAD has seen. A call this rare
can be a genotyping or sequencing error; confirm it before relying on it.
Allelix can only check variants with population counts in gnomAD, so no
warning is not a guarantee. The suspicious hover text and detail-panel row
say the same: at or below the limit for one copy, by default no carriers.
Classification is unchanged.

The JSON report is unchanged, with no schema change: every annotation keeps
`plausibility`, and a `not_assessed` one keeps `plausibility_reason`. The
HTML page's embedded data no longer carries the reason's display text,
which nothing shows. The multi-homozygote report warning is unchanged.

## Amendment 2.6.2 (2026-11-07): at or under the limit for one copy

The 2.6.2 amendment above words the explanation as "rarer in the population
(gnomAD) than Allelix's limit for that call". That is exact for two copies
only. The rules are unchanged: a homozygous call is implausible when its
allele's gnomAD frequency is strictly below `homozygous_max_af`, and a
heterozygous call is suspicious when the frequency is at or below
`heterozygous_max_af` (with the default 0, every suspicious call is at the
limit). The HTML "Reading This Report" section, the terminal's legend line
and the README now say the call looks doubtful because the variant's
frequency in gnomAD is at or under Allelix's limit for that call (under it
for two copies). Classification is unchanged.

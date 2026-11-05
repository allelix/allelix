# ADR-0060: The selected sample's genotype filter (FORMAT/FT) is a filter

**Status:** Accepted
**Date:** 2026-10-29

## Context

A VCF can reject a call in two places. The site's FILTER column rejects the
record for every sample. The sample's FORMAT/FT field rejects one sample's
genotype: it holds `PASS` or the semicolon-separated codes of the filters
that call failed (VCF 4.3, 1.6.2). Through 2.5.2 the reader honoured FILTER,
GQ and DP 0, but never read FT. A record with FILTER `PASS` and FT `LowAD`
was annotated as passing by every source, and `filtered_rows` stayed 0.

## Decision

- **One gate.** The selected sample's FT is part of the same gate as the
  site's FILTER (`_parse_record`). An FT that is not `PASS` and not `.`
  (or empty or absent) makes the call filtered.
- **Counting.** A filtered record is counted once in `filtered_rows`,
  whether FILTER, FT or both reject it. `filtered_by_label` counts it under
  each site label and each FT code. A code in both is counted once.
- **Kept or skipped.** By default a filtered call is never annotated. At a
  watched position it is a no-call with the reason `filtered`, labelled
  "Caller FILTER or sample FT rejected". `--keep-filtered` keeps calls
  rejected by FILTER and by FT alike, and they are still counted.
- **Only the selected sample.** Another sample's FT in a multi-sample file is
  not consulted.

## Consequences

The gate is shared, so every reader of `_parse_record` follows it:
- the analysis
- the panel and high-value observation
- reference-call watching
- sequential and indexed `extract` (a rejected record is a no-call)
- `compare`, `stats` and `export plink`, which read the same parser

`--vcf-out` stamps a record only with findings from that record's own call
(ADR-0059, `record_rows`). So a filtered call gets nothing unless it was
kept. `ALLELIX_FILTER` still repeats the site's FILTER only; a sample's FT
stays visible in its own column. The JSON schema does not change: the
existing counters cover it.

## Amendment 2.6.0 (2026-11-02): a record takes only its own call's findings

The Consequences above say `--vcf-out` stamps a record only with findings
from that record's own call (`record_rows`). `record_rows` gave eligibility
only: a filtered call was not stamped, but a passing duplicate with another
genotype could take a finding its own call did not produce. From this
amendment the writer also requires the finding's `record_index` to be the
record's own (ADR-0059, amendment 2026-11-02). The FT behaviour itself does
not change.

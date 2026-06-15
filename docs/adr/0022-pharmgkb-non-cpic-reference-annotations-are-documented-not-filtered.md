# ADR-0022: PharmGKB reference-genotype annotations on non-CPIC genes are documented, not filtered

- **Date:** 2026-05-14
- **Status:** Accepted (scope reduced by ADR-0023)

## Context

ADR-0020 makes the PharmGKB non-finding filter a structured join: each `(rsid, base)` is looked up in a per-allele function table sourced from CPIC's API. If every base in the user's genotype maps to `Normal function`, the row is suppressed.

CPIC publishes function classifications only for the gene-drug pairs it has curated guidelines for. For genes outside that scope — MTHFR, F2, F5, and others — CPIC has no allele function entries. The lookup returns nothing for those rsids, and the filter abstains. The result is that **reference-genotype rows on non-CPIC PharmGKB annotations emit as findings.**

Three such rows appeared in real-data testing:

- `rs1801133` (MTHFR) `GG` with methotrexate
- `rs1799963` (F2) `GG` with oral contraceptives
- `rs6025` (F5) `CC` with oral contraceptives

In each case the user is homozygous reference. The PharmGKB annotation row itself is real — PharmGKB publishes one row per `(annotation, genotype)` including the reference homozygote. CPIC's structured filter cannot suppress it because there's no per-allele function data for those genes.

The temptation is to add a hand-curated exclusion list ("always suppress MTHFR GG, F2 GG, F5 CC, ...") or a regex against "decreased risk" phrasing. Both would recreate the failure mode ADR-0016 and ADR-0020 were written to end: hand-maintained code chasing upstream data, growing without bound, silently breaking when phrasings change or when CPIC adds (or removes) a gene.

## Decision

**Allelix does NOT filter PharmGKB reference-genotype rows on non-CPIC genes. The rows emit. The README documents the limitation and tells the user what the leaked rows look like so they can recognize and discount them in reports.**

The README's "Supported Databases" section gains a paragraph:

> A small number of PharmGKB rows on non-CPIC genes (MTHFR, F2, F5, etc.) appear in reports even when the user is homozygous reference. PharmGKB publishes one annotation per genotype including the reference homozygote, and CPIC — Allelix's structured per-allele function source — has no curation data for these genes. The structural filter cannot suppress them. These rows are identifiable by a homozygous-reference genotype combined with "decreased risk" or "may have a typical response" language. They are an upstream data gap, not an Allelix bug; we surface them honestly rather than hide them behind a curated exclusion list.

## Why not a curated exclusion list

- **Exclusion lists are the same maintenance trap as regex.** They chase upstream data with manual code edits. CPIC adds a gene → the list is stale. PharmGKB removes a row → the list still suppresses something nonexistent. The user pays for our maintenance lapse.
- **Exclusion lists silently hide signal.** A row we suppressed today because we judged it noise might be the most important row in tomorrow's report. A user can read prose and decide; static suppression cannot.
- **Documentation is honest.** "This row is here because CPIC doesn't cover the gene; you can ignore it" is more useful than "we hid this row for you," because the next person doesn't know what was hidden or why.
- **The PharmGKB filter's correctness story stays clean.** ADR-0020's join is the entire filter. No carve-outs, no exceptions, no special cases. When CPIC adds a gene, the filter starts working for it automatically with zero code change.

## Consequences

- **The reference-homozygous rows on MTHFR / F2 / F5 stay in the output.** They are correctly emitted given the data sources Allelix has.
- **The README explicitly names the limitation** and gives the user the recognition pattern (homozygous reference + comparative-risk language).
- **No exclusion list, no regex against PharmGKB phrasing, no per-gene carve-out shipped in code.** Future PRs adding any of these are rejected unless they also rewrite this ADR.
- **The path to broader coverage is upstream:** push CPIC to curate more genes, or integrate an additional structured source that publishes per-allele function for MTHFR / F2 / F5. Either improves the filter without code maintenance.
- **If a future contributor proposes "just a small filter for MTHFR" six months from now**, this ADR is the answer: we considered it, and we deliberately chose not to. The maintenance cost is the failure mode, not the leaked rows.
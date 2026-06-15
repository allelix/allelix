# ADR-0020: CPIC API is the structured per-allele function source

- **Date:** 2026-05-13
- **Status:** Accepted (supersedes ADR-0017 and ADR-0018; superseded as primary filter by ADR-0023)

## Context

Earlier PharmGKB filter approaches (prose regex, structured column, hybrid patterns, template extraction) each failed because they tried to derive per-allele function from within PharmGKB's own data rather than from the authoritative external source.

The architectural directive: **the PharmGKB filter is not a text classification problem. It is a table join.** Each allele at each position has a discrete function value. Look up both alleles in the table; if all map to `Normal function`, suppress.

CPIC (Clinical Pharmacogenetics Implementation Consortium) publishes this table at https://api.cpicpgx.org/v1/. The join across three CPIC tables — `sequence_location`, `allele_location_value`, `allele` — yields `(rsid, base) → clinicalfunctionalstatus`. CPIC's enumeration:

| Value | Meaning |
|---|---|
| `Normal function` | reference / wild-type metabolizer |
| `Decreased function` | reduced enzyme activity |
| `No function` | enzyme inactive |
| `Increased function` | hyperactive metabolizer |
| `Uncertain function` | classification pending |
| (gene-specific tags, e.g. `Malignant Hyperthermia associated`) | per-gene clinical labels |

The three v0.7.0 leakers all have explicit CPIC entries:

| rsid | C | T |
|---|---|---|
| `rs1800559` (CACNA1S) | Normal function (Reference) | Malignant Hyperthermia associated |
| `rs116855232` (NUDT15) | Normal function (*1) | No function (*3) |

Real CPIC data, two API calls away from being correct.

## Decision

**The `pharmgkb_allele_function` table is populated from CPIC's API. The filter is `is_nonfinding_by_allele_lookup(rsid, genotype, table)`. No prose parsing anywhere in production code.**

### The filter

```python
def is_nonfinding_by_allele_lookup(rsid, genotype, lookup):
    if rsid not in {k[0] for k in lookup}:
        return None  # CPIC has no opinion → caller emits row
    for base in set(genotype):
        if lookup.get((rsid, base)) != FUNCTION_CLASS_NORMAL:
            return False  # at least one non-Normal allele → finding → emit
    return True  # all Normal → non-finding → suppress
```

### The source

`fetch_cpic_allele_functions()` in `allelix.databases.cpic_loader` performs three HTTP GETs against CPIC's PostgREST API:

| CPIC table | Fields used |
|---|---|
| `sequence_location` | `id`, `dbsnpid` |
| `allele_location_value` | `alleledefinitionid`, `locationid`, `variantallele` |
| `allele` | `definitionid`, `clinicalfunctionalstatus` |

Joined client-side, filtered to single-base `variantallele` values (A/C/G/T), normalized to the project's `FUNCTION_CLASS_*` enum. Statuses outside the recognized set are skipped — never silently coerced to `Normal` (the failure mode v0.5–v0.7 kept producing).

### Coverage

CPIC publishes function classifications for the genes it has guidelines for: CYP2D6, CYP2C19, CYP2C9, NUDT15, TPMT, DPYD, SLCO1B1, CACNA1S, RYR1, IFNL3, HLA-A, HLA-B, and a few others. Approximately 1,400 `(rsid, base)` entries.

For PharmGKB rsids in non-CPIC genes (MTHFR, COMT, etc.), CPIC has no data; the lookup returns no entry; the filter abstains; the row emits. Acceptable: the alternative (silently suppressing rows we have no structured evidence about) is what produced four releases of leakers.

### Schema impact

None. `pharmgkb_allele_function` keeps its v0.7.0 schema. Only the data source changes; the `source` column now reads `cpic_api` instead of `cpic_template`. Existing v0.7.0 caches are schema-current; users must run `db update --force` to re-ingest with the CPIC API as the source.

## Consequences

- **v0.7.0 production leakers fixed.** Verified against the real PharmGKB July-2025 dump joined against the live CPIC API: all three leakers (`rs1800559 CC` ×7 anesthetics, `rs116855232 CC` ×3 thiopurines) classify as `is_nonfinding=1`. The v0.6.1 DPYD cluster (`rs115232898 TT`, `rs1801266 GG`, `rs3918290 CC`) also classifies correctly via the same path.
- **Zero over-suppression of known carriers.** `rs1801133 AG/AA`, `rs1799853 CT`, `rs4149056 CT`, `rs4244285 AG` all classified `is_nonfinding=0`.
- **ADR-0017 (prose fallback) and ADR-0018 (CPIC template extraction) are superseded.** Both relied on parsing description text. Both are gone from production code.
- **`is_nonfinding_for_row` is now a two-step pipeline:** structured `Allele Function` column (rarely populated on SNVs, but authoritative when it is), then the CPIC API lookup. If neither has data, the row emits. The filter never silently suppresses without structured evidence.
- **Mock fixture no longer ships rows whose sole purpose was exercising regex paths.** PA-010 (somatic-prose synthetic) and PA-011 (CPIC-template synthetic) deleted. PA-008 retained, simplified to exercise the structured-lookup path.
- **The annotator's `setup()` adds one HTTP call** to fetch the CPIC lookup. Small (~5 MB of JSON). Failure is loud — the loader is never invoked, so the existing cache stays intact.
- **The data classification principle (ADR-0016) is finally enforced end-to-end.** Filter input = structured table. Description text serves only its real purpose: display to the user in the report.
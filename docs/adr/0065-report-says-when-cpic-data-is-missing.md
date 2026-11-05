# ADR-0065: A report made without CPIC data says so

**Status:** Accepted
**Date:** 2026-10-30

## Context

ClinPGx rows for alleles that CPIC rates as normal function are not
findings, so they are hidden (ADR-0020). This needs CPIC's allele-function
data, fetched at `db update` time. When CPIC's API is down or empty, the
ClinPGx cache is still built, without that data, and its stamp records
`cpic:unavailable`. Those normal-function rows are then shown.

`db status` and the line printed after `db update` said so. A report made
without an update (`--no-update`, or any run that skips it) did not: not in
the console, the HTML or the JSON. A reader could not tell that the filter
was off.

## Decision

- **When:** a report says that CPIC data is missing when ClinPGx ran on a
  cache with no CPIC allele-function data. That is when either holds:
  - the stamp says `cpic:unavailable`, including a 2.5.0 stamp that claims
    CPIC data while the table is empty;
  - the CPIC table is empty, with or without a stamp.
- **Where:** the notice is a `report_warnings` entry with code
  `cpic_unavailable`. It shows in the console, as a `Warning:` line, and in
  the HTML and JSON reports, like the other report warnings.
- **What it says:** CPIC data is missing, ClinPGx rows for normal-function
  alleles may be shown, and `allelix db update` tries CPIC again.
- **No notice** for a healthy cache, or when ClinPGx did not run (switched
  off in the config, or not installed).
- **No other change:** the check is read-only: it reads the cache's stamp and checks whether
  the CPIC table has any row.
  Analysis stays offline and read-only. Findings are unchanged: the notice
  reclassifies nothing. `count` is null. No JSON field is added, so the
  schema does not change.

## Consequences

A reader can see when ClinPGx findings were made without CPIC's
normal-function filter. `db status` keeps its existing check, which reads
the cache signal; for a legacy stamp that claims CPIC, that signal check
also looks for CPIC rows. It marks "(no CPIC)" when the signal says CPIC
is unavailable, and this ADR does not change it.

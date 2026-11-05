# ADR-0067: A switched-off source is fully off

**Status:** Accepted
**Date:** 2026-10-30

## Context

A user can switch a source off: in the config (`sources.X = false`, or a license the user has not accepted), with a flag (`--no-gnomad`, `--no-alphamissense`, `--no-cadd`), or through a focused report's exclusions. Until 2.6 such a source was left out of the annotation loop, but it stayed wired into the others:
- ClinVar's reference allele and gnomAD's major allele still decided the normal allele for ClinPGx and SNPedia.
- The identity table still gave placements and strand readings.
- `stats`, `extract` and `export plink` still resolved calls through all three.
- Its cache was still prepared and refreshed before the run.

The report said "skipped" while the source changed results. On a real file (user1190), switching off ClinVar, gnomAD or the identity table gave different reports depending only on whether that source's cache was present.

## Decision

Owner decision D3 (2026-10-29), option A: a switched-off source is fully off. It is never annotated and never consulted as a dependency. Dependents take their existing unavailable path, and say so.

- **Choices come first.** `analyze` and the focused reports compute the user's choices before any source is prepared, refreshed or wired. A switched-off source is not prepared (`prepare_cache`), not refreshed, and not passed to another source.
- **Stand-ins equal absence.** Every other source is wired to a stand-in (`switched_off_stand_in`) for a switched-off ClinVar, gnomAD or identity table.
  - The stand-in is bound to `SWITCHED_OFF_DIR`, a path under the null device where no cache can exist.
  - Its connection answers as an absent cache's does, without opening anything.
  - Dependents therefore run the same code as when the cache is not installed. Switched off with the cache present gives the same report, byte for byte, as switched off with the cache absent.
- **Other commands.** `stats`, `extract` (sequential and indexed) and `export plink` read ClinVar, gnomAD and the identity table through `configured_source`. A source switched off in the config is a stand-in there too.
- **Disclosure.** When ClinVar, gnomAD or the identity table is switched off, the report carries a `report_warnings` entry `source_dependency_unavailable`, in the console, HTML and JSON. It names the source and what other sources ran without: the normal allele for ClinPGx and SNPedia, strand readings, rsID resolution. Provenance keeps the source as skipped with its reason.
- **Startup.** In `analyze` and the focused reports, the automatic update step leaves a switched-off source's cache files untouched: no WAL conversion, no preparation. A switched-off cache also never decides whether the run can start, so the outcome is the same with that cache installed or absent.
- **With every source on** (the default), nothing changes: the reviewer's default digests are byte-identical to dev before V6.

`db update` is cache management, not annotation. Its explicit flags (`--only`, `--no-gnomad`, `--no-alphamissense`, `--cadd`) are unchanged. Its WAL conversion still covers every cache in the data directory (ADR-0048); only the automatic step of a report command leaves a switched-off source's files alone.

## Consequences

A switched-off source can no longer change a report through another source. The cost is that results which used it are made without it, and the report now says so. ADR-0058's allowance for a skipped ClinVar provider to be used through another source's callback is reversed (see the amendment there).

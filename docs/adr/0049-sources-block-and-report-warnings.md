# ADR-0049: Every report names its sources; report-level warnings

**Status:** Accepted
**Date:** 2026-10-15
**Builds on:** ADR-0048 (one readiness state per cache). **Amends:** ADR-0036
(the annotated VCF's provenance block gains one line).

## Context

A report listed only the sources that ran (`annotators`). A source that
did not run left no trace beyond its absence. After an upgrade that
invalidates the ClinVar cache, `analyze` skipped ClinVar with one console
line and wrote a report that read as complete: empty summary cards, and a
`clinvar` block of zeros that looked like "ClinVar ran and found nothing".
A user's own choices (a source turned off in config, a `--no-` flag, the
license ladder) were not recorded either, and a source turned off and
never installed was reported as "not installed".

With only one ClinVar build installed (`db update --build grch37`),
`analyze` skipped ClinVar for every input. ClinVar's readiness needed
every build the annotator managed, and analyze managed both.

## Decision

- **One `sources` list per run**, built once by the CLI after the analysis,
  from the cache states (ADR-0048) and the user's choices. The renderers
  only read it. Each registered source, in registry order, has
  `{name, display_name, state, reason, version}`:
  - `used`: it ran, and it is in `annotators`, which keeps that meaning;
  - `ready_not_used`: it was ready and given to the run, but nothing
    queried it (for example gnomAD when no annotation needed a frequency);
  - `skipped`, with a reason. A failure is the cache's state: `missing`,
    `outdated`, `unreadable`, `corrupt`, `wal`, `readonly_wal`. A user's
    choice is `disabled` (config), `license` (the permission ladder),
    `excluded` (`--exclude-*`), `flag` (`--no-gnomad` and the others) or
    `not_applicable_build` (a source that has nothing for the input's
    build).
  - A choice wins over a failure: a source the user turned off is
    reported as turned off, whether or not its cache is there, and it
    gets no "Skipping" console line.
- **`report_warnings`**: one list of `{code, message, count}` for
  report-level warnings. JSON carries it as is; the HTML shows each as a
  `notice-warn` at the top, beside the build warnings; the console prints
  each. The first code is `clinvar_skipped`: ClinVar did not run, for any
  reason, so the report has no ClinVar findings. The report is still
  written.
- **HTML.** The failures form one notice, each with the console's reason
  and fix. The user's choices form one quiet line. The note on
  source-specific magnitude thresholds names only sources that ran. The
  subtitle already lists only the sources that ran, and the regulatory
  notice is unchanged.
  *Amended (2.5):* in the HTML, a short disclaimer and the warnings that
  change what to do stay at the top. The full regulatory notice, sources
  off by choice, input quality and the ClinVar and strand counts sit in a
  collapsed "Technical details" section. JSON and the console are
  unchanged.
- **Annotated VCF.** `##ALLELIX_SkippedSources=` lists each skipped
  source with its reason (`clinvar:outdated,cadd:disabled`), or `none`.
  `VCF_SCHEMA_VERSION` goes from 0.3.0 to 0.4.0 under ADR-0036's bump rule.
- **ClinVar with one build installed.** Analyze manages the ClinVar
  builds whose caches are ready, or both when neither is. With both ready
  nothing changes: the other build's cache still gives the normal allele
  for reference-carried variants. Once the input's build is decided, a
  ClinVar that does not manage it is skipped with that build's cache state
  (`Annotator.serves_build`), named on the console with the input's build,
  and left out of `annotators`.

## Consequences

- A report without ClinVar says so in every format, so it can no longer
  be read as complete.
- `sources` is the place for later per-source states; `report_warnings`
  is the place for later report-level warnings (the multi-homozygote
  warning is next).
- The JSON schema stays 9: both blocks are 2.5 additions under the same
  bump.

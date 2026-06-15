# ADR-0033: Schema Version Bump Policy

**Status:** Accepted
**Date:** 2026-06-11
**Supersedes:** The schema-version carve-out in ADR-0032 (Context section,
"JSON schema version stays at 3").

## Context

ADR-0032 stated that additive optional fields (specifically `cadd_phred`)
do not require a schema version bump. v1.6.1 adds `zygosity`, also an
additive optional field, and bumps the schema to version 4. The two
decisions are inconsistent.

## Decision

Every new field emitted in JSON report output bumps the schema version,
even if the field is optional and additive. This gives consumers a
reliable signal that the output shape changed.

The `cadd_phred` addition in v1.6.0 should have bumped the schema from
2 to 3 under this policy. It did not (ADR-0032 carved it out). That
ship has sailed — v1.6.0 emits schema version 3 with `cadd_phred`
present, and changing it retroactively would break existing consumers.
Going forward, all new fields bump.

## Consequences

- Any commit that adds, removes, or renames a field in the JSON report
  output must increment `SCHEMA_VERSION` and update
  `_SUPPORTED_SCHEMA_VERSIONS` in the diff loader.
- The diff loader's supported-versions set grows monotonically.
- ADR-0032's body is not edited (ADR immutability). This ADR supersedes
  its schema-version stance.

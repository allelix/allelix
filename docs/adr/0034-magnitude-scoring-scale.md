# ADR-0034: Magnitude Scoring Scale and Ceiling

**Status:** Accepted
**Date:** 2026-06-11

## Context

Allelix uses a 0-10 magnitude scale inspired by SNPedia to rank variant
significance for report sorting. Each annotator independently assigns a
score within this range. The composite score shown in reports is the
maximum across all sources.

No source's scoring logic produces 10. Per-source caps (verified
against code): ClinVar caps at 9.0 (Pathogenic), PharmGKB at 9.0
(Level 1A), GWAS at 9.0 (explicit `min(..., 9.0)` cap), SNPedia
passes through wiki-assigned magnitudes (no cap in code; nothing
scores 10 in the wiki data in practice).

## Decision

The 0-10 scale is retained. The practical ceiling of 9 is intentional,
not a bug. No upstream database asserts absolute certainty about any
variant, and the scoring system reflects that. 10 is reserved headroom
that no current evidence tier reaches.

Max-across-sources is the composite rule. This means the least
conservative source's score dominates. This is a known tradeoff accepted
for v1.x. The failure mode (score inflation from one generous source)
is mitigated by source attribution in reports — every score is labeled
with its origin.

### Future direction (v2.0.0)

Per-source scores will be surfaced alongside the composite in all report
outputs. The composite max rule stays, but the individual contributions
become visible so users can judge the evidence basis themselves. This
separates the three dimensions the current composite conflates: clinical
significance, evidence quality, and actionability.

## Consequences

- 10 is valid but unreachable by design. Do not add scoring logic that
  produces 10 without a new ADR.
- Report renderers continue to show a single composite score in v1.x.
- JSON schema does not change in v1.x. The per-source breakdown is a
  v2.0.0 schema change.

# ADR-0038: GWAS Magnitude Modifiers Require Odds-Ratio Evidence

**Status:** Accepted
**Date:** 2026-07-31
**Supersedes:** ADR-0024's effect-size modifier. Its p-value tiers, magnitude
ceiling, unknown-risk-allele cap, filtering, and rollup decisions remain in
force.

## Context

The GWAS Catalog exposes odds ratios and beta coefficients through the same
numeric field. ADR-0024 applied odds-ratio thresholds to either value. That is
unsound for beta coefficients: their magnitude depends on the phenotype's
unit, scale, and transformation. A beta of 2.0 centimetres and a beta of 2.0
standard deviations do not represent comparable effects.

The Catalog's confidence-interval text supplies a conservative structured
signal. Odds ratios use a bare positive numeric interval such as
`[0.95-1.05]`; beta rows generally carry units or direction text. Missing,
`NR`, descriptive, or unit-bearing text cannot establish that the numeric
effect is an odds ratio.

## Decision

The existing `+0.5` and `+1.0` effect-size modifiers apply only when the
confidence-interval text positively identifies a bare, positive ratio
interval. The numeric value must also be present and positive.

All ambiguous cases abstain from the modifier, including:

- missing or `NR` confidence-interval text;
- intervals with units, direction, or other descriptive text;
- negative interval bounds; and
- beta coefficients of any magnitude.

The p-value-derived base magnitude is still emitted. Allelix will not add
beta-specific thresholds until the data model carries enough unit and
phenotype-scale context to compare beta coefficients honestly.

## Consequences

- Quantitative-trait betas no longer receive arbitrary odds-ratio bonuses.
- Rows with insufficient effect-type context may lose a modifier that a true
  odds ratio deserved; conservative under-scoring is preferred to a false
  high-magnitude claim.
- A rare unitless beta with a bare interval remains a known blind spot because
  the current Catalog fields cannot distinguish it from an odds ratio.
- No cache-schema bump is required: confidence-interval text was already
  stored. Both direct and batch annotation queries now consume it.
- Tests pin the accepted interval grammar and the ambiguous cases that must
  abstain.

# ADR-0037: Build Detection Requires a Three-Seed, 80% Majority

**Status:** Accepted
**Date:** 2026-07-31
**Supersedes:** ADR-0021's unspecified confidence threshold. ADR-0021's
central decision that position evidence outranks file headers remains in
force.

## Context

ADR-0021 established position-based genome-build detection because provider
headers can be wrong, but described confidence only as "all (or all but a
noise tolerance)" of the seed matches. The implementation required exact
unanimity and at least three matches. One stale or repositioned seed could
therefore discard a strong position majority and fall back to the header.

The correctness review proposed two relaxations: accept a sufficiently strong majority,
and allow a unanimous sample smaller than three seeds to override the header.
The first repairs the missing noise tolerance. The second is unsafe: two
non-reference WGS calls can coincidentally be the only seed positions observed,
so a 2/2 result is too small a sample to override another build signal.

## Decision

A detected build is confident only when both conditions hold:

1. At least three inspected seed positions match the winning build.
2. The winner accounts for at least 80% of all inspected seeds.

The boundaries are intentional:

- 4/5 is confident and may override the header.
- 3/4 is not confident because 75% is below the ratio floor.
- 2/2 is not confident because it is below the match-count floor.

Detection below either threshold remains diagnostic evidence, not the
effective build. The pipeline reports the tentative position result and uses
the existing fallback order: explicit override, header, chromosome-prefix
inference, then GRCh37. ADR-0025's GRCh36 fail-safe remains unchanged.

## Consequences

- One noisy seed no longer defeats a well-supported 4/5 position result.
- A two-seed panel cannot silently override a header or chromosome-prefix
  signal.
- CLI diagnostics distinguish a confident detected build from a tentative
  position majority and report the build actually used.
- Tests pin 4/5, 3/4, and 2/2 as separate confidence boundaries.

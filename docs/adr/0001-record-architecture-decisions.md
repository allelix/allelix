# ADR-0001: Record architecture decisions

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Allelix has non-obvious design decisions that span privacy, regulatory posture, licensing, and the plugin architecture. Without an in-tree record, future contributors (and future-us) will see the code without the reasoning and either re-litigate settled questions or quietly violate constraints whose purpose isn't visible.

## Decision

Use Architecture Decision Records (ADRs) in `docs/adr/` to document non-obvious architectural and design decisions. Format: a short Markdown file per decision, numbered sequentially.

Each ADR captures: context (the forcing problem), decision (what we chose), consequences (what follows). ADRs are immutable — superseding decisions get new ADRs that reference the originals.

What does NOT belong in ADRs:
- Code-level conventions (lint config, naming).
- Reversible product decisions (CLI flag names, output formatting).
- Anything that can be derived from reading the current code.

## Consequences

- New non-trivial decisions require an ADR before merge.
- ADR review is part of code review.
- The ADR index in `docs/adr/README.md` is the entry point for understanding "why is it this way?"
- Existing decisions are backfilled as ADRs when they come up — not all at once.

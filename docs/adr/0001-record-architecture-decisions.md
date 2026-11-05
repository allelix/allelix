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

## Amendment 2026-10-30 (2.6.0): released ADRs are append-only, and a test enforces it

Owner order.

- A released ADR's lines are not changed, moved or deleted.
- A clarification that keeps the decision is appended at the end of the file as a dated amendment section.
- A changed decision gets a new ADR that references the old one.
- The one released line that may change is the Status line, and only to a "Superseded by ADR-NNNN" or "Superseded in part by ADR-NNNN" form naming an ADR that exists.
- `tests/data/adr_released.json` holds, per released ADR, the number of lines released and their SHA-256. `tests/test_adr_released_lock.py` fails when a locked line is changed or deleted. ADRs added since the last release are free until it ships.
- At each release `scripts/adr_lock.py` refreshes the lock: it adds the new ADRs, grows each entry to the full file, and refuses when a locked line changed. `scripts/adr_lock.py --covers` fails when the lock does not cover the released tree.
- The index, `docs/adr/README.md`, is not locked.

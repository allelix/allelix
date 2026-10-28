# ADR-0058: Own analysis providers through their final use

Date: 2026-10-21
Status: Accepted (2026-10-22)

## Context

Providers acquire SQLite and tabix handles lazily. Analysis used to close its
providers after streaming, then reopen them during enrichment, scoring and
source bookkeeping. A skipped ClinVar provider can also be used through a
different source's reference-filter callback. Being owned and being an active
annotation source are separate concerns.

A ZIP wrapper can fail before its member or archive exists. Cleanup that
assumes both exist can obscure the intended refusal during collection, or
leave the archive open when member cleanup raises.

## Decision

`run_analysis` owns the union of active annotators, optional providers and
explicit `resource_owners`, deduplicated by object identity. Ownership starts
before metadata/build work and ends after result construction. Each object
closes once per run; distinct instances with the same source name remain
distinct owners. Lazy acquisition, source selection, skipped-source reporting
and annotation/filtering semantics stay the same.

The CLI supplies its complete existing registry, including callback-only
providers. Direct API callers must supply callback-only providers explicitly:

```python
result = run_analysis(
    input_path,
    parser,
    active_annotators,
    resource_owners=registry,
)
```

The function does not recursively inspect arbitrary callback closures.
Submitted resources, including already-open ones, close on early failures,
interrupts and successful return. Every registered cleanup is attempted. A
primary computation error survives secondary cleanup errors; cleanup failures
after successful computation are reported. Independently short-lived cache
status/version connections retain their own context-managed cleanup.

The command-level resolved stream protects identity placement lookup with its
ownership scope. It explicitly owns the parser, build, resolution and output
iterators it creates, closing wrappers before providers on early context exit
even if the returned stream remains retained. It does not exhaust an incomplete
stream or close arbitrary caller-owned inputs to the shared iterator helpers.
Resolution, export and biological policy are unchanged.

Watch-reset callbacks are registered before fallible watch setup in both the
analysis and command paths. All registered resets and provider/iterator closes
are attempted even if another cleanup raises; an existing computation error or
interrupt remains primary. Without an existing error a cleanup failure remains
visible. A deliberately failing reset is known to have been attempted, not
necessarily completed.

`_ZipMember` initializes nullable owned fields before acquisition, detaches them
when closing, and attempts member, archive and base cleanup. Failed construction
preserves the original exception. Explicit close of a constructed wrapper
reports the first cleanup failure once, while still attempting remaining
cleanup. Repeated close and collection do not retry disposed owners. A resource
whose own close fails before releasing itself is only known to have received a
cleanup attempt; no stronger release guarantee is claimed.

Python 3.13 joins the existing 3.11/3.12 validation matrix and package classifiers.
Warnings remain errors and the coverage floor and publication workflow remain
unchanged. Routine private work uses local regression/static checks and
independent review. Hosted Linux validation is reserved for the finished
candidate before main/public release, together with the owner-run FULL Mac
battery on the exact final RC. Windows execution is not a release gate.

## Verification and limits

Offline regressions inspect every actual acquired SQLite handle before any
test-owned cleanup, including enrichment, both-build ClinVar rescoring,
callback-only providers, repeated runs and early/late errors. Tiny ZIPs cover
construction refusals, natural collection, member acquisition and cleanup
failures, stored/deflated reading and read-error unwinding. Synthetic tabix
handles exercise CADD's backend ownership without downloading source data.

No output schema, cache version, source freshness or gzip diagnostic changes
follow from this lifetime repair. Targeted checks do not establish whole-genome,
live-source, complete platform or final release readiness.

# ADR-0030: Data Source Update Cadence Drives Refresh and Integrity

**Status:** Accepted
**Date:** 2026-06-09

## Context

Allelix downloads from six external data sources. Each source has
different release cadences, different upstream infrastructure, and
different integrity guarantees. Before this ADR, the refresh and
integrity strategy was decided source-by-source without a named model,
making it unclear how a new source should behave. ADR-0028's
`local_version_tag` signal-loop bug was partially caused by this: the
SNPedia loader copied the gnomAD pattern without recognizing they belong
to different cadence tiers.

## Decision

Every data source belongs to one of two tiers. A new source picks a
tier, not a neighbor to copy.

| Tier | Sources | Freshness | Integrity |
|---|---|---|---|
| **Server-driven** | ClinVar, GWAS Catalog, PharmGKB | Re-check remote signal (ETag, MD5 sidecar) every `db update` | ClinVar: NCBI sidecar MD5. GWAS/PharmGKB: none possible — no upstream checksum (ADR-0029) |
| **Code-driven** | gnomAD, AlphaMissense, SNPedia | Commit-pinned URL + pinned SHA256; updates only via code change in an allelix release | Pinned SHA256 verified after download |

### Server-driven (Frequent)

The upstream publishes new data on its own schedule (ClinVar: monthly;
GWAS: ~weekly; PharmGKB: irregular). `db update` probes for a freshness
signal change. If the signal differs from the cached one, the annotator
re-downloads. `_maybe_refresh_databases` performs the same check during
`allelix analyze` for stale caches.

### Code-driven (Rare/Fixed)

The upstream data is a large pre-built cache hosted on HuggingFace. The
download URL is pinned to a specific commit SHA and the file content is
pinned to a specific SHA256. Neither can change without a code change in
allelix. Updates are explicit: bump the commit SHA, bump the expected
SHA256, ship a new release. Runtime freshness checks against the server
are unnecessary — the URL is immutable.

## Consequences

- New data sources must declare their tier in the loader module and
  follow the corresponding refresh/integrity pattern.
- Code-driven sources do not need runtime freshness probes (HEAD
  requests) — their `local_version_tag` + pinned SHA256 are sufficient.
- Server-driven sources without upstream checksums (GWAS, PharmGKB) are
  documented gaps in ADR-0029, not bugs to fix.

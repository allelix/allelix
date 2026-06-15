# ADR-0012: Freshness detection via per-annotator remote signal

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

`allelix db update` shipped in v0.2.0 with two failure modes:

1. **v0.2.0–v0.4.0**: re-downloaded every annotator on every invocation. Hundreds of MB of waste.
2. **v0.4.1**: skipped if `is_ready()`. Avoided waste, but had no way to detect that the remote source had changed — users had to remember to pass `--force` to get fresh data, with no signal that fresh data was even available.

The user-facing requirement: *"if remote has even a single change it needs to be updated."* A skip-if-cached policy that ignores remote state is unfit for a tool whose value depends on annotation currency (a new ClinVar release with a new pathogenic variant the user happens to carry would be silently missed).

## Decision

Each annotator implements `fetch_remote_signal() -> str | None`, which fetches a small remote artifact identifying the current published version. The signal is captured at download time and stored in `database_versions.remote_signal`. On the next `db update`, the annotator fetches the current signal and compares to the stored one.

### Decision tree (no `--force`)

| Cache state | Remote fetch | Cached signal | Action |
|---|---|---|---|
| Missing | — | — | Download |
| Present | Failed (network/None) | — | **Skip with notice** ("can't verify; pass --force") |
| Present | Succeeded | Equal | **Skip** ("already current") |
| Present | Succeeded | Different | **Refresh** |
| Present | Succeeded | None (legacy v0.4.1) | **Refresh** ("legacy cache; no stored signal") |

`--force` short-circuits to download regardless.

### Per-annotator signal mechanisms

- **ClinVar**: NCBI publishes `clinvar.vcf.gz.md5` next to the data file. Fetch it (~50 bytes), take the first whitespace-separated token, store as `"md5:<hex>"`. MD5 changes if any byte of the data file changes — exact-byte freshness for free.
- **PharmGKB**: HEAD request to the clinical-annotations URL. Prefer `ETag` (content-derived on most CDNs); fall back to `Last-Modified`. Store as `"etag:<value>"` or `"lm:<value>"`.

### Type-prefixed signal values

Stored signals carry a type prefix (`md5:`, `etag:`, `lm:`). If a server changes its signal type — e.g., PharmGKB stops sending ETags and only sends Last-Modified — the prefixed comparison fails (`etag:abc` ≠ `lm:Mon, 01 Jan…`) and we refresh, rather than silently skipping because two different-meaning strings happened to match.

### Schema migration

`database_versions` gains a nullable `remote_signal TEXT` column. New caches populate it; legacy v0.4.1 caches don't have the column at all. `get_database_info` tries the 5-column SELECT first and falls back to the 4-column SELECT on `OperationalError`, returning `remote_signal=None` for legacy caches. The decision tree treats `None ≠ remote` as "refresh," so legacy caches auto-upgrade on first v0.4.2 `db update`.

### Failure-safe network

Every fetch helper (`fetch_remote_text`, `head_request_headers`) returns `None` on any `OSError`/`ValueError` and never raises. A flaky network can't crash `db update` — at worst we see "can't be verified" and the user retries or `--force`s.

## Consequences

- **Bandwidth**: nominal `db update` invocations now do one tiny fetch (50-200 bytes) per annotator instead of downloading the full data files. Real refreshes still pay full bandwidth, but only when needed.
- **Currency vs. trust**: ClinVar's MD5 mechanism is exact (any byte change → new MD5). PharmGKB's Last-Modified mechanism is server-implementation-dependent: if PharmGKB regenerates the zip on a daily cron even with unchanged content, we'll refresh daily. Acceptable trade-off; documented in the docstring.
- **Offline operation**: nothing here changes the offline-first principle (ADR-0004). Since v1.2.0, freshness checks run by default before analysis (`--no-update` to skip). If the network is unreachable, a warning is printed and analysis proceeds against the local cache. No network access occurs during annotation itself.
- **Test surface**: the freshness logic is fully covered by the CLI test suite — match, differ, unverifiable, legacy-cache, and `--force` paths all have direct tests. Mocked HEAD/GET handlers in `tests/databases/test_remote_signal.py` exercise the network helpers without hitting live endpoints.
- **Future signal sources**: SNPedia and GWAS Catalog (now implemented) plug into the same protocol — `fetch_remote_signal` returns whichever opaque, type-prefixed string the source publishes. No CLI changes needed for new annotators.

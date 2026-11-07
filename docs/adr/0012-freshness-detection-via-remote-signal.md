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

## Amendment — 2.5.1

ClinVar's TSV freshness signal still comes from the submission-summary
checksum. ADR-0053 replaces first-token parsing with bounded, filename-aware
MD5 validation and requires variant-summary integrity before refresh
publication. Current equal-signal caches retain the small-probe update path.

## Amendment 2026-10-25: `--force` keeps ClinVar's integrity checks

`--force` skips the freshness comparison, not verification. ClinVar never
refreshes from a download it cannot verify (ADR-0053): if NCBI's checksum
for `submission_summary` or `variant_summary` cannot be fetched, the
refresh stops, the existing caches stay, and `db update` reports ClinVar
as failed, with or without `--force`.

## Amendment 2.6.0 (2026-10-30): header case, both fields, like-for-like comparison

Supervisor ruling, 2026-10-30. This replaces "Type-prefixed signal values"
above for HTTP signals.

- **Header case.** HTTP header names are case-insensitive (RFC 9110).
  `head_request_headers` lowercases them, and `http_signal` matches them in
  any case. Before, the code read `ETag` or `Etag` only. The ClinPGx
  download, served by a CDN that sends `etag`, therefore always fell back
  to Last-Modified.
- **Both fields stored.** An HTTP signal keeps the ETag and the
  Last-Modified when both are present: `etag:…|lm:…`, prefixed `pgkb:`
  and followed by `|cpic:…` for ClinPGx. A quoted ETag (RFC 9110 8.8.3)
  is read whole up to its closing quote, so a `|` or `:` inside it never
  splits it; an unquoted, non-standard ETag is read up to the next field.
  It stays in `remote_signal`, so a
  stamp of either older form (`etag:…`, `lm:…`) stays valid. No download is
  forced by the change.
- **Like with like** (`compare_signals`):
  - A cached and a remote signal are compared on a field both have: the
    ETag first, else Last-Modified.
  - Any other field (ClinPGx's CPIC date) must also be equal.
  - Signals that share no HTTP field cannot be compared. `db update` keeps
    the cache and says the change cannot be verified (`--force` refreshes),
    and the automatic update leaves it. A change only in which fields a
    server sends is never a change by itself.
  - Non-HTTP signals (ClinVar's checksum) are still compared whole.
- **Unchanged means untouched.** The unchanged path writes nothing; a newly
  seen header is recorded at the next real refresh.
- **Content.** A real change still downloads. For ClinPGx, ADR-0068 then
  skips the rebuild when the content is the same.

The GWAS Catalog and ClinPGx are the only sources whose signal comes from
HTTP headers.

## Amendment 2.6.1 (2026-11-05): a one-build ClinVar install is checked too

The pre-analysis freshness check and `db status` built ClinVar for both
builds. With only one build installed (`db update --build grch37` or
`grch38`), that annotator was not ready, so the check skipped ClinVar (no
probe, no refresh, no "N days old" notice) and `db status` showed it not
ready, although `analyze` used it. Both now build ClinVar from the builds
whose caches are ready, as the analysis already did
(`installed_clinvar_builds`), and fall back to both builds when none is:
a stale one-build cache is probed and refreshed for its own build only.

## Amendment 2.6.1 (2026-11-06): db status shows the builds that are there

`db status` builds ClinVar from the builds whose cache file exists, not only
the ready ones, so a damaged or outdated second build shows ClinVar not
ready instead of being left out. The freshness check keeps using the ready
builds, as above.

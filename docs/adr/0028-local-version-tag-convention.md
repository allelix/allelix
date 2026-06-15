# ADR-0028: Local Version Tag Convention

**Status:** Accepted
**Date:** 2026-06-09

## Context

Every annotator's SQLite cache has a `database_versions` row that tracks
freshness. Before v1.5.0, two kinds of information were mixed into a
single `remote_signal` column:

1. **Remote freshness signal** — what the server reported (ETag,
   Last-Modified, MD5). Used to decide whether to re-download.
2. **Local processing version** — which version of allelix's local
   processing logic (interpreter, parser, categorizer, schema) was used
   to build or process the cache. Used to decide whether to re-ingest or
   re-download even when the remote data hasn't changed.

These were concatenated with a `|` delimiter (e.g.
`etag:"abc"|iv:1`, `lm:Thu, 01 Jan 2026|pv:6`). This caused three
bugs:

- SNPedia's `cached_remote_signal()` returned the full string including
  the parser version suffix, which never matched the remote ETag,
  triggering a re-download loop on every `db update`.
- Signal comparison logic had to strip suffixes before comparing,
  spreading `split("|")` calls across multiple files.
- Adding a new annotator required discovering this undocumented
  convention and reproducing the suffix pattern correctly.

## Decision

Split `remote_signal` into two columns:

| Column | Purpose | Example |
|---|---|---|
| `remote_signal` | What the remote server reports | `etag:"abc123"`, `lm:Thu, 01 Jan 2026`, `md5:a1b2c3` |
| `local_version_tag` | What local processing version built the cache | `iv:1`, `pv:6`, `cv:3`, `sv:1` |

Every annotator stamps `local_version_tag` during cache construction.
The tag prefix indicates the kind of local processing:

| Prefix | Meaning | Used by |
|---|---|---|
| `iv:N` | Interpreter version — logic that classifies raw data into structured annotations | ClinVar, PharmGKB |
| `pv:N` | Parser version — logic that parses raw markup into structured rows | SNPedia |
| `cv:N` | Categorizer version — logic that classifies traits into categories | GWAS Catalog |
| `sv:N` | Schema version — pre-built cache schema that triggers re-download on change | gnomAD, AlphaMissense |

The prefix is a documentation convention, not a parsed field. Code
compares the full string (`local_version_tag == f"cv:{_CATEGORIZER_VERSION}"`),
never splits on `:`.

### Rules for all annotators

1. `remote_signal` holds only the remote freshness signal. Never append
   local state to it.
2. `local_version_tag` holds exactly one tag string. It is set during
   cache construction (`load_*_tsv`, `install_prebuilt_cache`, or
   equivalent).
3. `is_ready()` checks `local_version_tag` against the current version
   constant. A missing or stale tag means the cache needs rebuilding.
4. Backward compatibility for NULL or empty `local_version_tag`
   depends on the annotator family:
   - **Pure-prebuilt** (gnomAD, AlphaMissense — decompress-only, no
     local transform): `is_ready()` treats an empty tag as legacy-ready.
     Pattern: `return tag == f"sv:{VERSION}" or not tag`.
   - **Active-transform** (ClinVar, PharmGKB, GWAS, SNPedia — local
     parsing, ingest, or classification): an empty tag triggers a
     one-shot `_stamp_existing_*_cache()` migration that splits the
     legacy `|tag:N` suffix out of `remote_signal` into
     `local_version_tag`. `is_ready()` returns the migration's success
     status, not unconditionally True.

   In both families, a non-empty tag that doesn't match the current
   version means "stale, rebuild."
5. Version constants live in `allelix/databases/_versions.py`. Bump the
   constant when local processing logic changes in a way that would
   produce different cache contents from the same input data.

### Rules for new annotators

Any new annotator added to allelix must:

1. Define a version constant in `_versions.py`.
2. Choose a prefix (`iv:`, `pv:`, `cv:`, `sv:`, or a new one if none
   fits — document it in this ADR).
3. Stamp `local_version_tag` during cache construction.
4. Check `local_version_tag` in `is_ready()`.
5. Implement `fetch_remote_signal()` and `cached_remote_signal()` using
   only `remote_signal` — no suffix parsing.

### Migration

The `local_version_tag` column is added lazily:

- `get_database_info()` tries the full-schema SELECT first. On
  `OperationalError` (pre-v1.5.0 cache missing `local_version_tag` or
  `remote_signal`), it falls back to progressively simpler SELECTs and
  calls `_ensure_local_version_tag_column()` to add the column in
  place.
- For active-transform annotators, `is_ready()` calls the annotator's
  self-healing migration function (`_stamp_existing_*_cache`) when it
  finds an empty `local_version_tag`. The migration strips any
  `|tag:N` suffix from `remote_signal` and writes the tag to
  `local_version_tag`. This runs once per cache on first access after
  upgrade. (See Rule 4 above for how this connects to `is_ready()`
  behavior.)

## Consequences

- `remote_signal` is clean and directly comparable against the remote
  server's response. No more suffix stripping.
- Cache invalidation on local logic changes is explicit and uniform
  across all six annotators.
- The SNPedia signal-loop bug class is structurally eliminated.
- Adding a new annotator has a documented checklist (this ADR) instead
  of requiring reverse-engineering of the suffix convention.
- Pre-v1.5.0 caches self-heal without re-download. The migration is
  invisible to users.

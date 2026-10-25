# ADR-0048: Caches are read-only at analyze time

**Status:** Accepted
**Date:** 2026-10-15
**Amends:** ADR-0023 (ClinVar one-shot migration), ADR-0024 (GWAS
auto-reingest, 2026-05-22 amendment), ADR-0028 (migration in `is_ready()`).

## Context

`analyze` wrote to the data directory in several ways, and nothing said so:

- `get_database_info` added a `local_version_tag` column to a legacy cache.
- `is_ready()` stamped legacy ClinVar, ClinPGx and GWAS caches, re-ingested
  GWAS from the cached TSV after a categorizer bump, and parsed SNPedia's raw
  pages.
- ClinVar and gnomAD ran `CREATE INDEX IF NOT EXISTS` when they connected.
- `load_config` created `config.toml`.
- The gnomAD, AlphaMissense and CADD caches are published in WAL mode, and
  the installer kept that mode. A WAL file cannot be read, even with a
  plain `SELECT`, in a directory the user cannot write unless its `-wal` and
  `-shm` files already exist. In a writable directory, reading it creates
  those files and leaves them.

So a read-only install (a shared data directory, a container image) failed.
Every failure looked the same: a missing, out-of-date, unreadable or
corrupt cache and a read-only WAL cache were all "not ready", with the
advice "run `allelix db update`". That advice cannot fix an unreadable file
or a read-only directory.

## Decision

- **With `--no-update`, `analyze` writes nothing to the data directory.**
  No index, no `ALTER`, no stamp, no re-ingest, no parse, no config file and
  no WAL side files. A cache that would need any of these is skipped, and
  its reason is printed.
- **Without `--no-update`, the update step may write, as before**, but
  `is_ready()` never does. The update step runs before readiness is
  checked, in `db update` and in `analyze`:
  1. It converts every WAL cache in the data directory to rollback-journal
     mode.
  2. It runs `Annotator.prepare_cache()` for each outdated cache. That holds
     the local, offline work that used to live in `is_ready()`: the legacy
     stamps, the GWAS re-ingest from the cached TSV and the SNPedia parse,
     each unchanged.
- **Every cache opens read-only**, through a `mode=ro` URI with the path
  quoted (`connect_readonly`). This covers annotator connections,
  readiness checks, `get_database_info` and the loaders' schema checks.
  The URI resolves the physical installed file first, so SQLite finds its
  WAL/SHM at the same location on Windows and POSIX. Update-time conversion
  and subsequent sidecar cleanup use that same physical target. Unreadable,
  dangling or cyclic links remain conversion no-ops. This path repair does
  not change cache contents, interpretation or version tags and requires
  no rebuild or download.
- **A WAL cache is never opened**, not even for its version. The WAL check
  reads header bytes 18 and 19 without opening SQLite. `get_database_info`
  reports a WAL file as not present, so AlphaMissense's check of gnomAD's
  version follows its path for an unavailable gnomAD. `db status` shows
  "WAL mode from an older Allelix; run `allelix db update` to convert it"
  in place of a version. `immutable=1` is never used, because it ignores
  WAL pages that have not been checkpointed (a `-wal` file beside the real
  file of a symlinked cache, for example).
- **In `db update` the conversion runs first**, with write access, and the
  loaders' schema checks then open the converted file normally. A cache
  still in WAL mode after that step (no write access, or another process
  has it open) is reported and not downloaded again.
- **No runtime index.** `idx_clinvar_position` and `idx_gnomad_position`
  are in the loader schemas. Position lookups use the indexes those
  schemas build.
- **The loader owns the journal mode.** `install_prebuilt_gz_cache` sets
  `PRAGMA journal_mode=DELETE` on the staged file before moving it into
  place. Afterwards it removes any `<name>-wal` and `<name>-shm` left from the
  old cache.
- **One readiness truth.** `cache_state(path, name, ...)` returns a typed
  state:
  - `ready`;
  - `missing`;
  - `outdated`: the version tag, the version row or a required column;
  - `unreadable`: permission;
  - `corrupt`: not an SQLite database, a required table is absent, or
    the main table has no rows (such a cache used to pass as ready);
  - `wal`: a WAL file in a directory the user can write. Reading it
    read-only would leave `-wal` and `-shm` beside it, so it is skipped
    until the update step converts it;
  - `readonly_wal`: a WAL file in a directory the user cannot write.

  Each annotator's `cache_status()` is built on it, and `is_ready()` is
  `cache_status().ready`.
- **One console line per skipped source**, naming the state and the fix:
  - missing or outdated: run `allelix db update`;
  - unreadable: check the permissions on the file;
  - corrupt: run `allelix db update`, which replaces only the damaged source
    (`--force` would download every source again);
  - wal: run `allelix db update` once to convert it; the conversion
    downloads nothing;
  - readonly_wal: run `allelix db update` once with write access to the
    data directory.

  Exit codes are unchanged: 0 when a source is skipped, 1 when none is
  ready. In the second case the error lists each source's reason. The
  states are kept on the analysis result for the reports to show.
- A ClinVar cache that cannot be read is not opened, even as ClinPGx's REF
  source. An outdated one is still read for that, as before.

## Consequences

- Read-only installs work once their caches are in rollback-journal mode,
  and the console says why a source was skipped.
- A gnomAD, AlphaMissense or CADD cache installed by an earlier release is
  in WAL mode, so `analyze --no-update` skips it until one `db update`, or
  one `analyze` without `--no-update`, converts it. This release already
  needs a `db update` for ClinVar (iv:9).
- A legacy cache that only needs a stamp, and a GWAS cache after a
  categorizer bump, are brought up to date by the update step. Under
  `--no-update` they are skipped as outdated until the next run without it.
- The published WAL files stay as they are; the installer converts them on
  install.

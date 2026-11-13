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

## Amendment 2.6.0 (2026-11-01): switched-off sources

In `analyze` and the focused reports, the automatic update step leaves the
cache files of a source switched off for the run untouched: no WAL
conversion, no `prepare_cache`, no refresh, and its cache never decides
readiness. `db update` still converts and prepares every cache (ADR-0067).

## Amendment 2.6.1 (2026-11-05): SNPedia parses a staged copy

The SNPedia parse ran on the published cache: `setup()` installed the raw
pages over `snpedia.sqlite`, then `prepare_cache()` parsed that file in
place, in one long transaction, and returned False on an SQLite error. A
Ctrl-C or SQLite error during the parse left the cache unparsed; two
refreshes at once could leave it unparsed while both reported success; and
a parse killed outright left a hot `-journal` beside the published file,
which a read-only open cannot recover and the next read-write open would
replay into a newly installed file.

Both now work in an `owned_directory` of their own (ADR-0056, 2.6.0
amendment): `setup()` installs the download there, and `prepare_cache()`
copies the outdated cache there through SQLite's backup API. The parse,
or the legacy `|pv:` stamp, runs on that copy, and the copy is moved onto
`snpedia.sqlite` only once it reads as ready. Until then the cache you had
is the one published. A parse that leaves no usable cache raises in
`setup()`, so `db update` reports a failure, not a success. In
`prepare_cache()` the failure is printed and the outdated cache is left as
it was.

Every prebuilt install, and SNPedia's publish, goes through
`publish_cache`, with one rule for every side file the old cache leaves.
A hot `-journal` is what the old file needs to get back its committed rows,
and a `-wal` holds committed pages the old file has not absorbed yet; a
reader of a new file at the same name would apply either to the new file.
So, before anything is replaced, the old database is opened read-write:
SQLite rolls back a hot journal and checkpoints a WAL
(`PRAGMA wal_checkpoint(TRUNCATE)`), and closing it removes the WAL and its
`-shm`. Then no `-journal`, `-wal` or `-shm` may remain beside it. If
recovery fails (the old file is locked or cannot be opened) or any side file
remains, the install stops with `UnrecoveredJournalError`, which
names `db clean`, before `os.replace`: the old file keeps its committed rows
and its side files. Only then is the new file moved in, and its name's side
files removed. Side files whose old cache file is gone are orphans that hold
nothing committed of anything: they are removed before the new file moves
in, as the installer did before 2.6.1.

## Amendment 2.6.1 (2026-11-06): which installs settle, and an empty journal

The amendment above holds for the installs that go through `publish_cache`:
the prebuilt caches (gnomAD, AlphaMissense, CADD, the rsID identity table)
and SNPedia. ClinVar, GWAS Catalog and ClinPGx still publish their cache with
a bare `os.replace`, as in 2.6.0; moving them onto `publish_cache` is left
for a later release.

A `-journal` that holds nothing (empty, or without SQLite's journal header)
is not hot: SQLite ignores it and leaves it in place when the old file is
opened, and refusing over it blocked every later install. It is also what a
live writer's open transaction has, so `_settle` does not judge it by its
bytes. After the open, it runs a no-op write in an exclusive transaction
(`BEGIN EXCLUSIVE` with no busy timeout, `user_version` rewritten to its
own value, `COMMIT`): SQLite takes the lock and removes a stale journal on
commit. If the lock is busy, the old cache is in use and nothing is
replaced. A journal still left after the commit refuses, as above. No
journal is removed by hand while the database could be in use.

## Amendment 2.7.0 (2026-11-10): every staged cache replacement goes through `publish_cache`

The installs the 2.6.1 amendment left on a bare `os.replace` now go through
`publish_cache` too: the ClinVar (both builds), GWAS Catalog and ClinPGx
loaders, the rsID identity table's writer, and the offline builders
`scripts/build_gnomad_cache.py`, `build_alphamissense_cache.py` and
`build_cadd_cache.py`. The rule above applies to each of them. An old hot
journal or committed WAL is settled into the old file before the new one
moves in, or the install stops with `UnrecoveredJournalError` and the old
cache keeps its committed rows. ClinVar settles every selected build's old
cache before it replaces any (`settle_for_publication`), a recovery
preflight: a refusal found then changes nothing. Each build's move after it
is atomic on its own, not the two together; a later failure can leave an
earlier build updated, and the error names the builds already updated.
Downloaded archives and the retained GWAS table are not SQLite caches and
are moved as before. The in-place writers (`stamp_existing_*`,
`update_remote_signal`, `convert_to_rollback_journal`) are unchanged.

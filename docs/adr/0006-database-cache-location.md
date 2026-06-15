# ADR-0006: Database cache location and resolution

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

ADR-0004 established that reference databases live on the user's machine, not bundled in distribution. That leaves the question: *where* on the user's machine?

Three forces:
1. **Default must be sane on a fresh install.** Users running `pip install allelix && allelix db update` should not have to think about paths.
2. **Devs need to override.** Working from a project checkout, having `data/` next to the source tree is convenient.
3. **CI / tests must isolate.** Each test run gets its own scratch directory; nothing leaks into a shared cache.

## Decision

Resolve the data directory with this precedence, highest first:

1. CLI `--data-dir PATH` flag on commands that touch the cache.
2. `ALLELIX_DATA_DIR` environment variable.
3. `XDG_DATA_HOME/allelix` if `XDG_DATA_HOME` is set.
4. `~/.local/share/allelix` (XDG default fallback).

`resolve_data_dir(override)` creates the directory if missing and returns the resolved Path.

Stdlib only — no `platformdirs` dependency. Linux/Mac get the right behavior; Windows `%APPDATA%` support is deferred until demand materializes.

## Consequences

- One function (`allelix.databases.resolve_data_dir`) is the single entry point for cache location. Annotators never compute paths themselves; they receive `data_dir` at construction.
- Tests pass `tmp_path` as `--data-dir` and get a fresh, empty cache per test.
- Developers can `export ALLELIX_DATA_DIR=$PWD/data` to keep cache local to a checkout.
- Users on a stock install get `~/.local/share/allelix/` and don't have to know.
- Windows users on v0.2.0 get `~/.local/share/allelix/` too (works, but unidiomatic). PyPI release adds proper Windows support.

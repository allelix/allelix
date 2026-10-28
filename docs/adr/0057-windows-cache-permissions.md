# ADR-0057: Check Windows cache directory permissions without writing

## Status

Accepted

## Context

Windows `os.access` can report a directory writable despite an ACL that denies
file creation. Opening a directory with backup semantics also bypasses that
denial when the runner has backup privileges. Analysis must not write probe
files or open a WAL database just to classify it.

## Decision

Use Windows AccessCheck for FILE_ADD_FILE against the physical directory's
owner, group and DACL. Check the effective thread token, falling back to the
process token only when the thread has no token, and duplicate it for the check.
Bind pointer-sized handles explicitly, resize descriptor and privilege buffers
with bounded retries, and close every acquired token even after failures.
Keep API success separate from access granted; unavailable checks conservatively
report no write access with a diagnostic. Do not adjust privileges.

Use this check for WAL readiness, update-time journal conversion and explicit
`db path --check`. Preserve existing POSIX access masks and original diagnostic
paths. Resolve symlinked files before locating SQLite sidecars and their parent.
No cache format, payload or version change is needed.

## Consequences

Permission checks are snapshots; later I/O can still fail. Analysis remains
read-only, and only database update performs journal conversion. Tests preserve
POSIX permission assertions and use verified real ACL denial on Windows,
restoring the original descriptor in finally. Native call-contract controls
also run offline; they do not establish actual Windows integration success.

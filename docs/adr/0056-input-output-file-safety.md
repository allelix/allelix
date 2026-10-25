# ADR-0056: Protect input files and own staged output resources

## Status

Proposed

## Context

A report destination equal to its input replaces that input. Fixed report and
PLINK temp names can overwrite or delete unrelated files, including an input
named like a former temp. Direct VCF output can leave a partial destination on
a streaming failure. CLI safety checks must precede cache and output mutation.

## Decision

Validate every requested report/VCF destination and all three actual PLINK
paths before cache-directory setup, refresh or output allocation. Use resolved
identity for lexical/new paths and samefile for existing hard links; real
stat/resolve errors propagate. Refuse input aliases of former fixed report and
PLINK temps as well, without opening or deleting those old temp paths. Direct
report/VCF APIs validate their known input path. The iterator-only PLINK API
cannot infer an input origin; its CLI supplies that boundary check. The CLI's
existing VCF no-overwrite policy stays intact.

Allocate unique same-directory temps with mkstemp, prefix `.` plus destination
name plus `.`, suffix `.tmp`. Consume the returned descriptors, finish writing
and close successfully before os.replace. Stage all three PLINK files before
first publication. Keep new staging files private and apply only existing
ordinary permission bits (0777) to complete content before publication, never
copying setuid/setgid/sticky bits. Supported replacement retains those bits;
permission failures preserve the previous destination.

Clean only allocated descriptors and still-owned unpublished temp paths on
ordinary failures and interrupts. Successful publication ends ownership of the
old staging pathname; cleanup must leave any later unrelated file there alone.
On Windows, a still-owned temp made readonly for publication may need its own
readonly attribute cleared after replacement fails. Never change the final
file's permissions to make cleanup succeed, and preserve the primary error.

## Consequences

Prepublication failure preserves all prior outputs and creates no partial
newly named output. PLINK publication is atomic per file: later rename failure
can leave complete already-published files alongside unchanged prior files.
Cleanup never deletes those published outputs as purported rollback.

JSON/HTML serialization, streaming VCF UTF-8/LF, PLINK content/order/counts and
schema versions are unchanged. Filename/build/allele changes and output/output
collisions remain separate work. There is no claim of preserving ownership,
ACLs, xattrs or link topology, crash durability, a three-file transaction, or
protection from hostile concurrent path swaps. An uncatchable kill may leave
an owned temp. Cross-platform behavior requires platform-specific verification.

# ADR-0053: ClinVar refresh integrity and publication

Status: Accepted

## Decision

ClinVar refresh verifies compressed `variant_summary.txt.gz` against its
published MD5 before ingest. Submission-summary verification is retained
when its freshness checksum is available. Variant integrity is required
even with `--force`. A hash mismatch, including a gzip cut at a valid
member boundary, refuses the candidate. Gzip parsing and count validation
remain required after hashing.

Sidecars are read with a 4,096-byte body limit and the existing signal
timeout. Each nonblank line must contain a 32-digit hexadecimal digest,
optionally followed by the expected source basename, an absolute or
relative path ending in it, or the standard binary marker `*` before the
name. Digest-only bodies are supported. Digests are normalized to
lowercase; conflicting digests, incorrect filenames, malformed lines
and oversized bodies are rejected. MIME type is not used to parse them.

Each build requires at least 5,000,000 installed `clinvar_variants` rows
and `10 * new_count >= 9 * previous_count`. The baseline comes from that
build's existing metadata before any publication. Both per-build names
and the legacy `clinvar` name are supported. Missing or invalid counts
on an existing cache cause a clear refusal, rather than treating it as
absent or scanning the full old table. An absent cache uses only the
fixed floor. Historical metadata counts are not independently certified
against old table contents.

The fixed floors leave headroom below the retained 2026-10-17 `iv:9`
snapshot (6,342,163 GRCh37 and 6,341,897 GRCh38 rows). These are snapshot
counts, not current upstream measurements. The relative comparison uses
exact integer arithmetic; for example, an old count of 11 requires 10.
A legitimate larger drop needs investigation rather than an automatic
lowering of the guard.

All selected caches are built in one temporary directory on the cache
filesystem. Each candidate receives its final per-build name, date,
interpreter tag, installed count and freshness metadata before publication.
The loader writes this metadata even if the freshness signal is null;
null does not assert successful checksum verification. The broader CLI
missing-freshness override is a separate change. *Amended (2.5.1):* there
is no such override. ClinVar's freshness signal is the submission_summary
checksum, which also verifies that download, so a refresh without it would
be unverified. Setup retries the fetch once; if it still fails, it stops
without advising `--force` (which reaches the same step), keeping the
existing caches.

Only after every candidate passes validation are caches replaced, one
atomic file replacement per selected build. Prepublication failures leave
all selected old caches intact. If publication fails, the failing build
keeps its old file, and any earlier successful replacement remains a
complete, valid new cache; the error identifies that partial result.
There is no cross-file crash transaction or backup-copy mechanism.
The temporary directory owns downloads, partial files and candidates,
including cleanup if the second download fails.

The interpreter stamp stays at `iv:9`. These checks run when a cache is
built and do not change what the loader stores: a cache built with them
holds the same tables, rows and metadata fields as one built by 2.5.0 from
the same files. Under ADR-0023 that is no reason to invalidate existing
caches, so a current `iv:9` cache stays ready and the checks apply at its
next build or refresh. Analysis remains read-only. A healthy current cache
with equal freshness uses the existing remote-probe-only update path.

## Validation and limits

Offline synthetic regressions cover checksums, concatenated gzip cuts,
both build policies and boundaries, metadata, cleanup and publication
failures. The live integration ingest asserts the production fixed floor.

Live NCBI gzip member structure remains unverified. During the next
authorized normal refresh, inspect the already-required downloaded
archives with a bounded streaming decoder, recording complete member
counts, compressed sizes, digests and observation date before cleanup.
No additional archive download is needed for that observation. Synthetic
multi-member tests demonstrate the defense without establishing the live
archive structure or current upstream row counts.

## Amendment 2026-10-25 (2.5.2): build diagnostics in the log

Each ClinVar build writes two lines, whatever the outcome, to the
persistent build log `clinvar-build.log` in the data directory. They are
appended with a UTC timestamp, so `db update` and the automatic refresh
keep them with no logging setup. The same lines also go to the `allelix`
logger at INFO. The outcomes covered: a success, a fixed-floor or 90%-floor refusal, or a failure
after insertion that rolls the candidate back.

1. **Sources.** Each source's file name, the MD5 of its bytes on disk (for
   a gzip file, the compressed bytes) and its size in bytes. These are
   written before the build starts.
2. **Counts.** These are written in a `finally`, so they exist after any
   exception, connection close or cleanup of downloads and staging:
   - variant_summary rows accepted for the build;
   - submission_summary rows kept (not flagged, not a placeholder
     classification);
   - join hits and misses (the VariationID groups looked up in
     variant_summary, found or not);
   - records yielded;
   - rows SQLite stored (measured from the connection's change count, so
     a storage loss shows up as fewer rows inserted than yielded). It is
     kept after every insert batch, and when a batch itself fails, so a
     source or SQL failure partway through still logs the rows stored
     before the rollback;
   - the candidate's final row count.

   A stage that was not reached reads "not reached".

A refused build, such as a candidate far below the floor, can then be
traced to the stage that lost the rows. Nothing about the floors, the stamp or the
stored cache changes.

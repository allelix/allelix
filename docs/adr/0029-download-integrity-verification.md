# ADR-0029: Download Integrity Verification

**Status:** Accepted
**Date:** 2026-06-09

## Context

`db update` downloads database files from four upstream sources:

| Source | Provider | Integrity before v1.5.1 |
|---|---|---|
| ClinVar | NCBI | Content-Length only |
| gnomAD, AlphaMissense, SNPedia | HuggingFace | Content-Length only |
| GWAS Catalog | EBI/NHGRI | Content-Length only |
| PharmGKB | api.pharmgkb.org | Content-Length only |

Before v1.5.1, the download pipeline verified Content-Length (truncation
guard) but never verified the content of the downloaded bytes against a
known-good hash. A corrupted or tampered download would be silently
accepted and used for analysis.

Two separate integrity gaps among the four sources that **can** be
verified:

1. **ClinVar.** NCBI publishes an MD5 sidecar file (`clinvar.vcf.gz.md5`)
   for every VCF. `_fetch_remote_signal_for()` already fetches this MD5
   and stores it as the freshness signal, but the downloaded VCF was
   never verified against it.

2. **HuggingFace caches.** The three pre-built SQLite caches (gnomAD,
   AlphaMissense, SNPedia) were fetched from mutable `resolve/main/`
   URLs. A force-push to the HuggingFace repo would silently change the
   file content served at the same URL. No hash verification of any kind.

## Decision

### ClinVar: MD5 verification

After downloading the VCF, compute its MD5 and compare against the hash
already fetched from the `.md5` sidecar. On mismatch, delete the
downloaded file and raise. This uses infrastructure that already existed
(the signal fetch) — the only missing step was the comparison.

### HuggingFace: commit-pinned URLs + SHA256 verification

1. Pin each cache URL to a specific HuggingFace commit SHA
   (`resolve/<commit>/filename` instead of `resolve/main/filename`).
   This makes the download deterministic — a force-push to the repo's
   `main` branch does not change what allelix fetches.

2. Store the expected SHA256 (the Git LFS object ID) as a constant in
   each loader module. After downloading the `.sqlite.gz`, verify its
   SHA256 against the constant. On mismatch, delete and raise.

3. Updating a cache requires a code change: bump both the commit SHA in
   the URL and the expected SHA256 constant. This makes cache updates
   explicit and code-reviewable.

### GWAS Catalog and PharmGKB: no hash verification (documented gap)

Neither EBI (GWAS Catalog) nor PharmGKB publishes a checksum sidecar or
digest header for their download endpoints. Both serve mutable content
that changes on their release cadence, so there is no stable hash to pin
against from either direction.

These two sources rely on **TLS transport integrity** and the existing
**Content-Length truncation guard** only. This is a known, accepted gap.

**Revisit trigger:** if EBI or PharmGKB begins publishing a checksum
(sidecar file, `Digest` header, or API field), adopt ClinVar's sidecar
pattern via the existing `verify_file_hash()` — no other change needed.

### Shared `verify_file_hash()`

A single utility function in `manager.py` handles both cases:

```python
verify_file_hash(path, algorithm, expected_hex)
```

Streams the file in chunks (no multi-GB memory load), compares hex
digests, and on mismatch deletes the file before raising `OSError`.

## Consequences

- **4 of 6** sources are content-verified (ClinVar + 3 HF caches).
  GWAS Catalog and PharmGKB are not — no upstream checksum exists.
- HuggingFace cache URLs are immutable. A force-push to the repo does
  not affect existing allelix releases.
- Cache updates require a code change (new commit SHA + new SHA256),
  making them auditable in git history.
- The MD5 sidecar that ClinVar publishes is now actually used for its
  intended purpose, not just as a freshness signal.

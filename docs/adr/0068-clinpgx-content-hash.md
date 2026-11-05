# ADR-0068: ClinPGx keeps its cache when a download's content is unchanged

**Status:** Accepted
**Date:** 2026-10-30

## Context

`db update` refreshes ClinPGx when its freshness signal changes
(ADR-0012, ADR-0020). The signal is the archive's ETag (or Last-Modified)
together with CPIC's latest change date. The archive is an S3 object, and
its ETag is the file's MD5, so an ETag change means the bytes changed.

New bytes are not always new content. A re-upload can repack the same
files with new timestamps, another compression or another member order,
or add a file the build does not read. Each of these rebuilt the cache
and replaced it on disk. The owner's rule is no extra downloads and no
disk waste.

## Decision

Supervisor ruling, 2026-10-30, option A.

- **Unchanged signal:** only the lightweight freshness probes, as before:
  the ClinPGx HEAD request and CPIC's change-date check. No archive
  download, no CPIC allele-function data download, no write.
- **Changed signal:** the archive is downloaded and CPIC's allele-function
  data fetched, as before. Their content hash is then computed:
  - It is SHA-256 over the two files the loader reads,
    `clinical_annotations.tsv` and `clinical_ann_alleles.tsv` (name, length
    and bytes), and the CPIC lookup in sorted, canonical form.
  - The archive's packaging and the lookup's order do not count.
- **Same hash as the last build, and the cache current:** the cache is
  kept, with no rebuild and no copy. Only the new signal is stored, so the
  next update does not download again. `db update` prints "content
  unchanged; cache kept".
- **A real change:** a different hash, from either TSV or from CPIC data
  alone, rebuilds the cache, as does a cache that is not current (an older
  interpreter version). Equal bytes never keep an outdated interpretation.
- **Storage:** the hash is stored in a new optional `content_hash` column
  of `database_versions`, added only to ClinPGx caches.
  - A cache built before 2.6.0 has no such column. It stays current, is
    not rebuilt for that reason, and records the hash on its next real
    refresh.
  - The shared stamp schema of the other sources does not change.

## Consequences

A re-upload of the same content costs one download (about 1.2 MB) and no
rebuild. A change that moves only the ETag is still caught, because the
download is compared by content. The other sources are unchanged.

## Amendment 2.6.0 (2026-11-03): the kept signal lands only on the checked content

The unchanged-content shortcut compared the download's hash with the
installed cache, then stored the new signal on whatever cache the path
named by then. A refresh that installed different content between the two
had its own signal replaced, so its rows carried a signal that described
other content. The signal is now stored only while the cache still holds
the hash that was checked, in one statement (`update_remote_signal` with
`content_hash`). When another refresh has won, its cache and signal are
kept, a warning says so, and "content unchanged" is not reported.

The `content_hash` field is no longer ClinPGx's alone: a GWAS Catalog cache
records one too, for the table it keeps for re-ingest (ADR-0024, amendment
of this date). The stamp schema of the other sources does not change.

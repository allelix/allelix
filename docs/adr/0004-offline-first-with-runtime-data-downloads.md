# ADR-0004: Offline-first with runtime data downloads

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Two forces push the same direction:

1. **Privacy.** Genotype data is the most identifying data that exists. It cannot be un-leaked. Any architecture that puts genotype bytes on a network — telemetry, cloud lookups, "anonymous" usage stats — leaks data the user can never recall.
2. **Licensing.** Reference databases have incompatible licenses. ClinVar, dbSNP, and GWAS Catalog are public domain (NCBI). PharmGKB is CC BY-SA 4.0. SNPedia is CC BY-NC-SA 3.0 — non-commercial. Bundling SNPedia content with an MIT-licensed Allelix distribution would create a license conflict; bundling PharmGKB without attribution would violate ShareAlike.

A naive design would query reference databases over the network at analysis time. That violates privacy (the rsIDs you look up reveal what variants you carry) AND creates network/rate-limit dependencies for offline analysis.

## Decision

Allelix ships with **zero third-party data**. Reference databases are downloaded by the user via `allelix db update` and cached locally. Analysis runs entirely against the local cache — no network access required after setup.

Each downloaded database retains its original license on the user's machine. The user, not Allelix, is the licensee. Allelix's role is plumbing.

For SNPedia specifically: content is downloaded and cached the same way as ClinVar (the source's API permits research-purpose retrieval), but the CLI provides `--exclude-snpedia` so commercial users can run analysis without SNPedia-derived annotations. README and reports must attribute SNPedia content prominently.

No telemetry. No analytics. No "optional" crash reports that might contain file paths or rsIDs. Privacy is non-negotiable.

## Consequences

- First-run UX: users must run `allelix db update` before `allelix analyze` works. Setup time is real (databases are large).
- The cache lives in `data/` (gitignored). Disk cost is the user's responsibility.
- Allelix does not need a backend, an API key, or a hosted service. Forever.
- Database staleness is the user's problem; `allelix db status` reports versions and freshness so the user can decide when to refresh.
- README requires a "Data Sources & Licensing" section listing every database, its license, its URL, and any usage restrictions. This ships in v0.1.0 even though annotators don't land until v0.2.0 — contributors and users need to see the licensing posture from day one.
- Architectural rule: any code that calls out to the network during analysis is a bug. Network is acceptable only inside the `databases/` module's update path.

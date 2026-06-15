# Architecture Decision Records

This directory holds the design decisions that shape Allelix.

Each ADR captures one decision: the context that forced it, what we chose, and the consequences. ADRs are immutable — if a decision changes, write a new ADR that supersedes the old one rather than editing history.

## Index

- [ADR-0001: Record architecture decisions](0001-record-architecture-decisions.md)
- [ADR-0002: Plugin-based parsers and annotators](0002-plugin-based-parsers-and-annotators.md)
- [ADR-0003: Source-attributed annotations (regulatory posture)](0003-source-attributed-annotations.md)
- [ADR-0004: Offline-first with runtime data downloads](0004-offline-first-with-runtime-data-downloads.md)
- [ADR-0005: SNP count is parse-derived, not metadata](0005-snp-count-is-parse-derived.md)
- [ADR-0006: Database cache location and resolution](0006-database-cache-location.md)
- [ADR-0007: Genotype matching requires the user to carry the ALT allele](0007-genotype-matching-requires-alt-allele.md)
- [ADR-0008: Magnitude scoring is Allelix-derived from source classifications](0008-magnitude-scoring-from-clinsig.md)
- [ADR-0009: PharmGKB matches the user's exact normalized diploid call](0009-pharmgkb-genotype-matching.md)
- [ADR-0010: Strand-flip helpers ship in v0.4.0; liftover is deferred](0010-strand-handling-and-liftover.md)
- [ADR-0011: Indel-anchor protection for array-based parsers](0011-indel-anchor-protection-for-array-parsers.md)
- [ADR-0012: Freshness detection via per-annotator remote signal](0012-freshness-detection-via-remote-signal.md)
- [ADR-0013: PharmGKB non-finding suppression](0013-pharmgkb-non-finding-suppression.md) — superseded by ADR-0023
- [ADR-0014: Somatic-variant suppression for germline parsers](0014-somatic-variant-suppression-on-germline.md)
- [ADR-0015: Mock data generators are the contract](0015-mock-data-as-spec.md)
- [ADR-0016: Data Classification Principle (structured fields only; regex forbidden in production)](0016-data-classification-principle.md)
- [ADR-0017: PharmGKB SNV row-level prose fallback](0017-pharmgkb-snv-prose-fallback.md) — superseded by ADR-0020
- [ADR-0018: PharmGKB per-allele function via CPIC template extraction](0018-cpic-template-allele-function-extraction.md) — superseded by ADR-0020
- ADR-0019 was reserved during drafting and not used.
- [ADR-0020: CPIC API is the structured per-allele function source](0020-cpic-api-allele-function.md) — superseded as primary filter by ADR-0023
- [ADR-0021: Genome build is detected from position data, not file headers](0021-build-auto-detection.md)
- [ADR-0022: PharmGKB reference-genotype annotations on non-CPIC genes are documented, not filtered](0022-pharmgkb-non-cpic-reference-annotations-are-documented-not-filtered.md) — scope reduced by ADR-0023
- [ADR-0023: ClinVar REF allele is the primary PharmGKB non-finding filter](0023-clinvar-ref-as-primary-pharmgkb-filter.md)
- [ADR-0024: GWAS Catalog magnitude scoring from p-value + effect size](0024-gwas-magnitude-scoring.md)
- [ADR-0025: GRCh36 graceful degradation (no full build support)](0025-grch36-graceful-degradation.md)
- [ADR-0026: gnomAD exome frequency cache](0026-gnomad-exome-frequency-cache.md)
- [ADR-0027: AlphaMissense pathogenicity cache](0027-alphamissense-pathogenicity-cache.md)
- [ADR-0028: Local version tag convention (`local_version_tag`)](0028-local-version-tag-convention.md)
- [ADR-0029: Download integrity verification](0029-download-integrity-verification.md)
- [ADR-0030: Data source update cadence drives refresh and integrity](0030-data-source-update-cadence.md)
- [ADR-0031: Centralized license descriptors on annotator base class](0031-centralized-license-descriptors.md)
- [ADR-0032: CADD v1.7 integration](0032-cadd-integration.md)
- [ADR-0033: Schema version bump policy](0033-schema-version-bump-policy.md) — supersedes ADR-0032's schema-version stance
- [ADR-0034: Magnitude scoring scale and ceiling](0034-magnitude-scoring-scale.md)

## Writing a new ADR

1. Copy the latest ADR as a starting structure.
2. Pick the next number in sequence.
3. Status is `Accepted` when merged. Mark as `Superseded by ADR-NNNN` if a later ADR replaces the decision.
4. Add it to the index above.
5. If it supersedes an earlier ADR, update both the old ADR's status and its index entry.

Keep them short. ADRs are reading material for contributors, not novels.

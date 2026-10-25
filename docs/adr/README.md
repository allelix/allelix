# Architecture Decision Records

This directory holds the design decisions that shape Allelix.

Each ADR captures one decision: the context that forced it, what we chose, and the consequences. An ADR's decision is not rewritten. A later change is a new ADR, or a dated amendment note in the old one; the index says which.

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
- [ADR-0011: Indel-anchor protection for array-based parsers](0011-indel-anchor-protection-for-array-parsers.md); its premise that arrays do not call indels is superseded by ADR-0051
- [ADR-0012: Freshness detection via per-annotator remote signal](0012-freshness-detection-via-remote-signal.md)
- [ADR-0013: PharmGKB non-finding suppression](0013-pharmgkb-non-finding-suppression.md) — superseded by ADR-0023
- [ADR-0014: Somatic-variant suppression for germline parsers](0014-somatic-variant-suppression-on-germline.md)
- [ADR-0015: Mock data generators are the contract](0015-mock-data-as-spec.md)
- [ADR-0016: Data Classification Principle (structured fields only; regex forbidden in production)](0016-data-classification-principle.md)
- [ADR-0017: PharmGKB SNV row-level prose fallback](0017-pharmgkb-snv-prose-fallback.md) — superseded by ADR-0020
- [ADR-0018: PharmGKB per-allele function via CPIC template extraction](0018-cpic-template-allele-function-extraction.md) — superseded by ADR-0020
- ADR-0019 was reserved during drafting and not used.
- [ADR-0020: CPIC API is the structured per-allele function source](0020-cpic-api-allele-function.md) — superseded as primary filter by ADR-0023
- [ADR-0021: Genome build is detected from position data, not file headers](0021-build-auto-detection.md) — superseded in part by ADR-0037 (confidence threshold) and ADR-0044 (decision rule)
- [ADR-0022: PharmGKB reference-genotype annotations on non-CPIC genes are documented, not filtered](0022-pharmgkb-non-cpic-reference-annotations-are-documented-not-filtered.md) — scope reduced by ADR-0023
- [ADR-0023: ClinVar REF allele is the primary PharmGKB non-finding filter](0023-clinvar-ref-as-primary-pharmgkb-filter.md): one-shot migration moved to the update step by ADR-0048
- [ADR-0024: GWAS Catalog magnitude scoring from p-value + effect size](0024-gwas-magnitude-scoring.md) — effect-size modifier superseded in part by ADR-0038; auto-reingest moved to the update step by ADR-0048
- [ADR-0025: GRCh36 graceful degradation (no full build support)](0025-grch36-graceful-degradation.md)
- [ADR-0026: gnomAD exome frequency cache](0026-gnomad-exome-frequency-cache.md)
- [ADR-0027: AlphaMissense pathogenicity cache](0027-alphamissense-pathogenicity-cache.md)
- [ADR-0028: Local version tag convention (`local_version_tag`)](0028-local-version-tag-convention.md): migration moved from `is_ready()` to the update step by ADR-0048
- [ADR-0029: Download integrity verification](0029-download-integrity-verification.md)
- [ADR-0030: Data source update cadence drives refresh and integrity](0030-data-source-update-cadence.md)
- [ADR-0031: Centralized license descriptors on annotator base class](0031-centralized-license-descriptors.md)
- [ADR-0032: CADD v1.7 integration](0032-cadd-integration.md)
- [ADR-0033: Schema version bump policy](0033-schema-version-bump-policy.md) — supersedes ADR-0032's schema-version stance
- [ADR-0034: Magnitude scoring scale and ceiling](0034-magnitude-scoring-scale.md)
- [ADR-0035: `Variant.ref`, per-Annotation alt threading, and structured GWAS fields](0035-variant-ref-and-coordinated-annotation-model.md) — partially supersedes ADR-0033's strict-per-field schema-version bump policy via the cluster carve-out
- [ADR-0036: Annotated VCF output contract](0036-annotated-vcf-output-contract.md); provenance block gains `##ALLELIX_SkippedSources` (ADR-0049)
- [ADR-0037: Build detection requires a three-seed, 80% majority](0037-build-detection-confidence-threshold.md) — supersedes ADR-0021's confidence threshold
- [ADR-0038: GWAS magnitude modifiers require odds-ratio evidence](0038-gwas-odds-ratio-magnitude-modifier.md) — supersedes ADR-0024's effect-size modifier
- [ADR-0039: Locus-keyed ClinVar selection and Annotation allele-identity invariants](0039-locus-keyed-clinvar-selection-and-annotation-invariants.md) — completes ADR-0035
- [ADR-0040: Cross-build rsID identity table](0040-cross-build-rsid-identity-table.md) — restores GRCh38 enrichment for GRCh37 input through explicit translation
- [ADR-0041: Plausibility flag from zygosity and population frequency](0041-plausibility-flag.md) — marks very rare homozygous and carrier-less heterozygous calls for review; not assessed named, haploid calls and the multi-homozygote warning added by ADR-0050
- [ADR-0042: VCF quality gate](0042-vcf-quality-gate.md) — FILTER and GQ honoured, low-quality gVCF blocks are no-calls, spanning deletions symbolic, every exclusion counted
- [ADR-0043: rsID-less resolution as a shared step](0043-shared-rsid-resolution.md) — one resolver chain for analyze, stats, extract and export plink; build detection votes by position
- [ADR-0044: The genome build is decided from the file's evidence before annotation](0044-build-decision-evidence-before-labels.md) — contig lengths and whole-file seeds before labels; contradictions stop the run; one decision for every command
- [ADR-0045: ClinVar rows are scored against ClinVar's own classification of the variant](0045-clinvar-aggregate-aware-scoring.md): capped by a Conflicting or Benign aggregate, submitter split shown, outranked P/LP kept by ClinVar's tiered rule where the aggregate is not a disease call; amended by ADR-0052
- [ADR-0046: ClinVar variants without an rsID are matched for VCF input only](0046-clinvar-variants-without-rsid.md): kept with an empty rsID, matched only where REF came from the input, named by locus in every report
- [ADR-0047: VCF calls that were dropped silently are matched or counted](0047-vcf-calls-dropped-silently.md): untrimmed alleles trimmed per carried allele, half calls assessed by carrier-rule sources, VCF homozygous deletions to one base matched, dropped calls counted with their reason
- [ADR-0048: Caches are read-only at analyze time](0048-caches-read-only-at-analyze-time.md): `analyze --no-update` writes nothing to the data directory, caches open read-only and install in rollback-journal mode, one readiness state per cache with its reason; amends ADR-0023, ADR-0024 and ADR-0028
- [ADR-0049: Every report names its sources; report-level warnings](0049-sources-block-and-report-warnings.md): `sources` gives every registered source's part in the run with a reason, `report_warnings` carries report-level warnings (ClinVar skipped first), the annotated VCF names skipped sources; amends ADR-0036
- [ADR-0050: Zygosity codes, file-declared haploid calls, a named "not assessed", and the multi-homozygote warning](0050-zygosity-codes-haploid-calls-not-assessed.md): JSON zygosity and plausibility as codes, hemizygous for calls the file declares haploid, `not_assessed` with a reason, a report warning at three rare homozygous P/LP calls; amends ADR-0041
- [ADR-0051: Count the array calls ClinVar cannot match, and panel coverage in four states](0051-array-calls-clinvar-cannot-match.md): `input.probe_classes` (i-probes at P/LP loci, D/I calls at P/LP rsIDs), `input.array_quality`, `Annotator.run_counts()`, `panel_coverage` as `{requested, counts, variants}`; supersedes ADR-0011's premise
- [ADR-0052: A zero-star ClinVar P/LP submission on a common allele leaves the default report](0052-zero-star-plp-on-a-common-allele.md): a P/LP row without assertion criteria whose allele has a gnomAD frequency of 5% or more is scored 4.0, marked `common_allele_capped`, and counted; amends ADR-0045

- [ADR-0053: ClinVar refresh integrity and publication](0053-clinvar-refresh-integrity.md)
- [ADR-0054: Comparison identity and indexed extraction completeness](0054-compare-record-identity-and-indexed-extraction.md)
- [ADR-0055: Explicit VCF panel observations and constituent gene filters](0055-explicit-panel-observations-gene-filters.md): bounded coverage separate from findings, JSON panel reasons, exact source-specific gene tokens
- [ADR-0056: Protect input files and own staged output resources](0056-input-output-file-safety.md): input-alias preflight, unique owned temps and per-file publication
- [ADR-0057: Check Windows cache directory permissions without writing](0057-windows-cache-permissions.md): effective-token AccessCheck, physical paths and conservative failures
- [ADR-0058: Own analysis providers through their final use](0058-analysis-and-zip-resource-ownership.md): explicit callback owners, identity deduplication and partial ZIP cleanup
- [ADR-0059: Annotated VCF refresh and explicit reference calls](0059-annotated-vcf-refresh-and-reference-calls.md): replace old Allelix INFO on re-annotation; stamp reference-carried ClinVar findings on `0/0` records

## Writing a new ADR

1. Copy the latest ADR as a starting structure.
2. Pick the next number in sequence.
3. Status is `Accepted` when merged. Mark as `Superseded by ADR-NNNN` if a later ADR replaces the decision.
4. Add it to the index above.
5. If it supersedes an earlier ADR, update both the old ADR's status and its index entry.

Keep them short. ADRs are reading material for contributors, not novels.

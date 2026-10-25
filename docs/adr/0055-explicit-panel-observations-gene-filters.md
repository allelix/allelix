# ADR-0055: Explicit VCF panel observations and constituent gene filters

**Status:** Proposed
**Date:** 2026-10-20

## Context

Panel coverage accumulated only variants sent to annotators. VCF reference
calls and caller-filtered rows are intentionally omitted from that stream,
so reports incorrectly described explicit panel records as absent. Gene
filters compared complete source strings, missing symbols in list fields.

## Decision

Observe explicit requested VCF records separately from annotation inputs
and high-value state. Retain at most one observation per requested rsID,
with state, genotype and reason from the same first record. A requested ID
can match directly or at an unambiguous placement supplied by the existing
ready identity cache's per-build `forward_sites` result. Disagreement about
the locus on that build abstains, independently of row order or the other
build. Multiple alleles at the same locus remain supported; direct IDs do
not depend on coordinate fallback. Only the record's start is consulted;
reference blocks do not establish interior panel genotypes or presence.

Use the parser's selected sample and existing quality decisions. An explicit
confident reference call displays REF/REF using the existing panel genotype
convention. Rejected/missing calls are no-calls. Existing MT/build exclusions
are present but unusable, with an exclusion reason. These observations never
become artificial annotation inputs or high-value warnings. Shown findings
retain precedence over underlying no-call coverage. Arrays retain their
existing states, genotype conventions and HTML/terminal layout.

JSON schema 10 adds nullable `panel_coverage.variants[].reason`; arrays,
absent records and usable calls have null. Codes identify caller no-call,
FILTER rejection, low quality/zero depth, spanning deletion, symbolic or
unusable record, and existing mitochondrial exclusions. HTML adds a Reason
column for VCF only; terminal VCF coverage shows genotypes and reason text.
The diff loader accepts schemas 1 through 10 and rejects future versions.
The annotated VCF contract does not change.

Gene matching tokenizes comma/semicolon lists, trims whitespace and compares
case-insensitive exact symbols. GWAS annotation fields and queries additionally
split spaced ` - ` and ` x ` separators. Other sources do not use these
extra separators. Hyphenated symbols such as HLA-DQA1 stay whole. Requested
composite fields use the same source rules as annotations, preserving exact
whole-field query compatibility. Source text remains unchanged for display.
Gene/rsID OR semantics, magnitude/category/source floors and ordering remain.

## Limits and validation

No substring, synonym, allele, biological duplicate or interior gVCF inference
is added. First occurrence controls coverage; annotation policy is unchanged.
The broader high-value exclusion state and reference-carried genotype/ploidy
rules are separate. Offline regressions cover reference/rejected/absent calls,
quality options, duplicates, selected samples, build-specific tiny identity
data, mitochondrial exclusions, parser reuse, output reasons, actual JSON
writer-to-diff consumption, legacy/future schemas and source-specific filters.
Integration with other unmerged changes requires affected checks and review
again. No whole-genome or live-source validation is claimed here.

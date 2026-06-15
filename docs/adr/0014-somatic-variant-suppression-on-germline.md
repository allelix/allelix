# ADR-0014: Somatic-variant suppression for germline parsers

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

A subset of PharmGKB clinical annotations describes somatic variants — mutations acquired in tumor tissue, not inherited. Consumer DNA testing samples germline DNA from saliva or buccal cells. Emitting somatic annotations on germline data is a categorical type error: at best confusing, at worst suggesting the user has a tumor mutation that germline sampling cannot detect.

## Decision

PharmGKB exposes no structured germline/somatic flag at the per-row level. Per the architectural principle in ADR-0016, classification must use structured fields — regex on prose is not permitted. The common case (a germline reference user matching a somatic-context row whose description says they carry the normal allele) is correctly filtered by the non-finding rule in ADR-0023. The rare residual case — a user whose germline genotype matches a somatic-context annotation with a non-Normal allele function — surfaces in the report for manual review. If PharmGKB adds a structured germline/somatic field in a future schema revision, automated suppression can be reinstated on that field.

## Consequences

- Most somatic-context rows are filtered by the non-finding mechanism (ADR-0023) because the reference genotype on a germline file typically matches the normal-function allele.
- Rare residual somatic annotations surface for human review rather than being silently approximated by regex.
- A future VCF parser for tumor-sequencing data can set a parser-level capability flag to pass through somatic annotations in that context.

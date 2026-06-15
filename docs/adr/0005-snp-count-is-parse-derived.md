# ADR-0005: SNP count is parse-derived, not metadata

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Initial drafts of the `GenotypeParser` interface included `snp_count` in the `GenotypeMetadata` dict returned by `get_metadata()`. The MyHappyGenes implementation counted lines in a separate header scan, while `parse()` did stricter validation (skipping malformed column counts and non-integer positions).

This created a silent divergence: on a file with even one truncated line, `get_metadata()['snp_count']` and `len(list(parse()))` would disagree. The CLI displayed one number while downstream callers using `get_metadata` got another — for the same file. This is exactly the kind of inconsistency that erodes trust in a tool that handles sensitive data.

## Decision

`GenotypeMetadata` contains only header-derivable fields: `format`, `sample_id`, `build`. SNP count is removed from the metadata contract. Callers that need a count compute it from `parse()`:

```python
snp_count = sum(1 for _ in parser.parse(file_path))
```

`parse()` is the single source of truth for "how many valid variants are in this file." Anything else would have to duplicate `parse()`'s validation logic, which is the divergence we're avoiding.

The CLI's `stats` command computes its total from the parse loop directly, so the displayed `Total SNPs` always agrees with `len(list(parse()))`.

## Consequences

- `get_metadata()` stays cheap — a quick header scan, no full file pass.
- The metadata contract is honest: every field can be derived from the file's first few lines.
- Callers that want a count pay the cost of a full parse, but they get the truth.
- Tests pin this: `tests/parsers/test_myhappygenes.py::TestEdgeCases::test_metadata_count_matches_parse_count_with_malformed_lines` verifies that a file with malformed lines yields the parse-derived count, not an inflated header-line count.
- Annotators iterating `parse()` get the same count the CLI reports, so progress bars, percentages, and totals are consistent across the pipeline.

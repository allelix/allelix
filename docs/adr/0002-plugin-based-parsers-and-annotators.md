# ADR-0002: Plugin-based parsers and annotators

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Consumer DNA testing produces many incompatible file formats (MyHappyGenes/Tempus, 23andMe, AncestryDNA, LivingDNA, VCF from WGS providers). Reference annotation databases are similarly fragmented (ClinVar, dbSNP, PharmGKB, GWAS Catalog, SNPedia). New formats and databases appear regularly.

A monolithic design where every format and database is referenced from core logic would force every parser/annotator change to touch shared code, making contribution painful and review high-risk.

## Decision

Each input format is one file under `allelix/parsers/` implementing the `GenotypeParser` ABC (`name`, `display_name`, `file_extensions`, `url`, `can_parse`, `parse`, `get_metadata`). Each reference database is one file under `allelix/annotators/` implementing the `Annotator` ABC.

Parser auto-detection iterates `PARSERS` calling `can_parse` — first match wins. The user can force a parser with `--format`. The list is hand-maintained in `allelix/parsers/__init__.py`; entry-point-based plugin discovery is rejected for v0.1 because there are no third-party parsers yet.

Parsers normalize to a single internal `Variant(rsid, chromosome, position, allele1, allele2, build)`. Downstream code (annotators, reports) only sees `Variant`, never raw file formats.

## Consequences

- Adding a format = one new file + one line in `PARSERS`. The CONTRIBUTING tutorial can be a one-page walkthrough.
- The `Variant` schema is load-bearing: changes to it ripple through every annotator. Schema changes need an ADR.
- Tests live alongside the plugin (`tests/parsers/test_<name>.py`) — adding a parser without tests is visible at PR time.
- Auto-detection is order-sensitive: parsers with stricter signatures should come earlier in `PARSERS` so they match before more permissive ones.
- No premature plugin discovery framework. If/when the parser count exceeds ~10 or third parties want to ship out-of-tree parsers, revisit with a new ADR.

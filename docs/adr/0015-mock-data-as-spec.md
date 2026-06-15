# ADR-0015: Mock data generators are the contract

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Across v0.4.x and v0.5.0, four categorical clinical-safety bugs shipped despite 250+ tests and 6 rounds of code review:

- **Indel-anchor false positives** (ADR-0011): ClinVar's anchor-base indels matched single-base array readouts, emitting hundreds of false-positive "Pathogenic" calls in cancer-predisposition genes.
- **PharmGKB non-findings as findings** (ADR-0013): wild-type rows ("do not have a copy of") emitted at the same magnitude as real findings.
- **PharmGKB somatic on germline** (ADR-0014): tumor-only variants emitted on consumer DNA tests.
- **Fixture-shape regressions**: early hand-authored ClinVar/PharmGKB fixtures matched the buggy code's expectations rather than real source data shapes.

Every bug was caught the moment real source data hit the tool. Not one was caught by the test suite. The unifying failure was structural:

> **Test fixtures derived from the developer's mental model of how the system worked will pass by construction. Real data didn't match those assumptions, so real data exposed the bugs the fixtures couldn't.**

The MHG mock generator (`tests/generate_mock_data.py`) existed from v0.1.0 as the project's canonical reference for what MyHappyGenes data looks like. But ClinVar and PharmGKB used ad-hoc fixtures that matched the code's assumptions — and even the MHG generator itself contained one entry (`rs113993960 CTT/C`) that violated the MyHappyGenes format spec (single-base alleles only).

## Decision

**The mock data generators are the spec.** Three concrete rules:

1. **One generator per source.** Every input format and reference database has exactly one mock generator under `tests/`. New annotators and parsers cannot ship without an updated generator. Hand-authored ad-hoc fixtures are forbidden.

2. **Generators model real source data, not the code's expectations.** A generator entry that violates the source format (e.g. multi-base alleles in an array-based parser, single-allele ALT for an indel in ClinVar VCF, all-finding text in PharmGKB) is a generator bug. Fix the generator first, then fix any code that broke as a result.

3. **Structural invariants on generator output are tested every run.** `tests/test_mock_data_invariants.py` asserts properties of the generated output (every MHG genotype is single-base or no-call; ClinVar has both SNVs and indels; PharmGKB has carrier-findings, non-findings, and somatic rows). An invariant failure means the generator drifted from the source contract.

4. **End-to-end snapshot tests run every release.** `tests/test_end_to_end.py` runs the full `analyze` pipeline against the generators and pins specific output: counts, presence of known carriers, absence of known anti-patterns (CFTR indel-anchor leak, TP53 wild-type leak, non-finding leak, somatic leak). Snapshots can be updated, but only with documented reason in the commit message and a corresponding CHANGELOG entry.

## Consequences

- The generators become the single load-bearing source of truth. Reviewers reading a PR can verify "does this PR update the generator? does it update the snapshot? do those changes match the CHANGELOG narrative?" without re-deriving real-world facts from scratch.

- Ad-hoc fixtures used for narrow unit tests (e.g. `tmp_path.write_text("...")` in a parser test) are still allowed for testing edge cases not appropriate for the canonical generator — malformed lines, truncated files, missing columns. The rule is: such fixtures must be obviously synthetic and clearly scoped to the edge case under test.

- When a generator entry is changed because real-world behavior diverged from it (as happened with rs113993960), the change must be accompanied by:
  - A comment in the generator explaining why (with link to the incident or ADR).
  - A new test or assertion that pins the corrected behavior.
  - A CHANGELOG entry if the snapshot count changed.

- Future annotator additions (SNPedia, GWAS Catalog) must extend the relevant generator with real-shape data — including the edge cases (anchor-base indels, non-finding language, etc.) that the existing generators have learned to model.

- The end-to-end snapshot test (`test_annotation_count_snapshot`) is the gate that would have caught v0.4.2 and v0.5.0's regressions. Going forward, every annotator or parser change must pass it. A new false positive shows up as a count increase; a regression that drops a real carrier shows up as a count decrease. Neither can be "the test was wrong; let me update it" — both demand investigation.

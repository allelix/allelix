# ADR-0011: Indel-anchor protection for array-based parsers

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

ClinVar's GRCh37 VCF uses the standard VCF anchor-base convention for indels: a deletion of TT at position P is encoded as `REF=CTT ALT=C` (the C is the unchanged anchor base immediately preceding the deleted bases). Insertions work the same way in reverse: `REF=A ALT=AG`.

Array-based parsers — MyHappyGenes (Tempus), 23andMe, AncestryDNA — produce single-base genotypes at probe positions. They do not call indels; they just read whichever single base sits at the array probe.

The ClinVar annotator's carrier rule (ADR-0007) is a string-equality check: `alt in {allele1, allele2}`. For a SNV this is correct. For an indel against an array readout, it produces a categorical false-positive: ClinVar's single-character anchor base equals the array's single-character readout, the rule fires, and the user is told they carry a frameshift deletion they do not have.

In production data this affects ~95%+ of high-magnitude cancer-gene hits (MSH6, MLH1, MSH2, PTEN, RB1, BRCA1, BRCA2, TP53, APC, …). Hundreds of false-positive "Pathogenic" calls per typical array-based input. This is a clinical-safety failure that pre-v0.4.2 reports must not be trusted.

## Decision

Indel ClinVar rows (REF or ALT longer than one base) do not fire on parser inputs whose genotype is single-base only. The annotator detects this at row evaluation time:

```python
clinvar_is_indel = len(ref) > 1 or len(alt) > 1
user_is_multibase = len(allele1) > 1 or len(allele2) > 1
if clinvar_is_indel and not user_is_multibase:
    continue
```

This runs *before* the carrier rule from ADR-0007. The carrier rule still applies after this check; an indel-calling parser that reports multi-base genotypes (e.g., a future VCF parser) continues to receive indel annotations exactly as before.

This is the conservative default. We'd rather miss a hypothetical real array call of an indel — vanishingly rare, since arrays don't reliably call them — than emit a false "you carry a pathogenic frameshift" claim that almost certainly isn't true.

## Consequences

- Array-parsed genotypes no longer produce indel false positives. The cancer-predisposition gene clusters (MSH6, MLH1, MSH2, PTEN, RB1, BRCA1, BRCA2, TP53, APC) drop from hundreds of hits to whatever's actually a SNV match.
- A real array carrier of an indel is conservatively dropped. Trade-off is correct: false negatives on rare unreliable array indel calls > false positives on common syntactic coincidences.
- The regulatory posture in ADR-0003 is reinforced. We do not assert carrier status on the basis of an anchor-base coincidence.
- A future VCF parser that calls indels will deliver multi-base genotypes like CTT/C. The `user_is_multibase` check passes and indel matching works as designed (validated by the existing `test_indel_carrier_triggers` test on the CFTR ΔF508 fixture row).
- The mock ClinVar fixture's CFTR row (rs113993960 REF=CTT ALT=C) and the MyHappyGenes mock's matching CTT/C genotype already validate the positive path. Three new regression tests in `tests/annotators/test_clinvar.py::TestIndelAnchorProtection` validate the negative paths (single-base array readout, single-base homozygous reference, multi-base homozygous deletion).

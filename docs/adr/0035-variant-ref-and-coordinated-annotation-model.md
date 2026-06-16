# ADR-0035: Variant.ref, per-Annotation alt threading, and structured GWAS fields

**Status:** Accepted
**Date:** 2026-06-16
**Supersedes:** ADR-0033's strict-per-field schema-version bump, as applied to
coordinated field clusters. Per-field bumps remain in effect for non-clustered
field additions; clusters defined in an ADR share a single bump on first PR.

## Context

v2.0.1 shipped three "suppress-half" fixes (GH allelix-dev #14, #23, #24) and
deferred the matching structural-halves to v2.1. The four open v2.1 issues —
#14 (carrier match case + strand), #23 (per-Annotation alt threading for
non-ClinVar sources), #24 (structured GWAS trait / p-value / PheCode fields),
and #50 (R-1 annotator-level strand-aware carrier matching) — share one
underlying model deficiency:

- **No reference allele on `Variant`.** Strand-aware carrier matching requires
  knowing which side of the user's pair is the risk allele. Today the carrier
  rule (ADR-0007) tests set membership against `{allele1, allele2}` and the
  pipeline never carries REF through the model.
- **Per-Annotation `alt` carried for ClinVar only.** v1.4.1's #25 added the
  field but only ClinVar populates it. v2.0.1 #23 suppress-half left a gap:
  rsID-bearing GWAS / SNPedia / ClinPGx rows skip allele-specific gnomAD /
  AlphaMissense / CADD enrichment because their alt is unknown.
- **GWAS trait / p-value / PheCode live in `description` prose.** v2.0.1 #24's
  suppress-half promoted `_PHECODE_DELIM` to a shared constant; the deeper fix
  is to carry the values as structured fields and stop re-parsing prose for
  rollup / diff / display.

Landing the four issues piecemeal would mean four PRs each touching `models.py`,
the renderer set, the diff loader, and the JSON schema version. One coordinated
design ADR + four sequential implementation PRs is the cleaner shape.

## Decision

### 1. `Variant.ref` field

Add `ref: str | None = None` to the `Variant` dataclass. Carries the reference
allele on the forward strand when available:

- **VCF inputs** populate from the VCF REF column.
- **Array parsers** (MHG, 23andMe, AncestryDNA, FTDNA, MyHeritage, LivingDNA)
  populate via the existing `resolve_strand` path against gnomAD ref/alt during
  ingestion (v1.6.0 strand-normalization machinery already produces this).
- **Unknown / pre-resolution** stays `None`. Downstream consumers (strand-aware
  carrier matching, plausibility flagging, ACMG engine in v2.4) treat `None` as
  "cannot resolve" and degrade gracefully.

`Variant.ref` is optional and additive. Existing callers that don't read it
continue to behave identically.

### 2. Per-Annotation `alt` threading

`Annotation.alt` already exists (v1.4.1 #25). Coverage extends to non-ClinVar
sources:

- **GWAS Catalog loader** carries the EFO-mapped effect allele as `alt` per row
  when the catalog publishes it. Rows with no published risk-allele continue to
  emit empty `alt` and stay capped at magnitude 3.0 per ADR-0024.
- **SNPedia / ClinPGx loaders** carry per-genotype alt where the source row
  identifies a specific allele.
- **Pipeline enrichment** (`_pipeline.run_analysis` gnomAD / AlphaMissense /
  CADD attachment) uses per-Annotation `alt` for direct `(ref, alt)` lookup on
  rsID-bearing inputs. The position-fallback path (resolved-rsID on rsID-less
  inputs) keeps the v2.0.1 #23 suppress-half exactly as documented — that path
  is allele-specific via the user's actual carried alt at the resolved position.

This restores allele-specific enrichment on rsID-bearing input VCFs without
re-opening the wrong-allele safety hole #23 closed.

### 3. Structured GWAS fields on `Annotation`

Add three fields on `Annotation`:

- `trait: str = ""` — EFO trait label (e.g., "Type 2 diabetes")
- `p_value: float | None = None` — raw GWAS Catalog p-value
- `phecode: str | None = None` — PheCode rollup identifier

The GWAS loader writes these directly from upstream columns instead of
formatting them into `description` and re-parsing later. `description` keeps
its rendered prose for HTML / terminal display but is built from the structured
fields, not the source of truth.

Existing rollup readers (`_gwas_base_trait`, `_gwas_phecode_parent`) consume the
structured fields directly. The `_PHECODE_DELIM` shared constant introduced in
v2.0.1 #24 becomes dead code at the third PR and is removed.

### 4. Schema-version bump: cluster amendment to ADR-0033

ADR-0033 established: *every new field emitted in JSON report output bumps the
schema version, even if the field is optional and additive.*

This ADR adds a **cluster carve-out**: fields introduced as a documented
cluster share one schema-version bump. A cluster is a set of fields:

- Defined in advance in an ADR ("manifest"),
- All consumed by a single design decision (no field stands alone semantically),
- Landed within one minor release cycle.

The bump fires on the **first PR** in the cluster. Subsequent cluster PRs add
their fields within the bumped schema version. The ADR's manifest section is
the consumer contract — anyone reading the ADR sees every cluster field before
the first one arrives in code.

**Cluster manifest for ADR-0035** (v4 → v5):

| Field | Type | Source | Lands in PR |
|---|---|---|---|
| `variant.ref` | `str \| None` | Variant model | 1 |
| `annotation.trait` | `str` | Annotation model | 3 |
| `annotation.p_value` | `float \| None` | Annotation model | 3 |
| `annotation.phecode` | `str \| None` | Annotation model | 3 |

`annotation.alt` is unchanged (already present since v1.4.1). PR 2 extends its
population coverage, not the schema.

ADR-0033's spirit ("don't silently change the schema") is preserved: the change
is documented in advance, consumers see one named bump, and the cluster
manifest names every field that will appear within that bump.

Non-clustered field additions continue under ADR-0033's per-field policy.

## Migration plan — four sequential PRs

Each PR is independently shippable and reverts cleanly without entangling the
others.

### PR 1 — `Variant.ref` field + schema bump v4 → v5

- Add `Variant.ref: str | None = None` to `models.py`.
- Bump `SCHEMA_VERSION = "5"` in `json_report.py`.
- Extend `_SUPPORTED_SCHEMA_VERSIONS` in the diff loader to include both "4"
  and "5".
- Populate `ref` in the VCF parser (REF column) and array parsers (via existing
  `resolve_strand` against gnomAD ref/alt).
- Emit `ref` in the JSON output's embedded variant data (HTML report uses the
  same blob; field is optional so v4 consumers ignore it).
- Document the cluster manifest in the bump's docstring referencing this ADR.
- No behavior change at downstream consumers — they opt in field-by-field as
  PRs 2-4 land.

**Closes part of #49** (this ADR + PR 1 = #49 done). Tests add fixtures
covering `ref`-populated paths.

### PR 2 — Per-Annotation `alt` threading

- GWAS / SNPedia / ClinPGx loaders populate `alt` from upstream data where
  available.
- Pipeline enrichment uses per-Annotation `alt` for direct `(ref, alt)` lookup
  on rsID-bearing inputs.
- Position-fallback path for resolved rsIDs unchanged.
- Wrong-allele safety regression tests pin both paths.

**Closes #23 structural-half.**

### PR 3 — Structured GWAS fields

- Add `Annotation.trait`, `Annotation.p_value`, `Annotation.phecode`.
- GWAS loader writes structured fields directly.
- Rollup readers consume structured fields.
- `description` is rendered from the structured fields, not the source of
  truth.
- `_PHECODE_DELIM` removed.

**Closes #24 structural-half.**

### PR 4 — Strand-aware carrier matching

- `ClinVarAnnotator.annotate()` / `PharmGKBAnnotator.annotate()` / `compare`
  command consult `complement(user_allele)` against `{ref, alt}` using the
  `Variant.ref` field PR 1 introduced.
- Palindromic A/T and C/G SNPs return "ambiguous" (already handled in
  `compare` since v1.1; this PR extends the policy to annotator carrier rules).
- Mandatory test fixture: rs5742904 forward C/C and coding G/G produce the
  same set of annotations (per CLAUDE.md R-1 spec).

**Closes #14 strand-half and #50 (R-1).**

## Consequences

- Schema v4 consumers see v5 once PR 1 lands. The diff loader carries both v4
  and v5 baselines forward (additive only — every v4 field is unchanged).
- Existing JSON baselines stay diff-compatible (new fields are optional with
  `None` / empty defaults).
- Annotation tests gain new fixtures covering `ref`-populated paths, structured
  GWAS fields, and strand-aware carrier matching.
- ADR-0033's per-field policy stays in effect for non-clustered field
  additions. Future cluster ADRs (e.g., JSON v3 in v2.3 — see GH allelix-dev
  #53) follow the same pattern: cluster manifest in the ADR, single bump on
  first PR.
- `Annotator.__del__` removal (#36) is independent of this ADR and can land in
  parallel without merge conflicts.

## References

- ADR-0007 (genotype matching requires alt allele — strand-aware extension
  here)
- ADR-0010 (strand-flip helpers shipped v0.4.0; liftover deferred — unchanged
  by this ADR)
- ADR-0024 (GWAS magnitude scoring — interacts with PR 2: unknown-risk-allele
  cap at 3.0 still fires when `alt` is empty)
- ADR-0033 (schema version bump policy — cluster carve-out documented above)
- GH allelix-dev #14, #23, #24, #50 (the four implementation issues)
- CLAUDE.md "Feature R-1" (strand-aware carrier matching spec)

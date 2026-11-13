# ADR-0035: Variant.ref, per-Annotation alt threading, and structured GWAS fields

**Status:** Accepted
**Date:** 2026-06-16
**Supersedes:** ADR-0033's strict-per-field schema-version bump, as applied to
coordinated field clusters. Per-field bumps remain in effect for non-clustered
field additions; clusters defined in an ADR share a single bump on the first step.

## Context

v2.0.1 shipped three "suppress-half" fixes and
deferred the matching structural-halves to v2.1. The four open v2.1 problems —
carrier-match case + strand, per-Annotation alt threading for
non-ClinVar sources, structured GWAS trait / p-value / PheCode fields,
and annotator-level strand-aware carrier matching — share one
underlying model deficiency:

- **No reference allele on `Variant`.** Strand-aware carrier matching requires
  knowing which side of the user's pair is the risk allele. Today the carrier
  rule (ADR-0007) tests set membership against `{allele1, allele2}` and the
  pipeline never carries REF through the model.
- **Per-Annotation `alt` carried for ClinVar only.** v1.4.1's fix added the
  field but only ClinVar populates it. v2.0.1 suppress-half left a gap:
  rsID-bearing GWAS / SNPedia / ClinPGx rows skip allele-specific gnomAD /
  AlphaMissense / CADD enrichment because their alt is unknown.
- **GWAS trait / p-value / PheCode live in `description` prose.** v2.0.1's
  suppress-half promoted `_PHECODE_DELIM` to a shared constant; the deeper fix
  is to carry the values as structured fields and stop re-parsing prose for
  rollup / diff / display.

Fixing the four problems piecemeal would mean four changes each touching `models.py`,
the renderer set, the diff loader, and the JSON schema version. One coordinated
design ADR + four sequential implementation steps is the cleaner shape.

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

`Annotation.alt` already exists (v1.4.1). Coverage extends to non-ClinVar
sources:

- **GWAS Catalog loader** carries the EFO-mapped effect allele as `alt` per row
  when the catalog publishes it. Rows with no published risk-allele continue to
  emit empty `alt` and stay capped at magnitude 3.0 per ADR-0024.
- **SNPedia / ClinPGx loaders** carry per-genotype alt where the source row
  identifies a specific allele.
- **Pipeline enrichment** (`_pipeline.run_analysis` gnomAD / AlphaMissense /
  CADD attachment) uses per-Annotation `alt` for direct `(ref, alt)` lookup on
  rsID-bearing inputs. The position-fallback path (resolved-rsID on rsID-less
  inputs) keeps the v2.0.1 suppress-half exactly as documented — that path
  is allele-specific via the user's actual carried alt at the resolved position.

This restores allele-specific enrichment on rsID-bearing input VCFs without
re-opening the wrong-allele safety hole the suppress-half closed.

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
v2.0.1 becomes dead code at the third step and is removed.

### 4. Schema-version bump: cluster amendment to ADR-0033

ADR-0033 established: *every new field emitted in JSON report output bumps the
schema version, even if the field is optional and additive.*

This ADR adds a **cluster carve-out**: fields introduced as a documented
cluster share one schema-version bump. A cluster is a set of fields:

- Defined in advance in an ADR ("manifest"),
- All consumed by a single design decision (no field stands alone semantically),
- Landed within one minor release cycle.

The bump fires on the **first step** in the cluster. Subsequent cluster steps add
their fields within the bumped schema version. The ADR's manifest section is
the consumer contract — anyone reading the ADR sees every cluster field before
the first one arrives in code.

**Cluster manifest for ADR-0035** (v4 → v5):

| Field | Type | Source | Lands in step |
|---|---|---|---|
| `variant.ref` | `str \| None` | Variant model | 1 |
| `annotation.trait` | `str` | Annotation model | 3 |
| `annotation.p_value` | `float \| None` | Annotation model | 3 |
| `annotation.phecode` | `str \| None` | Annotation model | 3 |

`annotation.alt` is unchanged (already present since v1.4.1). Step 2 extends its
population coverage, not the schema.

ADR-0033's spirit ("don't silently change the schema") is preserved: the change
is documented in advance, consumers see one named bump, and the cluster
manifest names every field that will appear within that bump.

Non-clustered field additions continue under ADR-0033's per-field policy.

## Migration plan — four sequential steps

Each step is independently shippable and reverts cleanly without entangling the
others.

### Step 1 — `Variant.ref` field + schema bump v4 → v5

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
  steps 2-4 land.

**Outcome:** the model carries the reference allele, the prerequisite for
the strand part of carrier matching in step 4. Tests add fixtures covering
`ref`-populated paths.

### Step 2 — Per-Annotation `alt` threading

- GWAS / SNPedia / ClinPGx loaders populate `alt` from upstream data where
  available.
- Pipeline enrichment uses per-Annotation `alt` for direct `(ref, alt)` lookup
  on rsID-bearing inputs.
- Position-fallback path for resolved rsIDs unchanged.
- Wrong-allele safety regression tests pin both paths.

**Outcome:** completes the alt-threading fix that v2.0.1 started by
removing the wrong-allele fallback for alt-less rows.

### Step 3 — Structured GWAS fields

- Add `Annotation.trait`, `Annotation.p_value`, `Annotation.phecode`.
- GWAS loader writes structured fields directly.
- Rollup readers consume structured fields.
- `description` is rendered from the structured fields, not the source of
  truth.
- `_PHECODE_DELIM` removed.

**Outcome:** completes the PheCode fix that v2.0.1 started by sharing
one delimiter constant.

### Step 4 — Strand-aware carrier matching

- `ClinVarAnnotator.annotate()` / `PharmGKBAnnotator.annotate()` / `compare`
  command consult `complement(user_allele)` against `{ref, alt}` using the
  `Variant.ref` field step 1 introduced.
- Palindromic A/T and C/G SNPs return "ambiguous" (already handled in
  `compare` since v1.1; this step extends the policy to annotator carrier rules).
- Mandatory test fixture: rs5742904 forward C/C and coding G/G produce the
  same set of annotations.

**Outcome:** completes carrier matching across strands; the case part
was fixed in v2.0.1.

## Consequences

- Schema v4 consumers see v5 once step 1 lands. The diff loader carries both v4
  and v5 baselines forward (additive only — every v4 field is unchanged).
- Existing JSON baselines stay diff-compatible (new fields are optional with
  `None` / empty defaults).
- Annotation tests gain new fixtures covering `ref`-populated paths, structured
  GWAS fields, and strand-aware carrier matching.
- ADR-0033's per-field policy stays in effect for non-clustered field
  additions. Future cluster ADRs (e.g., JSON v3 in v2.3) follow the same
  pattern: cluster manifest in the ADR, single bump on the first step.
- `Annotator.__del__` removal is independent of this ADR and can land in
  parallel without merge conflicts.

## Amendment 2026-10-16 (2.5): no complemented REF in ClinVar

ClinVar accepted a call whose REF was the complement of the row's REF, for
a "reverse-stranded" input. That never happens on a correct file:

- A VCF is forward strand by specification.
- An array call's REF comes from the rsID identity table at the same
  forward locus as ClinVar's REF, so the two are equal.

A complemented REF therefore meant mis-coordinated input, and accepting it
could attach a finding to the wrong variant. ClinVar now requires the
call's REF to equal the row's, and the carrier check is direct membership.
A one-base VCF REF that is none of ClinVar's one-base REFs at the position
is counted (`input.vcf_quality.clinvar_ref_mismatches`), and the report
asks the user to check the file's build and coordinates.
`strand_aware_carrier_match` remains for ClinPGx.

## References

- ADR-0007 (genotype matching requires alt allele — strand-aware extension
  here)
- ADR-0010 (strand-flip helpers shipped v0.4.0; liftover deferred — unchanged
  by this ADR)
- ADR-0024 (GWAS magnitude scoring — interacts with step 2: unknown-risk-allele
  cap at 3.0 still fires when `alt` is empty)
- ADR-0033 (schema version bump policy — cluster carve-out documented above)

## Amendment 2.6.0 (2026-11-03): ClinPGx's strand helper

ClinPGx does not use `strand_aware_carrier_match`. It matches the
diploid call with `strand_aware_genotype_match` (ADR-0009). No annotator
calls `strand_aware_carrier_match`; it stays in `allelix.utils.allele`
with its tests.

## Amendment 2.7.0 (2026-11-09): ClinPGx's complement match takes site evidence

`strand_aware_genotype_match` accepted a ClinPGx complement match only when
the call's REF was absent from the row's genotype and the REF's complement
was in it. A homozygous-ALT row on the reverse strand never holds the REF's
complement: at a GRCh37 G>A site with rows CC, CT and TT, a G/A call met CT
but an A/A call did not meet TT.

The complement reading is now also accepted on site evidence. All of the
rsID's ClinPGx genotypes, non-finding rows included, are read together
against the identity table's alleles at the call's locus on the call's
build (`forward_strand_reading`, ADR-0066). Where only their complement
fits, a single-base call there matches by its complement, homozygotes
included (`_reads_complement` in `allelix/annotators/pharmgkb.py`). A
single-base call is one whose two carried alleles are single bases and
whose REF is one base or not known. The site's identity alleles include
the letters of indel rows: some deletion ALTs are just an anchor base. An
indel call can likewise have a single-base allele, so a letter of the call
that fits the site's alleles need not identify the same sequence variant;
site evidence says nothing about it. The evidence is
decided per call: in a batch, a call at another locus does not take it from
a call of the same rsID that is at the site. A call without a position (an
unplaced rsID-list row, ADR-0069) is read against the rsID's one site and
keeps no locus.

- The exact match still comes first and is unchanged.
- No new complement match is made where the identity table is unavailable
  or switched off, does not hold the rsID, places it at another locus or at
  more than one, or flags it palindromic, or where the source alleles fit
  as written (`direct`), both ways (`palindromic`, `both_fit`) or neither
  way (`no_fit`), where a genotype is not two single bases, or for a call
  that is not single-base (an indel call, or a REF longer than one base).
  The REF rule above still applies there, unchanged.
- A non-finding row is still never reported; its alleles only count as
  evidence of the site's orientation.
- The non-finding filters (ADR-0020, ADR-0023) are unchanged.

## Amendment 2.7.0 RC2 (2026-11-10): one reading per ClinPGx call

The 2.7.0 amendment above kept "the exact match still comes first" beside
the new site evidence, so one call could meet two genotypes. At a
three-allele site whose rows read only as their complement (GRCh37 G>A,C
with rows CC, CT, TT, GG), a G/G call met CC by the complement and GG as
written, the genotype of the other homozygote. With no site evidence, a
G/A call with REF G met AG as written and CT by its REF.

A call now meets ClinPGx genotypes by one reading (`_genotypes_met` in
`allelix/annotators/pharmgkb.py`):

- Where the site decides the reading (`_site_reading`: all of the rsID's
  genotypes read against the identity alleles at the call's site), that is
  the one reading. Where they fit only as their complement (the amendment
  above), a row spelled as the call is not matched
  (`strand_aware_genotype_match` with `site_reads_complement`). Where they
  fit only as written, only the exact match counts: the call's REF does
  not add a complement match (a G>T call at a G>A,C site does not meet the
  AC genotype).
- Elsewhere (no identity table, the rsID not placed at the call or placed
  more than once, both readings or neither fitting, an indel call) the
  exact match and the REF rule apply as before, over every
  genotype of the rsID, non-finding rows included. Where one genotype
  meets the call as written and another as its complement, the call meets
  neither: nothing decides which reading is the person's, and preferring
  the exact spelling would guess. Each such call is counted once
  (`run_counts()["ambiguous_strand_calls"]`) and the report says so in
  `report_warnings` (`clinpgx_ambiguous_strand`).
- A palindromic site, an indel call and the non-finding filters are
  unchanged. Other sources keep their own readings (ADR-0066).

# ADR-0036: Annotated VCF Output Contract

**Status:** Accepted
**Date:** 2026-07-02

## Context

v2.3.0 headline feature (GH allelix-dev #150) adds an annotated VCF
output mode alongside the existing HTML and JSON reports. Bioinformatics
users work in VCF-in / VCF-out pipelines and won't rearchitect around
HTML or a bespoke JSON schema — an annotated VCF is what slots into
their existing tooling and positions Allelix as a lightweight,
deterministic alternative to VEP.

The design forces three architectural commitments the codebase
hasn't made before:

- **A third output format** with its own governance surface (INFO field
  names, provenance line format, attribution format, schema version).
- **A two-pass write** — the existing analyze pipeline is genotype-aware
  (parser only yields ALTs the sample carries; drops ref/ref sites)
  which is correct for DTC reporting but wrong for pipeline-integration
  users who need every input line preserved.
- **A distinct schema-version namespace** from the JSON report. ADR-0033
  governs `SCHEMA_VERSION` in `allelix.reports.json_report`; that
  policy does not extend to the VCF writer's contract. Coupling them
  would force a JSON bump for every VCF-only tweak (and vice versa),
  training reviewers to ignore schema bumps as noise.

## Decision

### Two-pass writer architecture

The pipeline runs pass 1 unchanged: existing parser + annotator chain
produces annotations for the ALTs the sample carries. The pipeline
groups those annotations into a `(chrom, pos, ref, alt) → [Annotation]`
lookup dict (plus an optional `(...) → recovered_rsid` map). The writer
performs pass 2: opens the source VCF, streams line by line, preserves
every header and blank line verbatim, injects the provenance /
attribution / `##INFO` block immediately before `#CHROM`, and stamps
INFO fields onto each data line using the lookup.

Consequences of this shape:

- **Parser is untouched.** DTC report semantics (genotype-aware
  yielding, symbolic-ALT filtering, gVCF handling) are unchanged.
- **Multi-sample VCFs pass through for free.** Sample columns and
  FORMAT fields are byte-for-byte preserved because the writer only
  rewrites columns 3 (ID) and 8 (INFO).
- **Non-carried ALTs are covered.** At a het site (`0/1` on `C→A,T`)
  the ALT the sample doesn't carry still appears in the output VCF
  with its ALT-column slot; annotations for both ALTs are looked up
  independently. The parser dropping non-carried ALTs from DTC
  reports is orthogonal.
- **gVCF reference blocks (`END=…`) pass through unchanged.** Symbolic
  ALTs (`<NON_REF>`, `<DEL>`, `<*>`, `<CNV>`) skip annotation lookup
  and preserve their column slot.
- **Cost model is honest.** Annotation cost is unchanged from analyze;
  the second pass adds a sequential file re-read plus dict lookups —
  I/O bound, not compute bound. No `<Nx current analyze time`
  performance commitment is offered.

### Per-allele `Number=A` semantics

Every Allelix-emitted INFO field is declared `Number=A` in the VCF
header. Multi-allelic rows emit one comma-joined value per ALT in
column order; empty slots render as `.`. This is the VCF-standard
per-allele shape (matches how gnomAD, VEP, and bcftools consume
allele-tagged annotations) and it's the failure mode we can't afford:
if the writer picked `Number=1` (joint) semantics, multi-allelic sites
would silently drop the second ALT's annotation and bcftools would
still validate the file.

### Provenance and attribution header block

Injected immediately before `#CHROM`, structured to be parseable by any
VCF-conformant reader:

```
##ALLELIX_Version=<allelix version string>
##ALLELIX_VCF_SchemaVersion=<VCF_SCHEMA_VERSION constant>
##ALLELIX_License=AGPL-3.0-or-later
##ALLELIX_RunDate=<ISO 8601 UTC timestamp>
##ALLELIX_Build=<GRCh37 | GRCh38>
##ALLELIX_Database=<Name=<annotator name>,Version=<pinned version>>
##ALLELIX_Attribution=<Name=<display name>,License=<SPDX>,URL=<source URL>>
```

`##ALLELIX_Database` is emitted per entry in `annotators_used`.
`##ALLELIX_Attribution` mirrors the `license_attributions` block from
`allelix.reports.html._license_attributions` and
`allelix.reports.json_report._license_attributions` — three output
formats, three matching lists. Extraction into a shared helper waits
for a fourth call site per project convention (Rule of Three plus one).

### rsID recovery: ID column stamp with `ALLELIX_ORIGINAL_ID` provenance

When the pipeline recovers an rsID for a variant that arrived with
`.` in the ID column (via ClinVar or gnomAD position lookup — see
GH allelix-dev #128, #129), the writer stamps the recovered rsID into
the output ID column and stashes the pre-recovery value under
`ALLELIX_ORIGINAL_ID` in INFO. `ALLELIX_ORIGINAL_ID` is declared only
when the caller passes a non-empty recovered-rsID map; a file with no
recovery activity has no declaration.

### `VCF_SCHEMA_VERSION`: distinct namespace, experimental 0.x

The VCF output contract's schema version is a semver string starting
at `0.1.0` (`allelix.reports.vcf.VCF_SCHEMA_VERSION`). It is **not**
coupled to the JSON report `SCHEMA_VERSION` (ADR-0033) — a JSON-only
change does not bump the VCF version and vice versa.

**Bump rule.** Any change to the following surfaces bumps
`VCF_SCHEMA_VERSION`:

- Adding, removing, or renaming an Allelix-emitted `##INFO=<>` field
- Changing the pipe-delimited value layout of an existing INFO field
- Changing the `##ALLELIX_*` provenance line format (keys, structure,
  value grammar)
- Changing the `##ALLELIX_Attribution=<>` field grammar
- Changing `ALLELIX_ORIGINAL_ID` semantics or its declaration gate

While the schema is 0.x it is explicitly experimental. Reviewers may
approve INFO-shape churn without ceremony because 0.x conveys that
churn is expected. The 1.0.0 graduation is a follow-up ADR after real
user feedback per the acceptance criteria on GH allelix-dev #150; it
will lock the contract and revoke the 0.x churn freedom.

### Malformed-row tolerance

Data rows with fewer than 8 columns or a non-integer POS pass through
unchanged rather than raising. The writer never fails a whole file on
one bad row; downstream `bcftools view` will flag the malformed row
if the user cares. This matches how the existing VCF parser handles
partial-line garbage in production files.

## Consequences

- **New governance surface: `VCF_SCHEMA_VERSION` bump discipline.**
  Contributors touching `allelix/reports/vcf.py`'s INFO / provenance /
  attribution shape must bump the constant and note the change in
  the PR body. Cross-referenced from `CONTRIBUTING.md` § Release-content
  discipline.
- **JSON schema (ADR-0033) and VCF schema evolve independently.** No
  fan-out obligation across formats; each contract is versioned on
  its own axis.
- **`_license_attributions` duplication reaches three call sites.**
  `allelix/reports/html.py`, `allelix/reports/json_report.py`,
  `allelix/reports/vcf.py`. Extraction into a shared helper is
  explicitly deferred until a fourth call site materializes (Rule of
  Three plus one — no premature abstraction).
- **Parser and `AnalysisResult` remain untouched.** The pipeline hook
  in the follow-up PR builds the annotation dict inside the existing
  `_flush()` path without model changes.
- **Multi-sample and gVCF users get correct behavior without extra
  code paths.** Two-pass writer sidesteps the genotype-aware parser
  entirely for the write.
- **`bcftools view` becomes a CI gate.** The follow-up PR wires
  `bcftools view` validation into the CI matrix; adds a new tooling
  dependency but is the standard round-trip check.

## Explicit non-goals for v2.3.0

- **BCF binary output** (via `pysam`/`htslib`). Deferred until real
  demand appears; text VCF is what pipelines default to.
- **`.vcf.gz` compressed output.** Same rationale as BCF — deferred.
- **Consequence annotation** (missense_variant, synonymous_variant,
  etc. per SO terms; GH allelix-dev #59). Real VEP-replacement pitch
  eventually includes this; v2.3.0 ships annotation without
  consequence prediction and pitches on that basis honestly.
- **Full VCF 4.3 §1.6 percent-encoding of INFO values.** The v0.1
  writer substitutes underscore for `;`, `,`, `=`, whitespace — safe
  minimum. Percent-encoding lands when a real user file surfaces text
  that round-trips poorly through underscore substitution.

## Related

- ADR-0033 — JSON schema version bump policy (parallel discipline,
  independent contract)
- ADR-0031 — Centralized license descriptors (the source of truth for
  `##ALLELIX_Attribution` values)
- ADR-0035 — Variant.ref + per-Annotation alt threading (the field
  the writer keys on for per-allele lookup)
- GH allelix-dev #150 — feature issue for annotated VCF output
- GH allelix-dev #128, #129 — rsID recovery from ClinVar / gnomAD
  positions (the source of the writer's recovered-rsID input)
- GH allelix-dev #59 — variant consequence annotation (explicit
  non-goal for v2.3.0)

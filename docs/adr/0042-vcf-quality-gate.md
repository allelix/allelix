# ADR-0042: VCF quality gate

**Status:** Accepted
**Date:** 2026-10-10
**Builds on:** ADR-0036 (annotated VCF output contract), ADR-0040 (identity
table placements).

## Context

The VCF parser read only CHROM, POS, ID, REF, ALT and the sample's GT.
Three consequences:

- Rows the variant caller had rejected (FILTER `LowQual`, a VQSR tranche,
  a hard filter) were reported exactly like passing calls. The annotated
  VCF also added Allelix INFO to them, which other tools would read as an
  Allelix-endorsed call.
- A gVCF reference block with no reads (GQ 0, DP 0) was treated like a
  confident reference call. The report could not say that the genome was
  blind at a high-value position, and `high_value_no_calls` was always
  empty on VCF input.
- The VCF 4.3 spanning-deletion allele `*` passed through as an allele,
  giving genotypes like `A/*` and a `*` in PLINK `.bim` files.

## Decision

### 1. One gate, one parse path

All VCF reading goes through `VcfParser.parse` / `_parse_record`. `analyze`,
`stats`, `extract` and `export plink` therefore apply the same rules.
`analyze` exposes the two settings as flags.

- **FILTER.** A row whose FILTER is not `PASS` or `.` is skipped and
  counted, once per label. `--keep-filtered` keeps such rows.
- **Genotype quality.** A call whose GQ is below `--min-gq` (default 20,
  GATK's low-confidence convention; 0 turns the GQ floor off), or whose DP
  is 0, is a no-call. The comparison is strict: GQ 19 is a no-call and
  GQ 20 is not. A missing GQ or DP (key absent, or `.`) never demotes a
  call, because benchmark VCFs often carry GT only. QUAL is not consulted:
  it is a site score that FILTER already acts on.
- **The same rule for reference blocks.** A `0/0` record or gVCF block
  below the GQ floor, or with DP 0, is a no-call for every position it
  covers. A block covers POS through its INFO `END`; any other record
  covers its REF.
- **Spanning deletion.** `*` is symbolic, like `<NON_REF>`. A call with
  `*` on either haplotype yields no Variant. Calls with `*` on exactly one
  haplotype are counted and reported as "not assessed". Allelix has no
  genotype shape for one called allele with the other haplotype deleted;
  that is follow-up work. PLINK export writes only A/C/G/T alleles.

### 2. Nothing excluded silently

`VcfParser.quality_stats` records the following, and `AnalysisResult.vcf_quality`
carries them to every report:

- filtered rows (by label);
- GQ-demoted calls;
- one-haplotype `*` calls;
- whether the file is a gVCF.

They are shown as CLI notes, as an "Input quality" notice in HTML, and as
`input.vcf_quality` / `input.gvcf` in JSON.

### 3. High-value positions on VCF input

Reference blocks carry no rsIDs, so the pipeline asks the parser to watch
the high-value rsIDs and their identity-table positions on both builds.
The effective build is not known until streaming ends. Each watched
position gets the strongest state among the records covering it:

- a variant call;
- a no-call;
- a confident reference call.

A no-call carries a reason:

- `no_call`: the caller's own `./.`;
- `filtered`;
- `low_gq_call`;
- `low_gq_block`;
- `spanning_deletion`;
- `symbolic_allele`;
- `no_coverage`: a gVCF position that no record covers.

On a plain VCF, a position without a record is unknown, not reference and
not a no-call. It stays off the list, and the report says once that
positions without a record are unknown.

Watching costs little. In a position-sorted file the records covering a
position start at or before it. The parser therefore holds only the rows
that can still cover the current position:

- the latest start;
- earlier variant records whose REF still reaches it.

It checks a held row's span only when the stream passes a watched
position. An unsorted file falls back to checking every row against every
watched position.

### 4. Output contracts

- **Annotated VCF.** A caller-rejected row is copied unchanged and never
  stamped, unless `--keep-filtered` is set. Then it is stamped like any row
  and carries `ALLELIX_FILTER=<label>`. The new INFO field moves
  `VCF_SCHEMA_VERSION` to 0.2.0 under ADR-0036.
- **JSON.** `input.gvcf`, `input.vcf_quality` and the `reason` on each
  `high_value_no_calls` entry ship in schema 7, together with the
  plausibility field from ADR-0041. Schema 7 has not been released, so one
  bump covers all of these additions (ADR-0033).

## Consequences

- Filtered whole-genome VCFs no longer turn the caller's own rejected
  calls into findings. The report shows how many rows were left out and
  why.
- A gVCF can answer the no-call question precisely at high-value
  positions. A plain VCF says honestly that it cannot.
- High-value watching needs the identity table to place rsIDs that the
  file does not name. Without it, only records that carry the rsID are
  classified.
- Parsing a gVCF while watching costs roughly a third more CPU than
  before. Parsing without watching is slightly faster than before, because
  reference blocks that nobody watches are skipped before their sample
  column is read.

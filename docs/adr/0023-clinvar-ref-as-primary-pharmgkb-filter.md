# ADR-0023: ClinVar REF allele is the primary PharmGKB non-finding filter

- **Date:** 2026-05-15
- **Status:** Accepted (supersedes ADR-0020's role as primary; reduces ADR-0022's scope)

## Context

Earlier PharmGKB non-finding filter iterations (prose regex, structured column, hybrid patterns, template extraction, CPIC API join) each required the upstream source to have classified the user's allele in a specific vocabulary. CPIC's coverage is structurally heterogeneous across genes; no extension to the join can make it complete. Genes where CPIC departs from the expected vocabulary leak through:

- **CFTR**: 104 alleles classified as `"ivacaftor responsive"` / `"ivacaftor non-responsive"`. Zero `Normal function` entries. The v0.7.1 filter has no suppression coverage for any CFTR rsID. Real-world output: ~30+ CFTR × ivacaftor reference-homozygote rows leaked into a v0.7.3 report at `--min-magnitude 5`.
- **MTHFR, F2, F5, and others**: CPIC has no allele table for the gene at all. ADR-0022 documented this as a "small upstream gap."

Each prior release fixed one shape of the leak and shipped without checking the rest of the PharmGKB output. The recurring failure mode is implicit: every fix required CPIC to have classified the user's allele as `Normal function`. CPIC's coverage is structurally heterogeneous; no extension to the join can make it complete.

The right structural source for the inclusion/exclusion decision is **the reference allele itself**, not its function classification. ClinVar's per-build VCFs publish REF for every variant they catalog. Coverage is broad, the vocabulary is uniform (a single A/C/G/T base for SNVs), and the data is already in Allelix's per-build SQLite caches from ADR-0021's dual-build work.

## Decision

**The PharmGKB non-finding filter's primary check is: does the user's genotype match ClinVar's REF allele homozygously? If yes, the row is a non-finding and is suppressed. CPIC's per-allele function classification is demoted to a secondary tier used only when ClinVar has no usable REF data for the rsid.**

### The filter, in order

1. **Primary — ClinVar REF carrier rule.** For the row's `(rsid, build)`, look up the single-base REF in `clinvar.<build>.sqlite`. If REF is found and BOTH of the user's alleles equal REF, the user is homozygous reference → suppress the row. This is universal: it applies to every gene ClinVar catalogs, regardless of CPIC vocabulary or coverage.

2. **Secondary fallback — CPIC join.** If ClinVar has no single-base REF for the rsid (rsid not in ClinVar, or only indel REFs cataloged), fall through to the cache's pre-computed `is_nonfinding` flag set at load time from CPIC's `(rsid, base) → function_class` lookup (ADR-0020). Tiny residual; covers the rare rsid PharmGKB knows about that ClinVar doesn't.

3. **No data either way.** Emit the row. The user sees the annotation and decides; we never silently suppress without structured evidence.

### Wiring

- `ClinVarAnnotator.reference_for(rsid, build) -> str | None` exposes the per-rsid REF lookup, lazily built from the per-build cache (`SELECT DISTINCT rsid, ref FROM clinvar_variants WHERE length(ref) = 1`).
- `PharmGKBAnnotator.__init__` accepts an optional `clinvar_ref_provider` callable.
- `get_annotators()` constructs ClinVarAnnotator first and passes `clinvar.reference_for` to PharmGKBAnnotator.
- The build dispatched on is `variant.build` — set by the auto-detection pipeline (ADR-0021) per file. ClinVar's per-build cache supplies the REF appropriate to the user's data.

### CPIC's residual role

CPIC stays useful for two things:

- **Tier 2 suppression fallback** for rsids ClinVar doesn't have.
- **Annotation enrichment** — when a real carrier emits, CPIC's `function_class` (Normal / Decreased / No / Increased) can be surfaced as supplementary information labeling *which* function class the variant the user carries. This is a separate UX feature, not a filter input.

CPIC is NOT the inclusion/exclusion gate any more. The decision of "does this row appear in the report" is the ClinVar REF check; CPIC's function classifications inform the explanation, not the filter.

### Genotype display

A second, smaller change ships in the same release: both annotators set `Annotation.genotype_match` to the user's actual diploid genotype (e.g., `"AG"`, `"GG"`, `"AA"`) instead of the matched ALT base (which is what ClinVar previously stored, producing single-letter displays like `"A"`). The diploid representation is sorted for SNVs and verbatim for indels (`"CTT/C"`). Reports now show a consistent genotype column across ClinVar and PharmGKB rows.

## Consequences

- **CFTR's ~30+ "do not have a copy of the variant" rows disappear** from real-data output. ClinVar has REF for every CFTR rsID we've seen in the v0.7.3 leak; the homozygous-reference check fires cleanly.
- **MTHFR, F2, F5 reference-genotype residuals also disappear** for rsids ClinVar knows about. ADR-0022's scope is now meaningfully smaller: it covers only rsids where neither ClinVar nor CPIC has data. We expect this to be a true tiny residual.
- **Known-carrier emissions are unchanged.** rs1801133 AG, rs1799853 CT, rs4244285 GG (MTHFR, CYP2C9 and CYP2C19 carriers) all still emit because the user is not homozygous reference per ClinVar.
- **The v0.7.1 CPIC join survives as the secondary tier.** Code is not removed; the cache's `is_nonfinding` flag is still pre-computed at load time. The primary tier supersedes it but the secondary path catches whatever ClinVar misses.
- **Reports show consistent genotype columns.** ClinVar rows previously showed `"A"` (the matched ALT); they now show `"AG"` (the user's diploid). PharmGKB rows already showed `"AG"`.
- **ADR-0022 is reduced in scope.** It still applies, but only to rsids where both ClinVar and CPIC lack data. Annotations on those rsids continue to surface honestly with the README note.
## ClinVar interpreter version stamp (amended 2026-05-22)

The ClinVar annotator's interpretation logic — `_CLNSIG_MAGNITUDE` map, carrier
rule, indel-anchor protection, benign suppression — can change between Allelix
releases while the upstream ClinVar VCF data stays the same. When this happens,
existing caches produce stale results: the data is correct but the
interpretation is wrong for the current code.

**Decision.** `CLINVAR_INTERPRETER_VERSION` (integer, in
`allelix/databases/_versions.py`) is stamped into each per-build cache's
`database_versions.remote_signal` as `|iv:N` during `load_clinvar_vcf()`.
`ClinVarAnnotator.is_ready()` checks that every managed build's cached signal
contains `|iv:{current}`. Mismatch returns False, causing `db update` to
rebuild the cache.

### One-shot migration

*Amended by ADR-0048:* the migration runs in the update step (`ClinVarAnnotator.prepare_cache()`), not in `is_ready()`, which never writes.

Pre-mechanism caches (before the stamp existed) lack any `|iv:` marker.
`stamp_existing_clinvar_cache(db_path)` in `manager.py` appends
`|iv:{current}` to the existing `remote_signal` value. Called from
`is_ready()`, it self-heals once without re-downloading. Caches stamped with
a stale version (e.g. `|iv:0`) are NOT self-healed — they require a full
`db update` because the interpretation logic has changed.

### Signal stripping

`cached_remote_signal()` strips the `|iv:N` suffix from each per-build
signal before composing the freshness comparison string. The internal stamp
is never compared against the remote signal — it is metadata about the
interpreter, not about the upstream data.

### When to increment

Bump `CLINVAR_INTERPRETER_VERSION` only when what the loader stores
changes: which rows are kept or dropped, how they are written, or a
column. The stamp says what the cache holds, and a bump makes every user
rebuild ClinVar (about 860 MB to download).

**Amended 2026-10-16 (2.5).** This list used to include runtime logic
(the `_CLNSIG_MAGNITUDE` map, the carrier rule in `annotate()`, indel-anchor
protection, benign suppression). Those are read from the code at run time,
so a change to them takes effect without a rebuild, and a bump would make
every user download ClinVar again for nothing. Do NOT bump for runtime
matching, scoring or display changes.

**Tests.** `TestInterpreterVersionStamp` in `test_clinvar.py` with three
tests: matching stamp accepted, missing stamp self-heals, old stamp rejected.

## Amendment 2026-10-16 (2.5): ClinVar's REF only at the call's locus

The hom-ref filter took ClinVar's REF by rsID and build, with no check that
ClinVar's row sits at the call's position, so it depended on how each build
places the rsID. ClinPGx and SNPedia now use
`ClinVarAnnotator.reference_at(rsid, build, locus)`. It gives ClinVar's REF
(or the normal allele, ADR-0039) only where ClinVar places the rsID at the
call's own `(chromosome, position)`. An rsID placed elsewhere, or a call
with no position, is not suppressed.

## Amendment 2026-10-25 (2.5.2): reference-carried sites without the input build's ClinVar

ClinVar's REF still decides wherever its cache for the call's build can be
read. Without that cache, a site where the two builds' single-base REFs
differ (the reference carries a variant allele on one of them, F5 rs6025
on GRCh37) gets its normal allele from:

1. the other build's ClinVar, when its row there names this build's REF as
   the variant allele (or its own REF == ALT row names the other build's
   REF as one);
2. otherwise gnomAD: the major allele at the site is normal, and a
   homozygote for it is a non-finding. The major allele is read from every
   allele gnomAD lists there, the REF's frequency being what the ALTs
   leave: one ALT's low frequency does not make the REF major. A site
   gnomAD does not list, an ALT whose frequency gnomAD does not give (it
   is not counted as zero), or a tie (within 1e-9, since the REF's share
   carries rounding error) decides nothing;
3. otherwise nothing: no suppression, and the row shows.

GRCh38's REF is never taken as normal by itself. Hiding a real
pharmacogenomic finding is worse than showing a normal genotype's row. The
normal allele also sets ClinPGx's and SNPedia's `alt` there
(`normal_allele_for`). Every other site is unchanged: without ClinVar,
nothing is suppressed. Implementation: `ReferenceCarriedFallback`
(`annotators/reference_fallback.py`), wired in `get_annotators`.

### Amendment 2026-10-29 (2.6.0): identity contradictions

A named call whose locus contradicts its rsID's identity-table placement
on the analysis build reaches no source, this one included (ADR-0063).

## Amendment 2.6.0 (2026-11-01): switched-off providers

A ClinVar cache, gnomAD or rsID identity table the user switched off is
never consulted for the normal allele, even when its cache can be read: the
filter reads it as not installed, and the report carries
`source_dependency_unavailable` (ADR-0067).

## Amendment 2.6.0 (2026-11-02): a normal allele from the other build is written on this build's strand

The outside review of RC3 found that the 2.5.2 fallback above, and the
other-build steps of ADR-0039, read "the two builds' REFs differ" as "the
reference allele changed". REFs also differ where dbSNP's reference strand
was inverted between the assemblies: the two letters can then be one allele
written on opposite strands. An allele from the other build (its REF, its
ClinVar row, gnomAD's GRCh38 major allele) was compared with the call
without translation. In the reviewer's case a GRCh37 TT call matched
GRCh38's major T, the letter of GRCh37's minor A, and its SNPedia finding
was suppressed. In the identity cache, 3,436 of the 5,875 rsIDs whose
single-base REFs differ are such inversions; 1,339 are flagged palindromic.

Every allele carried across builds is now written on the input build's
strand first:
- **Identity table** (the one-ClinVar-build path and the fallback): the
  strand relation comes from the rsID's allele sets on both builds, as the
  table's loader pairs them (same strand, or inverted). A row flagged
  palindromic, or allele sets that fit both readings or neither, give no
  relation, and nothing is carried across: the site abstains.
- **The other build's ClinVar** (both builds installed): the same
  VariationID's ALT gives the relation. It must be the carried allele (same
  strand) or its complement (inverted); any other ALT abstains.
- **gnomAD** is GRCh38's own reading: its major allele is translated for
  GRCh37 input, and needs no translation for GRCh38 input.
- Where the two REFs are one allele on opposite strands, the reference did
  not change: the other build's ClinVar decides nothing there, and only
  gnomAD can.

**F5 rs6025 on a one-ClinVar-build install (supervisor ruling, option a).**
The identity table flags rs6025 palindromic (dbSNP lists every base), so
its strand relation between the builds cannot be told. With only one
ClinVar build installed (`db update --build grch37` or `--build grch38`),
the result depends on the input's build and the cache installed:
- **GRCh37 call, only GRCh37's ClinVar** (the identity-table path): the
  normal allele is no longer taken from the identity table. The 2.5.0 fix
  that gave C there is undone for this setup. A T/T call keeps ClinVar's F5
  finding and its SNPedia and ClinPGx T/T rows. A C/C call has no ClinVar
  finding, as before, but is no longer suppressed: its SNPedia and ClinPGx
  C/C rows show (the ClinPGx hormonal-contraceptive row, level 2B, among
  them).
- **GRCh37 call, only GRCh38's ClinVar** (the reference fallback, with no
  ClinVar for the input's build): neither GRCh38's ClinVar nor gnomAD's
  major allele can be written on GRCh37's strand, so nothing is decided and
  the C/C rows show. A T/T call has its SNPedia and ClinPGx rows but no
  ClinVar finding, as before: ClinVar is matched on the input's own build.
- **GRCh38 call:** with GRCh38's ClinVar, its own REF (C) suppresses C/C and
  T/T has the ClinVar finding, as before. With only GRCh37's ClinVar, a
  decisive gnomAD major allele, GRCh38's own reading, still suppresses C/C
  without any translation; without such gnomAD evidence the call abstains
  and its C/C rows show.

Abstaining here hides no finding; taking the letter as written could.

With both ClinVar builds installed (the default), ClinVar's own rows decide
rs6025 (C on both builds, same strand) and nothing changes.

## Amendment 2.6.0 (2026-11-03): the interpreter stamp is `local_version_tag`

The interpreter version is no longer stamped into `remote_signal`. Since
ADR-0028 it is written as `iv:N` to `database_versions.local_version_tag`
when `load_clinvar_tsv()` builds a cache (the VCF loader is gone; ClinVar
is loaded from NCBI's TSVs, ADR-0053), and readiness checks that tag.
`stamp_existing_clinvar_cache()` moves a legacy `|iv:N` marker out of
`remote_signal` into the tag, only when N is the current
`CLINVAR_INTERPRETER_VERSION`; it runs in the update step (ADR-0048).

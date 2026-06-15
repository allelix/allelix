# Changelog

All notable changes are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2026-06-15

### Known limitations

- **ClinVar per-(classification, condition) pairing (#42, tracked for
  v2.1).** ClinVar variants with multiple SCV submissions display the
  primary aggregate classification and the union of all associated
  conditions (semicolon-joined). The exact 1:1 pairing between a
  specific classification and the specific condition it was submitted
  for is not preserved — the VCF source doesn't carry that mapping.
  Full per-(classification, condition) rows via ClinVar's
  `submission_summary.txt.gz` structured export are planned for v2.1.
  Affects ~0.01% of cached ClinVar records — exclusively
  multi-submission variants.

### Changed

- **HuggingFace asset URLs moved to `allelix/` org (#37).** The gnomAD
  and AlphaMissense prebuilt-cache download URLs were previously
  hardcoded to `huggingface.co/datasets/dial481/...`; both datasets now
  live under the `allelix/` HF org to match the GitHub migration. HF
  serves redirects from the old paths, so existing pinned client
  installs continue to work, but new installs pull from the canonical
  org URL. Commit SHAs in the resolve paths are unchanged — the
  underlying asset bytes are identical and the SHA256 integrity checks
  pass without modification.
- **Terminal report: bare-min column set (#9).** The terminal analyze
  table previously had 9–12 columns (rsID, Gene, Source, Significance,
  Review Status, Magnitude, Genotype, Zygosity, Freq, AM, CADD,
  Condition), which Rich auto-squeezed to hairline-zero widths on a
  typical 100–120 col terminal — rsIDs got cut off mid-string, the
  middle columns (Significance/Review Status/Magnitude) rendered as
  empty hairline bands, and the "all `—`" columns (Freq in a
  GWAS-only filter, etc.) wasted space. The terminal is the
  quick-extract view; deep enrichment belongs in HTML/JSON. New
  layout: `rsID | Gene? | Source | Significance | Mag | GT |
  Condition?` (Gene and Condition dropped when no row has data).
  Source uses a compact display (`GWAS Catalog` → `GWAS`),
  significance drops the redundant `source_` prefix (`clinvar_pathogenic`
  → `pathogenic`), Genotype renamed to `GT`. **Terminal output only**
  — every dropped field still surfaces in HTML and JSON reports
  unchanged.

### Fixed

- **`Variant`: alleles normalized to uppercase at construction (#14,
  case half).** Reference databases (ClinVar, gnomAD, ClinPGx, ...)
  all ship uppercase alleles, and carrier matching is raw set
  membership against `{variant.allele1, variant.allele2}` — a
  lowercase user allele would silently miss every clinical annotator
  and a real carrier would get zero annotations with no warning.
  Production parsers all emit uppercase today, but a user-supplied
  filter file (custom panel) or future format variant could leak
  lowercase through. Normalization now happens in
  `Variant.__post_init__` so every emission point (parsers, filters,
  synthetic variants in tests, and any future code) is covered at the
  single chokepoint. No-call marker preserved as-is; multi-base
  indels uppercased in place. Strand-naive half of #14 (using
  `resolve_strand` at carrier-match time) is deferred to v2.1+ per
  ADR-0010.
- **SNPedia parser: allowlist guard on `raw_table` interpolation (#12).**
  The `raw_table` name from `detect_raw_table` flows into three SQL
  queries via f-string interpolation (SQLite doesn't accept parameter
  binding for identifiers). Today the function only ever returns one
  of two literals or None — safe — but the signature is `str | None`
  with no constraint, so any future drift (config-driven name,
  scraped metadata) could open an injection path. Added a
  `_VALID_RAW_TABLES = frozenset({"_raw_pages", "pages"})` guard at
  the top of `_parse_raw_pages_inner`; an unexpected name raises
  `ValueError`. No behavior change on current paths.
- **GWAS rollup: consistent PheCode delimiter (#24, partial).**
  `_gwas_base_trait` was splitting on `" (PheCode "` (leading space)
  while `_gwas_phecode_parent` was searching for `"(PheCode "` (no
  leading space) — the two readers reaching into the same formatted
  description string with different delimiters. Promoted the delimiter
  to a module-level `_PHECODE_DELIM` constant used by both. The deeper
  structural fix (carrying trait / p-value / PheCode as structured
  fields on `Annotation` instead of re-parsing the formatted
  description, per ADR-0024 followup) is deferred to v2.1+.
- **`split_csv_line` docstring corrected (#27).** The previous docstring
  claimed it handles single, double, and double-double-quoted fields
  per real CSV quoting rules. The implementation is
  ``line.split(",")`` + ``strip('"')``, which mishandles any quoted
  field containing a literal comma (silently dropped by callers'
  column-count guard). Docstring now states the limitation and the
  intended consumer set (FTDNA / MyHeritage / Living DNA, none of
  which embed commas in fields). No behavior change.
- **FTDNA `can_parse` rejects MyHeritage files (#26).** FTDNA and
  MyHeritage exports share the same `RSID,CHROMOSOME,POSITION,RESULT`
  header and a byte-identical data section. Both `can_parse`
  implementations were accepting MyHeritage files, with routing masked
  only by the registry's listing order (MyHeritage first). A reorder
  or `--format ftdna` on a MyHeritage file would silently mislabel the
  source format and display name. FTDNA now scans the first 50 lines
  for a `MyHeritage` substring and returns `False` if found; detection
  is mutually exclusive regardless of registry order.
- **GWAS `_magnitude`: `<=` at p-value boundaries (#17).** Strict `<`
  on the magnitude thresholds meant the canonical genome-wide-
  significance value `p = 5e-8` (which appears verbatim in the GWAS
  Catalog) landed in the suggestive bucket (4.0) rather than the
  significant bucket (6.0) — one full magnitude below a barely-
  significant 4.9e-8 hit. All boundary comparisons now use `<=`
  (`5e-100`, `5e-20`, `5e-8`, `5e-6`, `5e-4`). Practical surface
  area is small given the default `--gwas-min-magnitude 9.0`, but
  the scoring was simply wrong at the most meaningful threshold.
- **`diff`: validate baseline entries, magnitude tolerance, no implicit
  zero (#25).** Two bugs in the report-diff path. (1)
  `load_previous_report` only validated the top-level shape; a baseline
  with `"annotations": {}` or entries missing `source`/`rsid` would
  flow into `compute_diff` and raise `KeyError`/`TypeError` rather than
  the documented `ValueError`. Per-entry validation is now done at load
  with a precise error message naming the missing key. (2) "Changed"
  detection used exact float `!=` (so JSON round-trip representation
  noise — `7.5` vs `7.499999999999999` — flagged as changed) and
  substituted `0.0` for a missing previous magnitude (so a baseline
  entry with no magnitude key looked like a change against any nonzero
  current value). Now uses `math.isclose(abs_tol=1e-9)` and carries
  `previous_magnitude: float | None`; absent baseline magnitude is
  passed through as `None` and rendered as `—` in the terminal table.
- **ClinVar `condition` joins all CLNDN entries (#42, suppress-half).**
  The loader was treating `|` in `CLNSIG` and `CLNDN` as parallel arrays
  index-paired with the VCF's ALT, then picking `[0]` from each. That's
  wrong: the `|` separator is per-SCV (per submission), not per-ALT, and
  the two lists are not guaranteed to be in matching order. The
  symptomatic case was a "Frankenstein pairing": one SCV's classification
  shown next to another SCV's condition (real-world example: rs1063192
  surfaced as `Likely_pathogenic + Three_Vessel_Coronary_Disease` when the
  Likely_pathogenic submission was actually for Breast cancer; the
  Three_Vessel_CD submission was Protective). The minimal v2.0.1 fix
  joins the full `CLNDN` list into a single `condition` string per
  record. `clinical_significance` keeps the primary aggregate
  (`CLNSIG[0]`) — that value is correct as a variant-level claim, only
  the condition-pairing was misleading. Honest framing: "this variant
  has classification X; here are all the conditions it has been
  submitted for." Per-(classification, condition) pairing requires the
  `submission_summary.txt.gz` source and is tracked for v2.1.
- **Enrichment: drop MAX-fallback for alt-less annotations (#23, suppress-half).**
  GWAS rows set `alt=""` because the catalog is rsID-keyed. The previous
  enrichment path took a `MAX(af) / MAX(am_pathogenicity) / MAX(phred)
  GROUP BY rsid` fallback for those rows, which at a multi-allelic site
  stamped the highest-frequency / highest-pathogenicity / highest-CADD
  alt's value next to the annotation as if it described the user's
  variant. Same wrong-allele hazard as the strand path fixed in #18.
  Symmetric fix: skip the enrichment number rather than show one that
  may belong to a different allele. The position-fallback path
  (used for rsIDs resolved via ClinVar `bulk_resolve_rsids` from GH #8)
  is allele-specific and stays — it pairs the alt-less annotation's
  rsID with the user's actual carried alt. CADD enrichment grew the
  same position-fallback path (it previously only had the MAX-fallback,
  which the suppress-half removed); for symmetry, alt-less rows whose
  rsID resolved on the fly now get a correct allele-specific CADD score
  rather than nothing. The full fix — carrying the user's alt onto
  every Annotation so GWAS rows can take the exact-alt path even when
  the input VCF already had rsIDs — is architectural and tracked for
  v2.1.
- **`resolve_strand`: drop complement-fallback (#18).** At multi-allelic
  sites, the complement of the user's true forward allele can
  coincidentally equal a different alt at the same position. The
  previous fallback (when the user allele didn't match `{ref, alt}`
  directly, try `complement(user_allele) ∈ {ref, alt}`) would return
  that coincidental complement, and the CADD enrichment path then
  stamped a wrong-allele score onto an annotation that read as
  describing the user's variant. Audit reproduction:
  `resolve_strand('A', 'C', 'T') -> 'T'` (the user's `A` is neither
  REF nor ALT but `complement('A') = 'T' = ALT`). The function now
  returns `None` whenever the user allele isn't directly in
  `{ref, alt}`. ADR-0010 already deferred proper strand handling; this
  makes the deferral explicit at the resolution layer. The CADD
  enrichment helper `_lookup_user_allele` is simplified — the
  complement loop is removed. Bounded blast radius (enrichment-only,
  not the carrier decision), but the previous behavior could read as
  describing the user's variant when it didn't.
- **`normalize_build_label`: anchored tokens, ambiguity rejection (#16).**
  The previous implementation was first-substring-wins over a bare
  `"3X" in s` chain. Bugs caught by the audit: `2038-01-01` registered
  as GRCh38, `hg19/hg38` picked GRCh37 (the first substring hit),
  `38 (37 liftover)` picked GRCh37. The normalizer is now driven by
  anchored regex tokens (`\bgrch3[678]\b`, `\bhg(18|19|38)\b`,
  `\b(36|37|38)\b`, etc.) and collects all distinct builds mentioned:
  exactly one → that build, zero or multiple → `None`. Word boundaries
  reject digits embedded in dates/versions; ambiguity rejection is the
  safer answer when a label genuinely names two builds. Direct callers
  of `normalize_build_label` (compare, export plink, the header-vs-
  detected mismatch warning) now receive `None` instead of a
  silently-wrong build label.
- **LivingDNA parser: only inspect build-marker comment lines (#16).**
  `get_metadata` previously fed every comment line into
  `normalize_build_label` and let the *last* match win, so a download-
  date comment could retag the file. It now skips comment lines that
  don't mention ``build`` / ``reference`` / ``genome``, and breaks on
  the first valid build-marker line — matching the format spec which
  puts ``# Human Genome Reference Build 37 (GRCh37.p13).`` near the
  top of the header.
- **gnomAD `is_ready()`: tagless caches fail (#22).** The previous check
  was `tag == f"sv:{GNOMAD_SCHEMA_VERSION}" or not tag`, so any cache
  written before the schema-version mechanism existed was accepted as
  ready. That defeats invalidation: if `GNOMAD_SCHEMA_VERSION` is ever
  bumped, every legacy cache silently passes as the new version. The
  escape is dropped — `is_ready()` now requires the exact current tag.
  Production caches built via `install_prebuilt_cache` already stamp
  the current tag, so the only caches affected are pre-tag legacy ones,
  which now correctly require `db update`.
- **ClinVar `.md5` body validated as 32 hex chars (#21).** When the
  ClinVar CDN returns a transient HTML error page, the previous parser
  treated the first whitespace-separated token (`<!DOCTYPE`) as the
  expected hash. That token then flowed into `verify_file_hash`, which
  on mismatch deletes the just-downloaded ~400 MB VCF and raises. So a
  five-second 503 upstream wiped the freshly fetched cache. The
  signal parser now requires the token to match `^[0-9a-fA-F]{32}$`
  before treating it as an md5 — non-hex, short hex, or HTML is logged
  and returns `None`, which `db update` handles as "freshness unknown,
  skip" rather than nuking the VCF.
- **`db update`: tagless caches re-download instead of being stamped (#20).**
  A cache with no stored freshness signal almost always predates the signal
  mechanism, i.e. is old. The previous behavior was to stamp the live
  remote signal onto it on the next `db update` and call it current, which
  permanently marked stale data as fresh — every subsequent run saw
  `cached == remote` and skipped, so only `--force` could escape.
  `cached is None` now triggers `setup()` like a remote-signal mismatch
  does, with a distinct message (`cache predates the freshness signal;
  re-downloading…`) so the user knows why. The dead helper
  `_stamp_remote_signal` is removed.
- **GWAS / ClinPGx loaders: post-load row-count floor (#19).** Closes the
  silent-truncation trap. Both upstream sources stream over chunked
  transfer with no `Content-Length`; a connection drop mid-stream produces
  a parseable-but-short TSV that previously committed to the cache and
  shipped fewer associations / clinical annotations than the user
  expected. Loaders now accept a `min_rows=` kwarg and raise `OSError`
  before `os.replace` if the ingested count falls under it. Production
  call sites pass `GWAS_MIN_ROWS=100_000` and `PHARMGKB_MIN_ROWS=5_000`;
  on failure the on-disk cache stays untouched and the user is told to
  `allelix db update --force` to retry. Test fixtures bypass the floor
  via an autouse conftest patch; loader-level `TestMinRowsFloor` covers
  the gate by passing explicit `min_rows=`.
- **Build detection: warn when inconclusive despite rsIDs (#15).** Closes
  the silent-coords trap. When position-detection inspected known-rsID
  positions but couldn't pick a build (votes tied across builds or no row
  matched any), the pipeline silently routed through `header_build` — so a
  GRCh36 file with a GRCh37-mislabeled header used to silently get the
  GRCh37 ClinVar cache. Now fires a yellow warning naming the
  rsID-checks-attempted count, what the header claims, and how to force
  `--build`. The CLI-layer warning leaves `BuildDetectionResult`
  unchanged.

### Changed

- `scripts/fetch_testdata.sh` now downloads from
  `allelix/allelix/releases/download/v2.0.0/test_data.tar.gz` (the new
  ~1.27 GB asset that includes the v2.0.0 VCF/gVCF additions: real
  GIAB HG002 GRCh38 + GRCh37 benchmarks, a real GATK-HC gVCF, a 1000
  Genomes chr22 multi-sample, plus the synthetic mocks). Previous URL
  pointed at `dial481/allelix/releases/download/v1.1.1/test_data.tar.gz`,
  which lacked all the VCF samples.

## [2.0.0] - 2026-06-15

Major version bump. New architectural anchor: VCF + gVCF support with a
batched annotation pipeline tuned for WGS scale, plus position-keyed rsID
resolution that makes the flagship feature actually work on real-world
caller output. Internal CLI restructure into a multi-file package. Second
FTDNA file shape (Illumina raw). A CI sentinel for ClinVar CLNSIG
vocabulary drift.

### Added

- **Position-keyed rsID resolution** (`ClinVarAnnotator.bulk_resolve_rsids`
  in `allelix/annotators/clinvar.py`, wired through
  `_pipeline.run_analysis`). Real-world VCFs from variant callers
  (GATK HaplotypeCaller, DeepVariant) emit `ID=.` — no rsID. Files
  from 1000 Genomes-style pipelines emit positional synthetics like
  `22:10519265:CA:C`. Without rsIDs the rsID-keyed annotation
  pipeline would return zero hits. The resolver queries ClinVar by
  `(chromosome, position)` per build, picks the row whose `ALT`
  matches the user's carrier allele (ADR-0007), and stamps the rsID
  back onto the Variant in-place. ClinVar's own batch_annotate plus
  ClinPGx, GWAS Catalog, and SNPedia all benefit unchanged. New
  `idx_clinvar_position` index on `clinvar_variants(chromosome,
  position)`; existing caches self-heal via a `CREATE INDEX IF NOT
  EXISTS` on first connection. ClinVar coverage is the
  clinically-actionable subset, which is what's in scope for this
  release; the load-bearing test on real GIAB HG002 GRCh38 (4M
  variants) drove annotation count from 0 → 4,867. GH #8.
- **Position-keyed enrichment fallback for gnomAD + AlphaMissense**
  (`bulk_lookup_by_position(set[(chrom, pos, ref, alt)])` on each).
  After the existing rsID-keyed enrichment pass, the pipeline does
  a second pass for annotations whose rsID was resolved on the fly
  AND whose gnomAD/AM rsID index didn't include that rsID. PK
  lookup recovers the allele frequency / pathogenicity score.
- **High-value rsID match runs post-resolution.** The
  `high_value_rsids` set is checked at the tail of each batch flush,
  after `bulk_resolve_rsids` has settled `v.rsid` for every variant
  in the batch. A variant entering the pipeline as `ID=.` that
  resolves to a high-value rsID is now collected into
  `AnalysisResult.hv_variants` (previously the match ran against
  the pre-resolution empty rsid and missed). GH #11.
- **VCF + gVCF parser** (`allelix/parsers/vcf.py`). Plain VCF, gVCF
  (GATK reference-confidence), and `.vcf.gz` all handled via the
  standard library — no `pysam` dependency on the base path. gVCF
  reference blocks (lines with `<NON_REF>` ALT and `END=` INFO) are
  skipped entirely: they match nothing in any annotation database
  and would push the pipeline through millions of zero-hit queries.
  Streams variants one at a time; never materializes the file.
- **Multi-sample VCF handling.** `--sample <ID>` selects the sample
  column. Multi-sample VCF without `--sample` fails loudly at the
  parser layer with the list of available samples — surfaced via a
  pre-flight check in `_resolve_parser`, before the analysis
  pipeline kicks off. Wired into `analyze`, `methylation`,
  `pharmacogenomics`, `stats`, `extract`, and `export plink`.
  `compare` accepts `--sample1` / `--sample2` for two VCFs.
- **Optional `[vcf-index]` extras** for `pip install allelix`. Pulls
  `pysam>=0.22` and enables a tabix-backed fast path for
  `extract --snps` on huge VCFs. Falls back to sequential scan
  without the extras — correct for any VCF, just O(file size).
- **FTDNA Illumina raw parser** (`allelix/parsers/ftdna_illumina.py`).
  Tab-delimited variant of the FTDNA export, distinct from the
  comma-delimited CSV variant that already shipped. Detection via
  the canonical `RSID/CHROMOSOME/POSITION/RESULT` tab-delimited
  header.
- **Chromosome normalization helper** (`normalize_chromosome` in
  `parsers/_helpers.py`). VCF/UCSC `chr`-prefixed and lowercased
  chromosomes are canonicalized to bare form (`1`, `X`, `MT`) at
  the parser boundary. Used by both VCF and FTDNA Illumina raw
  parsers; available for any future format that needs it.
- **R-4: ClinVar CLNSIG vocabulary drift CI test**
  (`tests/databases/test_clinvar_clnsig_drift.py`). Pins the
  current ClinVar CLNSIG vocabulary as a snapshot
  (`allelix/data/clinvar_clnsig_snapshot.yaml`); the fast tier
  validates code consistency in default CI, the
  `pytest -m slow` tier queries the installed ClinVar SQLite
  cache and fails on any unknown term. Surfaces drift introduced
  by fresh ClinVar releases. Snapshot covers 47 canonical terms
  including 2024+ vocabulary (penetrance modifiers, VUS
  subdivisions, risk-allele tiers, slash-joined aggregates).

### Changed

- **`analyze`'s annotation pipeline restructured to two-phase
  batched flow** (`allelix/reports/_pipeline.py`). Phase 1: build
  detection prelude (unchanged, ADR-0021). Phase 2: variants
  accumulate into a buffer of `_BATCH_SIZE = 5000` and are handed
  to each annotator's `batch_annotate` on every fill. EOF flushes
  the final partial batch. At WGS scale (4–6 M variants) this
  replaces millions of per-variant SQLite roundtrips with a few
  thousand chunked `WHERE rsid IN (...)` queries — the difference
  between annotation completing in minutes and in hours.
  High-value SNP collection folded into the same single pass; the
  second `parser.parse(file_path)` iteration the v1.x pipeline did
  is gone, which matters for huge VCFs where re-streaming
  gigabytes is non-trivial.
- **Annotator interface gains `batch_annotate(variants)`.** Base
  class default loops over `annotate(variant)` so existing
  annotators work unchanged via the fallback. ClinVar, ClinPGx,
  GWAS, and SNPedia each implement an SQL-optimized override with
  chunked `IN`-clause queries. ClinVar groups by build (ADR-0021)
  before dispatching. Parity tests in each annotator's test file
  pin the contract: `batch_annotate(vs)` produces the same
  annotation set as `[a for v in vs for a in annotate(v)]`.
- **PharmGKB rebranded to ClinPGx** (dial481/allelix#2). Upstream
  rebrand: `pharmgkb.org` 301-redirects to `clinpgx.org`; both API
  endpoints 303 to the same `s3.pgkb.org` archive. User-facing
  surfaces — `display_name`, `attribution`, terminal/HTML/JSON
  report strings, annotation description prefix, license
  attribution text, README and CHANGELOG narrative — now say
  ClinPGx. Internal identifiers (`Annotation.source == "pharmgkb"`,
  `pharmgkb_annotations` / `pharmgkb_allele_function` table names,
  `pharmgkb.sqlite` cache filename, `PharmGKBAnnotator` class) are
  unchanged: existing user caches keep working unmigrated, and
  downstream tools that key off `source == "pharmgkb"` still match.
  The upstream URL flipped from `api.pharmgkb.org` to
  `api.clinpgx.org`. ADRs that reference PharmGKB stay as-is —
  they describe decisions made when that was the name.
- **CLI module restructured from `allelix/cli.py` (1596 lines) to
  the `allelix/cli/` package** (8 files, 36–516 lines each).
  External API unchanged (`from allelix.cli import main` still
  works via `cli/__init__.py`). Subcommands moved to dedicated
  modules: `analyze.py`, `db.py`, `utility.py`, `focused.py`,
  `config.py`. Shared option decorator constants live in
  `cli/_options.py`; orchestration helpers and the Rich console
  in `cli/_helpers.py`. The `_helpers.func(...)` namespace-call
  pattern (rather than destructured imports at use sites) gives
  the test suite a uniform `monkeypatch.setattr` surface.

### Fixed

- **rsID resolver: indel-anchor disambiguation (ship-blocker, found in
  external review).** The position-keyed resolver matched ClinVar rows
  by ALT only. At positions where two co-located ClinVar entries share
  an ALT but differ in REF anchor length (e.g. `CA>C` and `CAG>C` at
  the same position — common in clustered indel regions like CFTR /
  BRCA), sort order picked the wrong row and stamped the user's
  variant with the neighbor's rsid. Downstream annotators bound to
  the wrong rsid; in the failure direction where Benign was suppressed
  by default, the user's actual Pathogenic variant silently vanished.
  Fix: the resolver now requires `user_alleles ⊆ {ref, alt}` so the
  row's REF must also be consistent with the user's diploid call, not
  just the ALT. The 0/1 het case (user carries `{REF, ALT}`) resolves
  correctly. The 1/2 multi-allelic case and hom-alt-at-shared-anchor
  case remain ambiguous from the Variant alone (no `Variant.ref` field
  yet — v2.1+) and now skip resolution rather than silently
  mis-stamping. Real-data measurement on GIAB HG002 GRCh38: pre-fix
  4,867 surviving annotations dropped to 4,843 after the fix — two
  distinct rsIDs (`rs34079299` at chr12:68158714 IFNG, `rs59439148` at
  chr10:45374099 ALOX5; both anchor-base indels) were silently
  misrouted to a co-located neighbor's rsid before this fix.
- **VCFs from variant callers now annotated** (GH #8). Pre-fix, real
  GIAB HG002 GRCh38 benchmark VCF (4,048,342 variants, `ID=.` per
  GATK convention) produced **0 annotations** against every
  rsID-keyed source. The pipeline filter was `if not v.rsid`; pivoted
  to `if not v.rsid.startswith("rs")` so positional synthetic IDs
  (`22:10519265:CA:C` from 1000 Genomes pipelines) and COSMIC IDs are
  also routed through resolution. Post-fix, the same file: 4,867
  surviving annotations from all four rsID-keyed sources, including
  ClinPGx at the bottom of the resolution chain.
- **`SampleNotFoundError` and `MultiSampleError` messages truncate
  the sample list.** On a 1000 Genomes multi-sample VCF (3,202
  samples), the error message previously dumped every sample ID
  verbatim — many screens of unreadable terminal output. New
  `format_sample_list(samples, max_shown=10)` shows the first 10 and
  appends "... and N more".
- **Build auto-detect failure now warns and recommends `--build`.**
  When the input has no rsIDs AND no `##contig assembly` tag (the
  common shape for VCFs from variant callers), build detection
  silently fell back to GRCh37 default. Now also prints a yellow
  warning telling the user to pass `--build grch37` or
  `--build grch38` explicitly — using the wrong build silently misses
  every annotation.
- **VCF `_canonical_rsid` now lowercases the `rs` prefix on
  return.** Original-case `RS1801133` would pass parser detection
  but silently miss every annotation lookup (downstream SQL all
  assumes lowercase). Non-`rs` IDs (e.g., COSMIC `COSV12345`)
  preserve case — their own conventions apply.
- **FTDNA Illumina parser applies `normalize_chromosome`
  defensively.** Real FTDNA Illumina exports observed in the
  wild use plain numeric chromosomes today, but a non-canonical
  variant (chr-prefix, lowercase) would otherwise silently miss
  all rsID lookups downstream.

### Known limitations

- **rsID resolver abstains on two cases that can't be disambiguated
  from the Variant model alone.** The `Variant` dataclass carries
  `(allele1, allele2)` per ADR-style genotype representation, not the
  VCF's REF column. When the resolver finds multiple co-located
  ClinVar rows that all pass the
  ``user_alleles ⊆ {ref, alt}`` consistency check, it has no way to
  pick the right one and abstains — the variant flows through without
  rsID-keyed annotations rather than risk a silent mis-stamp. Two
  shapes hit this:
  1. **1/2 multi-allelic** at multi-rsid positions. User carries two
     different ALTs (e.g. VCF `REF=A, ALT=G,T, GT=1/2` → user alleles
     `{G, T}`). Each ClinVar row at the position has `REF=A` and one of
     the ALTs; the user's alleles aren't a subset of either single row's
     `{REF, ALT}`, so all rows fail the subset check → zero matches →
     abstain.
  2. **1/1 hom-alt at shared-anchor indel positions.** User carries one
     ALT homozygously (e.g. `(C, C)`) at a position where two ClinVar
     rows share that ALT but differ in REF (e.g. `CA>C` and `CAG>C`).
     `{C}` is a subset of both `{CA, C}` and `{CAG, C}`, so multiple
     matches → abstain.

  A future `Variant.ref` field (v2.1+) carrying the VCF's REF column
  would let the resolver pick the unique correct row in both cases.
  Pinned by `tests/annotators/test_clinvar.py::TestBulkResolveRsids::
  test_multi_allelic_1_2_genotype_no_resolution` and
  `test_indel_anchor_homalt_abstains_when_ambiguous`.

## [1.9.0]

### Added
- **`--filter-file` option on `analyze`.** Accepts a plain text file with
  rsIDs and gene names (one entry per line) for custom-panel filtering of
  the analyze report. Lines matching `^rs\d+$` (case-insensitive) are
  treated as rsIDs; everything else is a gene name. Comments (`#`) and
  blank lines are ignored. Gene and rsID matches combine with OR — an
  annotation passes if it matches either set. An empty filter file (or
  one with only comments) produces an empty report. Not added to
  `methylation` or `pharmacogenomics`, which already operate on curated
  panels.

## [1.8.4]

### Added
- **`--no-cadd` flag on analyze, methylation, and pharmacogenomics commands.**
  Per-invocation CADD enrichment exclusion. Required for commercial
  deployments without a CADD license from UW CoMotion (CADD is
  `commercial_ok=False`). Mirrors the existing `--exclude-snpedia`
  licensing-exclusion pattern.

## [1.8.3]

### Changed
- **README quickstart.** Lead with `pip install allelix` for end users.
  Development setup moved to the Development section.
- **PyPI publish workflow hardened.** Pinned `pypa/gh-action-pypi-publish`
  to full commit SHA (v1.14.0) instead of moving branch ref.
- **README links.** Relative markdown links replaced with absolute GitHub
  URLs so they resolve on PyPI.
- **Security policy scope.** Narrowed supported versions to current minor
  only (1.8.x). v1.x feature line is frozen at v1.8.3.

## [1.8.2]

### Fixed
- **HTML report link colors.** Links in dark mode were browser-default
  neon blue (#0000EE), unreadable against the dark background. Dark mode
  now uses #93c5fd; light mode uses #1976d2 (matching the existing accent).
- **Favicon SVG rendering.** The report's inline SVG favicon used a
  `linearGradient` with an internal IRI reference (`url(#g)`) that doesn't
  resolve inside `data:` URIs, rendering the icon invisible. Replaced with
  a solid fill and URL-encoded angle brackets.
- **Genotype column header.** Abbreviated "Genotype" to "GT" to prevent
  column header overlap on narrow/mobile screens.

### Changed
- **PyPI metadata.** Homepage now points to allelix.io. Added Source and
  Changelog links.
- **Automated PyPI publishing.** GitHub Actions workflow publishes to PyPI
  via Trusted Publishing on every GitHub Release.

## [1.8.1]

### Changed
- **Updated sample reports.** Regenerated `examples/sample_reports/` with
  v1.8.0 report format (5-column table, detail sidebar, dark mode).

### Fixed
- **Test fixture genotype format.** Corrected the `test_html.py` fixture's
  `genotype_match` default from `A/G` to the concatenated `AG` form that
  production emits for SNVs.

## [1.8.0]

### Changed
- **HTML report redesign.** Replaced the 12-column scrollable table with
  a compact 5-column layout (Magnitude, Gene, Genotype, Repute, Summary).
  Annotations from multiple sources for the same variant are grouped into
  a single row. Clicking a row opens a sliding detail sidebar showing all
  source annotations vertically — genotype, zygosity, significance, review
  status, condition, description, frequency, AlphaMissense, CADD PHRED,
  and references.
- **Dark / light mode.** Follows system preference (`prefers-color-scheme`)
  by default. A toggle button lets the user override. All component CSS
  uses custom properties — no hardcoded colors outside semantic badges and
  the accent color.
- **CADD and AlphaMissense legend.** The "Understanding Magnitude Scores"
  section now includes CADD PHRED tier thresholds (≥30 top 0.1%, ≥20 top
  1%, ≥10 top 10%) and AlphaMissense classification bands (≥0.564
  likely pathogenic, 0.340–0.564 ambiguous, <0.340 likely benign).
- **CADD tier context in sidebar.** CADD scores display the PHRED
  percentile tier inline (e.g. "38.0 (top 0.1% most deleterious)")
  instead of a bare number.
- **Embedded variant JSON uses numeric types.** `allele_frequency`,
  `am_pathogenicity`, and `cadd` are now floats in the `<script
  id="variant-data">` blob, matching the v4 JSON schema. Field names
  aligned: `am_pathogenicity`, `am_class`, `allele_frequency` (was
  `amScore`, `amClass`, `frequency`).

### Fixed
- **PLINK split-chromosome error.** MHG exports can have straggler
  autosomal variants appended after the Y chromosome section, producing
  non-contiguous chromosome blocks in the .bim. PLINK1.9 rejects these.
  `export plink` now sorts variants by chromosome then position before
  writing. The exporter itself remains a single-pass writer — the sort
  lives in the CLI layer.

## [1.7.0]

### Added
- **PLINK export (#29).** `allelix export plink` converts any supported
  genotype format to PLINK1 binary (.bed/.bim/.fam) for downstream
  tools (plink2 PCA, ADMIXTURE, PRSice). Single-sample, SNP-major
  encoding. Uses gnomAD ref/alt for allele coding when available.
  No-call variants skipped.

### Fixed
- **Multi-allelic strand collision in PLINK coord selection.** At sites
  where one alt is the complement of another (e.g. ref=G, alts=A,T),
  single-pass coord selection picked the complement match over the
  forward match — same bug class as CADD #45. Fixed with a two-pass
  loop that prefers forward allele matches.

### Documentation
- **ADR-0034: Magnitude scoring scale and ceiling (#23).** Formalizes the
  0-10 scale with practical ceiling of 9. Documents max-across-sources
  composite rule and reserves per-source scoring for v2.0.0.

## [1.6.1]

### Added
- **Zygosity column in all report outputs (#41).** Every annotation row
  now shows `Heterozygous`, `Homozygous`, or `No Call` — derived from
  the genotype call. Appears in HTML, terminal, and JSON reports.
  Functional-medicine `+/−` notation deferred until the risk allele
  field lands (v2.1+).
- **CADD PHRED styling in HTML reports (#42).** CADD scores are now
  color-coded by deleteriousness tier: red (≥30, top 0.1%), orange
  (≥20, top 1%), gray (10–20, top 10%), gray (<10, no tooltip).
  Tooltip shows the percentile tier on hover for scores ≥10.
- **Config file path in `config get` / `config set` output (#43).**
  Both commands now print the resolved config file path, so users know
  which file is being read or written.

### Changed
- **JSON schema version bumped to 4.** New `zygosity` field on every
  annotation. Diff between v1–v4 reports still works.
- **Methylation panel expanded to 34 genes (#31).** Added ACAT1, DHFR,
  GNMT, MAOA, NOS3, SUOX, VDR. Sorted alphabetically.

### Fixed
- **CADD multi-allelic scoring bug (#45).** At multi-allelic positions,
  `_enrich_cadd` max-reduced CADD PHRED across all alts, stamping the
  highest score regardless of which allele the user carries. Now looks
  up the score for the user's specific allele (direct match preferred
  over complement). Biallelic sites were unaffected.
- **Test protocol AM column name.** `FULL_TEST_PROTOCOL.md` referenced
  "AM Score" but the actual HTML header is "AM".

## [1.6.0]

### Added
- **CADD v1.7 variant deleteriousness scores (ADR-0032).** PHRED-scaled
  scores ranking how deleterious any single-nucleotide variant is, using
  100+ annotation tracks. Enrichment-only annotator following the
  gnomAD/AlphaMissense pattern. Two modes: cache (pre-built SQLite from
  HuggingFace, ~5 GB, ~120M variant keys) and full (81 GB tabix file via pysam, GRCh38
  only). CADD column appears in HTML, terminal, and JSON reports when
  scores are present.
- **Non-commercial source opt-in pattern.** CADD is the first source
  with `commercial_ok=False`. Disabled by default (`sources.cadd =
  false`). Users opt in via `allelix config set sources.cadd true` or
  `allelix db update --cadd`. First download shows a license
  confirmation prompt.
- **Strand normalization for array data.** `resolve_strand()` maps
  array-reported alleles to reference-forward orientation using gnomAD
  ref/alt as ground truth. Palindromic SNPs (A/T, C/G) return None
  rather than guessing.
- **CADD cache build script.** `scripts/build_cadd_cache.py` filters
  the full CADD SNV and indel files to positions present in gnomAD,
  AlphaMissense, and ClinVar (GRCh38). Uses int64 packing for SNV
  keys to fit the ~120M position set (117M SNV + 3M indel).
- **`options.cadd_full` config key.** Enables full CADD mode (tabix
  queries against the complete CADD file). Requires `pip install
  allelix[cadd]` for pysam.

## [1.5.3]

### Fixed
- **Added missing GWAS Catalog attribution to HTML and JSON reports.**
  GWAS Catalog was registered as an annotator but absent from the license
  attribution maps in both renderers.

### Added
- **SECURITY.md.** Vulnerability reporting policy (GitHub private
  vulnerability reporting), supported versions, and scope definition.

### Changed
- **License metadata centralized on annotator base class via
  LicenseDescriptor (ADR-0031).** Each annotator now declares its
  license as a required `license` ClassVar. Non-commercial gating
  derived from SPDX identifier instead of a hand-maintained set.
  Report attribution text generated from the descriptor at render
  time. The `NON_COMMERCIAL_SOURCES` frozenset in `config.py` is
  deleted.
- **JSON schema version bumped to 3.** The `license_attributions`
  block now carries `source_url` (source website) and `license_url`
  (license deed) as separate keys. `license` field uses SPDX
  identifiers. Diff between v2 and v3 reports still works.

## [1.5.2]

### Fixed
- **Coverage gate enforced by pytest again.** v1.5.1 moved the threshold
  to `[tool.coverage.report]`, but pytest only fails on low coverage when
  `--cov-fail-under` is set (CI runs `pytest`, not `coverage report`), and
  the config key is honored version-dependently on unpinned pytest-cov.
  Restored `--cov-fail-under=92` to `addopts`; `precision = 2` retained,
  so the v1.5.1 rounding fix stands. Pinned `pytest-cov>=7,<8` to prevent
  future drift.

## [1.5.1]

### Fixed
- **Download integrity verification.** ClinVar downloads now verify
  their md5 checksum against the NCBI sidecar file. HuggingFace
  downloads (gnomAD, AlphaMissense, SNPedia) pin to specific commit
  SHAs and verify SHA256 after download. Mismatches delete the corrupt
  file and raise. GWAS Catalog and PharmGKB are documented gaps — no
  upstream checksum exists. See ADR-0029.
- **Coverage gate rounding.** pytest-cov with default `precision=0`
  rounded 91.91% to 92%, silently passing `--cov-fail-under=92`. Set
  `precision=2` and `fail_under=92.00` in `[tool.coverage.report]`.
- **Circular import (#36).** `_versions.py` moved from
  `allelix/annotators/` to `allelix/databases/` to break the
  `gnomad_loader → annotators → alphamissense → gnomad_loader` cycle.

### Changed
- **Extracted shared loader utilities.** `install_prebuilt_cache()` was
  duplicated across gnomad_loader, alphamissense_loader, and
  snpedia_loader. Extracted to `install_prebuilt_gz_cache()` in
  `databases/loader_utils.py`.
- **Two-tier data source model (ADR-0030).** Server-driven sources
  (ClinVar, GWAS Catalog, PharmGKB) probe for freshness at runtime.
  Code-driven sources (gnomAD, AlphaMissense, SNPedia) use
  commit-pinned HuggingFace URLs — no HEAD requests, no signal
  stamping, refresh only via `--force` or code bump of the pinned
  commit SHA. Vestigial `probe_http_signal()` and all redundant test
  monkeypatches removed.

## [1.5.0]

### Changed
- **Version tag consolidation across all six annotators.** Local
  processing stamps are now stored in a dedicated `local_version_tag`
  column in `database_versions` instead of being appended to
  `remote_signal` as `|iv:N` / `|pv:N` suffixes. This eliminates the
  fragile suffix-parsing pattern that caused the SNPedia signal-loop
  bug: `remote_signal` now holds only the remote ETag/Last-Modified,
  `local_version_tag` holds the local processing state. All six
  annotators use the same dual-version mechanism:
  - ClinVar: `iv:N` (interpreter version)
  - PharmGKB: `iv:N` (interpreter version)
  - SNPedia: `pv:N` (parser version)
  - GWAS Catalog: `cv:N` (categorizer version)
  - gnomAD: `sv:N` (schema version)
  - AlphaMissense: `sv:N` (schema version)

  The `sv:` tag is new for gnomAD and AlphaMissense — pre-built caches
  now stamp their schema version so a future schema change forces
  re-download, matching the cache-invalidation behavior the other four
  annotators already had. Existing caches self-heal on first run — no
  re-download required. `get_database_info()` lazily adds the
  `local_version_tag` column when reading pre-v1.5.0 caches.

### Fixed
- **Multi-allelic enrichment accuracy (#25).** gnomAD and AlphaMissense
  enrichment now uses exact alt-allele matching instead of `MAX()`
  aggregation at multi-allelic sites. Added `alt` field to the
  Annotation model; `bulk_lookup_by_alt()` on both enrichment
  annotators; pipeline splits exact-match and MAX-fallback paths.
- **Disk preflight multiplier (#27).** Bumped from 5x to 6x for both
  gnomAD and AlphaMissense loaders. AlphaMissense compresses at 4.4x,
  so peak disk (gz + decompressed) is 5.4x gz — the old 5x check
  would greenlight a disk that ENOSPC'd at ~90% decompression.
- **Test suite disk usage.** db update tests were downloading real
  databases (7.8 GB AlphaMissense, 678 MB GWAS TSV) into pytest
  tmp_path because they lacked monkeypatches for all annotators. All
  db update tests now stub every annotator. Real-data GWAS tests
  delete the extracted TSV after SQLite load. Total pytest tmp_path
  reduced from ~3.4 GB to ~376 MB.
- **SQLite variable limit portability (#33).** `bulk_lookup_by_alt()`
  batched at 900 keys (1800 bound variables) — over the 999 limit on
  SQLite < 3.32. Now batches at 450 keys (900 variables).
- **GWAS enrichment regression.** GWAS annotations set `alt` to the
  risk allele, but GWAS risk alleles are not VCF-normalized ALT. This
  caused exact-match lookups to miss, skipping the MAX fallback and
  losing gnomAD/AM enrichment on GWAS rows. Fixed by not setting `alt`
  on GWAS annotations (risk alleles are conceptually different from
  VCF ALT alleles).
- **SNPedia `db update` crash and re-download loop.** Three related
  bugs in the SNPedia download flow: (1) `install_prebuilt_cache`
  crashed with `no such table: database_versions` because the
  third-party HuggingFace cache doesn't include that table — fixed by
  creating the table on demand before stamping, consistently across all
  three pre-built cache loaders and the CLI signal-stamp fallback path.
  (2) `parse_raw_pages` overwrote the ETag remote signal with only the
  parser version tag, causing every subsequent `db update` to see a
  signal mismatch and re-download — fixed by preserving the existing
  ETag when appending the parser version. (3) `cached_remote_signal`
  returned the raw signal including the parser version suffix, which
  never matched the remote ETag — root cause eliminated by the version
  tag consolidation above.

### Added
- **ADR-0028: Local version tag convention.** Documents the
  `local_version_tag` mechanism, the tag prefix vocabulary per
  annotator, the new-annotator checklist, and the lazy migration
  strategy.
- **AlphaMissense gnomAD version stamping (#28).** Build script stamps
  which gnomAD version provided the rsID mapping. Runtime warning on
  version mismatch: "AlphaMissense cache was built against gnomAD X
  but installed gnomAD is Y."
- **SNPedia HuggingFace download (#30).** `db update` now downloads
  the SNPedia cache from HuggingFace automatically — same pattern as
  gnomAD and AlphaMissense. The manual scraper scripts remain as a
  rebuild-from-source option.
- `allelix/databases/snpedia_loader.py` — pre-built cache download
  and decompression.
- `test_data/FULL_TEST_PROTOCOL.md` — external reviewer checklist for
  full real-data verification.

## [1.4.0]

### Added
- **AlphaMissense variant pathogenicity enrichment.** New
  `AlphaMissenseAnnotator` enriches annotations with missense variant
  pathogenicity scores from DeepMind's AlphaMissense (71M variants,
  CC BY 4.0). Pre-built SQLite cache downloaded from HuggingFace via
  `db update`. AM Score column in terminal, HTML, and JSON reports.
  PharmGKB rows show AM scores as neutral with caveat (protein
  structure impact only — tooltip in HTML, dimmed `*` footnote in
  terminal, `am_caveat` field in JSON). `--no-alphamissense` flag to
  skip.
- **Config file system.** `config.toml` with per-source on/off toggles
  and `license.commercial = true` safety switch that auto-disables
  non-commercial sources (SNPedia). `allelix config show/set/reset`
  CLI commands. CLI flags override config per-invocation.
- `scripts/build_alphamissense_cache.py` — AlphaMissense cache build
  script with Zenodo HTTPS streaming (default) and local TSV modes.
  Joins against gnomAD cache for coordinate-to-rsID mapping.
- AlphaMissense CC BY 4.0 attribution in HTML and JSON reports.
- Magnitude scoring legend in HTML report (collapsible, per-source
  scoring tables for ClinVar, PharmGKB, GWAS, SNPedia).
- Source floor note in HTML report when per-source magnitude minimums
  are active.
- Repute row background tints in HTML report (red for pathogenic/risk,
  green for protective/benign) derived from existing significance
  field.
- Sortable columns in HTML report (magnitude, gene, source, AM score)
  via inline JavaScript.
- ADR-0027 documenting the AlphaMissense enrichment cache architecture.
- `scripts/run-tests.sh` — detached background test runner with log
  rotation.

### Fixed
- HTML report table overflows viewport, columns clipped on left (#20).
  Added `overflow-x: auto` container, sticky rsID column,
  `max-width` on description cells, refs collapsed into `<details>`
  toggle, conditional Review Status column (hidden when all empty),
  stat card `flex-wrap`.
- AlphaMissense build script has zero unit-test coverage (#24). Added
  25 tests covering TSV parsing, gnomAD rsID join, chr prefix
  normalization, `--no-gnomad` NULL-rsid path, multi-allelic composite
  PK, batched insert, and end-to-end integration.
- Download integrity: Content-Length check after downloads catches
  truncated files.
- Disk space preflight before decompressing `.sqlite.gz` caches uses
  5x gz size (accounts for gz + decompressed tmp on disk
  simultaneously).
- `_connection()` guards on gnomAD and AlphaMissense annotators raise
  `FileNotFoundError` with actionable message when cache is missing.
- Dead `cache_exists()` removed from gnomAD and AlphaMissense loaders.
- Legacy caches stamp remote signal instead of re-downloading on
  `db update`.
- README database sizes updated to match actual on-disk measurements.

### Changed
- `db update` display includes gnomAD and AlphaMissense in "Analyzing
  against" annotator list.
- Both build scripts (`build_gnomad_cache.py`,
  `build_alphamissense_cache.py`) run `VACUUM` for smaller output
  files.

## [1.3.1]

### Fixed
- Test suite downloaded real ~6 GB gnomAD cache on every run, filling CI
  runner disk. All `db update` tests now use a 792-byte mock fixture via
  `file://` URL -- same pattern as ClinVar, PharmGKB, and GWAS. No
  production code changes.

### Changed
- CI: job timeout (20 min), pytest step timeout (15 min),
  `workflow_dispatch` trigger, verbose output (`pytest -v --tb=short`)
- Ship tooling: `scripts/tag-release.sh` derives tag from pyproject.toml
  (single source of truth)
- Git hooks: raw `.githooks/pre-push` replaces pre-commit framework shim,
  blocks tag pushes where version doesn't match
- CONTRIBUTING.md: corrected slow-test documentation (CI skips them, not
  runs them), added "Run the full suite locally" section emphasizing that
  developers must run the full suite with real-data fixtures locally
  before pushing
- Documentation: fixed stale hook instructions, added missing changelog
  comparison links and ADR index entries
- Removed dead code: `scripts/check_version_tag.sh`
- Removed `version-tag-match` entry from `.pre-commit-config.yaml`

## [1.3.0]

### Added
- **gnomAD population allele frequencies (R-6).** New `GnomadAnnotator`
  enriches report annotations with population frequency context from
  gnomAD v4.1 exomes (~16M rsIDs). Pre-built cache downloaded from
  HuggingFace via `db update`. Frequency column appears in terminal,
  HTML, and JSON reports when gnomAD data is available. `--no-gnomad`
  flag on `analyze`, `methylation`, `pharmacogenomics`, and `db update`
  to skip.
- **CPIC fallback for PharmGKB (R-5).** `db update` succeeds when the
  CPIC API is unreachable — PharmGKB downloads complete and the
  non-finding filter degrades gracefully. Signal carries
  `cpic:unavailable` so recovery auto-triggers a refresh.
- `scripts/build_gnomad_cache.py` — streaming VCF build script for
  the gnomAD frequency cache. Downloads ~185GB over HTTPS, never saves
  VCFs to disk, outputs ~6GB SQLite.
- `scripts/extract_array_manifest.py` — extracts rsID superset from
  genotype files for filtered gnomAD cache builds.
- gnomAD ODbL v1.0 attribution in HTML and JSON reports.
- JSON report `schema_version` bumped to `"2"` (added `allele_frequency`
  field on annotations). Diff engine accepts both v1 and v2 baselines.
- `db update` now handles individual annotator failures gracefully —
  prints error and continues to remaining annotators instead of aborting.

### Fixed
- Offline claim in README and ADR-0012 corrected: analysis runs offline
  by default with an opt-out freshness check (`--no-update`), not
  opt-in network access (#10).
- `.gitignore` updated for GWAS Catalog test data (#12).
- `scripts/fetch_testdata.sh` downloads GWAS Catalog associations
  from EBI FTP (#12).

### Changed
- Pre-push hook reduced to version-tag check only; pytest removed
  (CI runs the full suite on every PR, pre-push pytest caused SSH timeouts).
- CI version-tag guard job added to `.github/workflows/ci.yml` (#11).

## [1.2.0] — 2026-06-07

### Fixed
- `pyproject.toml` version corrected to match release (was `1.1.0` on
  the v1.1.1 release).
- **GRCh36 fallback bug.** Non-confident GRCh36 detection (e.g., 3/4
  probe SNPs matched) was falling back to GRCh37 as the effective build,
  silently bypassing the ClinVar safety guard and annotating GRCh36
  positions against GRCh37 coordinates. Fixed in both the end-of-stream
  `flush()` path AND the buffer-limit path (large files where probe SNPs
  appear past the 100K-variant buffer cap). The pipeline now uses GRCh36
  as the effective build whenever detection points to GRCh36, even
  non-confidently or with a single probe SNP match.

### Added
- **Auto-refresh stale databases.** `analyze`, `methylation`, and
  `pharmacogenomics` now check database file ages before running. If any
  database is older than 7 days and the remote signal (MD5/ETag) has
  changed, the database is refreshed automatically. If the network is
  unreachable, a warning is printed and analysis continues with the stale
  cache. SNPedia is excluded (no remote download). Use `--no-update` to
  skip the freshness check entirely.
- Corpas family exome VCF attribution in `test_data/edge_cases/README.md`
  with paper DOI (Corpas et al., *BMC Genomics* 2015,
  doi:10.1186/s12864-015-1973-7). Licensing table in `test_data/README.md`
  updated. Every genotype fixture in the repo now has documented
  provenance and license.
- **Version-tag drift guard.** Pre-push hook
  (`scripts/check_version_tag.sh`) asserts any pushed `v*` tag matches
  the version in `pyproject.toml`. Prevents the class of bug where a
  release ships with a stale version string.

## [1.1.1] — 2026-06-06

### Changed
- Relocated real genotype test data (`test_data/real/` and
  `test_data/transcoded/`) to GitHub release assets. Fresh clone size
  reduced from ~650 MB to ~150 MB. Tests skip gracefully when data is
  absent; `scripts/fetch_testdata.sh` restores it.
- Clarified `.gitignore` and `test_data/README.md`: the "never commit"
  rule applies to private genetic data, not CC0 public-domain openSNP
  fixtures hosted as release assets.

### Fixed
- Orphaned `[Unreleased]` changelog sections assigned proper version
  numbers (`[0.7.2]` and `[0.8.0]`) matching their chronological
  position in the development history.
- Duplicate `[0.7.1]` changelog header consolidated into a single entry.
- Dead compare links for internal pre-release versions removed (0.x tags
  were never pushed to the public repository).

## [1.1.0] — 2026-06-06

### Added
- **`allelix compare` command** with strand-aware concordance classification:
  concordant, strand-flip match, discordant, strand-ambiguous, no-call. Build
  detection via `detect_build()` with `get_metadata()` fallback. Per-chromosome
  breakdown. Build rows in Coverage Summary table.
- **High-value SNP no-call flagging.** YAML data file with 12 clinically
  important SNPs (APOE, BRCA1, MTHFR, CYP2D6, etc.). No-call warnings surface
  in `stats`, `analyze`, and all report formats. Cluster-incomplete detection
  (e.g., APOE genotype cannot be determined). Loader supports merging
  user-provided YAML overrides with error handling for malformed input.
- **ClinVar review status column.** CLNREVSTAT surfaced in terminal, HTML, and
  JSON renderers including all diff tables (new, changed, removed). Users can
  distinguish expert-panel-reviewed from single-submitter pathogenic calls.
- **GRCh36 position-based build detection.** All 11 probe SNPs now have GRCh36
  positions. 3-way voting across builds. Headerless files (FTDNA, MyHeritage)
  with GRCh36 positions now detected correctly.
- **CONTRIBUTING.md** with "How to add a parser" and "How to add an annotator"
  tutorials, development setup instructions, and coding standards summary.

### Changed
- **23andMe parser detection tightened.** Anchored to canonical first-line
  header `# This data file generated by 23andMe`. Bare substring matches
  rejected. Fallback loop for user-prepended comments before the canonical line.

### Fixed
- `is_must_include` internal field no longer leaks into public JSON output
  (filtered from `annotations`, `diff.new`, and `diff.changed` paths).
- Build detection docstrings updated for GRCh36 (module docstring and
  `BuildDetectionResult`).
- `_ready_annotators` return type annotation corrected.
- Compare command now uses `detect_build()` instead of reading the parser's
  default build from `variants[0].build`.
- README status updated from "Pre-release" to "Production".

## [1.0.0] — 2026-06-05

> Six parsers, four annotators, three report formats, report diffing,
> 794 tests, 94% coverage. All array-based consumer DNA formats supported.

### Added
- **23andMe parser** (`parsers/twentythreeandme.py`). Four-column tab-delimited
  format with concatenated genotype. Handles I-prefixed probe IDs, haploid
  MT/Y calls, no-calls (`--`), CRLF line endings. Detection by "23andMe" in
  early comment lines. Build from header comments (supports build 36/37/38).
- **AncestryDNA parser** (`parsers/ancestrydna.py`). Five-column tab-delimited
  format with separate allele columns. Chromosome mapping: 23→X, 24→Y, 25→X
  (PAR), 26→MT. No-calls as `0`. Detection by `#AncestryDNA` first-line
  signature. Default build 37.
- **FTDNA parser** (`parsers/ftdna.py`). CSV format with double-quoted fields,
  concatenated genotype in RESULT column. Handles quoted/unquoted headers,
  haploid MT/Y calls. Detection by `RSID,CHROMOSOME,POSITION,RESULT` header
  pattern. Default build 37.
- **MyHeritage parser** (`parsers/myheritage.py`). CSV format, structurally
  identical to FTDNA. Detected by "MyHeritage" in first comment line. Handles
  double-double-quoted field variant. Shares `_helpers.py` with FTDNA.
- **Living DNA parser** (`parsers/livingdna.py`). Tab-delimited despite `.csv`
  extension. Handles AX-, AFFX-prefixed probe IDs and CHR:POS positional
  notation. Build detection from header comments via `normalize_build_label`.
- **Shared parser helpers** (`parsers/_helpers.py`). `split_csv_line` and
  `split_genotype` extracted from FTDNA to share across CSV-family parsers.
- **GRCh36/hg18 build detection.** `normalize_build_label` recognizes "36",
  "hg18", "build 36", "GRCh36". 23andMe parser detects "build 36"/"hg18" in
  header comments. CLI emits warning that GRCh36 positions won't match modern
  references.
- **`BUILD_GRCH36`** constant in `build_detect.py`.
- **Build mismatch warning in HTML report** (R-9). When the file header claims
  a different build than position data indicates, the HTML report renders a
  visible warning banner matching the CLI warning.
- **"Reading This Report" education section** (R-7). Static HTML block after
  the regulatory notice covering pseudogene cross-hybridization, ClinVar
  aggregation, carrier vs. affected, and confirmatory testing.
- **Signal guard on all annotator `setup()` methods.** ClinVar, PharmGKB, and
  GWAS annotators abort `setup()` if `fetch_remote_signal()` returns None,
  preventing persistence of incomplete cache stamps.
- **Diff key collision fix.** Diff key extended from `(source, rsid, condition)`
  to `(source, rsid, condition, description)` — prevents silent data loss when
  PharmGKB has multiple annotations for the same rsid/condition with different
  drugs.
- **Terminal diff test coverage.** Four tests covering new-only, changed-only,
  removed-only, and no-changes branches of `render_terminal_diff`.
- **`--exclude-snpedia` flag** on `analyze`, `methylation`, and
  `pharmacogenomics` commands. Suppresses SNPedia annotations at the CLI
  level — required for commercial use (CC BY-NC-SA 3.0). Wired through
  existing `exclude_sources` plumbing.
- 39 FTDNA tests, 37 23andMe tests, 31 AncestryDNA tests, 34 MyHeritage
  tests, 34 Living DNA tests, 8 new HTML tests. Full suite: 790.

### Changed
- **`--diff` CLI help** reframed as a dev/QA tool for version-to-version
  validation, not monitoring.
- **CLI build banner** now shows "header (no position confirmation)" when probe
  SNPs don't match but a header build is present, instead of the misleading
  "fallback (no known SNPs matched)".
- Parser registry includes all six parsers in detection order: MyHappyGenes,
  23andMe, AncestryDNA, Living DNA, MyHeritage, FTDNA.
- Removed all internal phase numbering from user-facing documentation.

### Known limitations
- **GRCh36 build detection is incomplete on headerless formats.** The build
  detector's probe table (`KNOWN_SNP_POSITIONS`) only has GRCh37/GRCh38
  positions. FTDNA files with GRCh36 coordinates are silently labeled as
  GRCh37 (no warning). PharmGKB, GWAS, and SNPedia continue to fire via
  rsID-only lookups, but ClinVar would query against the GRCh37 cache with
  GRCh36 positions and miss any variants whose coordinates shifted between
  builds. Tracked as R-12.
- **No ClinVar GRCh36 cache.** `CLINVAR_SUPPORTED_BUILDS` is GRCh37/GRCh38
  only. When a file is correctly detected as GRCh36 (e.g. 23andMe with
  "build 36" in the header), `analyze` produces zero ClinVar annotations —
  PharmGKB, GWAS, and SNPedia still fire normally. The CLI prints a warning.
  Full GRCh36 annotation support requires either an NCBI GRCh36 VCF or
  liftover. Tracked as R-12.

## [0.9.2] — 2026-05-22

### Fixed
- SNPedia parser now drops stale unique index before recreating with
  COALESCE, fixing silent no-op on caches that already had the old
  `idx_snpedia_genotype_dedup` definition.
- Backfill dedupe removes pre-existing NULL-summary duplicate rows
  (2 G6PD-family entries) on parser re-run.
- SNPedia annotator now auto-reparses when parser version changes.
  `_PARSER_VERSION` stamped into `database_versions.remote_signal`;
  `is_ready()` rejects stale caches and triggers re-parse automatically.
  Eliminates the recurring "code fix shipped but cache not rebuilt" failure.
- SNPedia parser extracts alleles from the page title when the
  `{{Genotype}}` template omits `allele1`/`allele2`. 22 genotype pages
  (e.g. `Rs104894073(A;G)`) carried alleles only in the title. These
  were silently dropped before; now parsed correctly. All 79
  originally-dropped pages verified against live SNPedia API — zero
  content differences between our scrape and source.
- SNPedia parser and annotator now handle I-prefixed 23andMe probe IDs
  (`I3000043`, `I5006212`, etc.). 1,402 genotype pages (1,401 rows
  after 1 allele-order dedupe on I4000178) and 2,851 SNP pages use
  23andMe internal probe IDs instead of rs-numbers. Gene mapping reads
  `Gene_s` from `{{23andMe SNP}}` templates. 392 I-probe SNP pages
  carry gene mappings. Total structured rows: 104,720 (101,328 with
  gene, 3,392 without). Prepares for 23andMe parser I-probe annotation.

### Added
- **ClinVar interpreter version stamp (ADR-0023).** `CLINVAR_INTERPRETER_VERSION`
  in `annotators/_versions.py` is stamped as `|iv:N` into each per-build
  cache's `remote_signal` during ingest. `is_ready()` rejects caches with a
  stale or missing stamp. One-shot migration self-heals pre-mechanism caches
  without re-downloading. Eliminates the "annotator logic changed but cache
  wasn't rebuilt" failure mode for ClinVar.
- **GWAS auto-reingest from cached TSV (ADR-0024).** `setup()` now retains
  the raw GWAS Catalog TSV (~200 MB) after ingest. When `_CATEGORIZER_VERSION`
  bumps and the cached TSV is present, `is_ready()` auto-reingests from
  the local file without re-downloading. Users who already have the data
  get updated classification automatically on first run.
- **`TestInterpreterVersionStamp`** (3 tests) pinning ClinVar interpreter
  version stamp behavior: matching stamp accepted, missing stamp self-heals,
  old stamp rejected.
- **`TestAutoReingest`** (2 tests) pinning GWAS auto-reingest: categorizer
  bump with TSV present triggers reingest, missing TSV returns False.

### Fixed
- **`cached_remote_signal()` now strips internal stamps.** ClinVar's
  `cached_remote_signal()` was returning the raw stored signal including
  `|iv:N`, causing freshness comparisons to always show a mismatch and
  triggering unnecessary re-downloads on every `db update`. Internal stamps
  are now stripped before composing the comparison string.

### Docs
- ADR-0023 amended with ClinVar interpreter version stamp mechanism.
- ADR-0024 amended with GWAS auto-reingest from cached TSV.
- Added roadmap feature R-11: supplemental genotype file merging
  (custom panels, Sanger confirmations).

## [0.9.1] — 2026-05-21

### Fixed
- SNPedia parser dedupes source-level genotype duplicates. Five rsIDs
  (rs4950928 and four G6PD-family entries) had paired identical rows
  in snpedia_genotypes from SNPedia source pages whose titles differed
  only in whitespace.

### Changed
- GWAS rollup: collapse MTAG twins and PheCode hierarchical
  sub-classifications before rendering. rs10455872 (LPA) drops from
  8 to 5 distinct findings. Must-include rows exempt. See ADR-0024.

### Docs
- ADR-0008 corrected SNPedia/ClinVar overlap claim (~11% complementary,
  not 0%).

## [0.8.2] — 2026-05-20

### Fixed
- **GWAS categorizer now classifies against `DISEASE/TRAIT` in addition to
  `MAPPED_TRAIT`.** UKB data-field rows with empty `MAPPED_TRAIT` (e.g.
  impedance traits) were leaking into reports at mag-9 in v0.8.1.
  `classify_gwas_trait()` now concatenates both fields for keyword matching.

### Changed
- **`_CATEGORIZER_VERSION` marker stamped into
  `database_versions.remote_signal`.** `schema_is_current()` rejects caches
  built under a different categorizer version, forcing auto-rebuild on
  `db update`. Eliminates the stale-cache failure mode from v0.8.0/v0.8.1.

## [0.8.1] — 2026-05-20

### Fixed
- **GWAS categorizer keyword leaks.** 11 UKB body-composition and
  metabolite-ratio rows leaked into the "other" category at mag-9.
  Three structural noise patterns addressed:
  - UKB bioimpedance traits ("whole body water mass", "impedance of arm",
    etc.) — 7 keywords added to `_BODY_MEASUREMENT_KW`.
  - NMR metabolite ratios ("cholesterol-to-phospholipid ratio") — new
    `_is_metabolite_ratio()` helper catches the `-to-…ratio` pattern.
  - Uncharacterized analytes ("X-12345 level") — new
    `_is_uncharacterized_analyte()` helper catches the `x-` prefix pattern.
  ADR-0024 amended with step 1.5 structural noise detection.

### Added
- **Must-include rsID allowlist (ADR-0024).** `_MUST_INCLUDE_RSIDS`
  frozenset in `gwas.py` lists clinically significant GWAS associations
  that bypass `--gwas-min-magnitude` floor: rs10737680 (CFH/AMD),
  rs11209026 (IL23R/IBD), rs9271366 (HLA-DRB1/MS). Global
  `--min-magnitude`, trait-category filter, and carrier rule still apply.
- **`Annotation.is_must_include` field.** Boolean flag set by the GWAS
  annotator; `AnalysisResult.filter()` exempts flagged annotations from
  per-source magnitude floors.
- **`TestStructuralNoiseDetection`** (18 parametrized tests) covering
  UKB body-composition, metabolite ratio, and uncharacterized analyte
  classification plus disease non-misrouting assertions.
- **`TestMustInclude`** (5 tests) covering constant shape, carrier flag,
  source floor bypass, global min-magnitude enforcement, and trait
  filter enforcement.

## [0.8.0] — 2026-05-13

> **ClinVar REF allele is the primary PharmGKB non-finding filter (ADR-0023).**
> Five prior releases (v0.5.x–v0.7.1) iterated on a CPIC-based filter that
> required CPIC to classify the user's allele as `Normal function`. CPIC's
> vocabulary and coverage are heterogeneous across genes — CFTR uses
> `"ivacaftor responsive"`, MTHFR/F2/F5 have no entries at all. Real-world
> v0.7.3 output leaked ~30+ CFTR × ivacaftor reference-homozygote rows.
>
> The fix is structural: the inclusion/exclusion question is "does the
> user carry the variant allele?" — and ClinVar publishes REF universally
> for every variant it catalogs. The new primary filter is one check:
> if ClinVar's REF is single-base and matches both of the user's alleles,
> the row is a non-finding. CPIC is demoted to a secondary tier for the
> rare rsid ClinVar doesn't know about.

### Fixed
- **CFTR-class reference-homozygote leak.** Real-world report at
  `--min-magnitude 5` previously surfaced ~30+ CFTR × ivacaftor rows
  saying "do not have a copy of the variant." All of these are now
  suppressed by the ClinVar REF check.
- **Inconsistent genotype display across annotators.** ClinVar rows
  used to show single-letter `genotype_match` (the matched ALT base);
  PharmGKB rows showed the user's diploid. Both now show the user's
  sorted diploid (`"AG"`, `"GG"`, etc.). Indel passthrough is verbatim
  (`"CTT/C"`).

### Added
- **`--include-benign` flag.** ClinVar Benign/Likely_benign annotations are
  now suppressed by default at the annotator level. Pass `--include-benign`
  on `analyze`, `methylation`, or `pharmacogenomics` to restore them.
  ADR-0008 amended to document the policy.
- **`--gwas-min-magnitude` flag (default 9.0).** Per-source magnitude floor
  for GWAS Catalog annotations. On real data the GWAS Catalog produces
  ~88,000 associations; the 9.0 floor keeps only hyper-significant + large
  effect size signals. Floor raised from 7.0 to 9.0 after real-data
  testing showed 30,000+ mag-7 common-trait rows (Height, BMI, blood
  counts). ADR-0024 amended.
- **`--include-gwas` flag on `methylation` and `pharmacogenomics`.** Focused
  reports exclude GWAS Catalog annotations by default — methylation biology
  is interpreted from ClinVar + PharmGKB, not GWAS trait associations.
  Pass `--include-gwas` for completeness.
- **GWAS trait-category filtering (ADR-0024).** Each GWAS Catalog row is
  classified into a trait category (disease, cancer, drug_response, immune,
  cardiovascular, metabolic, neurological, body_measurement, lipid_measurement,
  hematological_measurement, other_measurement, behavioral, other) using EFO
  ontology labels. Default excludes measurement and behavioral categories.
  On real data: 605 GWAS rows pass at mag-9 (down from 3,314 pre-filter).
  `--gwas-all` disables trait-category filtering for the full unfiltered dump.
- **`TestBenignSuppressionEndToEnd`** and **`TestDefaultReportSanity`** in
  `test_end_to_end.py` pinning that benign suppression works and that
  default filters produce a tractable annotation count (≤ 20).
- **`TestMethylationSanity`** in `test_end_to_end.py` pinning that
  methylation output (ClinVar + PharmGKB, no GWAS) stays under 20 rows.
- **`TestRealDataGwasSanity`** (`@pytest.mark.slow`) in
  `test_end_to_end.py` — runs against the real GWAS Catalog (test_data/,
  gitignored) and pins that default filters keep output bounded. Skips
  when real data hasn't been downloaded.
- **`test_benign_suppressed_by_default`**, **`test_include_benign_flag`**,
  **`test_gwas_min_magnitude_default`**, and
  **`test_gwas_min_magnitude_lowered`** in `test_cli.py` covering the
  new CLI flags.
- **`test_gwas_excluded_by_default`** and **`test_include_gwas_flag`** on
  `TestMethylationCommand`; **`test_gwas_excluded_by_default`** on
  `TestPharmacogenomicsCommand` — pin GWAS exclusion from focused reports.
- **GWAS Catalog annotator.** Trait–SNP associations from EBI/NHGRI.
  Carrier rule (ADR-0007): only fires when the user carries the risk
  allele. P-value magnitude scoring (ADR-0024): six tiers from 2.0
  (weak) to 8.0 (hyper-significant) with an OR/beta modifier (+0.5
  or +1.0, capped at 9.0). Unknown-risk-allele entries fire on rsID
  match alone but are capped at magnitude 3.0 so they don't pass
  typical `--min-magnitude 5` thresholds.
- **ADR-0024** documenting GWAS Catalog magnitude scoring from p-value
  and effect size, unknown-risk-allele cap rationale, and the
  forward-strand assumption (deferred to R-1).
- **`TestGwasMockInvariants`** (four tests) pinning the GWAS fixture
  shape per ADR-0015: single-base risk allele, unknown risk allele,
  p-value tier coverage, and haplotype-skip row.
- **`ClinVarAnnotator.reference_for(rsid, build) -> str | None`.**
  Per-build lazy-built `(rsid -> single-base REF)` lookup. Indel REFs
  are excluded by the loader's SQL filter (`WHERE length(ref) = 1`).
- **`PharmGKBAnnotator(..., clinvar_ref_provider=callable)`** accepts
  the REF provider. `get_annotators()` wires `clinvar.reference_for`
  into the PharmGKB instance.
- **ADR-0023** documenting the architectural shift.
- **`TestClinvarRefPrimaryFilter`** (six tests) pinning the
  REF-primary filter behavior: homozygous-reference suppression,
  heterozygous emission, CFTR-class leak fix, fallback to CPIC for
  rsIDs ClinVar doesn't know, multi-base REF (indel) fall-through,
  and consistent diploid genotype display.

### Changed
- **Default `--min-magnitude` raised from 0.0 to 5.0.** Previous default
  dumped every annotation regardless of importance. The new default surfaces
  only clinically meaningful findings (ClinVar Pathogenic/Likely_pathogenic,
  PharmGKB LoE 1–2, GWAS genome-wide significant). Pass `--min-magnitude 0`
  for a full dump.
- **GWAS Catalog download URL updated to FTP ZIP.** EBI deprecated the old
  API endpoint (Nov 2025). The loader now downloads the ZIP archive from
  `ftp.ebi.ac.uk` and extracts the TSV.
- **Focused reports exclude GWAS by default.** `methylation` and
  `pharmacogenomics` commands no longer run the GWAS annotator. On real
  data, 98.5% of methylation output was GWAS noise (342/347 rows, 228
  from FUT2 human-milk-oligosaccharide studies). `--include-gwas` opts in.
- **`PharmGKBAnnotator.annotate()` filter order.** Primary tier:
  ClinVar REF check (suppress hom-ref, emit carrier WITHOUT the cache's
  `is_nonfinding` filter so a CPIC-driven false-suppression can't hide
  a real carrier). Secondary tier (no ClinVar REF data for rsid):
  fall through to the cache's pre-computed `is_nonfinding` flag.
- **`Annotation.genotype_match` semantics.** Now always the user's
  sorted diploid for SNVs; passthrough for indels. ClinVar previously
  set this to the matched ALT base — that was lossy and inconsistent
  with PharmGKB's display.
- **ADR-0022 scope reduced.** Still applies, but only to rsIDs where
  BOTH ClinVar and CPIC lack data. The CFTR/MTHFR/F2/F5 cases that
  motivated ADR-0022 are now caught by the ClinVar REF check.

### Verification
- 347 tests pass (was 341; +6 new for the primary tier). Coverage
  steady. Lint and format clean.
- The recurring "PharmGKB filter still leaks on a new gene" failure
  mode that spanned v0.5.x–v0.7.1 is structurally resolved: the
  filter no longer depends on CPIC's per-gene vocabulary.

## [0.7.2] — 2026-05-13

> **Genome build auto-detection (ADR-0021 + ADR-0022) with Round 23
> audit follow-up.** A real-world MyHappyGenes/Tempus export was
> confirmed to ship GRCh38 positions while its header claims "build
> 37.1." Cross-build REF/ALT comparison against ClinVar's GRCh37 VCF
> produced a false-positive pathogenic call on NIPA1 `rs104894490`.
> The carrier check (ADR-0007) was correct — it was matching against
> the wrong build's REF/ALT. The fix is structural: detect build from
> position data, hold per-build ClinVar caches, dispatch per variant.
> ADR-0022 documents the deliberate decision NOT to filter PharmGKB
> reference-genotype rows on non-CPIC genes.
>
> An external audit (Round 23) subsequently found three fixture-layer
> defects and one `.gitignore` omission in the build-detection work.
> The new code was right; the test fixtures were lying about what they
> covered.

### Added
- **`allelix/utils/build_detect.py`** with a hardcoded ~11-entry
  `(rsid, build) → (chromosome, position)` table covering chromosomes
  1, 10, 11, 12, 17, 19, and 22. `detect_build()` streams variants and
  returns the matching build once any table entry's position confirms.
- **`allelix.databases.cpic_loader`'s API + `manager.CLINVAR_URL_BY_BUILD`**
  — per-build URL map containing both `CLINVAR_URL_GRCH37` and a new
  `CLINVAR_URL_GRCH38` (`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/`).
- **`ClinVarAnnotator(data_dir, builds=("GRCh37", "GRCh38"))`** holds
  one SQLite cache per build (`clinvar.GRCh37.sqlite` and
  `clinvar.GRCh38.sqlite`). `annotate(variant)` dispatches by
  `variant.build`. `is_ready()` requires every managed build to be
  cached. `fetch_remote_signal()` returns a composite signal
  `"GRCh37:md5:…|GRCh38:md5:…"` so either side's update triggers a
  refresh.
- **CLI flags.** `allelix db update --build {grch37, grch38, both}`
  (default `both`) and `allelix analyze --build {auto, grch37, grch38}`
  (default `auto`). `methylation` and `pharmacogenomics` accept the
  same `--build` flag.
- **Build mismatch warning.** When `analyze` detects a build different
  from the file header's claim, a yellow warning surfaces explaining
  the discrepancy and naming both builds.
- **ADR-0021** documenting the auto-detection feature and the policy
  of distrusting file headers in favor of position data.
- **ADR-0022** documenting the deliberate decision NOT to filter
  PharmGKB reference-genotype rows on non-CPIC genes. The rows emit;
  the README documents the limitation.
- **Mock data generator `--build` and `--header-build` flags.**
  Default emits GRCh38 (matching real MHG behavior). Three test
  fixtures committed: `mock_myhappygenes.txt` (clean GRCh38),
  `mock_myhappygenes_grch37.txt` (clean GRCh37), and
  `mock_myhappygenes_mislabeled.txt` (GRCh38 positions, GRCh37 header
  — replicates the real-world MyHappyGenes mislabel).
- **`tests/utils/test_build_detect.py`** with 18 tests pinning the
  detector's behavior on confident matches, single-match unambiguity,
  inconsistent files, table-shape invariants, and label normalization.
- **Three end-to-end CLI tests** in `test_cli.py`:
  `test_analyze_warns_on_build_header_data_mismatch`,
  `test_analyze_no_warning_on_clean_grch37`, and
  `test_analyze_build_override_skips_detection`.
- **README "Known PharmGKB limitation" section** documenting non-CPIC
  reference-genotype rows per ADR-0022.
- **`tests/fixtures/mock_clinvar_grch37.vcf` and `mock_clinvar_grch38.vcf`**
  with 11 single-allele records + 1 multi-allelic row. rs104894490
  NIPA1 included with build-specific REF/ALT (the strand-inverted
  regression case).
- **`tests/test_mock_data_invariants.py::TestClinvarFixturePositionInvariants`**
  pinning that every rsID in the build-detect table uses
  build-authoritative positions in the matching fixture, and that
  rs104894490 specifically has the strand-inverted REF/ALT pair
  preserved across the two fixtures. ADR-0015's "mock-as-spec"
  invariant now applies to the ClinVar fixture too.
- **rs104894490 NIPA1 entry in the MHG generator's known SNPs**
  (chr15, G/G genotype, both build positions). The mock MHG default
  fixture grew from 2,015 to 2,016 SNPs to carry the regression case.

### Fixed
- **Per-build ClinVar fixtures.** The pre-Round-23 mock_clinvar.vcf
  had mixed-build positions (some GRCh37, some GRCh38, rs80357906 at
  41245466 matching neither). This was the exact failure mode
  ADR-0021 was written to detect, hardcoded into the project's own
  test fixture. Replaced by two build-correct fixtures:
  `tests/fixtures/mock_clinvar_grch37.vcf` and `mock_clinvar_grch38.vcf`.
  Generator rewritten (`tests/generate_clinvar_fixture.py`) to emit
  per-build VCFs with positions verified against NCBI dbSNP. The old
  `mock_clinvar.vcf` is removed.
- **rs80357906 BRCA1 position corrected.** Old fixture had 41245466;
  authoritative is GRCh37 41209080 / GRCh38 43057063 (verified
  against NCBI Variation API and the live ClinVar GRCh37 VCF). The
  `build_detect.py` table was correct; the fixture was wrong.
- **Dual-cache dispatch test coverage is no longer fictional.** The
  conftest now loads `mock_clinvar_grch37.vcf` into the GRCh37 cache
  and `mock_clinvar_grch38.vcf` into the GRCh38 cache, so per-build
  dispatch produces DIFFERENT results across caches. Two new
  end-to-end tests pin the contract directly using rs104894490 NIPA1:
  - `test_nipa1_strand_inversion_no_emission_on_grch38_data` — the
    smoking-gun case from the user's report. MHG fixture has G/G at
    the GRCh38 NIPA1 position; auto-detection identifies GRCh38;
    dispatch queries the GRCh38 cache (REF=G ALT=A); zero A alleles
    → no annotation. The false positive is gone.
  - `test_nipa1_grch37_dispatch_reproduces_legacy_false_positive` —
    forces `--build grch37` on the same data; dispatch queries the
    GRCh37 cache (REF=C ALT=G); user's G matches ALT=G → annotation
    DOES emit. This is the OLD wrong behavior pinned so a future
    "default to GRCh37" regression visibly flips this assertion.

### Changed
- **Auto-detection is the analyze default.** Files without a CLI
  `--build` override get detected from position data; the
  parser-reported header build is treated as informational. A `[dim]
  Build:` banner prints to terminal showing the effective build, the
  source (detected / override / fallback), and how many known SNPs
  matched.
- **README** updated to note that real-world MyHappyGenes files ship
  GRCh38 despite their header claim. Developer documentation explicitly
  tells future contributors not to trust the header.
- **Conftest** populates both build caches via `clinvar_data_dir` so
  the rsID-based annotator queries work regardless of detected build.
- **Database versioning** uses per-build record names
  (`clinvar.GRCh37` / `clinvar.GRCh38`) in `database_versions`.
- **`mock_clinvar.vcf` fixture removed.** The single mixed-build file
  predated ADR-0021's per-build split. Callers that referenced the
  generic name (test_manager.py, test_cli.py legacy paths) now use a
  back-compat alias fixture that points to `mock_clinvar_grch37.vcf`.
  New tests should use the build-specific fixtures directly.
- **Count snapshots updated for the +1 NIPA1 row.** MHG 2,015 → 2,016;
  ClinVar 12 → 13 per build; composite status display 24 → 26.

### Tooling
- **`.gitignore` excludes `test_data/`** — exploratory scripts kept
  out of CI lint/format and out of the committed tree.
- **`pyproject.toml [tool.ruff]` adds `extend-exclude = ["test_data"]`**
  belt-and-braces: even if ruff is run with `--no-respect-gitignore`,
  exploratory scripts under `test_data/` are skipped.

### Migration
- v0.7.x users have a single `clinvar.sqlite` cache, which the new
  annotator ignores. **Run `allelix db update --force`** to populate
  the new per-build caches (`clinvar.GRCh37.sqlite` and
  `clinvar.GRCh38.sqlite`). The legacy `clinvar.sqlite` can be
  deleted; nothing references it anymore.

### Verification
- 341 tests pass (was 336; +5: two end-to-end NIPA1 regressions, three
  fixture invariants). Coverage 95.61%. Lint and format clean.
- The NIPA1 case from the user's real-world report has BOTH directions
  of the dispatch contract pinned: correct behavior on GRCh38 data,
  reproducible legacy false positive when GRCh37 is forced.

### Notes
- The user's NIPA1 false positive was diagnosed in this cycle: header
  mislabel + cross-build strand inversion + correct carrier rule
  matching the wrong build's REF/ALT = false pathogenic call. v0.7.1
  shipped a carrier check the user briefly thought was broken; this
  release vindicates that code and addresses the actual root cause.
- Coverage of the seven remaining ClinVar pathogenic hits (NCR3, IL10,
  GP9, ADRA2A, TMPRSS6, PKD1, FLG) requires running `db update --force`
  on real data and re-annotating; some may be the same cross-build
  artifact as NIPA1. The tool now uses the correct build automatically.

## [0.7.1] — 2026-05-13

> **Patch release for v0.7.0's PharmGKB regression.** v0.7.0 left three
> production rows leaking — CACNA1S rs1800559 CC (×7 anesthetic
> annotations) and NUDT15 rs116855232 CC (×3 thiopurine annotations).
> Root cause: v0.7.0 populated the per-allele function table from CPIC
> template sentences embedded in PharmGKB's annotation_text. That's
> regex on description text — and many CPIC genes (including CACNA1S
> and NUDT15 in real data) don't publish those sentences at all,
> leaving the lookup empty for those rsids and the prose-fallback tier
> failing to catch the carrier-vs-reference distinction.
>
> The filter is a table join, not a text classifier. The per-allele
> function table now comes from CPIC's structured API
> (`api.cpicpgx.org/v1`), where every CPIC-curated allele carries a
> discrete `clinicalfunctionalstatus` value. ADR-0020 documents this.

### Fixed
- **PharmGKB non-finding filter is now a pure structured-data lookup.**
  `is_nonfinding_by_allele_lookup(rsid, genotype, lookup)` is the
  filter: for the user's two bases, look up each in
  `(rsid, base) → function_class`; if both map to `Normal function`,
  suppress the row. No regex, no description parsing, no text
  classification.
- **v0.7.0 production leakers all classify correctly.** Verified
  end-to-end against the real PharmGKB July-2025 dump joined against
  the live CPIC API:
  - rs1800559 CC (CACNA1S, 7 anesthetic annotations) → all `is_nonfinding=1`.
  - rs116855232 CC (NUDT15, 3 thiopurine annotations) → all `is_nonfinding=1`.
  - v0.6.1 DPYD cluster (rs115232898 TT ×2, rs1801266 GG ×1, rs3918290 CC ×5)
    → all `is_nonfinding=1`.
  - Known carriers (rs1801133 AG/AA ×20, rs1799853 CT ×3, rs4149056 CT ×39,
    rs4244285 AG ×7) → all `is_nonfinding=0`. Zero over-suppression.
- **C-1: test suite no longer hits the real CPIC API.**
  `test_db_update_with_file_url` monkey-patches
  `fetch_cpic_allele_functions` on the annotator module to return the
  same `MOCK_CPIC_LOOKUP` the rest of the suite uses. Run time on that
  single test dropped from ~4 s (network-bound) to 0.6 s (offline).
- **C-2: README refreshed to v0.7.1.** Status banner, supported-
  databases table, and architecture-decision summary reflect ADR-0020
  (CPIC API as the structured per-allele function source) and drop
  obsolete references to ADR-0013 / ADR-0014.
- **M-3: CPIC conflict-resolution policy flipped.** When the same
  `(rsid, base)` appears under multiple allele definitions with a
  Normal-vs-non-Normal disagreement, the loader now picks the
  non-Normal classification. Suppressing happens only when EVERY base
  is Normal — biasing conflicts toward non-Normal ensures we never
  silently suppress a real variant just because one CPIC row flagged
  it Normal. Pinned by `test_conflict_prefers_non_normal`.

### Added
- **ADR-0020: CPIC API as the structured per-allele function source.**
  Documents the three-way client-side join across CPIC's
  `sequence_location` (rsid), `allele_location_value` (base), and
  `allele` (clinicalfunctionalstatus) tables.
- **New module `allelix.databases.cpic_loader`** with
  `fetch_cpic_allele_functions()` — returns `(rsid, base) → function_class`.
- **`load_pharmgkb_tsv` accepts `allele_function_lookup` directly.**
  Production fetches it from CPIC's API in the annotator's `setup()`;
  tests inject a synthetic dict.
- **`mock_cpic_lookup` pytest fixture** in `conftest.py` holds the
  canonical synthetic CPIC lookup the test suite joins against.
- **M-1: retry-with-backoff on CPIC fetches.** `_http_get_json`
  retries up to 3 times with `(1, 2, 4) s` backoff on `URLError`,
  `TimeoutError`, and `JSONDecodeError`. A single transient TCP RST
  during `db update` no longer aborts the refresh.
- **M-2: composite PharmGKB + CPIC freshness signal.**
  `PharmGKBAnnotator.fetch_remote_signal()` now returns
  `pgkb:<pgkb-signal>|cpic:<cpic-signal>` where the CPIC portion is
  the latest date from CPIC's `change_log` table. If CPIC publishes
  new allele functions while PharmGKB's zip is unchanged, the next
  `db update` detects it and refreshes. CPIC probe failure returns
  None (existing "can't verify, pass --force" UX).
- **M-4: mutation-gap tests pinning critical policies.**
  `_classify_cpic_status` returning None for unknown statuses, and
  `fetch_cpic_allele_functions` skipping multi-base / non-ACGT values,
  both have direct pins. A future refactor can't silently break
  either without a test failure.
- **m-2: zip cleanup wrapped in try/finally** in
  `PharmGKBAnnotator.setup()`. A failure between download and ingest
  no longer leaves the staged `clinicalAnnotations.zip` on disk.
- **m-6: `MOCK_CPIC_LOOKUP` shape invariants.** New
  `TestCpicLookupMockInvariants` checks every fixture entry is
  `(rsid, single-base)` keyed and every value is a recognized
  `function_class`. Mirrors what real `fetch_cpic_allele_functions`
  returns, so the same fixture-as-spec violation that bit v0.6.0 and
  v0.7.0 can't recur with CPIC data.
- **`tests/databases/test_cpic_loader.py`** with 16 tests covering
  the three-way join, multi-base filtering, unknown-status filtering,
  null-dbsnp filtering, conflict resolution, network-error
  propagation, retry-then-success, timeout retry, malformed-JSON
  retry, and the CPIC freshness probe's failure modes.

### Removed (regressions reverted)
- **All regex-on-prose tiers gone.** `_NONFINDING_PROSE_FALLBACK`,
  `_CPIC_ALLELE_FN_RE`, `extract_cpic_allele_function`, and
  `build_allele_function_lookup` (the CPIC-template extractor)
  deleted. ADR-0017 and ADR-0018 superseded.
- **Prose-fallback safety-net test deleted.** No prose tier to safety-
  net any more.
- **Mock fixture trimmed.** PA-010 (somatic-prose synthetic) and
  PA-011 (CPIC-template synthetic) removed. PA-008 retained but
  rewritten to exercise the structured-lookup path (reference
  homozygote suppression).

### Migration
- v0.7.0 caches have the same SQLite schema as v0.7.1 (no schema
  change), but the `pharmgkb_allele_function` table is empty/sparse
  on v0.7.0 caches built from PharmGKB's clinical_ann_alleles. **Users
  on v0.7.0 must run `allelix db update --force`** to re-ingest with
  the CPIC API as the lookup source.

## [0.7.0] — 2026-05-13

> **Architectural fix for PharmGKB SNV classification.** v0.6.1's hybrid
> classifier still leaked DPYD rows because its prose fallback couldn't
> catch every CPIC phrasing variant. The fix is not more regex — it is
> recognizing that PharmGKB publishes per-allele function for SNV
> alleles inside a **canonical CPIC template sentence** embedded in
> `annotation_text`. Parsing a bounded template to extract fielded
> data is structurally distinct from regex-on-prose intent
> classification. v0.7.0 extracts that field at load time and stores it
> in a new structured table.

### Fixed
- **PharmGKB SNV per-allele function now structurally extracted (ADR-0018).**
  The DPYD cluster the v0.6.1 reviewer flagged (rs115232898 TT,
  rs1801266 GG, rs3918290 CC and ~17 others of the same shape) now
  classifies correctly. Reference-allele homozygotes do not emit.

### Added
- **ADR-0018: PharmGKB per-allele function via CPIC template
  extraction.** Refines ADR-0016 (Data Classification Principle) by
  drawing a line: bounded canonical templates with enumerated fields →
  structured extraction (allowed); arbitrary prose intent inference →
  forbidden. The CPIC sentence
  `"The {A|C|G|T} allele of {rsid} is assigned [a] {normal|decreased|no|increased} function [allele] by CPIC."`
  is parsed at load time.
- **New SQLite table `pharmgkb_allele_function`.** Stores
  `(rsid, allele) -> function_class` with a `source` discriminator so
  future per-allele function sources (e.g. var_fa_ann aggregation or a
  new PharmGKB schema column) populate the same table.
- **Pre-pass / main-pass loader.** `load_pharmgkb_tsv` now does a
  pre-pass over `clinical_ann_alleles.tsv` to build the global
  `(rsid, allele) -> function_class` lookup, populates the new table
  atomically, then re-iterates to classify each row using the lookup.
- **Hybrid classifier priority (`is_nonfinding_for_row`):**
  (1) structured `Allele Function` column → (2) per-allele CPIC lookup
  → (3) prose fallback (ADR-0017 residual). Tier 2 is the primary path
  for SNV rows where CPIC has published guidance. Tier 3 becomes inert
  wherever tier 2 has entries.
- **Mock fixture row PA-011 / rs900000111.** Three allele rows (CC, CT,
  TT) each carrying
  `"The C allele of rs900000111 is assigned decreased function by CPIC."`
  Exercises the per-allele lookup path end-to-end. Fixture now ships
  10 annotations + 21 allele rows = 17 stored records.
- **Tests:** `TestExtractCpicAlleleFunction` (direct regex unit tests),
  `TestBuildAlleleFunctionLookup` (fixture-integration), and
  `TestIsNonfindingByAlleleLookup` (per-allele classifier including
  the reviewer's exact leaker shapes).

### Migration
- v0.6.x caches lack `pharmgkb_allele_function`, so
  `schema_is_current()` returns False and `db update` refreshes
  automatically. Run `allelix db update --force` to re-ingest with the
  new classifier.

## [0.6.1] — 2026-05-13

> **Recovers v0.5.2's filtering.** v0.6.0 was a regression on real data —
> the structured classifier keyed off `Allele Function`, which PharmGKB
> populates only on haplotype rows (the loader rejects those). 100% of
> in-scope SNV rows had empty `Allele Function` → ~13,500 rows emitted in
> production. v0.6.1 restores filtering via a documented hybrid:
> structured signal when present, prose pattern set as the row-level
> fallback per ADR-0017.

### Fixed
- **PharmGKB classifier is now hybrid (ADR-0017).**
  `is_nonfinding_for_row(allele_function, annotation_text)` checks the
  structured field first; if empty, falls back to a bounded named prose
  pattern set. Structured-first invariant preserved wherever PharmGKB
  publishes the structured signal. If PharmGKB adds `Allele Function`
  on SNV rows in a future revision, the fallback becomes inert
  automatically.

### Changed
- **Mock fixture now models real PharmGKB shape (ADR-0015 + ADR-0017).**
  SNV genotype rows have `Allele Function = ""` (matches real PharmGKB).
  Haplotype rows retain populated `Allele Function` (matches real
  PharmGKB and exercises the structured path through unit tests). The
  earlier inverted-shape fixture was the proximate cause of the v0.6.0
  regression; the new fixture-shape invariant test
  (`test_snv_rows_have_empty_allele_function`) fails loudly if anyone
  inverts the shape again.

### Added
- **ADR-0017: PharmGKB SNV row-level prose fallback.** Documents the
  narrow ADR-0016 exception: when a database publishes a structured
  classification field but leaves it empty on every row the consumer
  code path processes, a bounded prose pattern set may serve as the
  row-level fallback. Four guardrails enforced (structured-first;
  named bounded patterns; real-shape invariant test; ADR-scoped
  boundary).
- **Data Classification Principle amendment.** New subsection defining
  the narrow row-level exception. The general rule is unchanged: regex
  on prose is forbidden for classification. The exception is precisely
  scoped.
- **Test: `test_snv_rows_have_empty_allele_function`.** Pins the
  real-data shape on the mock fixture. If a future fixture revision
  re-populates Allele Function on SNV rows, the test fails with a
  message naming the v0.6.0 regression and pointing at ADR-0015 +
  ADR-0017.
- **Test: `TestIsNonfindingForRow`.** Exercises the hybrid classifier
  across structured-wins, prose-fallback, and edge cases.

### Migration
- v0.6.0 caches (built between this morning and now) lack the new
  hybrid classification — every row has `is_nonfinding=0`. The schema
  is the same as v0.6.1, so `schema_is_current()` returns True and the
  freshness check would skip the rebuild. **Users on v0.6.0 must run
  `allelix db update --force`** to re-ingest with the hybrid
  classifier.
- v0.5.x caches → automatic refresh (schema check fails on missing
  `function_class`).

### Owed to follow-up
- A redacted slice of real PharmGKB committed to the repo and run
  through the pipeline in CI nightly. Mock-shape invariant is necessary
  but not sufficient; real-source-data smoke is the next layer. Tracked
  in ADR-0017; not blocking v0.6.1.
- Audit of `variantAnnotations.zip` / `genes.zip` for a structured
  per-SNV classification field. If one exists, replace the prose
  fallback with it and ADR-0017 becomes inert.

## [0.6.0] — 2026-05-13 [SUPERSEDED by 0.6.1]

> **Functional regression.** ADR-0016's structured classifier consumed
> PharmGKB's `Allele Function` field, which is empty on every SNV row.
> The annotator (per ADR-0009) processes only SNV rows. Net effect:
> every PharmGKB row emitted in production. Use v0.6.1 or later. The
> ADR-0016 principle is correct and remains in force; ADR-0017
> documents the row-level fallback that v0.6.0 lacked.

> **Architectural correction.** Classification by regex against prose was
> the wrong mechanism. v0.6.0 enforces the Data Classification Principle
> (ADR-0016): structured database fields
> are the only classification input. Regex against description text is
> permitted only in test safety nets.

### Changed
- **PharmGKB non-finding suppression now uses the structured `Allele
  Function` column.** The previous regex-on-prose mechanism (eight
  patterns in v0.5.0, expanded to ten in v0.5.2) is deleted from
  production code. `pharmgkb_loader.classify_function()` maps the
  authoritative field to a normalized enum (`normal` / `decreased` /
  `no_function` / `increased` / `unknown`). A row is a non-finding iff
  `function_class == "normal"`. Stable against PharmGKB editorial drift.
- **`is_nonfinding` is now a derived structured signal**, not a regex
  hit count.
- **Annotator SELECT** filters by `is_nonfinding = 0` and no longer
  references `is_somatic`.

### Removed
- `_NONFINDING_PATTERNS`, `_is_nonfinding(annotation_text)`,
  `_SOMATIC_PATTERN`, `_is_somatic(annotation_text)` — all four
  regex-on-prose classifiers excised from production code per ADR-0016.
- `is_somatic` column from `pharmgkb_annotations`. PharmGKB has no
  structured germline/somatic flag, so per ADR-0016 the decision cannot
  be automated. The common case (somatic-context rows describing
  reference genotypes) is correctly caught by the non-finding filter via
  `Allele Function = Normal function`; rare residuals surface for
  manual review.

### Added
- **ADR-0016: Data Classification Principle** — structured fields only,
  regex forbidden in production, regex permitted only as a test safety
  net. Codifies the project's non-negotiable architectural directive.
- `function_class TEXT NOT NULL` column on `pharmgkb_annotations`,
  storing the normalized enum from `classify_function()`. Required by
  `schema_is_current()`, so v0.5.x caches refresh automatically on next
  `db update`.
- `tests/test_pharmgkb_safety_net.py` — the one place the regex now
  lives. Runs as a canary: regex match must agree with the structured
  `is_nonfinding` column. Disagreement = loader bug, not regex bug.

### Schema migration
- v0.5.x caches lack `function_class` →
  `PharmGKBAnnotator.is_ready()` returns False → next `db update`
  refreshes. No user action beyond running `allelix db update` after
  upgrading.

### Documentation
- **ADR-0013 amended**: mechanism revised to structured-field-based;
  user-facing contract unchanged.
- **ADR-0014 superseded**: somatic suppression removed; structured
  signal does not exist in PharmGKB; gap documented.
- Full repo audit (recorded in commit message): two regex violations
  identified, both in `pharmgkb_loader.py`, both removed. Remaining
  regexes (`_RSID_RE`, `_TWO_LETTER_GENOTYPE_RE`, VCF header parsers,
  test fixture validators) are structural-format checks, not prose
  classification.

## [0.5.2] — 2026-05-13

### Fixed
- **PharmGKB non-finding suppression: two more patterns.** Real-data
  review on the v0.5.1 report surfaced 30+ rows that still leaked
  through ADR-0013's filter. The DPYD cluster was the most visible
  (15+ entries at magnitude 9.0 saying "Both variants of rsX are
  assigned normal function by CPIC"). Pattern set extended to include
  `assigned normal function` and `may not have altered risk`. The
  expanded patterns are pinned by new tests under
  `TestNonFindingClassifier`. Mock PharmGKB fixture gains two PA-008
  allele rows demonstrating the new patterns (ADR-0015 contract).
- New regression tests pin that the classifier does NOT over-filter:
  protective findings ("decreased risk of neutropenia") and clinical
  dosing guidance ("require a decreased dose of warfarin") are
  preserved as findings.

### Added
- **`allelix extract --snps rs1,rs2,...`**: spot-check diploid
  genotypes at specific rsids without running a full analyze. Implements
  the CLI command that had been on the deferred list.
  Useful for verifying ClinVar/PharmGKB hits against the actual file
  before trusting them — particularly for the residual
  high-magnitude ClinVar hits flagged in the v0.5.1 review (NIPA1,
  NCR3, IL10, GP9, ADRA2A, TMPRSS6).

### Migration
- No schema change. v0.5.1 PharmGKB cache is compatible with v0.5.2,
  but the new patterns only apply to rows ingested after upgrade. Run
  `allelix db update --force` to re-ingest with the expanded
  classifier and drop the leaked non-finding rows.

## [0.5.1] — 2026-05-11

> **Process correction.** Closes the underlying failure mode that allowed
> the v0.4.2 and v0.5.0 clinical-safety bugs to ship. See ADR-0015.

### Fixed
- MHG mock generator wrote `CTT/C` for rs113993960 (CFTR ΔF508). Real
  MyHappyGenes arrays cannot call indels — every genotype is one of
  `A`/`T`/`G`/`C` or the no-call marker `-`. The buggy entry violated
  the MHG format spec and is the proximate reason the ClinVar
  indel-anchor incident (ADR-0011) wasn't caught by tests before
  shipping in v0.4.0. Corrected to `-/-`.

### Added
- **ADR-0015: Mock data generators are the contract.** The mock
  generators are the canonical model of real source data. Code that
  doesn't work against generator output is buggy; generators that don't
  model real data are buggy. Hand-authored ad-hoc fixtures are
  forbidden.
- `tests/test_mock_data_invariants.py`: structural assertions on every
  generator's output (MHG must produce single-base or no-call alleles
  only; ClinVar must include both SNVs and indels and multi-allelic
  rows; PharmGKB must include carrier-findings, non-findings, and
  somatic rows). Drift in a generator now fails the suite explicitly.
- `tests/test_end_to_end.py`: full pipeline against the canonical mock
  generators with snapshot assertions. Pins the absence of every
  categorical bug closed across v0.4.2/v0.5.0 (CFTR indel-anchor leak,
  TP53 wild-type leak, PharmGKB non-finding leak, PharmGKB somatic
  leak) plus a numeric snapshot. Any drift in the snapshot demands code
  review, not a blind update.

### Changed
- `test_indel_variant` renamed to `test_cftr_indel_position_is_no_call`
  and rewritten to assert the corrected fixture behavior. The old test
  asserted `CTT/C` was a valid MHG genotype — exactly the assumption
  that hid the bug.
- Stats count assertions updated for the corrected fixture: 102 no-call
  → 103, 564 het → 563. Hom unchanged.
- End-to-end snapshot: 7 ClinVar + 4 PharmGKB annotations against the
  corrected fixture. CFTR no longer fires (no-call short-circuits
  before the indel-anchor check).

## [0.5.0] — 2026-05-11

> **Clinical-safety release.** v0.5.0 closes three categorical false-positive
> bugs in the annotator layer. v0.4.2 already shipped the ClinVar indel
> fix (#1 below); v0.5.0 adds the two PharmGKB fixes (#2 and #3) and is
> the first version where reports against array-based parsers can be
> read with confidence. **Treat all v0.4.x and earlier reports as
> untrusted; regenerate against v0.5.0.**

### Fixed (clinical safety)

1. **ClinVar: indel-anchor false positives** *(originally shipped in
   v0.4.2; restated for users skipping straight from v0.4.1).* ClinVar
   encodes indels with an anchor base (e.g. `REF=AT ALT=A` for a single-T
   deletion). Array-based parsers report single-base genotypes. The
   carrier rule treated the anchor character as equivalent to the array
   readout and emitted hundreds of false-positive "Pathogenic" calls in
   cancer-predisposition genes (MSH6, APC, PTEN, MLH1, MSH2, RB1, BRCA1,
   BRCA2, TP53, …). Indel rows are now suppressed when the user's
   genotype is single-base only. See ADR-0011.

2. **PharmGKB: non-findings emitted as findings.** PharmGKB stores one
   row per genotype at each variant position, including the reference
   (wild-type) genotype, with text like "do not have a copy of the
   variant" or "decreased but not absent risk". These rows were emitted
   at the same magnitude as real variant calls. Real-data measurement:
   ≥12% of PharmGKB matches in a production array report were explicit
   non-findings. Rows whose annotation text matches the non-finding
   pattern set are now classified at load time and suppressed at query
   time. See ADR-0013.

3. **PharmGKB: somatic-tumor annotations on germline data.** A subset of
   PharmGKB rows (EGFR T790M / L858R, BRAF/KRAS/PIK3CA tumor markers)
   describe somatic variants only present in tumor tissue. Consumer DNA
   tests sample germline DNA. These rows are now suppressed on germline
   parsers. See ADR-0014.

### Schema

- `pharmgkb_annotations` gains three columns: `allele_function TEXT`,
  `is_nonfinding INTEGER NOT NULL`, `is_somatic INTEGER NOT NULL`.
- `PharmGKBAnnotator.is_ready()` now consults
  `schema_is_current(db_path)` (PRAGMA-based check). v0.4.x caches
  return False → `db update` sees "not ready" and refreshes into the
  v0.5.0 schema. No silent degradation.

### Migration

- Treat all v0.4.x reports as untrusted. Regenerate against v0.5.0.
- Run `allelix db update` to rebuild the PharmGKB cache (no `--force`
  needed; the schema check forces a refresh automatically).

## [0.4.2] — 2026-05-11

### Fixed (CRITICAL — clinical safety)

- **ClinVar annotator emitted hundreds of false-positive "Pathogenic" calls
  in cancer-predisposition genes (APC, MSH6, MLH1, MSH2, PTEN, RB1, BRCA1,
  BRCA2, TP53, …) when run against array-based genotype files
  (MyHappyGenes, 23andMe, AncestryDNA).** Root cause: ClinVar encodes indels
  with an anchor base (e.g. `REF=AT ALT=A` for a deletion of T). Array
  parsers report single-base genotypes at probe positions. The carrier
  rule's string equality treated ClinVar's anchor character as equivalent
  to the array's readout, producing pathogenic calls for users who carried
  only the wild-type sequence. Indel rows are now skipped when the user's
  genotype is single-base only. Multi-base genotypes (future VCF parsers
  that actually call indels) are unaffected. See ADR-0011.
- **Affects every report generated by v0.4.0 and v0.4.1 against array data.
  Treat all v0.4.x outputs as untrusted until regenerated against v0.4.2.**

### Added

- **Freshness detection** (ADR-0012). `allelix db update` now detects when
  the remote source has changed and refreshes only when needed, without
  requiring `--force`. Each annotator implements `fetch_remote_signal()`
  that fetches a small published signal (NCBI's `clinvar.vcf.gz.md5` for
  ClinVar; HEAD-request `ETag`/`Last-Modified` for PharmGKB). Stored
  signals are type-prefixed (`md5:…`, `etag:…`, `lm:…`) so a server
  switching signal types triggers a refresh rather than a silent miss.
  Network helpers `fetch_remote_text` and `head_request_headers` swallow
  all `OSError`/`ValueError` and return `None`, so a flaky network never
  crashes `db update` — at worst the user sees "freshness can't be
  verified" and can `--force` to refresh.
- Schema migration: `database_versions` gains a nullable `remote_signal`
  column. `get_database_info` falls back to a 4-column SELECT for
  pre-v0.4.2 caches and reports `remote_signal=None`. The decision tree
  treats `None ≠ remote` as "refresh" so v0.4.1 caches auto-upgrade on
  first v0.4.2 `db update`.
- ADR-0011 (indel-anchor protection) and ADR-0012 (freshness detection).
- Three regression tests under `TestIndelAnchorProtection` pin the
  negative cases the carrier rule used to get wrong; existing
  `test_indel_carrier_triggers` covers the multi-base positive path.
- Tests for `fetch_remote_text` / `head_request_headers` against a local
  HTTPServer fixture; CLI tests for every freshness decision-tree branch
  (match, differ, unverifiable, legacy-cache, --force).

## [0.4.1] — 2026-05-11

### Fixed
- PharmGKB loader looked for `clinical_ann.tsv`, but the real PharmGKB
  `clinicalAnnotations.zip` ships `clinical_annotations.tsv` (plural).
  Loader updated; mock fixture renamed; new schema-pinning test snapshots
  the real zip file listing so a future PharmGKB rename trips a test
  rather than a `db update` failure.
- `allelix db update` now skips annotators where `is_ready()` is True;
  pass `--force` to refresh anyway. Previously every invocation
  re-downloaded every annotator.

### Changed
- Round-11 minors (originally listed under Unreleased; folded into 0.4.1):
  - **m-1** Removed dead `render_annotations` from `reports/terminal.py`;
    `render_terminal` is the single entry point. Test file rewritten to
    assert against `render_terminal` + `AnalysisResult` fixtures.
  - **m-2** Dropped the unreachable `else` branch in
    `_run_analysis_command`; renderer dispatch is now a ternary.
  - **m-4** `analyze --category` help text no longer lists "methylation"
    (no annotator emits that category — use `allelix methylation` instead).
  - **m-5** HTML and JSON outputs are now atomic via a shared
    `atomic_write_text` helper in `allelix.reports` (`.tmp` then
    `os.replace`). A killed process mid-write leaves either the previous
    file or no file, never a half-written one.
  - **m-6** `REGULATORY_NOTICE` moved from `json_report.py` to
    `allelix.reports`. Both `json_report` and `html` import from there.

### Added
- **m-3 coverage**: PA-007 fixture row in the synthetic PharmGKB dump
  exercises the inner `_normalize_genotype` skip path
  (single-rsid annotation with only non-SNV per-allele rows).
- `tests/reports/test_init.py` covers `REGULATORY_NOTICE` presence and
  `atomic_write_text` (success path + simulated-rename-failure cleanup).

## [0.4.0] — 2026-05-11

### Added
- **Reports (Phase 6).** New `allelix/reports/_pipeline.py` builds an
  `AnalysisResult` once and hands it to format-specific renderers, so the
  CLI streams the file and queries each annotator exactly once per `analyze`
  invocation regardless of output format.
- **JSON report** (`reports/json_report.py`): versioned schema (v1),
  embedded regulatory notice (ADR-0003), per-annotation source attribution
  preserved, applied filters echoed in the payload.
- **HTML report** (`reports/html.py`): single self-contained file with inline
  CSS, no external resources, attribution column, informational-only banner.
  XSS-safe via `html.escape` on every user-supplied string.
- **`allelix analyze --output report.{html,json}`** dispatches by extension;
  `--report-format html|json` overrides explicitly.
- **`allelix methylation`** subcommand: focused report filtered to the
  curated methylation gene panel (MTHFR, MTR, MTRR, COMT, CBS, BHMT, …).
- **`allelix pharmacogenomics`** subcommand: focused report filtered to
  `category=pharma` annotations (PharmGKB-style drug response).
- **Strand-flip helpers** (Phase 7 polish, partial): `allelix.utils.allele`
  ships `complement`, `flip_genotype`, `is_strand_ambiguous` for parser /
  annotator authors. Not wired into matching logic yet — opt-in.
- **ADR-0010**: documents the strand-flip decision and explicitly defers
  liftover to v1.0.0 with the rationale (chain-file weight, no real GRCh38
  consumer yet).

### Changed
- Annotators are entered into a `contextlib.ExitStack` inside the pipeline
  module, replacing the inline ExitStack in the CLI's `analyze` command.
  Same deterministic-cleanup contract (C-1) applies across all report
  commands now.
- `allelix/cli.py` factored: `_resolve_parser`, `_ready_annotators`,
  `_run_analysis_command` are shared by `analyze`, `methylation`, and
  `pharmacogenomics`.

### Fixed
- **R-1**: `test_version_flag` no longer hardcodes a version string; it
  asserts the rendered output contains `__version__` and a new
  `test_pyproject_version_matches_metadata` test in `tests/test_version.py`
  pins that `pyproject.toml`'s version equals the installed package
  metadata. Catches the regression class where someone bumps
  pyproject.toml without reinstalling the editable package.

### Added
- **R-2**: PharmGKB loader gets the W-1-style batch-flush spy test
  (`test_batched_insert_flushes`) — `INSERT_BATCH_SIZE = 3` against the
  9-record fixture asserts `[3, 3, 3]`. A mutation flipping the flush
  condition to `if False:` now fails the suite.
- Asymmetric no-call test (`Variant("A", "-")` and `("-", "A")`) for
  both ClinVar and PharmGKB annotators (r-2). Verified that the ClinVar
  test catches `if variant.is_no_call:` → `if variant.allele1 == "-":`.
- `_safe_float` direct unit tests cover the empty-string and ValueError
  branches in `pharmgkb_loader` (r-1).
- `test_leftover_tmp_file_is_cleared` covers the stale-`.tmp` cleanup
  branch in `load_pharmgkb_tsv` (r-5).

### Changed
- `pharmgkb_loader._open_directory` now imports `Path` at module top
  instead of inline-aliasing inside the function (r-3). `extractall`
  now has an explicit code comment about Python 3.11+'s built-in
  `..`/absolute-path sanitization (r-4).

## [0.3.0] — 2026-05-11

### Added
- **PharmGKB annotator** (Phase 3, first cut). `allelix db update` now also
  downloads PharmGKB's `clinicalAnnotations.zip` into `pharmgkb.sqlite`;
  `allelix analyze` layers in pharmacogenomic annotations alongside ClinVar's
  clinical-significance calls.
- `allelix.databases.pharmgkb_loader` parses the joined
  `clinical_ann.tsv` + `clinical_ann_alleles.tsv` view, normalizes 2-letter
  SNV genotypes (sorted, uppercased), and skips star alleles, multi-rsid
  composites, and indel genotypes (Phase 7+).
- `PharmGKBAnnotator` matches the user's exact normalized diploid call (per
  ADR-0009), with magnitude derived from PharmGKB Level of Evidence
  (1A → 9.0, 1B → 8.0, …, 4 → 2.0; per ADR-0008).
- ADR-0009 documents the per-genotype matching rule that diverges from the
  ClinVar carrier rule (ADR-0007).
- Synthetic PharmGKB fixture (`tests/fixtures/mock_pharmgkb/`) with two TSVs
  whose genotypes line up with the MyHappyGenes mock fixture so end-to-end
  `analyze` tests fire on known carriers.
- 37 new tests (162 total): genotype normalization, single-rsid filtering,
  ZIP + directory ingestion, atomic SQLite load, end-to-end `db update` /
  `db status` / `analyze` against both annotators.

### Changed
- `db status` table now shows record counts for any annotator that exposes
  `record_count()` (no longer a ClinVar-specific switch in the CLI).
- `allelix/databases/schema.py` factors the shared `database_versions`
  table into a constant; both `CLINVAR_SCHEMA` and `PHARMGKB_SCHEMA`
  embed it so `get_database_info` works uniformly across annotators.

### Added
- `.pre-commit-config.yaml` and `pre-commit` in dev deps. Two hook stages
  installed via `pre-commit install --hook-type pre-commit --hook-type pre-push`:
  - **pre-commit** runs `ruff check` + `ruff format --check` so a commit
    that doesn't lint or format clean is blocked.
  - **pre-push** runs the version-tag check only (fast, no test suite).
    The full test suite runs in CI on every pull request.
- Annotator test fixture now closes connections on teardown (round-4 N-1).
- `pytest.PytestUnraisableExceptionWarning` is escalated to error in
  `filterwarnings`, so future fixtures without teardown fail CI loudly.
- Tests covering the INSERT batch flush path, the download timeout path,
  the `parse_clinvar_version` end-of-file fall-through, and friendly
  `db update` failure output.
- Round 5 (W-1, W-2): the batch-flush test now spies on `executemany`
  via a delegating connection proxy so a mutation that disables
  mid-iteration flushing fails the suite. The `db update` friendly-error
  test now asserts Click's `Error: clinvar:` prefix and that
  `result.exception` is not a `URLError`, so removing the
  `ClickException` wrap now fails the suite.
- More cheap coverage: VCF rows with <8 columns, `_chrom_sort_key`
  unknown-chromosome fallback, `_percent(0, 0)`, MyHappyGenes
  junk-before-header warning, and metadata blank-line skip.

### Changed
- `db update` now wraps `annotator.setup()` in a `ClickException`, so a
  network/disk/parse failure shows a one-line "clinvar: …" message
  instead of a Python traceback.

## [0.2.0] — 2026-05-11

### Added
- **First annotator: ClinVar.** Downloads the GRCh37 ClinVar VCF (`allelix db update`), parses it into a local SQLite cache, and annotates carriers — `allelix analyze` flags variants where the user carries the ALT allele (ADR-0007).
- **Multi-allelic VCF rows** are split into one record per ALT during parse, with parallel-indexed CLNSIG/CLNDN/ALLELEID.
- **Indel matching** for variants like CFTR ΔF508 (REF=CTT ALT=C).
- `allelix db update` and `allelix db status` commands. Status table shows version (from VCF `##fileDate`) and record count.
- `allelix analyze` with `--min-magnitude` and `--category` filters; reports versions used and an annotation count.
- Data-directory resolution: `--data-dir` > `$ALLELIX_DATA_DIR` > `$XDG_DATA_HOME/allelix` > `~/.local/share/allelix` (ADR-0006).
- Annotators implement the context-manager protocol; the CLI uses `contextlib.ExitStack` for deterministic SQLite-connection cleanup.
- `allelix/py.typed` marker (still present from v0.1.0; restated for completeness).
- ADRs 0006 (data-dir), 0007 (carrier rule), 0008 (Allelix-derived magnitude scoring from CLNSIG).
- `pytest-cov` with a 92% coverage floor; `filterwarnings = ["error::ResourceWarning"]` in pytest config so future leaks fail CI.
- README adds "Supported Databases" table and a regulatory posture summary.

### Changed
- `download()` is now atomic (`.part` then `os.replace`), times out after 60s, sends an `allelix/<version>` User-Agent, and `fsync`s before rename.
- `load_clinvar_vcf()` is now atomic (writes to a `.tmp` SQLite then `os.replace`s); a failed mid-parse leaves the previous cache intact.
- `GenotypeMetadata` no longer carries `snp_count`; the only authoritative source is `parse()` (ADR-0005, restated).
- Annotation `category` field is documented as a non-diagnostic filter bucket; never bare medical terms.

### Removed
- The old non-atomic `db_path.unlink()` step in `load_clinvar_vcf` that destroyed the previous cache before writing the new one.

## [0.1.0] — 2026-05-11

### Added
- Project scaffolding: `pyproject.toml`, package layout, MIT `LICENSE`, comprehensive `.gitignore`.
- Core models: `Variant` (with `is_heterozygous`, `is_no_call`, `genotype` properties) and `Annotation` (with required `attribution` field).
- Plugin parser architecture: `GenotypeParser` ABC, registry with auto-detection, `MyHappyGenesParser` (streaming, malformed-line tolerance, logged warnings).
- `allelix stats` CLI command with parser-warning surfacing and a clean logger lifecycle.
- Synthetic MyHappyGenes test fixture and deterministic generator script.
- README with regulatory posture, data sources & licensing table.
- ADRs 0001–0005 documenting the meta process, plugin architecture, source-attributed annotations, offline-first data model, and the parse-derived SNP count contract.
- GitHub Actions CI matrix on Python 3.11 and 3.12.


[2.0.0]: https://github.com/allelix/allelix/compare/v1.9.0...v2.0.0
[1.9.0]: https://github.com/dial481/allelix/compare/v1.8.4...v1.9.0
[1.8.4]: https://github.com/dial481/allelix/compare/v1.8.3...v1.8.4
[1.8.3]: https://github.com/dial481/allelix/compare/v1.8.2...v1.8.3
[1.8.2]: https://github.com/dial481/allelix/compare/v1.8.1...v1.8.2
[1.8.1]: https://github.com/dial481/allelix/compare/v1.8.0...v1.8.1
[1.8.0]: https://github.com/dial481/allelix/compare/v1.7.0...v1.8.0
[1.7.0]: https://github.com/dial481/allelix/compare/v1.6.1...v1.7.0
[1.6.1]: https://github.com/dial481/allelix/compare/v1.6.0...v1.6.1
[1.6.0]: https://github.com/dial481/allelix/compare/v1.5.3...v1.6.0
[1.5.3]: https://github.com/dial481/allelix/compare/v1.5.2...v1.5.3
[1.5.2]: https://github.com/dial481/allelix/compare/v1.5.1...v1.5.2
[1.5.1]: https://github.com/dial481/allelix/compare/v1.5.0...v1.5.1
[1.5.0]: https://github.com/dial481/allelix/compare/v1.4.0...v1.5.0
[1.4.0]: https://github.com/dial481/allelix/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/dial481/allelix/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/dial481/allelix/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/dial481/allelix/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/dial481/allelix/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/dial481/allelix/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/dial481/allelix/releases/tag/v1.0.0

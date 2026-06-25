<div align="center">

<img src="https://allelix.io/apple-touch-icon.png" width="180" alt="Allelix" />

**Open-source genotype analysis. Format-agnostic. Offline-first.**

[Website](https://allelix.io) · [Quickstart](#quickstart) · [Changelog](https://github.com/allelix/allelix/blob/main/CHANGELOG.md)

[![python](https://img.shields.io/pypi/pyversions/allelix.svg)](https://www.python.org/downloads/)
[![pypi](https://img.shields.io/pypi/v/allelix.svg)](https://pypi.org/project/allelix/)
[![license](https://img.shields.io/pypi/l/allelix.svg)](LICENSE)
[![CI](https://github.com/allelix/allelix/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/allelix/allelix/actions/workflows/ci.yml)
[![privacy: local-first](https://img.shields.io/badge/privacy-local--first-success)](#privacy)

</div>

# Allelix

> **Status:** Production. Eight parser formats (including VCF + gVCF), four annotators (ClinVar + ClinPGx + GWAS Catalog + SNPedia), three enrichment sources (gnomAD + AlphaMissense + CADD), licensable-source gating for commercial users, dual-build ClinVar caches (GRCh37 + GRCh38), HTML/JSON/terminal reports, methylation + pharmacogenomics focused commands, report diffing, persistent config with commercial-mode safety switch. Build auto-detection from position data (ADR-0021). No regex on prose anywhere in production. The [Changelog](https://github.com/allelix/allelix/blob/main/CHANGELOG.md) tracks every release.

## Quickstart

**Requires Python 3.11 or newer.** On Python 3.10 or older, `pip install allelix` succeeds but the first command fails with `ImportError: cannot import name 'UTC' from 'datetime'` — that's the symptom; the fix is to install on 3.11+. Check with `python --version`.

```bash
pip install allelix

# Download reference databases (~16 GB default install; ~22 GB if you opt into CADD).
# Use --no-gnomad / --no-alphamissense to skip the large ones.
# CADD is opt-in: allelix db update --cadd
allelix db update

# Analyze a genotype file
allelix analyze your_genotype_file.txt --output report.html

# VCF / gVCF input — same command, auto-detected
allelix analyze your_wgs.vcf.gz --output report.html

# Multi-sample VCF — pick which sample to analyze
allelix analyze trio.vcf.gz --sample HG002 --output report.html

# Filter to a custom panel (rsIDs + gene names, one per line; '#' comments and blank lines ignored)
allelix analyze your_genotype_file.txt --filter-file my_panel.txt --output report.html
```

See [Development](#development) for source installs and running tests, [Managing your data](#managing-your-data) for cache locations and cleanup, and [Troubleshooting](#troubleshooting) for common failure states.

## Supported Formats

| Format | Status | Notes |
|---|---|---|
| MyHappyGenes (Tempus) | ✓ | Tab-delimited, 5 columns. **Build is auto-detected** — real-world MHG exports mislabel the header as "build 37.1" while shipping GRCh38 coordinates. Allelix detects from position data and warns on header/data disagreement (ADR-0021). |
| 23andMe | ✓ | Tab-delimited, 4 columns, concatenated genotype. Supports build 36/37/38 from header. I-prefixed probe IDs passed through. |
| AncestryDNA | ✓ | Tab-delimited, 5 columns. Chromosome mapping: 23→X, 24→Y, 25→X (PAR), 26→MT. V1 and V2 chip layouts. |
| Family Tree DNA | ✓ | CSV, double-quoted fields, concatenated genotype. Build 37 default. |
| MyHeritage DNA | ✓ | CSV, same structure as FTDNA. Detected by "MyHeritage" in comment header. Handles double-double-quoted field variant. |
| Living DNA | ✓ | Tab-delimited despite `.csv` extension. Handles AX-, AFFX-prefixed and CHR:POS positional SNP IDs. |
| FTDNA Illumina raw | ✓ | Tab-delimited variant of the FTDNA export (distinct from the CSV format above). `RSID/CHROMOSOME/POSITION/RESULT` columns. Build 37 default. |
| FTDNA FamFinder | ✓ | Third FTDNA file shape: tab-delimited with **separate `ALLELE1` / `ALLELE2` columns** instead of a concatenated `RESULT`. Detected by the `famfinder` substring in the file header plus the 5-column canonical header. Build 37 default. |
| VCF / gVCF | ✓ | REF/ALT encoding, `0/1` genotype notation. Plain VCF: absence at a position means reference. gVCF: explicit reference blocks (lines with `<NON_REF>` ALT and `END=` INFO) are skipped — they match nothing in any annotation database. Multi-sample files require `--sample <ID>`. Streams via stdlib; `.vcf.gz` handled transparently. Optional `pip install allelix[vcf-index]` enables pysam-backed tabix random access for fast `extract --snps` on huge VCFs. |

Adding a new format means adding one file to `allelix/parsers/` and registering an instance in the `PARSERS` list in `allelix/parsers/__init__.py`.

### Roadmap

| Release | Theme |
|---|---|
| **v2.2** | Annotation model expansion (cross-source conflict surfacing, risk-allele display, structured MTAG flag), custom-panel UX improvements, methylation panel retune + citations, ClinVar source switch (VCF → `variant_summary.txt`), CLI cleanups (`db clean` / `db path`). |
| **v2.3** | Per-source magnitude decomposition, Good / Bad / Neutral annotation repute, JSON v3 (AI-legible output contract), plausibility flagging via zygosity × gnomAD MAF, ClinVar review-status weighting, functional-medicine +/− notation. |
| **v2.4** | ACMG/AMP automated classification engine, variant consequence annotation (VEP-equivalent). |
| **v2.5+** | CNVs and repeat-expansion support, dbSNP resolution for rsID-less VCFs, gnomAD genomes (non-coding population frequencies), supplemental genotype file merging. |
| **v3.0** | Ancestry estimation (PCA + 1KG/HGDP reference panels), Polygenic Risk Score (PRS) integration. |

Issues and milestones are tracked at [github.com/allelix/allelix/issues](https://github.com/allelix/allelix/issues). Annotator-level strand-aware carrier matching (R-1) shipped in v2.1.0's ADR-0035 cluster.

## Supported Databases

| Database | Status | Notes |
|---|---|---|
| ClinVar (GRCh37 + GRCh38) | ✓ | Public domain (NCBI). SNVs + indels + multi-allelic sites. **Both builds cached**; `analyze` dispatches by detected build (ADR-0021). Carrier rule (ADR-0007) requires the user to carry the ALT allele. Indel-anchor protection (ADR-0011) prevents single-base array readouts from matching anchor-base indels. |
| ClinPGx (formerly PharmGKB) | ✓ | CC BY-SA 4.0. Clinical annotations only — single-rsid SNVs; star alleles and haplotypes deferred (ADR-0009). **Primary non-finding filter is the ClinVar REF carrier rule (ADR-0023):** if ClinVar publishes a single-base REF for the rsid and the user is homozygous for it, the row is suppressed. CPIC's `(rsid, base) → function_class` join (ADR-0020) survives as a secondary tier for rsids ClinVar doesn't catalog. Earlier prose tiers (ADR-0013, ADR-0017, ADR-0018) are superseded. |
| CPIC (per-allele function table) | ✓ | Internal data source for the ClinPGx filter. Fetched from `api.cpicpgx.org` at `db update` time. Used to populate the `pharmgkb_allele_function` table — not surfaced to end users as its own annotator. |
| SNPedia | ✓ | CC BY-NC-SA 3.0 US. Pre-built cache downloaded via `db update` (~216K wiki pages, ~105K genotype rows). If the SNPedia database is absent, analysis runs without it. For commercial use, pass `--exclude-snpedia` — `analyze` runs using all other databases and omits SNPedia annotations. The cache can also be rebuilt from source via `scripts/scrape_snpedia.py` + `scripts/parse_snpedia.py`. |
| GWAS Catalog | ✓ | Public domain (EBI/NHGRI). Trait–SNP associations with p-values and effect sizes. Carrier rule (ADR-0007) requires the user to carry the risk allele. P-value magnitude scoring (ADR-0024) maps continuous p-values to the 0–10 scale; unknown-risk-allele entries fire on rsID match alone but are capped at 3.0. |
| gnomAD | ✓ | ODbL v1.0. **Enrichment annotator** — adds population allele frequency context to existing annotations. Shows how common each variant is in the general population (~16M exome variants from 730K individuals). A pathogenic variant that 35% of people carry reads very differently from one seen in 0.001%. Pre-built cache downloaded via `db update` (~6GB on disk). Use `--no-gnomad` to skip. |
| AlphaMissense | ✓ | CC BY 4.0. **Enrichment annotator** — adds DeepMind's protein-structure-based pathogenicity predictions to existing annotations. Scores 71M missense variants on a 0–1 scale: <0.34 = likely benign, >0.564 = likely pathogenic. Complements ClinVar's expert classifications with computational predictions — especially valuable for variants ClinVar hasn't reviewed yet. Pre-built cache downloaded via `db update` (~8GB on disk). Use `--no-alphamissense` to skip. |
| CADD | ✓ | LicenseRef-CADD (non-commercial). **Enrichment annotator** — adds PHRED-scaled deleteriousness scores from CADD v1.7. **Opt-in**, disabled by default; see [CADD modes](#cadd-modes) below for details. |

### CADD modes

CADD ranks how deleterious any single-nucleotide variant is using 100+ annotation tracks (coding, non-coding, regulatory). PHRED 10 = top 10% most deleterious, 20 = top 1%, 30 = top 0.1%.

**Opt-in**: enable via `allelix db update --cadd` or `allelix config set sources.cadd true`. Use `--no-cadd` to skip enrichment for a single run.

Two modes:

- **Cache mode** (default): pre-built ~5.8 GB SQLite cache (~120M variant keys). Covers nearly every position allelix can annotate from its other databases (gnomAD, AlphaMissense, ClinVar). For genotyping chip data (23andMe, AncestryDNA, MyHappyGenes, etc.), this mode is functionally complete — chip probes overwhelmingly target known, cataloged variants.
- **Full mode** (`options.cadd_full = true`): tabix queries against the complete CADD v1.7 file via pysam (GRCh38 only). Adds coverage for novel or private variants that appear only in WGS/WES data and are not in any pre-computed database. Requires `pip install allelix[cadd]` for the pysam dependency.

If your input is a genotyping chip file, cache mode is all you need.

### Build coverage asymmetry (GRCh37 vs GRCh38)

ClinVar dispatches per-build (ADR-0021) and ships with both GRCh37 and GRCh38 caches. The two caches are essentially equivalent in coverage: 2,896,063 rows / 2,645,206 distinct rsIDs in GRCh37 vs 2,896,102 / 2,645,243 in GRCh38 — a difference of 39 rows.

Despite that equivalence, the same person's WGS file produces noticeably more annotations as GRCh37 than as GRCh38. The mechanism is in the resolution step, not in upstream-data shape. Position-keyed rsID resolution requires exact `(chromosome, position, ref, alt)` alignment between the user's variant call and ClinVar's stored row. Lift-over between builds does not preserve that alignment perfectly: the `~0.4%` of the genome where the reference assembly was rebuilt has different REF alleles, multi-allelic sites split differently, and some benchmark VCF positions drop out entirely in the GRCh38 lift. Each misalignment loses one resolution, which in turn loses all the rsID-keyed downstream annotations that rsID would have driven (ClinVar's own carrier annotation, plus GWAS Catalog, SNPedia, and ClinPGx).

Real GIAB HG002 benchmark, surviving the default `--min-magnitude 5.0` filter: GRCh37 surfaces 520 distinct rsIDs across all sources, GRCh38 surfaces 341. The two sets overlap on 331 rsIDs; 189 are GRCh37-only and 10 are GRCh38-only — pure asymmetric loss in the GRCh38 lift, not different upstream coverage. The unfiltered totals (65,965 vs 4,867) magnify the same pattern at lower magnitudes, mostly via GWAS-Catalog weak-association rows.

If you have a choice of build for the input, GRCh37 surfaces more annotations today on rsID-less VCFs that flow through position-keyed resolution. GRCh38 still surfaces every ClinVar carrier hit it has an exact alignment for.

### Known ClinPGx limitation: reference-genotype rows where ClinVar and CPIC both lack data

ADR-0022 + ADR-0023: a tiny residual of ClinPGx rows may appear in reports even when the user is homozygous reference. ClinPGx publishes one annotation per genotype including the reference homozygote, and for the reference-homozygote row to be suppressed Allelix needs structured data on the variant from either:

- **ClinVar's REF allele** (the primary filter — see ADR-0023). Covers any rsID ClinVar catalogs.
- **CPIC's per-allele function table** (the secondary fallback — see ADR-0020). Covers rsIDs CPIC has classified.

For the rare rsID where ClinPGx has an annotation but *neither* ClinVar nor CPIC has data, the row emits. These are identifiable by a homozygous-reference genotype combined with "decreased risk," "may have a typical response," or similar comparative language. They are an upstream data gap, not an Allelix bug — we surface them honestly rather than hide them behind a curated exclusion list (which would recreate the maintenance trap the v0.5–v0.7 prose filters were trying to escape).

The CFTR × ivacaftor leak (~30+ rows on real data, pre-v0.7.3) is fixed by the ADR-0023 ClinVar REF check: CPIC's CFTR vocabulary (`"ivacaftor responsive"`) doesn't match the four-class enum the secondary tier expects, but ClinVar publishes REF for every CFTR rsID, so the primary tier catches them universally.

### Known ClinVar upstream data quality issues

Two ClinVar rows in real-world reports are known upstream artifacts, not Allelix bugs:

- **PKD1 rs199476100 GG (Pathogenic/Likely pathogenic, magnitude 8.5).** This is a stop-gained variant with a gnomAD frequency of 0.0005% (7 observations in 1.38 million chromosomes). Homozygosity for this variant is biologically implausible — PKD1 is autosomal dominant and the nonsense variant would be embryonic-lethal or devastating in homozygous state. The chip genotyping call is almost certainly a probe artifact. The code correctly reports what ClinVar says and what the chip reads; the error is upstream of Allelix. Future work: population-frequency filtering could flag ultra-rare variants where the chip call is likely unreliable.

- **IL10 rs1800896 CT (Pathogenic, magnitude 9.0).** This is a common polymorphism (MAF ~20–40%) in the IL-10 promoter. ClinVar's Pathogenic classification comes from a single submitter for hepatitis C susceptibility; a second submitter classifies the same allele as "Uncertain risk allele" for leprosy susceptibility. The ClinVar VCF aggregates across conditions, so the report may pair the Pathogenic classification with the wrong condition. Future work: ClinVar review-status weighting (number of submitters, star rating) could down-weight single-submitter classifications on common variants.

Neither issue affects Allelix's filter logic. Both are inherent to ClinVar's aggregation model and the limitations of array-based genotyping chips.

## Regulatory Posture

Allelix is an informational research tool. It reports classifications made by external databases. It does not independently classify variants, diagnose conditions, or make health recommendations. All variant significance is attributed to its source — Allelix says "ClinVar classifies this variant as pathogenic," never "this variant is pathogenic."

This is not a disclaimer afterthought. It is a design constraint that affects model naming, report wording, and category labeling throughout the codebase.

## Privacy

- No data leaves your machine. No telemetry. No uploads. No analytics.
- Reference databases are downloaded via `allelix db update` and cached locally.
- Analysis runs offline against local database caches. A brief freshness check runs before analysis by default (skipped with `--no-update`).

### Output files contain real annotations of your genome

The JSON / HTML / terminal output of `allelix analyze` and its
focused subcommands contains real annotations against your specific
variants — drug-response calls, carrier-status flags, hereditary-
disease findings. Wherever you write them via `--output <path>`,
that's where they sit until you delete them. Allelix doesn't
auto-clean and won't warn you when you write to `/tmp/` or any
other shared location. Treat the files as personal data: read them,
move them somewhere you control, or delete when you're done.
`allelix db clean` (shipped in v2.2.0) handles the cache side of
the lifecycle; output reports are still the user's responsibility
to manage.

## Configuration

Allelix stores persistent configuration in `config.toml` (in the data directory, default `~/.local/share/allelix/`). A default config is created on first run.

```bash
# View current config (annotated with license notes)
allelix config show

# Read a single key
allelix config get sources.cadd
allelix config get license.commercial

# Disable a source permanently
allelix config set sources.gnomad false

# Enable commercial mode (auto-disables non-commercial sources)
allelix config set license.commercial true

# Assert that you hold a commercial CADD license
allelix config set license.cadd true
```

CLI flags (`--no-gnomad`, `--no-alphamissense`, `--no-cadd`, `--exclude-snpedia`, `--cadd`) override the config for a single run. The config sets the baseline; flags override per-invocation.

### Database sizes and download times

Not all databases are equal in size. `allelix db update` downloads them all by default, but you can skip the large ones if disk space or bandwidth is a concern:

| Database | On disk | Download time | What it adds |
|---|---|---|---|
| ClinVar (GRCh37 + GRCh38) | ~900MB | 1–2 min | Core clinical variant classifications. Required. |
| ClinPGx + CPIC | ~6MB | seconds | Drug-gene interactions. |
| GWAS Catalog | ~200MB | 1–2 min | Trait-SNP associations from genome-wide studies. |
| gnomAD | ~6GB | 5–15 min | Population allele frequencies (how common is this variant?). |
| AlphaMissense | ~8GB | 5–15 min | Missense pathogenicity predictions (how likely to break protein function?). |
| CADD (opt-in) | ~5GB | 5–15 min | Variant deleteriousness scores (how damaging is this variant?). Enable with `--cadd`. |

gnomAD and AlphaMissense are the largest but add the most interpretive context. gnomAD answers "is this variant rare or common?" — a pathogenic variant carried by 35% of the population reads very differently from one seen in 3 people. AlphaMissense answers "does this missense change likely damage the protein?" — especially valuable for the thousands of variants ClinVar hasn't reviewed yet.

To skip either during download: `allelix db update --no-gnomad --no-alphamissense`. To disable permanently: `allelix config set sources.gnomad false`.

### Managing your data

Allelix stores reference database caches at `~/.local/share/allelix/` (Linux/macOS) or the XDG-equivalent on other platforms. Default install (no CADD) is ~16 GB on disk; with the CADD opt-in it's ~22 GB.

```
~/.local/share/allelix/
├── alphamissense.sqlite       ~7.8 GB
├── gnomad.sqlite              ~6.1 GB
├── cadd.sqlite                ~5.8 GB   (only if you opted into --cadd)
├── clinvar.GRCh37.sqlite      ~470 MB
├── clinvar.GRCh38.sqlite      ~480 MB
├── gwas_catalog_associations.tsv  ~700 MB
├── snpedia.sqlite             ~390 MB
├── gwas.sqlite                ~190 MB
└── pharmgkb.sqlite            ~14 MB
```

**The caches are disposable and re-downloadable.** Everything here was fetched by `allelix db update` and can be re-fetched at any time. There's no per-user state in the cache — your genotype files never enter this directory. Safe to delete to reclaim space (or to force a full refresh if you suspect a cache went stale):

```bash
# Reclaim all ~16-22 GB, with size report and confirmation prompt
allelix db clean

# Preserve the CADD cache (avoids re-downloading the slowest one, ~5.8 GB)
allelix db clean --keep-cadd

# Show what would be deleted without acting
allelix db clean --dry-run

# Scripted use — skip the confirmation prompt
allelix db clean --yes

# Re-populate from scratch next time you analyze a file
allelix db update
```

`db clean` refuses to act on a directory that doesn't contain any recognized allelix cache files — a typo'd `--data-dir ~/Documents` won't result in `rm -rf` against an unrelated tree, even with `--yes`. Pass `--force` if you've staged caches under a non-standard name and want to delete them anyway.

For scripting and backup integration, `allelix db path` prints the resolved data directory:

```bash
ALLELIX_DATA=$(allelix db path)
du -sh "$ALLELIX_DATA"

# --check additionally verifies the path is writable; exit non-zero if not
allelix db path --check
```

**The one thing worth backing up: `config.toml`.** It lives separately (XDG config dir, typically `~/.config/allelix/config.toml`) and stores your license assertions (commercial-mode toggle, CADD opt-in confirmation) and per-source enable/disable settings. It's the only non-reproducible state in the install. The caches will rebuild on the next `db update`; `config.toml` won't. `allelix db clean` preserves `config.toml` even when it lives inside the data directory (older installs).

To fully uninstall allelix and reclaim everything:

```bash
pip uninstall allelix
allelix db clean --yes            # reference database caches
rm -rf ~/.config/allelix/         # license + source-toggle config
```

## Data Sources & Licensing

Allelix source code is licensed under the **GNU Affero General Public License v3.0 or later** (AGPL-3.0-or-later). Allelix ships with **zero third-party data**. All reference databases are downloaded by the user at runtime via `allelix db update`. Each database retains its original license on the user's machine:

| Database | Source | License | Usage |
|---|---|---|---|
| ClinVar | NCBI | Public domain | No restrictions |
| GWAS Catalog | EBI/NHGRI | Public domain | No restrictions |
| ClinPGx (formerly PharmGKB) | clinpgx.org | CC BY-SA 4.0 | Attribution required |
| CPIC | cpicpgx.org | CC BY-SA 4.0 | Attribution required. Per-allele function data fetched from `api.cpicpgx.org` at `db update` time; used internally for the ClinPGx non-finding filter (ADR-0020), not surfaced as its own annotator. |
| SNPedia | snpedia.com | CC BY-NC-SA 3.0 US | Attribution required, **non-commercial only**. Use `--exclude-snpedia` to omit. |
| gnomAD | gnomad.broadinstitute.org | ODbL v1.0 | Attribution required. Population allele frequencies for context; not a clinical annotator. Use `--no-gnomad` to omit. |
| AlphaMissense | zenodo.org/records/10813168 | CC BY 4.0 | Attribution required. Cheng et al., Science 2023. Missense variant pathogenicity predictions. Use `--no-alphamissense` to omit. |
| CADD | cadd.gs.washington.edu | LicenseRef-CADD | Attribution required, **non-commercial by default**. Commercial licenses available from UW CoMotion. Opt-in via `allelix db update --cadd`. Use `--no-cadd` to omit. |

**Commercial users:** When `license.commercial = true`, non-commercial sources are gated by a three-state permission model. SNPedia is permanently blocked (no commercial license is available). CADD is blocked by default but can be unlocked — the University of Washington offers commercial licenses at `https://els2.comotion.uw.edu/product/cadd-scores`; after purchasing, assert your license with `allelix config set license.cadd true` to re-enable CADD in commercial mode. All other databases (ClinVar, ClinPGx, GWAS Catalog, gnomAD, AlphaMissense) are compatible with commercial use. `allelix config show` displays the permission state for each source.

### SNPedia data download

SNPedia data is downloaded automatically by `allelix db update` from a pre-built cache. If the SNPedia database is not present, `allelix analyze` runs normally using all other databases and prints a note that SNPedia data is not available.

To rebuild the cache from source (not normally needed):

```bash
python scripts/scrape_snpedia.py   # scrape 216K pages from bots.snpedia.com (1-4 hours)
python scripts/parse_snpedia.py    # parse raw wiki markup into structured genotype rows
```

### Known SNPedia source data quality notes

SNPedia appears frozen — no edits have been observed since mid-2023. The data below reflects the state of the wiki at scrape time (May 2026) and is unlikely to change.

Of the 104,806 genotype pages in the archive:

- **103 pages have empty or missing allele fields.** These are incomplete entries on the source wiki — the `{{Genotype}}` template was created but the `allele1`/`allele2` fields were never filled in (e.g., `Rs1131692198(;)` with `|allele1=\n|allele2=\n`). All 103 were verified against the live site on 2026-05-21; every one matches the source exactly. The annotator silently skips these — they cannot match any user genotype.

- **1 page has no `{{Genotype}}` template at all.** `Rs1799853(T)` is a malformed single-allele page (`{{is a|genotype}}` instead of a proper genotype template). Skipped by the parser.

- **2 pages have a space before the parenthesis in the title** (`Rs52820871 (G;G)` and `Rs52820871 (G;T)` instead of the standard `Rs52820871(G;G)` format). The annotator handles both title styles.

None of these are scraping errors. They are editorial inconsistencies on the source wiki. The annotator handles all of them correctly: incomplete entries are skipped, variant title formats are matched, and no false annotations are produced.

## Troubleshooting

### `ImportError: cannot import name 'UTC' from 'datetime'`

You're on Python 3.10 or older. Allelix requires 3.11+. `python --version` to check; upgrade and reinstall.

### `allelix db update` failed partway through

Run it again. The update is idempotent — completed downloads are skipped (signal-matched against the upstream freshness probe), so re-running picks up where it left off. If a specific source keeps failing, you can isolate it: `allelix db update --no-gnomad` (or `--no-alphamissense` etc.) to skip the failing source for now, then retry that source standalone later.

If the failure is a network timeout on the same source repeatedly, the issue is likely on the source's end (HuggingFace / NCBI / EBI / etc. transient outage). Wait an hour and retry.

### `db status` shows an annotator as "not ready" but the cache files exist on disk

Almost always a schema or interpreter version bump in a recent allelix release that invalidates the pre-bump cache. The fix is `allelix db update` — it'll re-ingest only the affected annotators (everything else signal-skips). See the relevant `CHANGELOG.md` entry for "Cache invalidation" notes per release.

### Build mismatch warning when analyzing a known-good file

Either: (a) the file's header claims one build and the position data is actually a different build (real-world MHG / Tempus exports do this — header says 37, positions are 38; ADR-0021 documents the auto-detection), or (b) the file is a GRCh36-era export (pre-2012), which allelix detects accurately but doesn't ship a ClinVar cache for (see [`docs/grch36-liftover.md`](https://github.com/allelix/allelix/blob/main/docs/grch36-liftover.md) for the external-liftover path).

If you trust the file's header more than allelix's detection, force the build: `--build grch37` or `--build grch38`. The build banner in the analyze output will reflect what was actually used.

### Is my database stale? How fresh is the data?

`allelix db status` shows the ingestion date of each cache. The freshness model is signal-driven: `allelix db update` queries an upstream freshness signal (md5 / etag) per source and re-ingests only when the signal changed. There's no time-based staleness — a 6-month-old cache against an unchanged source IS current.

If you want to force a refresh: `allelix db update --force` re-downloads every cache regardless of freshness signal.

### Output reports contain real genome annotations

Anything `allelix analyze` writes to disk contains real annotations of your genome. Treat the output files (HTML, JSON, terminal redirects) as sensitive — they're as private as the input genotype file. Allelix does not transmit them; the lifecycle is local.

### Where does allelix store its data?

See [Managing your data](#managing-your-data) above for cache locations, the disposability rule, what to back up, and the uninstall procedure.

## Architecture & Design Decisions

The "why" behind major design choices lives in [`docs/adr/`](https://github.com/allelix/allelix/blob/main/docs/adr/README.md) as Architecture Decision Records. Read these before proposing changes that touch the parser/annotator interfaces, the regulatory posture, or the data-handling model.

Notable load-bearing ADRs:

- **ADR-0016 — Data Classification Principle.** Classification reads structured fields only. Regex on prose is forbidden in production code.
- **ADR-0020 — CPIC API as the per-allele function source.** The ClinPGx non-finding filter is a table join keyed on `(rsid, base) → clinicalfunctionalstatus`, sourced from CPIC's structured API. Supersedes the prose-extraction tiers from earlier versions (ADR-0017, ADR-0018).
- **ADR-0007 — Genotype matching requires the user to carry the ALT allele.** Applies to ClinVar.
- **ADR-0009 — ClinPGx matches the user's exact normalized diploid call.**
- **ADR-0015 — Mock data generators are the contract.** Fixture shape must mirror real data shape; invariants tested.

Release history: see [`CHANGELOG.md`](https://github.com/allelix/allelix/blob/main/CHANGELOG.md).

## Development

```bash
source .venv/bin/activate
pip install -e ".[dev]"

# One-time: install pre-commit hooks
pre-commit install --hook-type pre-commit

ruff check .
ruff format --check .
pytest
```

The pre-commit hook enforces `ruff check` + `ruff format --check`. If a commit is blocked, fix the underlying problem rather than skipping the hook.

For the full real-data release-validation battery — every parser format, VCF/gVCF, PLINK round-trip, edge cases, the GIAB/1000G corpus — see [`test_data/FULL_TEST_PROTOCOL.md`](test_data/FULL_TEST_PROTOCOL.md). That document is the project's real-data testing story; running `pytest` covers the fast suite that gates every commit.

## License

GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later). See `LICENSE`.

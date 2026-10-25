<div align="center">

<img src="https://allelix.io/apple-touch-icon.png" width="180" alt="Allelix" />

**Open-source genotype analysis. Format-agnostic. Offline-first.**

[Website](https://allelix.io) · [Demo](https://analyze.allelix.io) · [Quickstart](#quickstart) · [Changelog](https://github.com/allelix/allelix/blob/main/CHANGELOG.md)

[![python](https://img.shields.io/pypi/pyversions/allelix.svg)](https://www.python.org/downloads/)
[![pypi](https://img.shields.io/pypi/v/allelix.svg)](https://pypi.org/project/allelix/)
[![license](https://img.shields.io/pypi/l/allelix.svg)](https://github.com/allelix/allelix/blob/main/LICENSE)
[![CI](https://github.com/allelix/allelix/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/allelix/allelix/actions/workflows/ci.yml)
[![privacy: local-first](https://img.shields.io/badge/privacy-local--first-success)](#privacy)

</div>

# Allelix

> **Status:** Production. Reads eleven file formats, including VCF and gVCF. Reports what ClinVar, ClinPGx, the GWAS Catalog and SNPedia say about your variants, with frequencies and scores from gnomAD, AlphaMissense and CADD. It decides the genome build from the file itself and stops rather than guess. HTML, JSON, terminal or annotated VCF. Everything runs on your machine. The [Changelog](https://github.com/allelix/allelix/blob/main/CHANGELOG.md) tracks every release.

## Quickstart

**Requires Python 3.11 or newer.** On Python 3.10 or older pip cannot install this release (it may fall back to an old one, which then fails with `ImportError: cannot import name 'UTC' from 'datetime'`); install on 3.11+. Check with `python --version`.

```bash
pip install allelix

# Download reference databases (~19 GB default install; ~25 GB if you opt into CADD).
# Use --no-gnomad / --no-alphamissense to skip the large ones.
# CADD is opt-in: download with `allelix db update --cadd`, then enable it
# with `allelix config set sources.cadd true`.
allelix db update

# Analyze a genotype file
allelix analyze your_genotype_file.txt --output report.html

# VCF / gVCF input — same command, auto-detected
allelix analyze your_wgs.vcf.gz --output report.html

# Multi-sample VCF — pick which sample to analyze
allelix analyze trio.vcf.gz --sample HG002 --output report.html

# VCF quality gate: rows the variant caller filtered (FILTER not PASS) are
# skipped and counted; calls with GQ below 20 or no reads become no-calls.
allelix analyze your_wgs.vcf.gz --keep-filtered --min-gq 10 --output report.html

# Filter to a custom panel (rsIDs and gene names, one per line; a # starts a comment, on its own line or after an entry; blank lines are ignored)
allelix analyze your_genotype_file.txt --filter-file my_panel.txt --output report.html
# Panel coverage puts each requested rsID in one state: annotated, typed (no annotation),
# no call, or not on chip ("not in file" for a VCF)
# VCF coverage includes explicit REF/REF calls and explains rejected/excluded records.
# Gene entries match constituent symbols in lists; hyphens inside symbols stay intact.

# Annotated VCF output — VCF-in / VCF-out for pipeline integration
# Writes a valid VCF with annotations as INFO fields + full provenance in the header.
# Missing ##contig declarations are synthesized; raw CHROM names are preserved.
# Refuses to overwrite; requires VCF or gVCF input.
allelix analyze your_wgs.vcf.gz --vcf-out annotated.vcf
```

See [Development](#development) for source installs and running tests, [Managing your data](#managing-your-data) for cache locations and cleanup, and [Troubleshooting](#troubleshooting) for common failure states.

## Supported Formats

| Format | Status | Notes |
|---|---|---|
| MyHappyGenes (Tempus) | ✓ | Tab-delimited, 5 columns. **Build is auto-detected** — real-world MHG exports mislabel the header as "build 37.1" while shipping GRCh38 coordinates. Allelix decides the build from position data and warns on header/data disagreement (ADR-0021, ADR-0044). See [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#myhappygenes-tempus). |
| 23andMe | ✓ | Tab-delimited, 4 columns, concatenated genotype. Reads build 36, 37 or 38; the header's build is used when the known-SNP positions agree. I-prefixed probe IDs passed through; they and `D`/`I` indel calls are not matched against ClinVar, and the report counts them (ADR-0051). |
| Genes for Good | ✓ | The 23andMe-style `.txt` from a Genes for Good download (23andMe's columns under a PLINK banner). The imputed text file (about 7.7 million rows) and the imputed VCF are labelled imputed in every report; the measured files are not. The download's VCFs are read as VCF; see [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#genes-for-good). |
| AncestryDNA | ✓ | Tab-delimited, 5 columns. Chromosome mapping: 23→X, 24→Y, 25→X (PAR), 26→MT. V1 and V2 chip layouts. |
| Family Tree DNA | ✓ | CSV, double-quoted fields, concatenated genotype. Build 37 default when known-SNP positions are present but do not decide; with none at all and no build line, the run stops (ADR-0044). |
| MyHeritage DNA | ✓ | CSV, same structure as FTDNA. Detected by "MyHeritage" in the file's preamble. Handles double-double-quoted fields and files re-saved through a spreadsheet (padded rows, signature not on line 1, tab-delimited copies). Two exports: the **array** export (pre-2026 kits) and the **low-pass WGS** export (kits from January 2026; genotypes imputed from about 2x coverage), which the report labels as imputed and whose `##reference=build37` line is used as the header build. See [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#myheritage). |
| Living DNA | ✓ | Tab-delimited despite `.csv` extension. Handles AX-, AFFX-prefixed and CHR:POS positional SNP IDs. |
| deCODEme | ✓ | `deCODEme_scan.csv` (deCODE genetics, service closed 2013): comma-separated, a per-row `Strand` column; minus-strand calls are complemented to the plus strand. Build 36: decided from known-SNP positions, then run with position lookups off (ADR-0025). See [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#decodeme). |
| FTDNA Illumina raw | ✓ | Tab-delimited variant of the FTDNA export (distinct from the CSV format above). `RSID/CHROMOSOME/POSITION/RESULT` columns. Build 37 default when known-SNP positions are present but do not decide. |
| FTDNA FamFinder | ✓ | Third FTDNA file shape: tab-delimited with **separate `ALLELE1` / `ALLELE2` columns** instead of a concatenated `RESULT`. Detected by the `famfinder` substring in the file header plus the 5-column canonical header. Build 37 default when known-SNP positions are present but do not decide. |
| VCF / gVCF | ✓ | REF/ALT encoding, `0/1` genotype notation. Plain VCF: absence at a position means reference. gVCF: explicit reference blocks (lines with `<NON_REF>` ALT and `END=` INFO) are not annotated — they match nothing in any annotation database — but they are read at high-value positions: a block with low genotype quality or no reads there makes the position a no-call, and so does a gVCF position no record covers. High-value no-calls list the reason (for example `low_gq_block`, `no_coverage`, `filtered`, `low_gq_call`). A plain VCF has no reference blocks, so a position without a record is reported as unknown rather than as reference or as a no-call. Multi-sample files require `--sample <ID>`. Rows the variant caller rejected (FILTER other than `PASS` or `.`) are skipped and counted in the report; `--keep-filtered` keeps them. A call with genotype quality below `--min-gq` (default 20) or with no reads (DP 0) is reported as a no-call. The spanning-deletion allele `*` is treated as symbolic: such calls are not assessed, and the report counts them and lists their positions. A half call (`./1`) is assessed on its called allele by ClinVar and the GWAS Catalog and shown as `T/.` with zygosity "Half call"; ClinPGx and SNPedia, which need both alleles, skip it. Malformed rows, rows without a genotype, genotypes naming an allele the record lacks and genotypes naming a symbolic allele such as `<NON_REF>` beside a real ALT are counted in the report; a no-call written `./.` stays a no-call (ADR-0047). A VCF that holds only structural-variant or copy-number records (such as a provider's separate `sv` / `cnv` file) is rejected with a message naming the file to use instead; see [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#sequencingcom-30x-whole-genome) for per-provider deliverables. A VCF whose header declares an imputation tool's fields (Minimac, Beagle, GLIMPSE) is labelled imputed in every report; see [docs/formats.md](https://github.com/allelix/allelix/blob/main/docs/formats.md#imputed-vcfs). Streams via stdlib; `.vcf.gz` handled transparently. Optional `pip install allelix[vcf-index]` enables pysam-backed tabix random access for fast `extract --snps` on huge VCFs. |

Files can be passed as delivered: gzipped (Family Tree DNA's `.csv.gz`, for example), or a zip holding one genotype file. A zip with several files, a password or an unusual compression stops with a message to unzip it first. Windows, Linux and old-Mac line endings are read; a leading byte order mark or blank line is ignored.

If Allelix does not recognise a file, use the original file from your testing company. `--format <name>` forces a reader at your own risk; every report then carries a warning (see [docs/formats.md, Forcing a format](https://github.com/allelix/allelix/blob/main/docs/formats.md#forcing-a-format)).

Adding a new format means adding one file to `allelix/parsers/` and registering an instance in the `PARSERS` list in `allelix/parsers/__init__.py`.

### Roadmap

| Release | Theme |
|---|---|
| **v2.2** | Custom-panel UX improvements, methylation panel retune + citations, ClinVar source switch (VCF → `variant_summary.txt`), CLI cleanups (`db clean` / `db path`). |
| **v2.3** | Annotated VCF output plus v2.3.1 correctness hardening: cross-build enrichment guards, GWAS effect-type disambiguation, conservative build confidence, SNPedia orientation handling, and strict slow-suite drift verification. |
| **v2.4** | Correctness hardening: cross-build rsID identity table with explicit allele translation, plausibility flag for genotyping-artifact detection, VCF quality gate honouring FILTER and genotype quality, unified rsID-recovery step across all commands, and parser/format hygiene. |
| **v2.5** | ClinVar's own classification caps each row; ClinVar variants without an rsID; VCF calls counted, not dropped; read-only caches and named sources in every report; zip input; Genes for Good and deCODEme; imputed input labelled. |
| **Later** | Variant consequences, CNVs and repeat expansions, gnomAD genomes, supplemental genotype-file merging, ancestry estimation, and PRS support. These are directional themes rather than release commitments. |

Issues and milestones are tracked at [github.com/allelix/allelix/issues](https://github.com/allelix/allelix/issues).

## Supported Databases

| Database | Status | Notes |
|---|---|---|
| ClinVar (GRCh37 + GRCh38) | ✓ | Public domain (NCBI). SNVs + indels + multi-allelic sites. **Both builds cached**; `analyze` dispatches by the decided build (ADR-0044). Rows are matched by chromosome, position, REF and ALT, never by rsID alone, and reference-allele records (REF = ALT) are never reported (ADR-0039). Carrier rule (ADR-0007) requires the user to carry the ALT allele. Indel-anchor protection (ADR-0011) prevents single-base array readouts from matching anchor-base indels. Array calls ClinVar cannot match are counted in the report: vendor probes without an rsID at a locus ClinVar classifies as Pathogenic or Likely pathogenic, and `D`/`I` calls (ADR-0051). ClinVar variants without an rsID are matched for VCF input only, where REF comes from the file (ADR-0046). A VCF allele written with bases it shares with REF (`GTT>GT,G`, `CA>TA`) is trimmed per carried allele before matching, and a VCF homozygous deletion to one base matches; indels are not left-aligned, so one the file places at a different position from ClinVar's still does not match (ADR-0047). |
| ClinPGx (formerly PharmGKB) | ✓ | CC BY-SA 4.0. Clinical annotations only — single-rsid SNVs; star alleles and haplotypes deferred (ADR-0009). **Primary non-finding filter is the ClinVar REF carrier rule (ADR-0023):** if ClinVar publishes a single-base REF for the rsid and the user is homozygous for it, the row is suppressed. CPIC's `(rsid, base) → function_class` join (ADR-0020) survives as a secondary tier for rsids ClinVar doesn't catalog. Earlier prose tiers (ADR-0013, ADR-0017, ADR-0018) are superseded. |
| CPIC (per-allele function table) | ✓ | Internal data source for the ClinPGx filter. Fetched from `api.cpicpgx.org` at `db update` time. Used to populate the `pharmgkb_allele_function` table — not surfaced to end users as its own annotator. |
| SNPedia | ✓ | CC BY-NC-SA 3.0 US. Pre-built cache downloaded via `db update` (~216K wiki pages, ~105K genotype rows). The cache carries SNPedia's orientation metadata: minus-strand alleles are complemented before matching. Single-base genotypes at rsIDs are then checked against the site's forward alleles in the rsID identity table: a page's stated orientation must put the genotype among them, an unknown orientation is read directly or else complemented, and A/T and C/G sites abstain. Sites that abstain are counted in the report. If the database is absent, analysis runs without it. For commercial use, pass `--exclude-snpedia`. The cache can also be rebuilt from source via `scripts/scrape_snpedia.py` + `scripts/parse_snpedia.py`. |
| GWAS Catalog | ✓ | Public domain (EBI/NHGRI). Trait–SNP associations with p-values and effect sizes. Carrier rule (ADR-0007) requires the user to carry the risk allele. The reported allele is the study's effect (risk) allele, which can be the reference allele (ADR-0039). P-value magnitude scoring maps continuous p-values to the 0–10 scale; only a confidence-interval field that positively identifies a bare odds ratio receives the odds-ratio modifier. Beta coefficients and ambiguous or unit-bearing effects receive no modifier (ADR-0038). Unknown-risk-allele entries fire on rsID match alone but are capped at 3.0. Risk alleles are read on the forward strand against the site's alleles in the rsID identity table: as written, else complemented; A/T and C/G sites and rsIDs the table does not place at the call abstain, and the report counts them. Every row says whether the effect allele raises or lowers the odds or the measured value (JSON `effect_direction`). GWAS rows must meet the larger of `--min-magnitude` and `--gwas-min-magnitude` (default 9.0); the must-include rsIDs skip the GWAS floor (ADR-0024). |
| gnomAD | ✓ | ODbL v1.0. **GRCh38-only enrichment annotator** — adds population allele frequency context to existing annotations. Shows how common each variant is in the general population (~57M exome variants with an rsID, from 730K individuals). A pathogenic variant that 35% of people carry reads very differently from one seen in 0.001%. Pre-built cache downloaded via `db update` (~7 GB on disk). GRCh37 input is enriched only through the rsID identity table's explicit allele translation (ADR-0040); gnomAD's own rsID resolution stays GRCh38-only, and GRCh36 input receives nothing. Use `--no-gnomad` to skip. |
| AlphaMissense | ✓ | CC BY 4.0. **GRCh38-only enrichment annotator** — adds DeepMind's protein-structure-based pathogenicity predictions to existing annotations. Scores 71M missense variants on a 0–1 scale: <0.34 = likely benign, >0.564 = likely pathogenic. Complements ClinVar's expert classifications with computational predictions — especially valuable for variants ClinVar hasn't reviewed yet. Pre-built cache downloaded via `db update` (~8GB on disk). GRCh37 input is enriched only through the rsID identity table (ADR-0040); GRCh36 input receives nothing. Use `--no-alphamissense` to skip. |
| CADD | ✓ | LicenseRef-CADD (non-commercial). **GRCh38-only enrichment annotator** — adds PHRED-scaled deleteriousness scores from CADD v1.7. **Opt-in**, disabled by default; GRCh37 input is scored only through the rsID identity table (ADR-0040) and GRCh36 input never; see [CADD modes](#cadd-modes) below for details. |
| dbSNP rsID identity | ✓ | Public domain (NCBI dbSNP). Not an annotator: a cross-build table of each annotator rsID's GRCh37 and GRCh38 placement and alleles (ADR-0040). Supplies the reference allele for GRCh37 and GRCh38 array data, recovers rsIDs for rsID-less VCFs on both builds, and translates a GRCh37 allele to its GRCh38 coordinates so gnomAD, AlphaMissense and CADD can be looked up. Strand-inverted loci are translated explicitly; palindromic sites (A/T, C/G, or a site whose dbSNP alleles include all four bases) and alleles that became the GRCh38 reference abstain. Pre-built cache downloaded via `db update` (~4.9M rows covering ~3.0M rsIDs, from dbSNP b157). |

### CADD modes

CADD ranks how deleterious any single-nucleotide variant is using 100+ annotation tracks (coding, non-coding, regulatory). PHRED 10 = top 10% most deleterious, 20 = top 1%, 30 = top 0.1%.

**Opt-in**: download the cache with `allelix db update --cadd`, then enable it with `allelix config set sources.cadd true`. Downloading the cache alone does not change persistent configuration. Use `--no-cadd` to skip enrichment for a single run when it is enabled in config.

Two modes:

- **Cache mode** (default): pre-built ~5.8 GB SQLite cache (~120M GRCh38 variant keys). Covers nearly every GRCh38 position allelix can annotate from its other databases (gnomAD, AlphaMissense, ClinVar). GRCh37 input is scored only after translation through the rsID identity table; GRCh36 input never is.
- **Full mode** (`options.cadd_full = true`): tabix queries against the complete CADD v1.7 file via pysam (GRCh38 only). Adds coverage for novel or private variants that appear only in WGS/WES data and are not in any pre-computed database. Requires `pip install allelix[cadd]` for the pysam dependency.

If your input is a genotyping chip file, cache mode is normally all you need. GRCh37 chip files are scored through the rsID identity table; GRCh36 files are not scored.

### Build coverage asymmetry (GRCh37 vs GRCh38)

ClinVar dispatches per build and ships with both GRCh37 and GRCh38 caches. gnomAD v4.1 exomes, AlphaMissense, and CADD v1.7 are GRCh38 datasets. Every enrichment lookup is made at GRCh38 coordinates by chromosome, position, reference and alternate allele, never by rsID. For GRCh37 input, the cross-build rsID identity table (ADR-0040) translates the allele the user carries at its GRCh37 position into the matching GRCh38 allele, written explicitly on the GRCh38 forward strand where the reference strand was inverted between assemblies. Palindromic sites (A/T, C/G, or all four bases in dbSNP, whose strand cannot be told from sequence) and alleles that are the GRCh38 reference receive no values rather than a guess. Without the identity table installed, GRCh37 analyses receive no gnomAD, AlphaMissense or CADD values; GRCh36 analyses never do.

Position-based rsID recovery for rsID-less VCFs uses ClinVar (on the selected build), the identity table (on GRCh37 and GRCh38), and gnomAD (GRCh38 only). The same recovery step runs for `analyze`, `stats`, `extract` and `export plink`, so high-value warnings, extracted genotypes and PLINK allele coding all see the recovered rsIDs (ADR-0043). The identity table covers the rsIDs the annotators hold, so an rsID-less GRCh37 VCF now recovers ClinPGx, GWAS Catalog and SNPedia rsIDs as well as ClinVar ones. If you convert between assemblies, re-normalize REF/ALT as part of the conversion and pass `--build` explicitly when the file does not carry a trustworthy build signal.

### Known ClinPGx limitation: reference-genotype rows where ClinVar and CPIC both lack data

ADR-0022 + ADR-0023: a tiny residual of ClinPGx rows may appear in reports even when the user is homozygous reference. ClinPGx publishes one annotation per genotype including the reference homozygote, and for the reference-homozygote row to be suppressed Allelix needs structured data on the variant from either:

- **ClinVar's REF allele** (the primary filter — see ADR-0023). Covers any rsID ClinVar catalogs.
- **CPIC's per-allele function table** (the secondary fallback — see ADR-0020). Covers rsIDs CPIC has classified.

For the rare rsID where ClinPGx has an annotation but *neither* ClinVar nor CPIC has data, the row emits. These are identifiable by a homozygous-reference genotype combined with "decreased risk," "may have a typical response," or similar comparative language. They are an upstream data gap, not an Allelix bug — we surface them honestly rather than hide them behind a curated exclusion list (which would recreate the maintenance trap the v0.5–v0.7 prose filters were trying to escape).

The CFTR × ivacaftor leak (~30+ rows on real data, pre-v0.7.3) is fixed by the ADR-0023 ClinVar REF check: CPIC's CFTR vocabulary (`"ivacaftor responsive"`) doesn't match the four-class enum the secondary tier expects, but ClinVar publishes REF for every CFTR rsID, so the primary tier catches them universally.

### Known ClinVar upstream data quality issues

One ClinVar row in real-world reports is a known upstream artifact, not an Allelix bug:

- **PKD1 rs199476100 GG (Pathogenic/Likely pathogenic, magnitude 8.5).** This is a stop-gained variant with a gnomAD frequency of 0.0005% (7 observations in 1.38 million chromosomes). Homozygosity for this variant is biologically implausible — PKD1 is autosomal dominant and the nonsense variant would be embryonic-lethal or devastating in homozygous state. The chip genotyping call is almost certainly a probe artifact. The code correctly reports what ClinVar says and what the chip reads; the error is upstream of Allelix. The row is now marked implausible (see [Plausibility flag](#plausibility-flag)) but is still shown.

It does not affect Allelix's filter logic. It comes from the limits of array-based genotyping chips.

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

Allelix stores persistent configuration in `config.toml` (in the data directory, default `~/.local/share/allelix/`). A default config is created on first run, except by `analyze --no-update`, which writes nothing to the data directory (ADR-0048).

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

# Tune the plausibility flag (allele frequencies between 0 and 1)
allelix config set plausibility.homozygous_max_af 0.0001
allelix config set plausibility.heterozygous_max_af 0
```

CLI flags (`--no-gnomad`, `--no-alphamissense`, `--no-cadd`, `--exclude-snpedia`) turn a source off for a single run. The config sets the baseline; flags override per-invocation.

### Plausibility flag

Genotyping arrays sometimes miscall very rare variants, most visibly as a homozygous call for an allele almost nobody carries. Allelix compares each call with the allele's gnomAD frequency and marks it for review. The row is never hidden (ADR-0041).

| Call | gnomAD frequency of the allele | Flag |
|---|---|---|
| Homozygous | below `plausibility.homozygous_max_af` (default 0.0001, i.e. 0.01%) | **implausible**: often a genotyping artifact |
| Heterozygous | at or below `plausibility.heterozygous_max_af` (default 0: gnomAD has the site but no carriers) | **suspicious** |
| Either | otherwise | plausible |

The HTML report shows a ⚠ icon next to a flagged genotype, with the reason on hover and in the detail panel, and a muted ○ with the reason when the call was not assessed. The terminal table marks implausible calls with ⚠, suspicious ones with ? and not-assessed ones with ○. JSON reports carry `plausibility` on every annotation (`plausible`, `suspicious`, `implausible` or `not_assessed`) and, when not assessed, `plausibility_reason` (`no_frequency`, `gnomad_not_run`, `chromosome`, `half_call` or `no_allele`; ADR-0050).

A call is not assessed when gnomAD has no frequency for the allele. Allelix's gnomAD cache covers exome sites with an rsID, so a missing value says nothing about rarity. Calls are also not assessed when the gnomAD database is skipped. Homozygous calls on X are not assessed because a male's single X reads as homozygous; heterozygous X calls are. A call the file itself declares as one copy (a one-letter array call such as 23andMe's male X calls, or a one-index VCF genotype) is shown as Hemizygous and assessed as one copy. Y and mitochondrial calls are not assessed. Confirm a flagged call with a clinical-grade test before relying on it.

When three or more variants on chromosomes 1 to 22 are called homozygous for an allele ClinVar classifies as Pathogenic or Likely pathogenic, and are rare in gnomAD or have no gnomAD frequency, the report warns at the top. This is counted before any display filter, so `--min-magnitude` cannot hide it (ADR-0050).

### Database sizes and download times

Not all databases are equal in size. `allelix db update` downloads them all by default, but you can skip the large ones if disk space or bandwidth is a concern:

| Database | On disk | Download time | What it adds |
|---|---|---|---|
| ClinVar (GRCh37 + GRCh38) | ~2.4 GB | 1–2 min, then up to about 20 min to build both caches | Core clinical variant classifications. Required. |
| ClinPGx + CPIC | ~6MB | seconds | Drug-gene interactions. |
| GWAS Catalog | ~950 MB | 1–2 min | Trait-SNP associations from genome-wide studies. |
| gnomAD | ~7 GB | 5–15 min | Population allele frequencies (how common is this variant?). |
| AlphaMissense | ~8GB | 5–15 min | Missense pathogenicity predictions (how likely to break protein function?). |
| dbSNP rsID identity | ~460 MB | ~1 min | Cross-build placements: reference allele for array data, rsID recovery on both builds, GRCh37 enrichment translation. |
| CADD (opt-in) | ~6 GB | 5–15 min | Variant deleteriousness scores (how damaging is this variant?). Download with `db update --cadd`; enable in config. |

gnomAD and AlphaMissense are the largest but add the most interpretive context. gnomAD answers "is this variant rare or common?" — a pathogenic variant carried by 35% of the population reads very differently from one seen in 3 people. AlphaMissense answers "does this missense change likely damage the protein?" — especially valuable for the thousands of variants ClinVar hasn't reviewed yet.

To skip either during download: `allelix db update --no-gnomad --no-alphamissense`. To disable permanently: `allelix config set sources.gnomad false`.

### Managing your data

Allelix keeps its caches in `~/.local/share/allelix/`, or `$XDG_DATA_HOME/allelix` when `XDG_DATA_HOME` is set. To use another folder, set `ALLELIX_DATA_DIR` or pass `--data-dir`. Default install (no CADD) is ~19 GB on disk; with the CADD opt-in it's ~25 GB.

```
~/.local/share/allelix/
├── alphamissense.sqlite       ~7.8 GB
├── gnomad.sqlite              ~7.1 GB
├── cadd.sqlite                ~5.8 GB   (only if you opted into --cadd)
├── clinvar.GRCh37.sqlite      ~1.2 GB
├── clinvar.GRCh38.sqlite      ~1.2 GB
├── gwas_catalog_associations.tsv  ~745 MB
├── rsid_identity.sqlite       ~457 MB
├── gwas.sqlite                ~202 MB
├── snpedia.sqlite             ~194 MB
└── pharmgkb.sqlite            ~6 MB
```

**The caches are disposable and re-downloadable.** Everything here was fetched by `allelix db update` and can be re-fetched at any time. There's no per-user state in the cache — your genotype files never enter this directory. Safe to delete to reclaim space (or to force a full refresh if you suspect a cache went stale):

```bash
# Reclaim all ~19-25 GB, with size report and confirmation prompt
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

**The one thing worth backing up: `config.toml`.** It lives inside the resolved data directory (typically `~/.local/share/allelix/config.toml`) and stores your license assertions (commercial-mode toggle, CADD opt-in confirmation) and per-source enable/disable settings. It's the only non-reproducible state in the install. The caches will rebuild on the next `db update`; `config.toml` won't. `allelix db clean` always preserves it.

To fully uninstall allelix and reclaim everything:

```bash
ALLELIX_DATA=$(allelix db path)
allelix db clean --yes                 # reference database caches
rm -f "$ALLELIX_DATA/config.toml"      # license + source-toggle config
pip uninstall allelix
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
| gnomAD | gnomad.broadinstitute.org | ODbL v1.0 | Attribution required. GRCh38 population allele frequencies for context (GRCh37 input through the identity table); not a clinical annotator. Use `--no-gnomad` to omit. |
| AlphaMissense | zenodo.org/records/10813168 | CC BY 4.0 | Attribution required. Cheng et al., Science 2023. GRCh38 missense variant pathogenicity predictions (GRCh37 input through the identity table). Use `--no-alphamissense` to omit. |
| dbSNP (rsID identity table) | ncbi.nlm.nih.gov/snp | Public domain | No restrictions. Cross-build rsID placements derived from dbSNP; not surfaced as an annotator. |
| CADD | cadd.gs.washington.edu | LicenseRef-CADD | Attribution required, **non-commercial by default**. Commercial licenses available from UW CoMotion. GRCh38 scores (GRCh37 input through the identity table); download with `allelix db update --cadd`, then enable in config (`allelix config set sources.cadd true`). Use `--no-cadd` to omit it for one run. |

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

The v2.3.1 cache format also preserves each page's SNPedia orientation. Known minus-strand rows are complemented to forward orientation before genotype matching. For palindromic A/T or C/G genotypes whose orientation is unknown, Allelix abstains because the source data cannot distinguish the strands safely. Re-run `allelix db update` after upgrading so a pre-v2.3.1 SNPedia cache is rebuilt.

## Troubleshooting

### `ImportError: cannot import name 'UTC' from 'datetime'`

You're on Python 3.10 or older, where pip cannot install this release and may fall back to an old one. Allelix requires 3.11+. `python --version` to check; upgrade and reinstall.

### `allelix db update` failed partway through

Run it again. The update is idempotent — completed downloads are skipped (signal-matched against the upstream freshness probe), so re-running picks up where it left off. If a specific optional source keeps failing, use its `--no-*` flag to finish the other sources, then remove the skip flag and re-run the full update later. `--only clinvar` exists for the ClinVar-only CI and recovery path; other sources do not currently have standalone `--only` selectors.

If the failure is a network timeout on the same source repeatedly, the issue is likely on the source's end (HuggingFace / NCBI / EBI / etc. transient outage). Wait an hour and retry.

### `db status` shows an annotator as "not ready" but the cache files exist on disk

`analyze` prints one line per skipped source with the reason and the fix (ADR-0048):

- **out of date:** almost always a schema or interpreter version bump in a recent allelix release that invalidates the pre-bump cache. The fix is `allelix db update`; it re-ingests only the affected annotators (everything else signal-skips). See the relevant `CHANGELOG.md` entry for "Cache invalidation" notes per release.
- **cannot read:** the file's permissions do not let you read it.
- **damaged:** the file is not an SQLite database, lacks a table, or its main table has no rows; `allelix db update` replaces that source only (`--force` would download every source again).
- **WAL mode from an older Allelix:** gnomAD, AlphaMissense and CADD caches installed by earlier releases are in WAL mode. Reading one leaves `-wal` and `-shm` files beside it, and in a directory you cannot write SQLite cannot read it at all, so `analyze --no-update` skips it. Run `allelix db update` once (with write access to the data directory); it converts them without downloading. `analyze` without `--no-update` converts them too.

The report records the same reasons: the JSON `sources` block lists every source with its part in the run, and the HTML lists the sources that could not be used. A report without ClinVar carries a warning at the top (ADR-0049).

With `--no-update`, `analyze` writes nothing to the data directory, so a data directory can be read-only.

### How the genome build is decided, and why a run can stop

Before anything is annotated, Allelix decides the file's genome build from the evidence the file itself carries (ADR-0044):

- **VCF contig lengths.** Every nuclear chromosome has a different length in GRCh37 and GRCh38, so the `##contig` lengths settle the build whatever the contig names look like (`chr1` or `1`). If chrM has the old hg19 length (16,571 bp), no mitochondrial record is annotated and a note says so, also when the build comes from `--build` or the file's assembly header. In a GRCh37 VCF, mitochondrial records named `chrM` or `M` (UCSC hg19 names) are left out, counted and noted, unless the file declares that name with the 16,569 bp reference length; `MT` and `chrMT` are read.
- **Known-SNP positions across the whole file.** A set of well-known SNPs sits at different positions in each build. Allelix reads the whole file for them before annotating (a few seconds on a full chip). At least three concordant positions and an 80% majority make the evidence confident (ADR-0037). A VCF record without an rsID counts when it sits at such a position.
- **The file's own label** (a build line in an array header, or a VCF `assembly=`) is used only when the evidence does not contradict it and at least one known-SNP position is there to check it.

The banner says where the build came from, for example `Build: GRCh38 (from the contig lengths: 24 contig length(s) match GRCh38)`. If confident positions contradict the label (real-world MHG / Tempus exports say 37 and ship 38), Allelix uses the positions and prints a **Build mismatch** notice. An unlabelled chip whose known-SNP positions do not decide uses its format's documented default, with a visible note, provided at least one of those positions matches a build and the build most of them match is the default (a few that disagree do not stop it); the note gives the counts. A tentative GRCh36 winner is the safety exception: Allelix uses GRCh36 so it cannot query a modern position cache (see [`docs/grch36-liftover.md`](https://github.com/allelix/allelix/blob/main/docs/grch36-liftover.md)).

**The run stops** instead of guessing when:

- a VCF's contig lengths contradict its `assembly=`, contradict each other, or match no supported assembly;
- weaker (tentative) position evidence contradicts the file's label, or points an unlabelled chip away from its default;
- the known-SNP positions are tied between builds and the build that would otherwise decide (the label, the format default, or the contig-length build) is not one of them. When it is one of them, Allelix proceeds on it and prints an "inconclusive" warning, also shown in the HTML report. Without contig lengths, a tie that includes GRCh36 instead takes the GRCh36 safety exception described above, with the same warning;
- a VCF has no contig lengths, no assembly and too few known-SNP positions;
- a chip states no build and carries none of the known-SNP positions at all (for example an X-chromosome-only export, since they are all on autosomes);
- a file states a build but carries none of the known-SNP positions to check it (for example an mtDNA-only export, or a small panel VCF without contig lengths);
- the file has known-SNP positions but none of them sits at any supported build's coordinate, whether or not it states a build (contig lengths, when a VCF declares them, still decide);
- contig lengths decided, but the known-SNP positions met while annotating point to the other build (for example a lifted-over file that kept its old header). The message says the annotation pass ran; nothing is written;
- `--build grch37` or `--build grch38` was given, but the file's contig lengths, or three or more known-SNP positions with an 80% majority, say the other build. The message names both builds. Contrary positions that fall short of that (fewer than three, or under an 80% majority) only print a **Build warning**, and so does a tie that leaves out the forced build. `--build grch36` is never checked, because it turns position-keyed lookups off.

The message names the evidence and asks you to say which build the file uses: `--build grch37` or `--build grch38`, or `--build grch36` for an old GRCh36 / hg18 file, which then runs with position-keyed lookups off (see [`docs/grch36-liftover.md`](https://github.com/allelix/allelix/blob/main/docs/grch36-liftover.md)). For `compare`, use `--build1` / `--build2`. A stop leaves no report, annotated VCF, diff, PLINK files or extract table behind. A forced build is recorded as such in every report.

### Keeping input and output files safe

Report `--output`, annotated `--vcf-out` and all three actual PLINK output paths
must differ from the input, including symlink and hard-link aliases. These
checks run before cache setup or output creation. Inputs that alias the former
fixed report or PLINK `.tmp` path are also refused; unrelated old temp files
are left untouched. Annotated VCF still refuses any existing output at the CLI.

Writers use unique temporary files in the destination directory and publish
only complete, closed content. Failure before publication preserves previous
outputs. PLINK stages all three files first, then publishes each atomically;
a later publication failure can leave complete new files next to unchanged old
ones. Existing ordinary permission bits are retained when replacement is
permitted; new staging files start private. Ownership, ACLs, extended attributes,
crash durability and concurrent path swaps are outside this guarantee.

### Is my database stale? How fresh is the data?

`allelix db status` shows the ingestion date of each cache. The freshness model is signal-driven: `allelix db update` queries an upstream freshness signal (md5 / etag) per source and re-ingests only when the signal changed. There's no time-based staleness — a 6-month-old cache against an unchanged source IS current.

If you want to force a refresh: `allelix db update --force` re-downloads every cache regardless of freshness signal. ClinVar is never refreshed unverified: its freshness signal is NCBI's checksum for the submission summary, which also verifies the download. If that checksum cannot be fetched (after one retry), ClinVar is not refreshed, with or without `--force`; the existing caches are kept, and a later `db update` tries again. The GWAS Catalog and ClinPGx publish no checksum, so `--force` refreshes them even when their freshness signal cannot be fetched; it stores no signal, and the next ordinary `db update` refreshes them once more when the signal is available again.

A ClinPGx cache built while CPIC's allele-function data could not be fetched is marked as such: `db status` shows "(no CPIC)", and the next `db update` tries CPIC again. A cache from 2.5.0 that claims CPIC data but holds none is repaired the same way by an ordinary `db update`.

### Output reports contain real genome annotations

Anything `allelix analyze` writes to disk contains real annotations of your genome. Treat the output files (HTML, JSON, terminal redirects) as sensitive — they're as private as the input genotype file. Allelix does not transmit them; the lifecycle is local.

### Where does allelix store its data?

See [Managing your data](#managing-your-data) above for cache locations, the disposability rule, what to back up, and the uninstall procedure.

## Architecture & Design Decisions

The "why" behind major design choices lives in [`docs/adr/`](https://github.com/allelix/allelix/blob/main/docs/adr/README.md) as Architecture Decision Records. Read these before proposing changes that touch the parser/annotator interfaces, the regulatory posture, or the data-handling model.

Notable load-bearing ADRs:

- **ADR-0007 — Genotype matching requires the user to carry the ALT allele.** Applies to ClinVar.
- **ADR-0009 — ClinPGx matches the user's exact normalized diploid call.**
- **ADR-0015 — Mock data generators are the contract.** Fixture shape must mirror real data shape; invariants tested.
- **ADR-0016 — Data Classification Principle.** Classification reads structured fields only. Regex on prose is forbidden in production code.
- **ADR-0020 — CPIC API as the per-allele function source.** The ClinPGx non-finding filter is a table join keyed on `(rsid, base) → clinicalfunctionalstatus`, sourced from CPIC's structured API. Supersedes the prose-extraction tiers from earlier versions (ADR-0017, ADR-0018).
- **ADR-0036 — Annotated VCF output contract.** JSON and annotated VCF schemas evolve independently.
- **ADR-0037 — Conservative position-based build confidence.** A winning build needs at least three matches and an 80% ratio; weaker signals remain diagnostic only.
- **ADR-0038 — GWAS effect-modifier disambiguation.** Only positively identified odds ratios receive an effect-size modifier; beta and ambiguous effects abstain.
- **ADR-0041 — Plausibility flag.** Very rare homozygous calls and heterozygous calls with no gnomAD carriers are flagged for review, never suppressed.
- **ADR-0042 — VCF quality gate.** The variant caller's FILTER and genotype quality are honoured, low-quality reference blocks are no-calls, and every exclusion is counted in the report.
- **ADR-0043 — One rsID-recovery step.** Every command that needs rsIDs recovers them for rsID-less input through the same shared step; build detection votes by position instead.
- **ADR-0044 — Evidence before labels.** The build is decided before annotation from contig lengths or whole-file known-SNP positions; evidence that contradicts the label, or no usable evidence, stops the run instead of producing a report on the wrong build.
- **ADR-0045: ClinVar's own classification.** A ClinVar row is capped by ClinVar's aggregate classification of the variant (Conflicting, Benign, Likely benign), with the submitter split shown; outranked P/LP submissions are kept by ClinVar's tiered rule where the aggregate is not a disease call.
- **ADR-0046: ClinVar variants without an rsID.** Kept in the cache and matched for VCF input only, by locus and the file's REF; reports name such a finding by `chrom:pos ref>alt`.
- **ADR-0047: VCF calls no longer dropped silently.** Untrimmed alleles are trimmed per carried allele, half calls are assessed by the sources that need only a carried allele, VCF homozygous deletions to one base match ClinVar, and dropped calls are counted with their reason (a genotype on a gVCF reference block's own `<NON_REF>` is not checked).
- **ADR-0048: Caches are read-only at analyze time.** With `--no-update`, `analyze` writes nothing to the data directory; every cache opens read-only, installs in rollback-journal mode, and a skipped source is named with its reason (missing, outdated, unreadable, damaged, or WAL mode from an older Allelix). Local migrations run in the update step.
- **ADR-0049: Every report names its sources.** JSON `sources` gives every registered source's part in the run (used, ready but not used, or skipped with a reason); `report_warnings` carries report-level warnings, starting with ClinVar not used; the annotated VCF names skipped sources.
- **ADR-0050: Zygosity codes and haploid calls.** JSON zygosity and plausibility are codes, a call the file declares haploid is shown as hemizygous, "not assessed" carries a reason, and three rare homozygous P/LP calls give a report warning.
- **ADR-0051: Array calls ClinVar cannot match are counted.** Vendor probes at P/LP loci and `D`/`I` calls are counted in the report, and panel coverage puts each requested rsID in one of four states.
- **ADR-0052: Zero-star P/LP on a common allele.** A ClinVar P/LP row without assertion criteria whose allele has a gnomAD frequency of 5% or more is scored down, marked and counted.

Release history: see [`CHANGELOG.md`](https://github.com/allelix/allelix/blob/main/CHANGELOG.md).

## Development

```bash
source .venv/bin/activate
pip install -e ".[dev,vcf-index]"   # the CI install: vcf-index (pysam) runs the indexed-VCF tests

# One-time: install pre-commit hooks
pre-commit install --hook-type pre-commit

ruff check .
ruff format --check .
pytest
```

The pre-commit hook enforces `ruff check` + `ruff format --check`. If a commit is blocked, fix the underlying problem rather than skipping the hook.

For the full real-data release-validation battery — every parser format, VCF/gVCF, PLINK round-trip, edge cases, the GIAB/1000G corpus — see [`test_data/FULL_TEST_PROTOCOL.md`](https://github.com/allelix/allelix/blob/main/test_data/FULL_TEST_PROTOCOL.md). That document is the project's real-data testing story; running `pytest` covers the fast suite that gates every commit.

## License

GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later). See `LICENSE`.

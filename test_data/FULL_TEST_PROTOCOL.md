# Full Test Protocol

External reviewer checklist for verifying an allelix release against real data.

**Requirements:** Fast machine, fast internet, ~50 GB free disk space.

## 1. Environment setup

```bash
git clone https://github.com/allelix/allelix.git
cd allelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,vcf-index]"
```

Verify Python 3.11+:

```bash
python --version
```

After any version change (pull, checkout, version bump) reinstall with
`pip install -e ".[dev,vcf-index]"`. The `vcf-index` extra installs pysam,
which the indexed-VCF tests need; without it they skip, and §2 allows no
skips. `allelix --version` must report the version in
this checkout's `pyproject.toml`.

## 2. Unit test suite (synthetic + auto-fetched real data)

**Run this section after §3**, in the same shell: the ClinVar drift
test reads the caches §3 builds. Run it as a non-root user, in the git
checkout from §1, with network access and with `plink2` and `bcftools`
installed. The real GWAS Catalog fixture (~65 MB) auto-fetches from EBI
on first run; subsequent runs are offline-fast against the local cache.

```bash
ALLELIX_CI_REQUIRE_CACHE=1 python -m pytest tests/ --tb=short -rs
```

No `-x`: one failure must not hide the tests after it. `-rs` prints the
reason for every skip.

**Expected:** every collected test passes, **0 skipped**, and the
collected count is at least the latest row of the floor table below.
Skipped tests are not a pass: a skip hides whether the code works.

Skips that can occur elsewhere, and why they cannot occur in this run:
- `bcftools` or `plink2` not installed (`tests/test_cli_vcf_out.py`,
  `tests/exporters/test_plink_export.py`): both are installed here.
- Running as root (`tests/test_cli.py`, the read-only directory check;
  root bypasses `chmod`): this run is as a non-root user.
- A real-data download failed (`tests/databases/test_clinvar_tsv_loader.py`):
  this run has network access; if a download fails, rerun.
- No ClinVar cache (`tests/databases/test_clinvar_clnsig_drift.py`): §3
  built it, and `ALLELIX_CI_REQUIRE_CACHE=1` turns this skip into a
  failure.

The fast tier for development is `pytest -m "not slow and not
integration"`; it is not the release run.

Test-count floor by release:
- v1.9.0: ~1,400
- v2.0.0: ~1,486 (VCF/gVCF parser, FTDNA Illumina raw, CLNSIG
  drift CI, rsID resolution, ClinPGx rename)
- v2.0.1: ~1,525 (correctness fixes,
  Variant case normalization, ClinVar CLNDN-join, terminal bare-min)
- v2.0.2: ~1,540 (auto-fetch GWAS fixture, chr-prefix build
  detection, pyproject version fallback, enrichment annotator
  stack-management, doc/process cleanup)
- v2.1.0: ~1,576 (ADR-0035: Variant.ref + per-Annotation
  alt threading + structured GWAS fields + strand-aware carrier
  matching; pipeline-level Variant.ref population for array data)
- v2.2.0: ~1,672 (per-SCV ClinVar TSV loader + dash-CLNSIG
  guard; streaming ClinVar ref-lookup + permission resolver
  dedup; FTDNA FamFinder parser; panel-coverage states +
  JSON schema_version 5→6; runtime nudges; branch
  coverage; methylation panel + magnitude floor;
  db clean / db path subcommands)
- v2.2.1: ~1,718 (CHANGELOG sanitize gate; per-SCV
  duplicate-row collapse; ClinVar curatorial-term
  filter expansion; VCF parser robustness on
  non-UTF-8 / binary input; §19 invariant-based protocol
  with HG002 ground-truth ride-along)
- v2.2.2: ~1,733 (rsID-less VCF coverage via gnomAD
  second-pass resolver; BOM + preamble tolerance;
  gnomAD position covering index + per-batch row-value IN
  collapse; stats Build em-dash + chr-prefix hint)
- v2.2.3: ~1,734 (23andMe v5 "is generated" banner-variant
  tuple; gnomAD position resolver chunk size to //2
  convention; README data-lifecycle docs; release-content
  discipline gate -- sanitize CI guard extended to flag
  drift-vulnerable counts in [Unreleased] without snapshot
  caveat)
- v2.3.1: 1,831 collected (v2.3.0 annotated-VCF coverage plus
  cross-build enrichment guards, GWAS effect-type disambiguation,
  conservative build confidence, SNPedia orientation, strict ClinVar
  drift execution, and complete VCF contig declarations)
- v2.4.0: 2,100 collected
- v2.4.1: 2,371 collected
- v2.5.0: at least 3,896 collected
- v2.5.1: at least 5,255 collected
- v2.5.2: at least 6,946 collected

Check lint:

```bash
ruff check . && ruff format --check .
```

**Expected:** All checks passed, 0 files reformatted.

## 3. Download all databases

Use an empty release-candidate data directory so this gate proves fresh
cache construction instead of reusing a developer cache. Keep this export
in the same shell for the rest of the protocol:

```bash
export ALLELIX_DATA_DIR="$(mktemp -d)"
echo "RC data directory: $ALLELIX_DATA_DIR"
allelix db update
```

This downloads ClinVar (GRCh37 + GRCh38), ClinPGx (formerly PharmGKB —
the upstream rebranded in 2026; `pharmgkb.org` 301-redirects to
`clinpgx.org`, but the cache file is still named `pharmgkb.sqlite`
internally for backward compatibility), GWAS Catalog, gnomAD
(~2.7 GB compressed download, ~7 GB on-disk SQLite cache),
AlphaMissense (~1.8 GB compressed, ~8 GB on-disk), and SNPedia from
HuggingFace. CADD is opt-in and not included here — step 7 installs
and enables it before the wrong-allele safety gate.
Total on-disk footprint after this step: roughly **18 GB**
(CADD adds another ~5.8 GB at step 7, bringing the
full install to ~24 GB).

**Expected:** All enabled annotators show green checkmarks. No errors.

Verify status:

```bash
allelix db status
```

**Expected:** Every source this step installed shows "yes" in the Ready
column with version strings and record counts. CADD is listed too and
shows "no" until step 7 installs it. SNPedia should show ~104K records.

The real `db update` above is the fresh-cache release gate for both GWAS
Catalog and SNPedia. Also run their isolated loader/parser regression
suites; every test receives a new `tmp_path` database:

```bash
python -m pytest tests/databases/test_gwas_loader.py \
  tests/databases/test_snpedia_parser.py --no-cov
```

**Expected:** All tests pass. The SNPedia cache version includes parser
format 7 and an `orientation` column; an older parser tag must be rejected
as stale. The GWAS loader must build from the current association archive
without reusing a pre-existing SQLite cache.

Now run the full test suite in §2, in this shell.

## 4. Fetch real test data

```bash
bash scripts/fetch_testdata.sh
```

This downloads the ~1.27 GB `test_data.tar.gz` release asset from the
v2.0.0 GitHub release. The bundle contains the full real-data
fixture set used by steps 5–6, 15, and 19:

- `real/23andme/` — 6 openSNP users
- `real/ancestrydna/` — 7 openSNP users (V1.0 + V2.0 arrays)
- `real/ftdna/` — 7 openSNP users (CSV + Illumina raw + gzip/zip variants)
- `real/livingdna/`, `real/myheritage/`, `real/mhg/` — user1190 transcoded into each format
- `real/vcf/` — GIAB HG002 GRCh37/38 benchmarks, HG002 chr22 slice, HG00187 GATK-HC gVCF, 1000G chr22 multi-sample, plus synthetic mocks
- `transcoded/` — user1190 also represented as AncestryDNA + FTDNA CSV (transcoded from the 23andMe source)

Full per-file index with sizes and provenance is in the "Tarball
contents (authoritative index)" subsection below.

As of v2.0.2, the GWAS Catalog zip is **auto-fetched by the test
suite itself** on first use, so this script's GWAS-fetch
is optional — left in place for users who prefer one-shot setup.

**Expected:** `test_data/real/` and `test_data/transcoded/` populated.
`test_data/gwas_catalog.zip` present (either from the script or
auto-fetched on first slow-test run).

### Tarball contents (authoritative index)

All DTC genotyping files sourced from [openSNP](https://opensnp.org/)
(CC0 public domain) unless noted otherwise. VCF / gVCF files from
NIST GIAB, 1000 Genomes, or synthetic.

#### `real/23andme/`

| File | Lines | Date | Notes |
|------|------:|------|-------|
| `user1_v1.txt` | 966,998 | 2011-05-03 | openSNP user 1. Early v1 chip (~967K SNPs). |
| `user10_file3_yearofbirth_1982_sex_XY.23andme.txt` | 966,998 | 2011-09-20 | openSNP user 10. v1 chip. |
| `user500_file2637_yearofbirth_1956_sex_XY.23andme.txt` | 960,629 | 2014-07-30 | openSNP user 500. v4/v5 era chip (~961K SNPs). |
| `user1190_v5.txt` | 960,628 | 2013-08-04 | openSNP user 1190. **Canonical cross-format test subject** — same individual transcoded to all other formats. |
| `user1500_file819_yearofbirth_1985_sex_XY.23andme.txt` | 960,628 | 2013-11-22 | openSNP user 1500. v5 chip. |
| `user3000_file1922_yearofbirth_unknown_sex_unknown.23andme.txt` | 574,533 | 2014-09-18 | openSNP user 3000. Older / smaller chip (~575K SNPs). |

#### `real/ancestrydna/`

| File | Array | Notes |
|------|-------|-------|
| `user1001.txt` | V1.0 | openSNP user 1001. |
| `user2393_file1486_yearofbirth_unknown_sex_unknown.ancestry.txt` | V1.0 | openSNP user 2393. |
| `user3672_file2417_yearofbirth_unknown_sex_unknown.ancestry.txt` | V1.0 | openSNP user 3672. |
| `user4440_file3043_yearofbirth_1954_sex_XX.ancestry.txt` | V1.0 | openSNP user 4440. Female. |
| `user4941_file3489_yearofbirth_1997_sex_XY.ancestry.txt` | V1.0 | openSNP user 4941. |
| `user5351_file3892_yearofbirth_1979_sex_XY.ancestry.txt` | V2.0 | openSNP user 5351. V2 array. |
| `user5715_file4190_yearofbirth_1981_sex_unknown.ancestry.txt` | V2.0 | openSNP user 5715. V2 array. |

#### `real/ftdna/`

| File | Type | Size | Notes |
|------|------|-----:|-------|
| `user288.csv` | CSV text | 23 MB | openSNP user 288. Standard FTDNA CSV format. |
| `user339_file150_yearofbirth_1950_sex_XY.ftdna-illumina.txt` | CSV text | 592 KB | openSNP user 339. An X-chromosome-only export (18,091 rows). It states no build and carries none of the autosomal known-SNP positions, so it stops without `--build`. Its X positions match neither GRCh37 nor GRCh38 (136 rsIDs checked against the identity table), and 25 of 25 sampled rsIDs sit at their GRCh36 (hg18) position in UCSC's snp130 track, so it is GRCh36: pass `--build grch36`. To check another such file, compare a few of its rsIDs' positions with dbSNP's on each build. |
| `user2503_file1534_yearofbirth_unknown_sex_unknown.ftdna-illumina.txt` | CSV text | 2.4 MB | openSNP user 2503. |
| `user3395_file2210_yearofbirth_1942_sex_XY.ftdna-illumina.txt` | Zip archive | 265 KB | openSNP user 3395. Compressed — contains VCF inside zip. |
| `user4706_file3310_yearofbirth_unknown_sex_unknown.ftdna-illumina.txt` | gzip compressed | 6.1 MB | openSNP user 4706. Gzipped raw data. |
| `user5404_file3917_yearofbirth_unknown_sex_unknown.ftdna-illumina.txt` | gzip compressed | 6.2 MB | openSNP user 5404. Gzipped raw data. |
| `user6056_file4561_yearofbirth_unknown_sex_unknown.ftdna-illumina.txt` | gzip compressed | 6.3 MB | openSNP user 6056. Gzipped raw data. |

#### `real/livingdna/`

| File | Size | Notes |
|------|-----:|-------|
| `user1190.csv` | 23 MB | **Transcoded** from user1190 23andMe v5. LivingDNA CSV format, GRCh37. |

#### `real/myheritage/`

| File | Size | Notes |
|------|-----:|-------|
| `user1190.csv` | 31 MB | **Transcoded** from user1190 23andMe v5. MyHeritage CSV format. |

#### `real/mhg/`

| File | Size | Notes |
|------|-----:|-------|
| `user1190.txt` | 24 MB | **Transcoded** from user1190 23andMe v5. MHG format. |

#### `real/vcf/`

**Synthetic fixtures (old copies; no step uses them)**

These are copies of the repository's `tests/fixtures/` files from before
the fixtures were made consistent. Their contig lengths contradict each
other (chr1 at its GRCh37 length; chr22, or chr12 in
`mock_vcf_rsidless.vcf`, at its GRCh38 length), so each of the four stops without `--build` ("its contig lengths contradict each
other"). The protocol uses the committed `tests/fixtures/` versions
instead. This index changes when the tarball is refreshed.

| File | Size | Notes |
|------|-----:|-------|
| `mock_vcf.vcf` | 1.2 KB | Plain single-sample VCF with rsIDs. Stops without `--build` (contig lengths contradict). |
| `mock_gvcf.g.vcf` | 823 B | gVCF with reference blocks (`<NON_REF>` ALT, `END=` INFO). Stops without `--build` (contig lengths contradict). |
| `mock_multisample.vcf` | 404 B | 3-sample VCF (`SAMPLE_A`, `SAMPLE_B`, `SAMPLE_C`). Stops without `--build` (contig lengths contradict). |
| `mock_vcf_rsidless.vcf` | 616 B | rsID-less VCF (`ID=.`). Stops without `--build` (contig lengths contradict). |

**GIAB HG002 benchmarks** — Source: [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) v4.2.1, Ashkenazi son (HG002 / NA24385).

| File | Size | Build | Notes |
|------|-----:|-------|-------|
| `HG002_GRCh38_benchmark.vcf.gz` | 150 MB | GRCh38 | Full WGS, ~4.05M variants. Ships without rsIDs. |
| `HG002_GRCh38_benchmark.vcf.gz.tbi` | 1.6 MB | | Tabix index. |
| `HG002_GRCh37_benchmark.vcf.gz` | 169 MB | GRCh37 | Full WGS, ~4.03M variants. Carries rsIDs in the ID column (~95% of rows). |
| `HG002_GRCh37_benchmark.vcf.gz.tbi` | 1.6 MB | | Tabix index. |
| `HG002_GRCh38_chr22.vcf.gz` | 1.9 MB | GRCh38 | chr22 subset, ~50K variants. Fast smoke tests. |
| `HG002_GRCh38_chr22.vcf.gz.tbi` | 21 KB | | Tabix index. |

**GATK HaplotypeCaller gVCF** — Source: Broad Institute GATK test-data bucket.

| File | Size | Build | Notes |
|------|-----:|-------|-------|
| `HG00187_gatkhc.g.vcf.gz` | 280 MB | GRCh37 | 1000 Genomes Finnish sample (HG00187). ~19.3M lines, ~99.8% reference blocks. Exercises gVCF parser at real scale. |
| `HG00187_gatkhc.g.vcf.gz.tbi` | 3.1 MB | | Tabix index. |

**1000 Genomes multi-sample** — Source: [EBI 1000 Genomes 30x](https://www.internationalgenome.org/).

| File | Size | Notes |
|------|-----:|-------|
| `thousandG_chr22.vcf.gz` | 496 MB | GRCh38. chr22, 3,202 samples (phased, duohmm). Exercises `--sample <ID>` selection. Its only contig line has no length and no assembly, so a bound run needs `--build grch38` (rs4680 sits at its GRCh38 position, 22:19963748). |
| `thousandG_chr22.vcf.gz.tbi` | 36 KB | Tabix index. |

`real/vcf/README.md` — in-tarball documentation for the VCF test set (v2.0.0 additions).

#### `transcoded/`

All transcoded from `real/23andme/user1190_v5.txt` using the project's transcoder. See `transcoded/README.md` (in-tarball) for methodology.

| File | Format | Size | Notes |
|------|--------|-----:|-------|
| `user1190_as_ancestrydna.txt` | AncestryDNA | 24 MB | Structural transcode. Chromosomes remapped to AncestryDNA convention (X→23, Y→24, MT→26). |
| `user1190_as_ftdna.csv` | FTDNA CSV | 31 MB | Structural transcode. Concatenated genotype in `RESULT` column. |

## 5. Analyze real genotype files

Run analysis on each format against the live-downloaded databases.
Allelix does not create the output directory, so create it first:

```bash
mkdir -p /tmp/allelix-review
```

### 5a. 23andMe

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/user1190_23andme.json
```

**Expected:** Exit code 0. JSON report written. This GRCh37 consumer
array contains ClinVar, ClinPGx, GWAS, and SNPedia findings. gnomAD and
AlphaMissense are GRCh38 datasets; on this GRCh37 input they are reached
through the rsID identity table, which translates each allele to its
GRCh38 coordinates, so they contribute values and provenance wherever
the table can translate (rows the table flags palindromic whose reference
differs between the builds, and alleles that became the GRCh38 reference,
abstain). CADD is not enabled until step 7. Check that
the annotation count is in the hundreds (varies by database version).

### 5b. MHG / Tempus

```bash
allelix analyze test_data/real/mhg/user1190.txt --output /tmp/allelix-review/user1190_mhg.json
```

**Expected:** Exit code 0. JSON report written. This file is a clean
GRCh37 transcode of user1190_v5.txt — no build mismatch expected.
(The mismatch fixture is `edge_cases/mhg_grch38_with_grch37_header.txt`,
tested in step 15.)

### 5c. AncestryDNA

```bash
# One named file, so every machine tests the same input (V1.0 array)
allelix analyze test_data/real/ancestrydna/user1001.txt \
  --output /tmp/allelix-review/ancestrydna.json
```

**Expected:** Exit code 0. JSON report written.

### 5d. FTDNA

```bash
# Pin user288.csv — standard FTDNA CSV format. Don't `find | head` here:
# real/ftdna/ also contains six `*.ftdna-illumina.txt` files (different
# parser), some of which are zip / gzip archives despite the .txt
# extension — picking one of those would either route to the wrong
# parser (5g exercises the Illumina-raw variant) or fail to parse.
allelix analyze test_data/real/ftdna/user288.csv \
  --output /tmp/allelix-review/ftdna.json
```

### 5e. Living DNA

```bash
allelix analyze test_data/real/livingdna/user1190.csv --output /tmp/allelix-review/user1190_livingdna.json
```

### 5f. MyHeritage

```bash
allelix analyze test_data/real/myheritage/user1190.csv --output /tmp/allelix-review/user1190_myheritage.json
```

### 5g. FTDNA Illumina raw (tab-delimited)

```bash
allelix analyze tests/fixtures/mock_ftdna_illumina.txt --output /tmp/allelix-review/ftdna_illumina.json
```

**Expected:** Exit code 0. JSON report written. This is the second
FTDNA file shape (tab-delimited, `RSID/CHROMOSOME/POSITION/RESULT`
header), distinct from the CSV variant tested in 5d.

### 5h. FTDNA FamFinder (separate allele columns)

```bash
allelix analyze tests/fixtures/mock_ftdna_famfinder.txt \
  --output /tmp/allelix-review/ftdna_famfinder.json
```

**Expected:** Exit code 0. JSON report written. This is the third
FTDNA file shape: tab-delimited with separate `ALLELE1` / `ALLELE2`
columns and a FamFinder header. It must select the FamFinder parser,
not the Illumina or CSV parser.

### 5i. VCF / gVCF (committed fixtures + bundled real-scale files)

```bash
# Plain single-sample synthetic VCF
allelix analyze tests/fixtures/mock_vcf.vcf --output /tmp/allelix-review/mock_vcf.json

# Synthetic gVCF (reference blocks present, must be skipped at parse time)
allelix analyze tests/fixtures/mock_gvcf.g.vcf --output /tmp/allelix-review/mock_gvcf.json

# Real-scale gVCF — HG00187 GATK-HC, ~19.3M lines, ~99.8% reference
# blocks. Exercises the gVCF parser end-to-end at WGS scale, not just
# the 823-byte synthetic. GRCh37 build (bare contigs). Two runs.

# Run 3a: DEFAULTS (no --min-gq flag). Exercises the GQ floor on a
# low-coverage (~4x) genome: all ~44K non-ref calls fail GQ >= 20,
# producing 0 annotations and 28 high-value no-calls.
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --output /tmp/allelix-review/hg00187_gatk_default.json \
  --report-format json

# Run 3b: --min-gq 0. Disables the GQ floor and exercises the
# annotation path on position-based resolution. The §19 ground-truth
# harness (`hg00187_gatkhc_gvcf`) checks this file's §19 run against
# its floors.
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --min-gq 0 \
  --output /tmp/allelix-review/hg00187_gatk.json \
  --report-format json
# Expected (3b): exit 0; non-zero total annotations (low per-sample because
# GATK-HC raw output is rsID-less — hits come via position-based
# ClinVar + rsID-less resolution, not the rsID fast-path). Exact
# count drifts with each ClinVar / GWAS / ClinPGx refresh; use the
# §19 ground-truth harness for the regression-detection signal.
# Build banner reads "GRCh37 (forced with --build)". Without --build
# the file's 24 GRCh37 contig lengths decide the same build (ADR-0044).

# Multi-sample VCF (3 samples) — must fail without --sample
allelix analyze tests/fixtures/mock_multisample.vcf 2>&1 | head -3
# Expected: MultiSampleError listing SAMPLE_A, SAMPLE_B, SAMPLE_C.
# (No tail truncation here — the file only has 3 samples.)

# Multi-sample VCF (3 samples) — succeeds with --sample
allelix analyze tests/fixtures/mock_multisample.vcf --sample SAMPLE_A \
  --output /tmp/allelix-review/mock_multisample_A.json

# Truncation-tail demo — 3,202-sample 1000 Genomes chr22 file
allelix analyze test_data/real/vcf/thousandG_chr22.vcf.gz 2>&1 | head -3
# Expected: MultiSampleError listing the first 10 samples followed by
# "... and 3192 more" (or similar count). This is the only file in
# the bundle large enough to exercise the truncation path. The sample
# check runs before the build decision, so this error comes first.
# The file is GRCh38 (rs4680 sits at its GRCh38 position 22:19963748)
# and declares no contig lengths or assembly, so the bound run needs
# --build grch38; without it the run stops and asks for the build.
allelix analyze test_data/real/vcf/thousandG_chr22.vcf.gz \
  --sample HG00096 --build grch38 \
  --output /tmp/allelix-review/thousandG_HG00096.json
```

**Expected:** Single-sample, gVCF, and the real-scale HG00187 run
all exit 0. Multi-sample without `--sample` raises MultiSampleError
in both cases; the 1000G file demonstrates the "... and N more"
truncation tail.

### 5j. rsID-less VCF resolution (the flagship feature)

Real VCFs from variant callers (GATK HaplotypeCaller, DeepVariant)
emit `ID=.` — no rsID. Without resolution, every rsID-keyed annotator
returns zero hits. v2.0.0 resolves rsIDs by `(chrom, pos)` through the
ClinVar cache.

```bash
allelix analyze tests/fixtures/mock_vcf_rsidless.vcf \
  --output /tmp/allelix-review/rsidless.json --build grch37
python3 -c "
import json
d = json.load(open('/tmp/allelix-review/rsidless.json'))
print(f'annotations: {len(d[\"annotations\"])}')
# Pre-fix this returned 0. Post-fix MTHFR rs1801133 (chr1:11856378 G→A)
# resolves through the test ClinVar cache; the annotation appears.
assert any(a['rsid'] == 'rs1801133' for a in d['annotations']), 'rsID resolution failed'
print('rs1801133 resolved from chr1:11856378 ID=.')
"
```

**Expected:** Non-zero annotation count. rs1801133 present in
annotations. Pre-v2.0.0 this would have produced 0 annotations.

### 5k. Build decision: evidence, stops and confidence (ADR-0044)

The build is decided before annotation from the file's own evidence.
Run every part.

**(a) No evidence: the run stops.** A VCF with no contig lengths, no
`assembly=` and no known-SNP positions.

```bash
allelix analyze tests/fixtures/mock_vcf_no_signal.vcf; echo "exit=$?"
```

Expected: exit 1, no report, and a message that starts
`Cannot use mock_vcf_no_signal.vcf: the file carries no contig lengths,
no assembly and not enough known-SNP positions.`, lists the evidence and
ends `... (command line: --build grch37 or --build grch38; for an old
GRCh36 / hg18 file, --build grch36, see docs/grch36-liftover.md), and run
it again.`

**(b) Contig lengths decide, whatever the contig names.**

```bash
allelix analyze tests/fixtures/mock_vcf_chr_prefix_grch38.vcf 2>&1 | grep -E "^Build:"
allelix analyze tests/fixtures/mock_vcf_rsidless.vcf 2>&1 | grep -E "^Build:"
```

Expected: `Build: GRCh38 (from the contig lengths: 5 contig length(s)
match GRCh38)` for the chr-named file, and `Build: GRCh37 (from the
contig lengths: 5 contig length(s) match GRCh37)` for the bare-named one.
There is no "inferred from chr-prefixed contig names" message any more.

**(c) Contradictions stop.** The deliberately inconsistent fixture
declares chr1 with the GRCh37 length and chr22 with the GRCh38 length.

```bash
allelix analyze tests/fixtures/mock_vcf_mixed_contigs.vcf; echo "exit=$?"
allelix analyze tests/fixtures/mock_vcf_mixed_contigs.vcf --build grch38 \
  --output /tmp/allelix-review/mixed_forced.json; echo "exit=$?"
```

Expected: the first run exits 1 with `its contig lengths contradict each
other`; the second exits 0 and its banner reads `Build: GRCh38 (forced with
--build)` (a `; the file says …` suffix appears only when the forced build
differs from the file's own label). This run proves only that `--build`
overrides the stop; the fixture is deliberately inconsistent, so its
annotations are not checked.

**(d) Full-size wrong-label and no-label inputs.** Real files are
sorted by chromosome, so most known SNPs sit far into the file; up to
2.4.0 only the first 100,000 rows were looked at. Run the generated
full-size regression tests (120,000 filler rows, seeds beyond them,
a mock ClinVar record on the wrong build at a coordinate the input
carries):

```bash
python -m pytest tests/test_build_decision_full_size.py -v --no-cov
```

Expected: 13 passed. Then the same on a real full-size chip, relabelled
and unlabelled:

```bash
sed "s/assembly build 37/assembly build 38/" test_data/real/23andme/user1190_v5.txt \
  > /tmp/allelix-review/user1190_label38.txt
grep -v 'assembly build 37' test_data/real/23andme/user1190_v5.txt \
  > /tmp/allelix-review/user1190_nolabel.txt
allelix analyze /tmp/allelix-review/user1190_label38.txt \
  --output /tmp/allelix-review/user1190_label38.json 2>&1 | grep -E "^Build|mismatch"
allelix analyze /tmp/allelix-review/user1190_nolabel.txt \
  --output /tmp/allelix-review/user1190_nolabel.json 2>&1 | grep -E "^Build"
allelix analyze test_data/real/23andme/user1190_v5.txt \
  --output /tmp/allelix-review/user1190_reference.json
python3 -c "
import json
def rows(f):
    d = json.load(open(f))
    return d['input']['build'], d['input']['build_source'], sorted(
        (a['rsid'], a['source'], a.get('significance')) for a in d['annotations'])
ref = rows('/tmp/allelix-review/user1190_reference.json')
for name in ('label38', 'nolabel'):
    got = rows(f'/tmp/allelix-review/user1190_{name}.json')
    assert got[:2] == ('GRCh37', 'seeds'), got[:2]
    assert got[2] == ref[2], name + ': annotations differ from the correctly labelled file'
    print(name, 'OK', got[:2], len(got[2]), 'annotations')
"
```

Expected: the relabelled file prints `Build: GRCh37 (from known-SNP
positions: 11/11 match GRCh37)` and a `Build mismatch: the file says GRCh38, but
its known-SNP positions say GRCh37` notice; the unlabelled file prints
the same build with no mismatch; both produce exactly the reference
file's annotations.

**(e) Confidence rule and decision table.** A modern-build winner needs
at least three matching seeds and at least 80% (ADR-0037); weaker
evidence that contradicts a label stops the run; a tentative GRCh36
winner stays on GRCh36 (ADR-0025).

```bash
python -m pytest tests/utils/test_build_decision.py \
  tests/utils/test_build_detect.py::TestDetectBuild \
  tests/reports/test_pipeline.py::TestBuildDecisionInPipeline -v --no-cov
```

**Expected:** all pass, including the seed table check against the
dbSNP snapshot in `tests/fixtures/dbsnp_seed_placements.json`.

### 5l. Formats new in 2.5

Genes for Good, deCODEme and zip downloads, on the committed synthetic
fixtures.

```bash
F=tests/fixtures
allelix analyze $F/mock_genesforgood_imputed.txt --output /tmp/allelix-review/genesforgood_imputed.json
allelix analyze $F/mock_decodeme.csv --output /tmp/allelix-review/decodeme.json
```

**Expected:**

- `genesforgood_imputed.json`: `genesforgood`, with `input.method`
  `array-imputed`; the measured fixture (`mock_genesforgood.txt`) has
  no method.
- `decodeme.json`: `decodeme`, GRCh36 from `seeds`.
- **Zip:** a zip holding `mock_23andme.txt` (and a `__MACOSX/` entry)
  gives the same annotations as the plain file. A zip holding two
  genotype files stops with "is a zip file holding 2 files".

### 5m. Builds that stop in 2.5

```bash
F=tests/fixtures
sed '1a\
##contig=<ID=1,assembly=GRCh37>' $F/mock_vcf_no_signal.vcf > /tmp/allelix-review/label37_only.vcf
allelix analyze /tmp/allelix-review/label37_only.vcf     # must stop
allelix analyze $F/mock_vcf_chr_prefix_grch38.vcf --build grch37   # must stop
```

**Expected:**

- **`label37_only.vcf` stops:** "the file says GRCh37 but carries no
  known-SNP position to check it". A label alone no longer decides the
  build.
- **The `--build grch37` run stops:** "--build grch37 was given but the
  file's contig lengths are GRCh38".
- No stop writes a report.
- An array file labelled GRCh36 (Build 36) with no known-SNP position is
  exempt from the first stop and runs on GRCh36, as `--build grch36`
  does (ADR-0025); `tests/parsers/test_twentythreeandme.py::TestGrch36Label`
  covers it.

### 5n. Report warnings, sources and strand abstentions

- **`clinvar_skipped`:** make a data directory holding links to every
  cache except `clinvar.*`. Run `allelix analyze
  tests/fixtures/mock_23andme.txt --data-dir <that directory>
  --no-update`. The JSON `report_warnings` has `clinvar_skipped`,
  `sources` shows clinvar `skipped` with reason `missing`, and no
  annotation comes from ClinVar.
- **`input_trailing_bytes`:** gzip `mock_23andme.txt` and append bytes
  that are not gzip data (a saved web page). The file is read, and
  `report_warnings` has `input_trailing_bytes`.
- **Every `report_warnings` entry** is `{code, message, count}`.
- **The `sources` block, in every report:** each registered source
  appears with `state` `used`, `ready_not_used` or `skipped`; a skipped
  source has a reason; every annotator that ran is `used`. With
  `--exclude-snpedia` (§6) SNPedia is `skipped:excluded`; in commercial
  mode (§10) it is `skipped:license`.
- **Strand abstentions:** the §7 GRCh38 report has a
  `strand_abstentions` block with integer counts for `gwas` and
  `snpedia`.
- **Feature tests:**

  ```bash
  python -m pytest tests/test_sources.py tests/test_haploid_and_homozygous_warning.py \
    tests/parsers/test_gzip_damage.py tests/test_gwas_snpedia_strand.py \
    tests/parsers/test_zip.py tests/parsers/test_unidentified_files_rejected.py \
    -q -rs --no-cov
  ```

  Expected: all pass, none skipped.

## 6. Cross-parser identity check

The user1190 genotype exists in 6 format representations. All must
produce identical annotation sets: the same `(source, rsid,
significance)` rows, each the same number of times. The exact number
depends on database versions; the sets must match across formats.

Each run writes its own file, named by format, so no run overwrites
another and no other section's output is compared:

```bash
mkdir -p /tmp/allelix-review
for pair in \
  23andme:test_data/real/23andme/user1190_v5.txt \
  mhg:test_data/real/mhg/user1190.txt \
  livingdna:test_data/real/livingdna/user1190.csv \
  myheritage:test_data/real/myheritage/user1190.csv \
  ancestrydna:test_data/transcoded/user1190_as_ancestrydna.txt \
  ftdna:test_data/transcoded/user1190_as_ftdna.csv; do
  fmt=${pair%%:*}; f=${pair#*:}
  echo "=== $fmt: $f ==="
  allelix analyze "$f" --exclude-snpedia \
    --output "/tmp/allelix-review/cross_${fmt}.json" || echo "FAILED: $f"
done
```

**Expected:** every run exits 0 (no `FAILED:` line).

Then compare the six annotation sets:

```bash
python3 -c "
import json, glob
from collections import Counter
sets = {f: Counter((a['source'], a['rsid'], a['significance'])
                   for a in json.load(open(f))['annotations'])
        for f in sorted(glob.glob('/tmp/allelix-review/cross_*.json'))}
assert len(sets) == 6, sorted(sets)
ref_name, ref = next(iter(sets.items()))
for name, rows in sets.items():
    missing, extra = ref - rows, rows - ref
    assert not missing and not extra, (name, 'vs', ref_name, dict(missing), dict(extra))
print(f'ALL 6 IDENTICAL: {sum(ref.values())} annotations each')
"
```

**Expected:** `ALL 6 IDENTICAL: N annotations each`. Any difference
is a parser or build-decision bug; the assertion names the file and
its missing and extra rows.

## 7. Wrong-allele safety — enrichment uses exact alt match

Verify that enrichment lookups use exact alt-allele matching, not
MAX-aggregated fallback and not complement-resolved coincidence. The
v2.0.1 ship closed three independent paths where a wrong-allele
number could attach to an annotation:

- `resolve_strand` no longer falls back to complement
  matching at multi-allelic sites. CADD scores attached to alt-set
  annotations must come from a direct `(ref, alt)` match.
- Alt-less annotations (raw GWAS rows) no longer take the
  `MAX(af) GROUP BY rsid` fallback. They get enrichment only via the
  safe position-fallback path (rsID resolved on-the-fly via ClinVar,
  resolved tuple carries the user's specific alt). CADD enrichment
  now has the same position-fallback for symmetry.
- ClinVar `CLNDN` list joined per record instead of
  index-paired with `CLNSIG`, eliminating Frankenstein pairings at
  multi-SCV variants.

```bash
# Download the optional cache and enable it persistently for the remainder
# of this isolated RC run. `db update --cadd` alone does not mutate config.
allelix db update --cadd
allelix config set sources.cadd true
allelix db status | grep -i cadd
# Expected: CADD ready, sources.cadd = true.

allelix analyze test_data/real/vcf/HG002_GRCh38_benchmark.vcf.gz \
  --build grch38 \
  --output /tmp/allelix-review/enrichment_check.json
allelix analyze test_data/real/23andme/user1190_v5.txt \
  --output /tmp/allelix-review/grch37_enrichment_gate.json

# Both reports: the enrichment sources took part and attached values.
python - <<'PY'
import json

for path, build, sources in (
    ('/tmp/allelix-review/enrichment_check.json', 'GRCh38',
     {'clinvar', 'pharmgkb', 'gwas', 'snpedia', 'gnomad', 'alphamissense', 'cadd'}),
    ('/tmp/allelix-review/grch37_enrichment_gate.json', 'GRCh37',
     {'gnomad', 'alphamissense', 'cadd', 'rsid_identity'}),
):
    report = json.load(open(path))
    assert report['input']['build'] == build, (path, report['input']['build'])
    annotators = {a['name'] for a in report['annotators']}
    assert sources <= annotators, (path, 'missing', sorted(sources - annotators))
    counts = {
        field: sum(a.get(field) is not None for a in report['annotations'])
        for field in ('allele_frequency', 'am_pathogenicity', 'cadd_phred')
    }
    assert all(counts.values()), (path, 'vacuous enrichment gate', counts)
    print(f'{build} enrichment counts: {counts}')
PY

# Every attached value belongs to the allele the person carries.
python scripts/check_wrong_allele.py /tmp/allelix-review/enrichment_check.json \
  --input test_data/real/vcf/HG002_GRCh38_benchmark.vcf.gz
python scripts/check_wrong_allele.py /tmp/allelix-review/grch37_enrichment_gate.json \
  --input test_data/real/23andme/user1190_v5.txt
```

**Expected:** All seven annotators participate in the GRCh38 report.
The real GRCh37 consumer report reaches the same three sources through
the rsID identity table: its provenance lists gnomAD, AlphaMissense,
CADD and the identity table. In both reports the gnomAD, AlphaMissense
and CADD counts are each non-zero, and `scripts/check_wrong_allele.py`
prints `N/N allele-direct, 0 non-direct, 0 mismatch, 0 unresolved` and
exits 0.

**What the harness checks.** For every enriched annotation, the
attached gnomAD AF, AlphaMissense score and CADD score must be the
database values for the allele the person carries, at that locus, with
no strand complement. It uses the product's enrichment key, never the
rsID: the GRCh38 `(chrom, pos, ref, allele)`, where the allele is the
annotation's alt or, for an alt-less row (a GWAS row with no reported
risk allele), the single non-reference allele the user carries. GRCh37
keys are translated through the identity table by locus and alleles
(at strand-inverted sites the GRCh38 allele is the complement, an
explicit translation; a row the table flags palindromic translates only
where its GRCh37 and GRCh38 reference bases are equal). An allele
equal to the build's reference base is checked the same way only where
that reference carries a ClinVar variant allele (a single-base REF == ALT
row in a current ClinVar cache, as Factor V Leiden on GRCh37); anywhere
else an enriched row whose allele is the reference base is a mismatch.
A value that
belongs to another allele at the locus is "non-direct"; a value that
matches nothing is a "mismatch"; an annotation no locus can be found
for is "unresolved". The loci come from the identity table, ClinVar and
gnomAD; `--input` (required) gives the analysed file, whose ID column
places rsIDs no cache can (for example retired rsIDs in the GRCh37
benchmark). Only the caches of sources in the report's provenance are
required, and a value from a source missing from the provenance is a
mismatch. Any non-direct, mismatch or unresolved row, or an empty
denominator, fails.

**ClinVar real-cache content gates.** After a `db update` against
the v2.2 per-SCV TSV loader. **These are ship-gates** in the same
sense §14's GWAS gates are: zero expected results mean
ship; any non-zero result means **do not tag**. The flagship
ClinVar loader gets the same gating posture as every other annotator.

**(a) Per-SCV row shape (supersedes the previous join pin).**
The v2.2 loader stores one row per (variant, SCV submission), not one
row per variant with semicolons. Expect multiple rows per multi-SCV
variant, each carrying a single SCV's condition exactly as the submitter
recorded it. How many SCVs one variant has follows NCBI's data and
changes between ClinVar releases (rs1063192 had 2 rows in one cache and
1 in a cache built ten days later), so the check counts alleles with
several rows instead of relying on any single variant:

```bash
sqlite3 "$ALLELIX_DATA_DIR/clinvar.GRCh38.sqlite" \
  "SELECT COUNT(*) FROM clinvar_variants WHERE rsid = 'rs1800896';"
# Expected: ≥1 row.
sqlite3 "$ALLELIX_DATA_DIR/clinvar.GRCh38.sqlite" \
  "SELECT COUNT(*) FROM (SELECT allele_id FROM clinvar_variants \
   WHERE allele_id IS NOT NULL GROUP BY allele_id HAVING COUNT(*) > 1);"
# Expected: at least 100,000 ClinVar alleles with more than one row
# (898,674 of 2,891,556 alleles in a GRCh38 cache built in October 2026).
# A per-variant cache (one combined row with a semicolon-joined
# condition, the pre-2.2 shape) has none.
```

The previous Frankenstein-pairing safety question is preserved
structurally: each per-SCV row carries one submitter's
(significance, condition) pair, so there is no cross-product
pairing risk anymore. The semicolon-join in the old pin was the
loader's coping mechanism for the VCF's CLNSIG|CLNDN parallel
arrays; with per-SCV rows that mechanism is unnecessary.

**(b) Significance-sentinel guard (ship gate).** ClinVar's
`submission_summary.txt.gz` carries several placeholder values that
mean "no classification recorded." The hazard isn't every
placeholder equally — it's specifically the placeholders that ISN'T
in `allelix.annotators.clinvar._CLNSIG_MAGNITUDE`. Anything in that
dict scores below the 5.0 analyze display floor and is safely
suppressed; anything missing falls through to the 5.0 default and
surfaces as a meaningless annotation on real rsIDs (rs137854557 →
Meningioma was the canary). The loader filters exactly that set at
ingest; the gate scans for the same set.

The authoritative set is `_CLINVAR_PLACEHOLDER_CLNSIGS` in
`allelix/databases/manager.py`. **Do not broaden the SQL below
without broadening the loader's set in the same PR** — they must
agree, or the gate will false-block on values the loader already
handles correctly (e.g. `"not provided"` maps to 2.0, harmless).

```bash
# Sentinel scan: zero rows expected. Any row → do not tag.
# Set MUST match _CLINVAR_PLACEHOLDER_CLNSIGS in databases/manager.py.
# LOWER() because the loader's skip is case-insensitive (matches
# _magnitude()'s _normalize_clnsig normalization); a case-sensitive
# scan would show false-green if ClinVar emits a placeholder with
# different casing than the canonical form.
sqlite3 "$ALLELIX_DATA_DIR/clinvar.GRCh38.sqlite" \
  "SELECT clinical_significance, COUNT(*) FROM clinvar_variants \
   WHERE LOWER(clinical_significance) IN \
     ('-', '', 'not specified', 'no classification provided', \
      'other', 'association', 'association not found') \
   GROUP BY clinical_significance;"
# Expected: 0 rows. v2.2.1 added the three non-classification
# curatorial terms — "other", "association", "association not found".
# Keep this list in sync with _CLINVAR_PLACEHOLDER_CLNSIGS in
# databases/manager.py.

# Spot-check a known-affected rsID:
sqlite3 "$ALLELIX_DATA_DIR/clinvar.GRCh38.sqlite" \
  "SELECT rsid, clinical_significance FROM clinvar_variants WHERE rsid = 'rs137854557';"
# Expected: real classifications (e.g. Pathogenic / Likely pathogenic).
# Never a placeholder.
```

Repeat the sentinel scan against `clinvar.GRCh37.sqlite`. Both
builds share the same TSV source pair, so both must be clean.

**Note on denylist vs allowlist.** This gate is a *denylist* — it
catches the placeholders we've observed. Any future ClinVar term
not in `_CLNSIG_MAGNITUDE` (a new real classification, a freshly-
introduced placeholder, an alternative casing) falls to the 5.0
default, and both the loader's skip and this denylist miss it
(false-green). The durable form is an *allowlist* drift-guard
asserting every distinct `clinical_significance` is in
`_CLNSIG_MAGNITUDE ∪ _CLINVAR_PLACEHOLDER_CLNSIGS`, else "do not
tag." That's the CLNSIG vocabulary drift test
(`tests/databases/test_clinvar_clnsig_drift.py` +
`allelix/data/clinvar_clnsig_snapshot.yaml`). When that test evolves
from snapshot-equality to "every value ∈ {known magnitude vocab,
placeholder set}", point it at the same caches this gate scans;
this §7b gate then graduates from "denylist confirmation" to a
spot-check of what the drift test enforces structurally.

### 7d. ClinVar matching on the forward strand (ADR-0035)

ClinVar matches a call only when the call's REF equals the row's REF at
the same forward-strand position, and the user carries the row's ALT
(ADR-0035, amended in 2.5). A complemented reading of the same site
(coding-strand C/T with REF C where ClinVar has G>A) does not match: a VCF
is forward strand by specification, and an array call's REF comes from the
rsID identity table at ClinVar's own locus, so a complemented REF means the
input is mis-coordinated. A VCF record whose one-base REF is none of
ClinVar's one-base REFs at the position is counted in
`input.vcf_quality.clinvar_ref_mismatches`.

This is separate from parser normalization: a reader that knows a row's
strand (deCODEme's minus-strand rows) complements the alleles to the plus
strand while reading the file, so those calls reach ClinVar on the forward
strand already.

The unit-test suite pins the contract on rs1801133 (mock ClinVar: REF=G
ALT=A Pathogenic MTHFR). Re-run the test class as a fast in-protocol
sanity check:

```bash
.venv/bin/python -m pytest \
  tests/annotators/test_clinvar.py::TestStrandAwareCarrierMatch -v --no-cov
```

**Expected:** three tests pass:

- `test_hom_ref_both_strands_no_annotation`: hom-ref forward G/G and
  coding C/C both produce zero ClinVar annotations (the user does not
  carry the ALT).
- `test_complemented_ref_does_not_match`: the het carrier forward G/A
  (REF G) produces one ClinVar annotation; the coding-strand reading C/T
  (REF C) produces none, and as a VCF record it adds one to
  `ref_mismatches`.
- `test_variant_ref_disagreement_abstains`: a `Variant.ref` that is
  neither the row's REF nor its complement produces no annotation.

For real-data validation against the Section 19 battery, see
`tests/utils/test_allele.py::TestStrandAwareCarrierMatch` (the
helper ClinPGx still uses; covers palindromic / orientation /
multi-allelic branches)
and the wrong-allele safety invariant check at the end of
Section 19 — v2.1.0 produces a larger denominator than v2.0.2
because per-Annotation alt threading (ADR-0035) made more rows
eligible for the exact-`(rsid, alt)` CADD lookup. The absolute
denominator drifts with each gnomAD / CADD refresh; the invariant is
**a non-zero, 100% allele-direct denominator on both builds, with GRCh37
alleles translated through the rsID identity table**. All matches must
remain direct with zero complement-resolved scores — that's the v2.0.1
invariant holding against the broader denominator, while the
identity table replaces the blanket GRCh37 gate with an explicit,
allele-level translation that abstains where the strand is undecidable.

## 8. Report formats

### 8a. HTML report

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/report.html
```

Open `/tmp/allelix-review/report.html` in a browser. Verify:

- Table renders without horizontal overflow
- rsID column is sticky when scrolling
- Columns are sortable (click headers)
- Review Status column appears for ClinVar rows
- Pop. Freq, AM and CADD columns are present and populated on this
  GRCh37 input wherever the rsID identity table translates the allele to
  GRCh38 (palindromic sites and alleles that became the GRCh38 reference
  stay empty)
- ClinPGx AM scores show dimmed caveat indicator
- Row borders are color-coded (red = pathogenic, green = benign)
- Zygosity column shows Heterozygous / Homozygous / Hemizygous / Half call for each row; JSON `zygosity` holds the codes (`heterozygous`, `homozygous`, `hemizygous`, `half_call`, `no_call`)
- On a GRCh38 enrichment report, CADD scores are color-coded (red ≥30,
  orange ≥20, gray <20) with tooltips
- "Reading This Report" section is present
- Regulatory notice is present

### 8b. Terminal report (bare-min layout)

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt 2>&1 | head -50
```

**Expected (v2.0.1+):** Rich-formatted table with the bare-min column
set: `rsID | Gene? | Source | Significance | Mag | GT | Condition?`.
Gene and Condition are dropped when no row carries data. Source
displays as `GWAS` (not `GWAS Catalog`); significance drops the
redundant `source_` prefix (`clinvar_pathogenic` → `pathogenic`);
Genotype column is `GT`.

**Intentionally NOT in terminal:** Review Status, Zygosity, Freq,
AM, CADD. The terminal is a quick-extract view; deep enrichment
belongs in HTML/JSON. Pre-v2.0.1 those columns were rendered but
Rich auto-squeezed them to hairline-zero widths on typical 100–120
col terminals.

### 8c. JSON report

```bash
python3 -c "
import json, sys
data = json.load(open('/tmp/allelix-review/enrichment_check.json'))
print(f\"Schema version: {data.get('schema_version')}\")
print(f\"Annotations: {len(data['annotations'])}\")
print(f\"Sources: {set(a['source'] for a in data['annotations'])}\")
has_af = sum(1 for a in data['annotations'] if a.get('allele_frequency') is not None)
has_am = sum(1 for a in data['annotations'] if a.get('am_pathogenicity') is not None)
print(f\"With gnomAD freq: {has_af}\")
print(f\"With AM: {has_am}\")
gwas = [a for a in data['annotations'] if a.get('source') == 'gwas']
print(f\"GWAS rows w/ structured trait: {sum(1 for a in gwas if a.get('trait'))}\")
print(f\"GWAS rows w/ p_value: {sum(1 for a in gwas if a.get('p_value') is not None)}\")
clinpgx = [x for x in data.get('license_attributions') or [] if x.get('source') == 'ClinPGx']
print(f\"ClinPGx attribution URL: {[x.get('source_url') for x in clinpgx]}\")
"
```

**Expected:** Schema version 11 for 2.5.2, which is also the
`SCHEMA_VERSION` constant in `allelix/reports/json_report.py` (ADR-0033
governs schema bump policy). This is the
GRCh38 `enrichment_check.json` produced in §7: multiple sources are
present and gnomAD and AM enrichment counts are both > 0. The paired
GRCh37 report in §7 reaches the same three sources through the rsID
identity table. `license_attributions` has a `ClinPGx` entry with
`source_url` `https://www.clinpgx.org`. On any file
with GWAS hits, the structured
`trait` and `p_value` counts must both be > 0 — that's the v5 contract
delivery check.

## 9. Stats, extract, and focused reports

```bash
allelix stats test_data/real/23andme/user1190_v5.txt
allelix extract --snps rs1801133,rs429358,rs7412 test_data/real/23andme/user1190_v5.txt
```

**Expected:** Stats shows SNP count, no-call rate, het rate. Extract
returns the requested SNPs with genotypes.

### 9a. Focused subcommands

```bash
allelix methylation test_data/real/23andme/user1190_v5.txt
```

**Expected:** Methylation pathway report with annotations from the
methylation gene panel. Non-zero annotation count.

```bash
allelix pharmacogenomics test_data/real/23andme/user1190_v5.txt
```

**Expected:** ClinPGx-focused report. Non-zero annotation count.

### 9b. Compare

```bash
allelix compare test_data/real/23andme/user1190_v5.txt test_data/real/myheritage/user1190.csv
```

**Expected:** Per-chromosome concordance table, coverage overlap stats, and
the following summary for these same-biology user1190 transcodes:

```text
Concordant                 941,911   98.05%
Strand-flip match                0    0.00%
Discordant                       0    0.00%
Strand-ambiguous                 0    0.00%
No-call (either file)       18,702    1.95%
```

The headline 98.05% uses all 960,613 overlapping variants as its denominator,
including the 18,702 variants where either file is a no-call. It therefore
tracks the input no-call rate as well as genotype agreement and is **not** the
pass criterion. The called-genotype denominator is `Concordant +
Strand-flip match + Discordant + Strand-ambiguous = 941,911`; the observed
concordance among called genotypes is `941,911 / 941,911 = 100.0000%`.

**Pass/fail rule:** `Discordant = 0` **and** `Strand-flip match = 0`. Any
non-zero value in either row fails this same-biology cross-format comparison,
regardless of the headline percentage. A lower headline can be caused by
no-calls or by genuine disagreement, so never infer pass or fail from that
percentage; read the `Discordant` and `Strand-flip match` rows directly.

## 10. Config system

Commercial mode drops SNPedia silently (nothing is printed), so the
check reads the report. The previous value of `license.commercial` is
saved and restored, whatever the analysis does:

```bash
allelix config show
prev=$(allelix config get license.commercial)
allelix config set license.commercial true
allelix analyze test_data/real/23andme/user1190_v5.txt --no-update \
  --output /tmp/allelix-review/commercial.json; rc=$?
allelix config set license.commercial "$prev"
echo "analyze exit $rc; license.commercial restored to $(allelix config get license.commercial)"

python3 -c "
import json
d = json.load(open('/tmp/allelix-review/commercial.json'))
assert 'snpedia' not in {a['name'] for a in d['annotators']}
assert not any(a['source'] == 'snpedia' for a in d['annotations'])
p = json.load(open('/tmp/allelix-review/user1190_23andme.json'))
assert 'snpedia' in {a['name'] for a in p['annotators']}
print('SNPedia excluded in commercial mode, present in personal mode')
"
```

**Expected:** `analyze exit 0`, `license.commercial` restored to its
previous value, and `SNPedia excluded in commercial mode, present in
personal mode`. The personal-mode report is §5a's.

## 11. CADD opt-in and commercial gating

CADD is disabled by default. Step 7 exercised the complete opt-in path:
the download flag produced the license prompt and installed the cache,
then `sources.cadd = true` enabled it for analysis. Verify that state and
the resulting enrichment:

```bash
allelix config show | grep cadd
# Expected: sources.cadd = true

allelix db status | grep -i cadd
# Expected: CADD shows "yes" in Ready column with version "v1.7"

# §7 already analysed the GRCh38 benchmark with CADD enabled; count its scores.
python3 -c "
import json
data = json.load(open('/tmp/allelix-review/enrichment_check.json'))
has_cadd = sum(1 for a in data['annotations'] if a.get('cadd_phred') is not None)
print(f'With CADD score: {has_cadd}')
"
# Expected: Non-zero CADD enrichment count on the GRCh38 benchmark.
```

Verify commercial mode gates CADD. A small GRCh38 fixture that gets CADD
scores is analysed once in each mode; `license.commercial` is saved and
restored:

```bash
F=tests/fixtures/mock_vcf_chr_prefix_grch38.vcf
allelix analyze "$F" --no-update --output /tmp/allelix-review/cadd_personal.json
prev=$(allelix config get license.commercial)
allelix config set license.commercial true
allelix analyze "$F" --no-update --output /tmp/allelix-review/cadd_commercial.json; rc=$?
allelix config set license.commercial "$prev"
echo "analyze exit $rc; license.commercial restored to $(allelix config get license.commercial)"

python3 -c "
import json
def load(name):
    d = json.load(open(f'/tmp/allelix-review/{name}.json'))
    return {a['name'] for a in d['annotators']}, sum(a.get('cadd_phred') is not None for a in d['annotations'])
names, scores = load('cadd_personal')
assert 'cadd' in names and scores > 0, (sorted(names), scores)
names, scores = load('cadd_commercial')
assert 'cadd' not in names and scores == 0, (sorted(names), scores)
print('CADD scores in personal mode, none in commercial mode')
"
# Expected: analyze exit 0; license.commercial restored; the last line printed.
```

If full mode is available (pysam installed + local tabix file). The
previous value of `options.cadd_full` is saved and restored:

```bash
prev_cadd_full=$(allelix config get options.cadd_full)
allelix config set options.cadd_full true
allelix analyze test_data/real/23andme/user1190_v5.txt --no-update 2>&1 | grep -i "cadd\|grch38"
# Expected: CADD scores on this GRCh37 input, looked up in the full file at the
# GRCh38 coordinates the rsID identity table translates each allele to. No build
# warning; untranslatable alleles (a row flagged palindromic whose REF differs
# between the builds, or an allele that is the GRCh38 reference) get no score.
allelix config set options.cadd_full "$prev_cadd_full"
echo "options.cadd_full restored to $(allelix config get options.cadd_full)"
```

## 12. Diff command

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/baseline.json
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/current.json --diff /tmp/allelix-review/baseline.json
```

**Expected:** Diff reports no changes (same input, same databases).

## 13. Database update signals

This section runs against a scratch copy of the databases, so it cannot
change the snapshot every other section uses: if a source changed
upstream, only the copy is refreshed. The copy needs as much free disk
as the caches (~24 GB with CADD; on macOS APFS, `cp -cR` clones it
without using space).

```bash
SCRATCH_DB="$(mktemp -d)/databases"
cp -R "$ALLELIX_DATA_DIR" "$SCRATCH_DB"
allelix db update --cadd --data-dir "$SCRATCH_DB"
```

**Expected:** Most annotators show "already current". Per-annotator
states:

- ClinVar, GWAS Catalog (server-driven): "already current" or "can't
  be verified" (ETag/sidecar-dependent)
- ClinPGx (server-driven, CPIC-API dependent): "already current" or
  "can't be verified"
- gnomAD, AlphaMissense, SNPedia, CADD (code-driven, ADR-0030): always
  "already current" — refresh only via `--force` or code bump of
  pinned commit SHA. CADD only appears if previously downloaded via
  `--cadd`.

No re-downloads.

```bash
allelix db update --cadd --force --data-dir "$SCRATCH_DB"
rm -rf "$SCRATCH_DB"
```

**Expected:** All annotators re-download and show green checkmarks.
One exception: ClinVar never refreshes from a download it cannot verify,
even with `--force` (ADR-0053). If NCBI's checksum for `variant_summary`
or `submission_summary` cannot be fetched, ClinVar reports "not refreshed"
or "refresh aborted", keeps its existing caches, and the command exits 1
listing `clinvar` as failed. That is a network failure, not a pass: rerun
when NCBI answers.
Note: `--force` semantics differ by tier. Server-driven sources
override a "signal matches" skip; code-driven sources have no
signal-match path to override — `--force` is the only way to
re-trigger their download because pinned URLs are deterministic.
See ADR-0030.

## 14. GWAS Catalog real-data sanity (slow tests, auto-fetch)

These tests load the real GWAS Catalog (auto-fetched from EBI on
first run as of v2.0.2 — no manual download step) and verify that
the magnitude scoring formula produces bounded output.

```bash
python -m pytest tests/test_end_to_end.py -k "TestRealDataGwasSanity" -v --no-cov
```

**Expected:** 2 tests **pass** (not skip). Default floor (9.0) keeps
output under 50 rows. Old floor (7.0) produces more output than new
floor. First run auto-fetches `test_data/gwas_catalog.zip` (~65 MB,
~30 s on a decent connection); subsequent runs are offline-fast.

A silent skip here is a ship-gate defect (policy): either the
auto-fetch failed and was logged as `OSError`, or the test was
deselected. Investigate, do not tag.

Run the v2.3.1 effect-type and SNPedia-orientation boundary tests
explicitly as named release gates:

```bash
python -m pytest tests/annotators/test_gwas.py::TestMagnitudeScoring \
  tests/annotators/test_snpedia.py::TestAnnotateGenotype::test_minus_orientation_normalizes_genotype_and_alt \
  tests/annotators/test_snpedia.py::TestAnnotateGenotype::test_unknown_orientation_palindrome_abstains_for_homozygote \
  tests/annotators/test_snpedia.py::TestAnnotateGenotype::test_known_minus_orientation_resolves_palindromic_homozygote \
  --no-cov -v
```

**Expected:** All tests pass. Only a bare positive confidence interval
identifies an odds ratio and permits its magnitude modifier; beta,
unit-bearing, missing, ambiguous, or non-positive effects abstain.
Known SNPedia minus orientation normalizes both genotype and ALT, while
unknown palindromic orientation abstains.

## 15. Edge case files

```bash
# Build mismatch detection (analyze and stats run the same build decision).
# The banner line must match exactly; a notice may wrap, so its first line is
# matched as a prefix.
allelix analyze test_data/edge_cases/mhg_grch38_with_grch37_header.txt \
  > /tmp/allelix-review/15-mismatch.log 2>&1
grep -Fx 'Build: GRCh38 (from known-SNP positions: 7/7 match GRCh38)' \
  /tmp/allelix-review/15-mismatch.log
grep -F 'Build mismatch: the file says GRCh37, but its known-SNP positions say GRCh38' \
  /tmp/allelix-review/15-mismatch.log
# Expected: both greps print their line (exit 0).

# 23andMe canonical-header check — this file should NOT be recognized as 23andMe
allelix stats test_data/edge_cases/23andme_lookalike_rejected_by_PA.txt 2>&1
# Expected: rejected, exit non-zero: "Allelix doesn't recognise the format of
# …" (before 2.5's message change: "No parser recognized"). A file whose
# source Allelix cannot identify is rejected, not read (2.5).

# Genes for Good — 23andMe-format export from a different service
allelix stats test_data/edge_cases/23andme_format_from_genes_for_good_service.txt
# Expected: Recognized as 23andMe format, stats displayed

# GRCh36 FTDNA file (no build line; positions decide)
allelix analyze test_data/edge_cases/ftdna_grch36_positions.csv \
  > /tmp/allelix-review/15-grch36.log 2>&1
grep -Fx 'Build: GRCh36 (from known-SNP positions: 4/4 match GRCh36)' \
  /tmp/allelix-review/15-grch36.log
grep -F 'Warning: GRCh36 (hg18) detected from the known-SNP positions.' \
  /tmp/allelix-review/15-grch36.log
# Expected: both greps print their line (exit 0). The warning goes on to say
# that ClinVar and the GRCh38 enrichment are skipped. (Before the
# chr10 GRCh36 seed positions were corrected this was a tentative 3/4,
# kept on GRCh36 by ADR-0025's fail-safe.)

# Despite its name, this is a SNP,Chr,Allele1,Allele2 file with no position
# column (not deCODEme's layout; deCODEme is §5l). 2.5 does not read files
# without positions.
allelix stats test_data/edge_cases/unsupported_decodeme.txt 2>&1
# Expected: rejected, as the lookalike above

# Previously-unsupported, now-supported: the 23andMe exome VCF is a valid
# VCFv4.1 file and is recognized by the VCF parser added after this fixture
# was named. Pin the regression — if the rejection ("doesn't recognise the
# format", or "No parser recognized" before 2.5) comes back, VCF detection
# broke.
allelix stats test_data/edge_cases/unsupported_23andme_exome_vcf.txt 2>&1
# Expected: parser auto-detected as VCF, stats displayed (the
# `unsupported_` prefix is historical — the fixture name predates VCF
# support and is retained for git-mtime continuity).
```

## 16. PLINK export

### 16a. Basic export

```bash
allelix export plink test_data/real/23andme/user1190_v5.txt -o /tmp/allelix-review/user1190 --build grch37
```

**Expected:** Exit code 0. Three files produced: `user1190.bed`,
`user1190.bim`, `user1190.fam`. Console shows variant count, no-call
skip count, and monomorphic marker count.

Verify file structure:

```bash
python3 -c "
data = open('/tmp/allelix-review/user1190.bed', 'rb').read()
assert data[:3] == bytes([0x6C, 0x1B, 0x01]), 'Bad BED magic'
print(f'BED: {len(data)} bytes, magic OK')
bim = open('/tmp/allelix-review/user1190.bim').readlines()
print(f'BIM: {len(bim)} variants')
fam = open('/tmp/allelix-review/user1190.fam').read().strip()
print(f'FAM: {fam}')
assert len(bim) == len(data) - 3, 'BIM/BED row count mismatch'
print('BIM/BED alignment OK')
"
```

**Expected:** BED magic bytes correct. BIM variant count matches
BED data bytes (one byte per variant in SNP-major, single-sample mode).
FAM has one sample line.

### 16b. GRCh37 ref/alt resolution

```bash
allelix export plink test_data/real/23andme/user1190_v5.txt -o /tmp/allelix-review/user1190_gnomad --build grch37
python3 -c "
lines = open('/tmp/allelix-review/user1190_gnomad.bim').readlines()
with_alt = sum(1 for l in lines if l.strip().split('\t')[5] != '0')
mono = len(lines) - with_alt
bed = open('/tmp/allelix-review/user1190_gnomad.bed', 'rb').read()[3:]  # one sample: one byte per marker
het = sum(1 for b in bed if (b & 0b11) == 0b10)
print(f'With ref/alt: {with_alt} (heterozygous: {het}, resolved homozygous: {with_alt - het})')
print(f'Monomorphic (A2=0): {mono}')
"
```

**Expected:** Heterozygous markers always carry both alleles. On GRCh37
a homozygous marker gets its ref/alt from the rsID identity table's
GRCh37 pair at its position, else from gnomAD's pair when that pair is
not strand-ambiguous (A/T, C/G) and the call fits it as written; the
console's "Allele coding:" line counts each source. Neither covers most
array sites (gnomAD's cache holds exome variants, ADR-0026), so most
homozygous markers stay monomorphic (user1190: about 59% of all markers
are monomorphic). Pass when the ref/alt count is greater than the
heterozygous count (some homozygous markers were resolved) and the
monomorphic count equals the CLI's "N markers (P%) exported as
monomorphic" line.

### 16c. Roundtrip with plink2 (optional)

If `plink2` is installed:

```bash
plink2 --bfile /tmp/allelix-review/user1190 --freq --out /tmp/allelix-review/freq_check
```

**Expected:** plink2 reads the files without error. Frequency report
produced.

## 17. Preserve RC artifacts until verification is complete

Do not clean the temporary reports or `$ALLELIX_DATA_DIR` yet. Sections
18–20 reuse them, and the final review needs the exact logs and outputs
from this RC commit. Cleanup commands appear after the pass criteria.

## 18. v2.3.0 → v2.3.1 cache-upgrade verification

SNPedia is the only v2.3.1 cache migration. Parser format 7 adds the
source-orientation column and tag. A v2.3.0 cache with raw pages must be
recognized as stale, reparsed in place, deduplicated, and stamped
`pv:7`. The GWAS odds-ratio/beta change is runtime scoring, so the GWAS
cache does not rebuild. ClinVar, ClinPGx, gnomAD, AlphaMissense, and CADD
also have no cache-schema bump.

Run the deterministic upgrade matrix against isolated temporary caches:

```bash
python -m pytest \
  tests/databases/test_snpedia_parser.py::TestParseRawPages::test_reparse_migrates_old_structured_schema \
  tests/databases/test_snpedia_parser.py::TestParseRawPages::test_reparse_clears_old_data \
  tests/databases/test_snpedia_parser.py::TestParserVersionStamp::test_parser_is_current_false_old_version \
  tests/databases/test_snpedia_parser.py::TestParserVersionStamp::test_parse_stamps_version \
  --no-cov -v
```

**Expected:** All four tests pass. The migrated table contains
orientation, old structured rows do not survive a reparse, an old tag
is rejected, and the rebuilt cache carries `pv:7`. The fresh real-cache
construction in step 3 separately proves the current download path.

## 19. Gold-standard real-data VCF battery (~1.5 GB, required)

The strongest end-to-end check is the **chain-of-trust** workflow
the v2.0.1 ship used: run an analysis against the NIST GIAB truth
set (curated, multi-platform-validated variants), then run the same
sample through a variant caller's raw output, and assert the
variant-caller run is a superset of the truth set at the production
filter. Confirms allelix is faithfully reflecting whatever VCF it's
handed, regardless of caller.

### Bundled VCF fixtures

See the **Section 4 "Tarball contents" → `real/vcf/`** subsection
above for the authoritative file list, sizes, and per-file notes.
The battery in this section exercises those bundled files at WGS
scale.

**Sample coverage caveat:** the tarball bundles GIAB HG002 benchmarks
(GRCh37 + GRCh38) — the curated truth set — plus a GATK HaplotypeCaller
gVCF on **HG00187** (1000 Genomes Finnish sample, not HG002). **No
real caller output for HG002 is bundled.** Same-sample
caller-vs-truth chain-of-trust (DeepVariant or GATK-HC against the
GIAB HG002 truth) requires the optional extras listed below.

### Canonical upstream URLs

All verified 2026-06-16 (HTTP 200).

- **GIAB v4.2.1 HG002 GRCh38 benchmark — VCF + tabix index + region BED**:
  ```
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
  ```
  Files at that path:
  - `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (the truth-set VCF — same data as our bundled `HG002_GRCh38_benchmark.vcf.gz`, canonical NIST filename)
  - `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi` (tabix index)
  - `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` (high-confidence regions, with the "noinconsistent" stricter-truth subset of `..._benchmark.bed`)
- **GIAB v4.2.1 HG002 GRCh37 benchmark**: same path, swap `GRCh38` → `GRCh37`.
- **GIAB raw HG002 sequencing** (for regenerating caller outputs from FASTQ):
  ```
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/
  ```
- **DeepVariant case-study pipeline** (regenerate HG002 gVCF):
  `https://github.com/google/deepvariant/tree/r1.10/docs` — case-study guide runs DeepVariant against GIAB raw FASTQ. Multi-hour compute on a GPU machine.
- **GATK HaplotypeCaller pipeline** (regenerate HG002 GATK-HC):
  GATK 4.x best-practices short-variant calling against `hs37d5` (GRCh37) or GRCh38 with the GIAB HG002 raw FASTQ.

### Battery commands (against the bundled fixtures)

```bash
mkdir -p /tmp/allelix-review

# 1. GIAB GRCh38 truth set
allelix analyze test_data/real/vcf/HG002_GRCh38_benchmark.vcf.gz \
  --build grch38 \
  --min-magnitude 0 \
  --output /tmp/allelix-review/giab_grch38.json \
  --report-format json
# Build: GRCh38 (forced with --build). See "Floor invariants" and
# "Spot-check invariants" below for pass/fail criteria — exact
# counts are NOT pinned (they legitimately drift on dedup, vocab,
# and ClinVar refresh changes; pinning them produced stale-pin
# regressions in §7, §15, and §19 historically).

# 2. GIAB GRCh37 truth set
allelix analyze test_data/real/vcf/HG002_GRCh37_benchmark.vcf.gz \
  --build grch37 \
  --min-magnitude 0 \
  --output /tmp/allelix-review/giab_grch37.json \
  --report-format json
# Per-build ClinVar cache dispatches correctly (ADR-0021). The
# upstream NIST GRCh37 benchmark VCF carries rsIDs in the ID column
# (~95% of rows); the GRCh38 one has none.

# 3. HG00187 GATK-HC gVCF — real caller at WGS scale (GRCh37 build)
# Two runs. Their outputs carry a _s19 suffix: §5i runs the same file at
# the default --min-magnitude under the plain names.
#
# 3a: DEFAULTS (no --min-gq flag). Exercises the GQ floor on a
# low-coverage (~4x) genome: all ~44K non-ref calls fail GQ >= 20,
# producing 0 annotations and 28 high-value no-calls.
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --min-magnitude 0 \
  --output /tmp/allelix-review/hg00187_gatk_default_s19.json \
  --report-format json
#
# 3b: --min-gq 0. Disables the GQ floor and exercises the annotation
# path on position-based resolution. Expected: exit 0 and non-zero
# annotations; the ground-truth harness below (`hg00187_gatkhc_gvcf`)
# checks this run against its floors.
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --min-magnitude 0 \
  --min-gq 0 \
  --output /tmp/allelix-review/hg00187_gatk_mingq0_s19.json \
  --report-format json
# This file is GATK-HC raw output (no rsIDs in ID column), so all
# hits come via position-based ClinVar + rsID-less resolution
# rather than the rsID fast-path used by Step 2. ~19.3M lines in
# the input, ~99.8% are reference blocks that the parser must skip,
# leaving ~44K non-ref variants. Exercises the gVCF parser
# (reference-block skip, alt,<NON_REF> handling, 0/0 hom-ref
# filter) at real WGS scale. Annotation totals here are small
# (tens, not thousands) — see floors below.

# 4. 1000 Genomes chr22 multisample — --sample binding
allelix analyze test_data/real/vcf/thousandG_chr22.vcf.gz 2>&1 | head -3
# Expected: MultiSampleError listing first 10 samples + tail count.
# Re-run with --sample to bind one of the 3,202 samples:
allelix analyze test_data/real/vcf/thousandG_chr22.vcf.gz \
  --sample HG00096 --build grch38 \
  --output /tmp/allelix-review/thousandG_HG00096.json
# Expected: chr22-only annotations for the chosen sample.

# 5. Synthetic edge-case fixtures (parser exercise, no DB lookups needed).
# Use the committed fixtures: their contig lengths, assembly and seed rows
# each name one build (GRCh38, GRCh38, GRCh37), so no --build is needed.
allelix analyze tests/fixtures/mock_vcf.vcf \
  --output /tmp/allelix-review/mock_vcf.json
allelix analyze tests/fixtures/mock_gvcf.g.vcf \
  --output /tmp/allelix-review/mock_gvcf.json
allelix analyze tests/fixtures/mock_vcf_rsidless.vcf \
  --output /tmp/allelix-review/mock_rsidless.json
```

### Ground-truth invariants — published HG002 spot checks + floors

Exact annotation counts are NOT pinned. Pinning them produced
stale-pin regressions three times (§7, §15, §19 historically) and the
last re-pin would have masked an over-filter regression like the one
found when this section was calibrated against v2.2.1. Instead the
battery's pass/fail is decided by **`test_data/check_ground_truth.py`**
against **`test_data/HG002_GROUND_TRUTH.yaml`**. The harness reads
every key in the YAML; a key it does not know fails the run (exit 2),
so the YAML cannot define a check that is never made. It asserts:

0. **Identity per file** (`floors:` in the YAML) — the report is of the
   entry's `file`, on its `build`, and its `input.build_source` is the
   entry's `build_source` (`override`: the battery passes `--build`).

1. **Floor invariants per file** (`floors:` in the YAML) — the written
   annotation count and the unique-key count must clear floors sized at ~80% of
   v2.2.1 (2026-06-18) observed values. A drop below floor
   is an over-filter or vocab-leak regression: classify the removed
   annotations with a v2.2.0→v2.2.x baseline `--diff` run before
   ratcheting the floor down. **"Plausible drift" is not a license
   to re-pin.**

2. **Spot-check invariants per rsID** (`analyze:` in the YAML) —
   concrete published HG002 ground truth at MTHFR rs1801133 and APOE
   rs7412 (both benchmarks) and Factor V rs6025 (GRCh37 benchmark): row counts,
   required sources (`{clinvar, pharmgkb}` for pharmacogenomic
   rsIDs), significance-pattern conformance, condition non-empty.
   The GIAB benchmark file IS the source of truth; the harness just
   verifies allelix reflects it faithfully. Unfalsifiable by drift.

3. **Vocabulary invariants** (`vocabulary:` in the YAML) — every
   emitted `significance` value must be in
   `clinvar_allowed ∪ pharmgkb_allowed`, with explicit forbidden
   placeholders the loader must filter. A novel value means the
   `_CLINVAR_PLACEHOLDER_CLNSIGS` filter regressed (a filter
   leak). The denylist form here is reactive; the CLNSIG drift test
   (`tests/databases/test_clinvar_clnsig_drift.py`) is the durable
   structural form.

4. **Universal invariants** (`universal:` in the YAML) — every
   annotation has the required non-null fields, `source` ∈
   {clinvar, pharmgkb, gwas, snpedia}, `attribution` matches the
   `source`-to-display-name map, `magnitude` ∈ [0, 10].

Run the harness against each of the four battery outputs:

```bash
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/giab_grch38.json giab_grch38_benchmark
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/giab_grch37.json giab_grch37_benchmark
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/hg00187_gatk_mingq0_s19.json hg00187_gatkhc_gvcf
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/hg00187_gatk_default_s19.json hg00187_gatkhc_gvcf_default
```

Each invocation prints a single `✓ <key>: N annotations, M unique
keys, all invariants hold` line on success and exits 0. Any failure
prints the failing assertions to stderr and exits 1.

**To add a new spot check** (e.g. a newly-published HG002 genotype,
or a regression-canary rsID for a specific fix): edit
`HG002_GROUND_TRUTH.yaml` — append to `analyze:` with `applies_to`
scoping the check to whichever battery files surface the rsID. The
harness picks up new entries automatically.

**To recalibrate floors after a legitimate dedup or vocab change**:
edit the `# observed N` comments in `floors:` and re-derive the
`_min` values at ~80% of new observed. Include the diff-classification
output (from the `analyze --diff` run that justified the recalibration)
in the commit message.

### Recommended extras (NOT bundled, downloadable separately)

For the strongest chain-of-trust check — **caller disagreement on
the same sample (HG002)** — these three additional files extend the
battery. None are in the tarball; they're downloadable from the URLs
above or regeneratable from the GIAB raw FASTQ.

| File | Size | Source | Why you want it |
|------|-----:|--------|-----------------|
| `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` | ~12 MB | GIAB v4.2.1 (NIST FTP, link above) | High-confidence region BED for restricting caller-output comparisons to GIAB's confident scope |
| `HG002.child.g.vcf.gz` | ~1.03 GB | DeepVariant 1.10.0 case-study pipeline on HG002 | Same-sample variant-caller output to compare against the GIAB truth set (DeepVariant vs GATK-HC vs truth) |
| `HG002_ALLCHROM_hs37d5_novoalign_Ilmn150bp300X_GATKHC.vcf.gz` | ~266 MB | GATK 3.5 HC on HG002 GRCh37 (PrecisionFDA Truth Challenge pipeline) | Same-sample GRCh37 caller output — useful for per-build chain-of-trust |

Drop them into `test_data/real/vcf/` (gitignored). The chain-of-trust
invariants below assume `HG002.child.g.vcf.gz` is present for the
same-sample superset check; the others are useful but optional.

### Chain-of-trust invariants

**Same-sample superset (with the optional HG002 DeepVariant gVCF
in place):** every annotation in the GIAB benchmark output at the
default filter (magnitude ≥ 5, `analyze`'s default) must also
appear in the DeepVariant gVCF output, which runs at that filter. The
benchmark run above used `--min-magnitude 0`, so its rows are filtered
to magnitude ≥ 5 first: both sides then use the same filter. The gVCF
should surface the benchmark's annotations PLUS additional ones
that fall outside the GIAB high-confidence regions. **0 missing**
at default filter is the pass condition.

```bash
# Requires HG002.child.g.vcf.gz (Optional extras section).
allelix analyze test_data/real/vcf/HG002.child.g.vcf.gz \
  --build grch38 \
  --output /tmp/allelix-review/hg002_gvcf.json \
  --report-format json

python3 -c "
import json
bench = [a for a in json.load(open('/tmp/allelix-review/giab_grch38.json'))['annotations']
         if a['magnitude'] >= 5]
gvcf = json.load(open('/tmp/allelix-review/hg002_gvcf.json'))['annotations']
key = lambda a: (a['source'], a['rsid'], a.get('condition',''), a.get('description',''))
bk = {key(a) for a in bench}
gk = {key(a) for a in gvcf}
missing = bk - gk
print(f'benchmark rows (magnitude >= 5): {len(bench)}  gvcf rows: {len(gvcf)}')
print(f'missing in gvcf at production filter: {len(missing)}')
assert not missing, f'SUPERSET BROKEN — {len(missing)} benchmark rows absent from gvcf'
print('SUPERSET HOLDS')
"
```

**Wrong-allele safety invariants:** every gnomAD AF, AlphaMissense
score and CADD score attached to a battery annotation must be the
database value for the allele the person carries, at that locus. The
same harness as §7 checks each report against the file it analysed
(`--input`); see §7 for what it checks. GRCh37 reports (`giab_grch37`
and both HG00187 runs) are translated through the rsID identity table
by locus and alleles. Denominators: each GRCh38 report and
`giab_grch37` must have enriched rows; the HG00187 runs are too small to
require any (the default run has no annotations at all), so they pass
`--allow-empty`. `hg002_gvcf` is checked when the optional gVCF was run.

```bash
V=test_data/real/vcf
python scripts/check_wrong_allele.py /tmp/allelix-review/giab_grch38.json \
  --input $V/HG002_GRCh38_benchmark.vcf.gz
python scripts/check_wrong_allele.py /tmp/allelix-review/giab_grch37.json \
  --input $V/HG002_GRCh37_benchmark.vcf.gz
python scripts/check_wrong_allele.py /tmp/allelix-review/hg00187_gatk_mingq0_s19.json \
  --input $V/HG00187_gatkhc.g.vcf.gz --allow-empty
python scripts/check_wrong_allele.py /tmp/allelix-review/hg00187_gatk_default_s19.json \
  --input $V/HG00187_gatkhc.g.vcf.gz --allow-empty
if [ -f /tmp/allelix-review/hg002_gvcf.json ]; then
  python scripts/check_wrong_allele.py /tmp/allelix-review/hg002_gvcf.json \
    --input $V/HG002.child.g.vcf.gz
fi
```

**Expected:** every call exits 0 with `0 non-direct, 0 mismatch,
0 unresolved`.

**ClinVar per-SCV row shape:** multi-SCV variants surface
one row per submission, not a single semicolon-joined row. The
Frankenstein-pair safety question raised against the old VCF
loader is now structurally impossible — each row carries one
submitter's (significance, condition) pair from
`submission_summary.txt.gz`, so no cross-product pairing exists.
The legacy semicolon-joined pin from the VCF era is removed.
See §7 for the post-update real-cache ship-gate (per-SCV row shape
+ significance-sentinel scan).

### When to run this battery

- **Before every ship-gate tag.** v2.0.1's load-bearing chain-of-trust check.
- **After any change to** `_pipeline.py`, the VCF parser, `iter_clinvar_records`, or strand / enrichment logic.
- **Not for routine dev work** — unit tests cover the same logic against synthetic fixtures.

## 20. Annotated VCF output (v2.3.0+; ADR-0036)

Verifies that `allelix analyze --vcf-out` produces a spec-conformant
annotated VCF: preserves every existing input header line verbatim
except Allelix's own (an input Allelix already annotated has its
`##ALLELIX_...` lines and `##INFO=<ID=ALLELIX_...` declarations replaced
by the current run's, ADR-0059),
synthesizes declarations for emitted contigs that the input omitted,
injects the provenance / attribution / `##INFO` block immediately before `#CHROM`,
stamps the per-source INFO fields with per-allele `Number=A` semantics
(the record-level `ALLELIX_CLINVAR_REF`, `ALLELIX_ORIGINAL_ID` and
`ALLELIX_FILTER` are `Number=1`), and stamps
recovered rsIDs into the ID column with the original value preserved
under `ALLELIX_ORIGINAL_ID`.

Requires `bcftools` on `$PATH` for §20.6. Skip §20.6 if unavailable
locally; CI installs bcftools unconditionally and runs the
round-trip test as part of the hosted release gate.

All §20 outputs go into one new directory per run: the writer refuses to
overwrite an existing file, so fixed paths would fail on a second run.
Create it first, in the shell used for the rest of this section:

    mkdir -p /tmp/allelix-review
    export VCF_OUT="$(mktemp -d /tmp/allelix-review/vcf-out.XXXXXX)"

### 20.1: Smoke test on a small VCF fixture

    allelix analyze tests/fixtures/mock_vcf.vcf --vcf-out $VCF_OUT/mock_annotated.vcf

Verify:

- [ ] Exit code 0
- [ ] `$VCF_OUT/mock_annotated.vcf` exists and is non-empty
- [ ] Terminal prints `Wrote annotated VCF with N stamped row(s) to $VCF_OUT/mock_annotated.vcf`
- [ ] `head -30 $VCF_OUT/mock_annotated.vcf` shows:
  - [ ] `##fileformat=VCFv4.2` on line 1 (preserved from input)
  - [ ] Every input `##` line preserved verbatim before `#CHROM`, except
    `##ALLELIX_...` lines and `##INFO=<ID=ALLELIX_...` declarations, which
    are replaced (the mock fixture has none)
  - [ ] Every CHROM value emitted by the fixture has a corresponding
    `##contig=<ID=...>` declaration; declarations missing from the input
    are synthesized without renaming bare or `chr`-prefixed contigs
  - [ ] Injected `##INFO=<ID=ALLELIX_*,Number=A,Type=...,Description="...">`
    declarations for every source that produced annotations; the
    record-level `ALLELIX_CLINVAR_REF`, `ALLELIX_ORIGINAL_ID` and
    `ALLELIX_FILTER`, when present, are `Number=1`
  - [ ] Provenance block: `##ALLELIX_Version=<current version>`,
    `##ALLELIX_VCF_SchemaVersion=0.6.0`,
    `##ALLELIX_SkippedSources=` (each skipped source and its reason, or `none`),
    `##ALLELIX_License=AGPL-3.0-or-later`,
    `##ALLELIX_RunDate=<ISO 8601 UTC>`,
    `##ALLELIX_Build=<ID=GRCh37|GRCh38,Source=override|contig_lengths|seeds|assembly_header|file_label|default>`
  - [ ] `##ALLELIX_Database=<Name=X,Version=Y>` line for every entry in
    the run's `annotators_used`
  - [ ] `##ALLELIX_Attribution=<Name=X,License=Y,URL=Z>` line for every
    licensed source

### 20.2: Multi-allelic site preserves per-allele semantics

Use the §20.3 output: `tests/fixtures/mock_gvcf.g.vcf` annotated with
`--vcf-out` stamps three two-ALT rows (`ALT=G,<NON_REF>` and
`ALT=A,<NON_REF>`). Run §20.3's command first, then check every stamped
multi-ALT row:

```bash
python3 - <<'PY'
import os, re

def rows(path):
    return [l.rstrip('\n').split('\t') for l in open(path) if not l.startswith('#')]

out = os.path.join(os.environ['VCF_OUT'], 'mock_gvcf_annotated.g.vcf')
header = [l for l in open(out) if l.startswith('##INFO=<ID=ALLELIX_')]
per_allele = {re.search(r'ID=(\w+)', h).group(1) for h in header if 'Number=A' in h}
inputs = {(r[0], r[1]): r for r in rows('tests/fixtures/mock_gvcf.g.vcf')}
stamped = 0
for r in rows(out):
    alts = r[4].split(',')
    assert r[4] == inputs[(r[0], r[1])][4], ('ALT rewritten', r[:5])
    fields = dict(kv.split('=', 1) for kv in r[7].split(';') if '=' in kv)
    mine = {k: v for k, v in fields.items() if k in per_allele}
    if len(alts) < 2 or not mine:
        continue
    stamped += 1
    for key, value in mine.items():
        assert len(value.split(',')) == len(alts), (r[:5], key, value)
assert stamped > 0, 'no stamped multi-ALT row'
print(f'{stamped} stamped multi-ALT rows; every Number=A field has one slot per ALT')
PY
```

- [ ] The `ALT` column is the input's ALT list verbatim — the writer
  never rewrites ALT
- [ ] For every stamped `ALLELIX_*` field on these rows, the number
  of comma-separated slots equals the number of ALTs (`Number=A`
  contract). Empty slots render as `.`, e.g.
  `ALLELIX_CLINVAR=.,clinvar_pathogenic|rs...|GENE` on a two-ALT
  row where only the second ALT has a ClinVar row. Commas inside a
  value are written as `_`, so counting commas is safe
- [ ] `bcftools query -f '%CHROM\t%POS\t%INFO/ALLELIX_CLINVAR\n'`
  on the row returns the comma-joined value verbatim (bcftools
  honors `Number=A` when the header declares it)

Rows where no source annotates any ALT get no per-allele `ALLELIX_*`
fields, which is expected passthrough, not a failure. Two cases still
add fields to such a row:

- A record whose call carries the reference allele (`0/1`, `0/.`, `0/0`,
  haploid `0`, or a single-position record with ALT `.`) at a site where
  the reference itself carries a ClinVar allele gets that finding in
  `ALLELIX_CLINVAR_REF` (`Number=1`, VCF schema 0.6.0, ADR-0059), never in
  an ALT slot. gVCF reference blocks (`END=`) are never stamped; their
  findings stay in the JSON and HTML reports.
- A row whose rsID was recovered gets the new ID in the ID column and
  `ALLELIX_ORIGINAL_ID`, whatever its annotations (§20.4).

### 20.3: gVCF reference blocks and `<NON_REF>` pass through

    allelix analyze tests/fixtures/mock_gvcf.g.vcf --vcf-out $VCF_OUT/mock_gvcf_annotated.g.vcf

Verify:

- [ ] Reference-block rows (`.  A  <NON_REF>  ...  END=<pos>`) are
  byte-identical to the input row — no `ALLELIX_*` INFO stamped
- [ ] Mixed rows (`ALT=A,<NON_REF>`) have per-allele INFO with `.` in
  the symbolic-ALT slot
- [ ] `##ALT=<ID=NON_REF,...>` header line preserved

### 20.4: rsID recovery visible in ID column with `ALLELIX_ORIGINAL_ID`

    allelix analyze tests/fixtures/mock_vcf_rsidless.vcf --vcf-out $VCF_OUT/mock_rsidless_annotated.vcf

For at least one row where the input ID column was `.` and the pipeline
recovered an rsID:

- [ ] Output ID column carries the recovered `rs...` value
- [ ] Output INFO column contains `ALLELIX_ORIGINAL_ID=.`
- [ ] `##INFO=<ID=ALLELIX_ORIGINAL_ID,Number=1,Type=String,...>` present
  in header

### 20.5: Multi-sample VCF passes sample columns through verbatim

    allelix analyze tests/fixtures/mock_multisample.vcf --sample SAMPLE_A --vcf-out $VCF_OUT/mock_multi_annotated.vcf
    diff <(grep -v '^##' tests/fixtures/mock_multisample.vcf | cut -f 9-) \
         <(grep -v '^##' $VCF_OUT/mock_multi_annotated.vcf | cut -f 9-)

Note: strip ``##`` header lines on **both** sides. The bare
``cut -f 9-`` prints tab-less lines whole (both GNU and BSD ``cut``
behaviour), so leaving the input's ``##`` lines in the left stream
produces a spurious diff against the output's injected provenance
block.

Verify:

- [ ] Sample columns (FORMAT + every sample) are byte-identical to input
- [ ] ``--sample SAMPLE_A`` restricts the annotation pass to that
  sample's genotypes only; the annotated output still carries every
  sample column (write pass is genotype-agnostic)

### 20.6: `bcftools view` round-trip validates the output

    bcftools view $VCF_OUT/mock_annotated.vcf > /dev/null

Verify:

- [ ] Exit code 0
- [ ] No warnings on stderr. The mock fixture declares only
  ``##contig=<ID=1,...>`` and ``##contig=<ID=22,...>`` but emits records
  on `1`, `22`, `19`, `17`, `chr1`, `chrX`, and `chrM`. The writer must
  synthesize declarations for all five missing raw names without
  normalizing or collapsing the mixed naming styles. Capture stderr:

  ```bash
  bcftools view $VCF_OUT/mock_annotated.vcf > /dev/null 2> $VCF_OUT/bcf_stderr.txt
  test ! -s $VCF_OUT/bcf_stderr.txt
  ```

  and require exit code 0
- [ ] `bcftools view -h $VCF_OUT/mock_annotated.vcf | grep '^##INFO=<ID=ALLELIX_'`
  lists every declared Allelix INFO field
- [ ] `bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/ALLELIX_CLINVAR\n' $VCF_OUT/mock_annotated.vcf`
  on a known-hit row (rs1801133 in mock_vcf) returns the value
  the writer stamped verbatim, e.g.
  ``1  11796321  rs1801133  clinvar_drug_response|rs1801133|MTHFR``

Skip if bcftools is unavailable locally. Hosted CI runs it on pushes
to the public `main`, on pull requests into it and on release tags, not
on development merges; there, the local full suite with bcftools
installed (§2) is the check.

### 20.7: Provenance headers

The `##ALLELIX_` lines record what produced the file: versions, build,
date, skipped sources and attribution. They are not enough to repeat the
run: they do not record the sample, filters, panel or quality options.
To reproduce a run, keep the exact command line and configuration, the
input file and the same database caches.

    grep '^##ALLELIX_' $VCF_OUT/mock_annotated.vcf

Verify:

- [ ] `##ALLELIX_Version` — Allelix version string
- [ ] `##ALLELIX_VCF_SchemaVersion` — the output-contract version (0.6.0)
- [ ] `##ALLELIX_SkippedSources` — each skipped source and its reason, or `none`
- [ ] One `##ALLELIX_Database=<Name=X,Version=Y>` per annotator used
- [ ] `##ALLELIX_Build=<ID=...,Source=...>`: `ID` matches `--build` and
  `Source=override` if supplied, else the decided build and its source
  (ADR-0044)
- [ ] `##ALLELIX_RunDate` in ISO 8601 UTC (`+00:00` suffix)

### 20.8 (negative case): Non-VCF input rejected

    allelix analyze tests/fixtures/mock_23andme.txt --vcf-out $VCF_OUT/should_not_exist.vcf

Verify:

- [ ] Non-zero exit code
- [ ] Error message names the format: `--vcf-out requires VCF or gVCF input`
- [ ] `$VCF_OUT/should_not_exist.vcf` does NOT exist (no partial file)

### 20.9 (negative case): Existing output not silently clobbered

    touch $VCF_OUT/mock_annotated.vcf
    allelix analyze tests/fixtures/mock_vcf.vcf --vcf-out $VCF_OUT/mock_annotated.vcf

Verify:

- [ ] Non-zero exit code
- [ ] Error message: `refuses to overwrite`
- [ ] File is untouched (mtime unchanged)

### 20.10: `chr`-prefixed contigs stamp annotations

Every real GRCh38 VCF from modern callers (GATK / DeepVariant /
DRAGEN) uses ``chr``-prefixed contigs (``chr22`` etc.). The parser
normalizes to bare form for the pipeline's annotation dict, so the
writer must symmetrically normalize when reading the raw input line
back — otherwise every stamp silently misses.

    allelix analyze test_data/real/vcf/HG002_GRCh38_chr22.vcf.gz \
        --vcf-out $VCF_OUT/hg002_chr22_annotated.vcf

Verify:

- [ ] Exit code 0
- [ ] Terminal reports NON-ZERO stamped rows (real chr22 will produce
  thousands — a "0 stamped rows" line here is the chr-prefix
  regression, not sparse annotation)
- [ ] `grep 'ALLELIX_CLINVAR' $VCF_OUT/hg002_chr22_annotated.vcf | grep -v '^##' | head`
  shows real chr22 variants with stamped annotations, e.g.
  ``chr22  <pos>  rs<N>  <REF>  <ALT>  ... ALLELIX_CLINVAR=...``
- [ ] Output CHROM column preserves the input's ``chr22`` — writer
  keys lookups via bare-form ``22`` but writes verbatim

### 20.11: Gzip input transparently handled

    allelix analyze tests/fixtures/mock_vcf.vcf.gz --vcf-out $VCF_OUT/mock_gz_annotated.vcf

Verify:

- [ ] Exit code 0
- [ ] Output is plain text VCF (compressed output is a deferred non-goal
  per ADR-0036)
- [ ] Provenance block present in output header

## Pass criteria

All of the following must be true:

- [ ] Unit test suite (§2, run after §3 with `ALLELIX_CI_REQUIRE_CACHE=1`,
      no `-x`): every collected test passes, **0 skipped**, and the
      collected count is at least the latest row of the §2 floor table
- [ ] Ruff lint + format: zero warnings (`ruff check .` and `ruff format --check .` both clean)
- [ ] `db update` downloads all enabled annotators without errors
- [ ] Fresh-cache gates pass for GWAS Catalog and SNPedia: real `db update`
      started from the empty RC data directory, both isolated loader/parser
      suites pass, and SNPedia parser format 7 stores orientation
- [ ] `db status` shows all annotators ready with version and record count; ClinPGx row labeled "ClinPGx" (not "PharmGKB")
- [ ] All 11 parser formats produce successful analysis (23andMe, Genes for Good, AncestryDNA, FTDNA CSV, FTDNA FamFinder, FTDNA Illumina raw, MyHeritage, MHG / Tempus, Living DNA, deCODEme, VCF/gVCF)
- [ ] Cross-parser identity: every user1190 representation exits 0 and all six have the same `(source, rsid, significance)` rows, each the same number of times (6 array-based formats; VCF doesn't have a user1190 representation)
- [ ] VCF flagship feature: rsID-less VCF produces non-zero annotations after build-scoped coordinate resolution (step 5j); multi-sample VCF requires `--sample` (step 5i)
- [ ] Build decision (step 5k): a VCF with no evidence and the mixed-contig fixture stop and ask for `--build`; contig lengths decide whatever the contig names; full-size relabelled and unlabelled chips resolve from positions with the reference file's annotations; 4/5 is confident while 3/4 and 2/2 remain tentative
- [ ] **Wrong-allele safety invariants**: `scripts/check_wrong_allele.py` exits 0 on the §7 GRCh38 and GRCh37 reports: every attached gnomAD AF, AlphaMissense score and CADD score is the database value for the carried allele at its GRCh38 key (GRCh37 keys translated through the identity table by locus and alleles), with 0 non-direct, 0 mismatch and 0 unresolved rows and a non-zero denominator; the real GRCh37 consumer report has non-zero gnomAD / AM / CADD enrichment with their provenance and the rsID identity table listed; ClinVar emits one row per (variant, SCV submission) — the per-SCV pairing structurally retires the Frankenstein-pair hazard of the old VCF loader (v2.2 per-SCV TSV loader)
- [ ] **ClinVar significance-sentinel ship-gate (§7b)**: post-`db update` sentinel scan returns 0 rows on both clinvar.GRCh37.sqlite and clinvar.GRCh38.sqlite — `LOWER(clinical_significance) IN ('-', '', 'not specified', 'no classification provided', 'other', 'association', 'association not found')` (case-insensitive matches loader's `_normalize_clnsig`) — and rs137854557 carries a real classification. Set MUST match `_CLINVAR_PLACEHOLDER_CLNSIGS` in `databases/manager.py` — never broaden one without broadening the other (`not provided` is intentionally excluded; it maps to 2.0 in `_CLNSIG_MAGNITUDE`, safe below the floor). v2.2.1 added the three non-classification curatorial terms. Denylist form is reactive; durable allowlist form lives in the CLNSIG drift test (`test_clinvar_clnsig_drift.py`)
- [ ] **Terminal report**: bare-min columns only (`rsID | Gene? | Source | Significance | Mag | GT | Condition?`); Review Status / Zygosity / Freq / AM / CADD intentionally absent (still present in HTML/JSON)
- [ ] HTML report renders correctly in a browser; "Annotators:" subtitle uses display names (ClinPGx, not pharmgkb); enrichment columns (Review Status, Zygosity, Freq, AM, CADD) are present, and populated on GRCh37 input wherever the rsID identity table translates the allele
- [ ] JSON report has schema version 11, with `input.build_source` and `input.build_evidence`; `sources` lists all eight registered sources with a state, and `report_warnings` has no `clinvar_skipped` entry when ClinVar ran (other entries, such as `input_trailing_bytes` or `forced_format`, are checked on their own merits, §5n); the GRCh38 `enrichment_check.json` and the paired GRCh37 consumer report both have non-zero gnomAD + AM + CADD enrichment with provenance (GRCh37 through the rsID identity table); structured GWAS fields are populated (`annotation.trait` non-empty on GWAS rows, `p_value` parses as float when published, `phecode` populated when the upstream trait carried one); `panel_coverage` appears when `--filter-file` is used, as `{requested, counts, variants}` with the four counts (`annotated`, `typed_no_annotation`, `no_call`, `not_on_chip`) summing to `requested`; array reports have `input.array_quality` and, when ClinVar ran, `input.probe_classes` (both null for VCF input); `license_attributions[].source` shows "ClinPGx" with `source_url` `https://www.clinpgx.org`
- [ ] **ClinVar forward-strand matching (ADR-0035)**: a call matches only when its REF equals ClinVar's at the same forward position and it carries the ALT; a complemented reading of the same het carrier gets no ClinVar annotation (counted as a REF mismatch on VCF input); on the step 7d fixture, whose site's reference does not carry a ClinVar allele, hom-ref on either strand produces zero annotations (at a reference-carried site such as F5 rs6025 on GRCh37, a `0/0` call is assessed for the allele the reference carries, docs/formats.md); a `variant.ref` that fits neither strand abstains (sanity check at step 7d)
- [ ] Config system gates SNPedia on `license.commercial`: the commercial-mode report has no SNPedia provenance or rows, §5a's personal-mode report has SNPedia, and `license.commercial` is back to its previous value
- [ ] GWAS effect-type gate: only positively identified bare odds ratios
      receive a modifier; beta, unit-bearing, missing, ambiguous, and
      non-positive effects abstain
- [ ] SNPedia orientation gate: known minus rows normalize genotype + ALT;
      unknown palindromic rows abstain; known palindromic minus rows resolve
- [ ] CADD opt-in: `db update --cadd` downloads the cache and shows the
      license prompt; `sources.cadd = true` in config enables
      enrichment; GRCh38 scores are non-zero
- [ ] CADD commercial gate: on `tests/fixtures/mock_vcf_chr_prefix_grch38.vcf`, personal mode has `cadd` provenance and CADD scores, commercial mode has neither, and `license.commercial` is back to its previous value
- [ ] Edge case files produce expected behavior
- [ ] `db update` (second run, §13, on a scratch copy of the databases) skips already-current databases
- [ ] **v2.3.0 upgrade path:** all four step-18 migration tests pass;
      SNPedia reparses to parser format 7 with orientation, while the
      other source caches require no schema rebuild
- [ ] GWAS Catalog slow tests pass (auto-fetch the fixture; a skip is a failure)
- [ ] Manual `workflow_dispatch` of `.github/workflows/slow.yml` is green
      against the exact RC commit under review. The `Assert drift sentinel
      executed` step parses the generated JUnit XML and confirms
      `test_no_new_clnsig_values_in_live_cache` executed and passed,
      never skipped or merely absent; require that step to be green.
- [ ] `methylation`, `pharmacogenomics`, and `compare` subcommands produce
      output (`pharmacogenomics --help` says "ClinPGx-style sources"); the
      user1190 cross-format comparison has zero `Discordant` and zero
      `Strand-flip match` rows. Its headline percentage includes no-calls and
      is not the pass criterion.
- [ ] PLINK export produces valid .bed/.bim/.fam with correct magic and alignment
- [ ] PLINK export resolves ref/alt from gnomAD when available (on GRCh37, from the rsID identity table first)
- [ ] `allelix --version` reports the actual pyproject version even from a bare source checkout
- [ ] **Gold-standard VCF battery (step 19)**: the bundled GIAB GRCh38 benchmark, GIAB GRCh37 benchmark, and HG00187 GATK-HC GRCh37 gVCF all analyze cleanly; `check_ground_truth.py` exits 0 on all four (`giab_grch38_benchmark`, `giab_grch37_benchmark`, `hg00187_gatkhc_gvcf`, `hg00187_gatkhc_gvcf_default`) — enforces each report's file, build and build source, the written and unique-key floors, published HG002 spot-checks at rs1801133 / rs7412 / rs6025, the significance vocabulary union, and per-row universal invariants, and fails on any YAML key it does not read; the optional same-sample HG002 DeepVariant gVCF is a strict superset of the benchmark rows at magnitude ≥ 5 (the same filter on both sides) when staged; `scripts/check_wrong_allele.py` exits 0 on every battery report (GIAB GRCh38 and GRCh37 with non-zero denominators); ClinVar emits per-SCV rows (multiple rows per multi-SCV rsID, each with its submitter's exact condition) — semicolon-joined rows are a regression to the previous shape

## Cleanup after verification is posted

Keep the RC artifacts until the reviewer has posted the protocol result
on the release PR. Then reclaim the isolated cache and report directory:

```bash
allelix db clean --yes
rm -f "$ALLELIX_DATA_DIR/config.toml"
rmdir "$ALLELIX_DATA_DIR"
rm -rf /tmp/allelix-review
```

`ALLELIX_DATA_DIR` was created by `mktemp -d` in step 3, so these commands
cannot target the user's normal Allelix cache. `db clean` preserves the
config file intentionally; the explicit file removal lets `rmdir` remove
the now-empty temporary directory. A default cache occupies ~18 GB, or
~24 GB after the optional CADD step.

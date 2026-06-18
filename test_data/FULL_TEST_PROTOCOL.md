# Full Test Protocol

External reviewer checklist for verifying an allelix release against real data.

**Requirements:** Fast machine, fast internet, ~50 GB free disk space.
Estimated wall-clock time: 30–45 minutes (most of it is database downloads).

## 1. Environment setup

```bash
git clone https://github.com/allelix/allelix.git
cd allelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Verify Python 3.11+:

```bash
python --version
```

## 2. Unit test suite (synthetic + auto-fetched real data)

Run the full test suite. The real GWAS Catalog fixture (~65 MB)
auto-fetches from EBI on first run; subsequent runs are offline-fast
against the local cache (GH #45).

```bash
python -m pytest tests/ -x --tb=short
```

**Expected for v2.2.0:** **~1,672 passed, 0 skipped** (fast tier:
`pytest -m "not slow and not integration"`) when `plink2` is
installed locally and the GWAS Catalog auto-fetch succeeds. The
"0 skipped" line is the goal — silent skips are forbidden as a ship-
gate signal (GH #45). If `plink2` isn't installed, expect 1 skip on
`TestRoundtripWithPlink` (`@pytest.mark.integration`, external binary
precondition); install `plink2` before tagging.

Test-count floor by release:
- v1.9.0: ~1,400
- v2.0.0: ~1,486 (VCF/gVCF parser, FTDNA Illumina raw, R-4 CLNSIG
  drift CI, rsID resolution, ClinPGx rename)
- v2.0.1: ~1,525 (audit-driven correctness fixes #16–#27 cluster,
  Variant case normalization, ClinVar CLNDN-join, terminal bare-min)
- v2.0.2: ~1,540 (auto-fetch GWAS fixture #45, chr-prefix build
  detection #38, pyproject version fallback #34, enrichment annotator
  stack-management #36, doc/process cleanup #43/#44/#46/#47/#48)
- v2.1.0: ~1,576 (ADR-0035 Cluster B: Variant.ref + per-Annotation
  alt threading + structured GWAS fields + strand-aware carrier
  matching; pipeline-level Variant.ref population for array data)
- v2.2.0: ~1,672 (#42 per-SCV ClinVar TSV loader + dash-CLNSIG
  guard; #28 streaming ClinVar ref-lookup + permission resolver
  dedup; #70 FTDNA FamFinder parser; #75 panel-coverage states +
  JSON schema_version 5→6; #90/#91 runtime nudges; #79 branch
  coverage; #51/#52 methylation panel + magnitude floor; #77
  db clean / db path subcommands)

Check lint:

```bash
ruff check . && ruff format --check .
```

**Expected:** All checks passed, 0 files reformatted.

## 3. Download all databases

```bash
allelix db update
```

This downloads ClinVar (GRCh37 + GRCh38), ClinPGx (formerly PharmGKB —
the upstream rebranded in 2026; `pharmgkb.org` 301-redirects to
`clinpgx.org`, but the cache file is still named `pharmgkb.sqlite`
internally for backward compatibility), GWAS Catalog, gnomAD
(~2.7 GB compressed download, ~6 GB on-disk SQLite cache),
AlphaMissense (~1.8 GB compressed, ~8 GB on-disk), and SNPedia from
HuggingFace. CADD is opt-in and not included here — see step 11.
Total on-disk footprint after this step: roughly **16 GB**
(CADD adds another ~5.8 GB if you opt in at step 11, bringing the
full install to ~22 GB).

**Expected:** All enabled annotators show green checkmarks. No errors.

Verify status:

```bash
allelix db status
```

**Expected:** All annotators show "yes" in the Ready column with
version strings and record counts. SNPedia should show ~104K records.

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
suite itself** on first use (GH #45), so this script's GWAS-fetch
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
| `user339_file150_yearofbirth_1950_sex_XY.ftdna-illumina.txt` | CSV text | 592 KB | openSNP user 339. Small file — likely partial export or early chip. |
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

**Synthetic fixtures**

| File | Size | Notes |
|------|-----:|-------|
| `mock_vcf.vcf` | 1.2 KB | Plain single-sample VCF with rsIDs. Deterministic test fixture. |
| `mock_gvcf.g.vcf` | 823 B | gVCF with reference blocks (`<NON_REF>` ALT, `END=` INFO). |
| `mock_multisample.vcf` | 404 B | 3-sample VCF (`SAMPLE_A`, `SAMPLE_B`, `SAMPLE_C`). |
| `mock_vcf_rsidless.vcf` | 616 B | rsID-less VCF (`ID=.`). Exercises position-keyed ClinVar resolver (GH #8). |

**GIAB HG002 benchmarks** — Source: [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) v4.2.1, Ashkenazi son (HG002 / NA24385).

| File | Size | Build | Notes |
|------|-----:|-------|-------|
| `HG002_GRCh38_benchmark.vcf.gz` | 150 MB | GRCh38 | Full WGS, ~4.05M variants. |
| `HG002_GRCh38_benchmark.vcf.gz.tbi` | 1.6 MB | | Tabix index. |
| `HG002_GRCh37_benchmark.vcf.gz` | 169 MB | GRCh37 | Full WGS, ~4.03M variants. Ships without rsIDs. |
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
| `thousandG_chr22.vcf.gz` | 496 MB | chr22, 3,202 samples (phased, duohmm). Exercises `--sample <ID>` selection. |
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

### 5a. 23andMe

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/user1190_23andme.json
```

**Expected:** Exit code 0. JSON report written. Should contain ClinVar,
ClinPGx, GWAS, SNPedia, gnomAD, and AlphaMissense annotations. Check
that annotation count is in the hundreds (varies by database version).

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
# Pick any one file from the directory
allelix analyze "$(find test_data/real/ancestrydna -maxdepth 1 -type f | head -1)" \
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

### 5h. VCF / gVCF (committed fixtures + bundled real-scale files)

```bash
# Plain single-sample synthetic VCF
allelix analyze tests/fixtures/mock_vcf.vcf --output /tmp/allelix-review/mock_vcf.json

# Synthetic gVCF (reference blocks present, must be skipped at parse time)
allelix analyze tests/fixtures/mock_gvcf.g.vcf --output /tmp/allelix-review/mock_gvcf.json

# Real-scale gVCF — HG00187 GATK-HC, ~19.3M lines, ~99.8% reference
# blocks. Exercises the gVCF parser end-to-end at WGS scale, not just
# the 823-byte synthetic. GRCh37 build (bare contigs).
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --output /tmp/allelix-review/hg00187_gatk.json \
  --report-format json
# Expected: exit 0; ~30-60 total annotations (low per-sample because
# GATK-HC raw output is rsID-less — hits come via position-based
# ClinVar + #8 rsID-less resolution, not the rsID fast-path). Build
# banner reads "GRCh37 (override; 0/0 known-SNP positions matched)"
# because GATK-HC raw output is ID=. and the override pins the build.

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
# the bundle large enough to exercise the truncation path.
allelix analyze test_data/real/vcf/thousandG_chr22.vcf.gz \
  --sample HG00096 --build grch37 \
  --output /tmp/allelix-review/thousandG_HG00096.json
```

**Expected:** Single-sample, gVCF, and the real-scale HG00187 run
all exit 0. Multi-sample without `--sample` raises MultiSampleError
in both cases; the 1000G file demonstrates the "... and N more"
truncation tail.

### 5i. rsID-less VCF resolution (GH #8 — the flagship feature)

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

### 5j. Build auto-detection: blind fallback + chr-prefix inference (GH #38)

Two paths in this step. Run both.

**(a) Blind-default warning** — fires when no rsIDs in input AND no
chr-prefixed contigs AND no `##contig assembly=` tag. The pipeline
has nothing to go on and recommends `--build` explicitly.

```bash
allelix analyze tests/fixtures/mock_vcf_no_signal.vcf 2>&1 | grep -iE "build:|auto-detect|--build"
```

Expected (yellow warning visible):

```
Build: GRCh37 (fallback (no known SNPs matched); 0/0 known-SNP positions matched)
Could not auto-detect genome build (no rsIDs in input, no ##contig
assembly tag, no chr-prefixed contigs). Defaulted to GRCh37. If the
file is the other build, pass --build grch37 or --build grch38 …
```

The fixture uses bare contigs (`ID=1`, `ID=22`, …) with **no**
`assembly=` tag and **no** rsIDs — all three auto-detect paths fail
and the GRCh37 fallback fires loudly. (The similarly-named
`mock_vcf_rsidless.vcf` carries `assembly=GRCh37` and is used
elsewhere to test header-based resolution; don't confuse the two.)

**(b) chr-prefix inference** — GH #38 path. Fires when no rsIDs and
no `##contig assembly=` tag BUT the contigs are `chr`-prefixed
(modern variant-caller convention for GRCh38).

```bash
allelix analyze tests/fixtures/mock_vcf_chr_prefix_grch38.vcf 2>&1 | grep -iE "build:|inferred|chr-prefix"
```

Expected (positive info message, no yellow warning):

```
Build: GRCh38 (inferred from chr-prefixed contig names; 0/0 known-SNP positions matched)
Inferred GRCh38 from chr-prefixed contig names (GRCh38 convention).
Pass --build grch37 if this file is UCSC hg19 with chr-prefixed
contigs instead.
```

The fixture (`mock_vcf_chr_prefix_grch38.vcf`) is built specifically
to exercise this path: chr-prefixed contigs, no `assembly=` tag, no
rsIDs. Both halves of the matrix verified.

**Why not the bundled GIAB benchmarks?** They carry
`##contig=<ID=chr1,...,assembly=human_GRCh38_no_alt_analysis_set.fasta>`
— the `assembly=` field wins via `header_build` before the
chr-prefix tertiary signal ever gets a vote. Same for the chr22
slice (it's a strict subset of the same VCF). The synthetic fixture
above is the cleanest demo.

**Pre-v2.0.0 baseline:** the pipeline silently defaulted to GRCh37
for both shapes, which would mis-annotate a GRCh38 file. v2.0.2 #38
closes the chr-prefix case; the blind-default path stays as the last
resort when no signal is available.

## 6. Cross-parser identity check

The user1190 genotype exists in 6 format representations. All should
produce identical annotation sets (same rsIDs, same significance, same
sources). The exact annotation count depends on database versions, but
the counts must match across formats.

```bash
mkdir -p /tmp/allelix-review
for f in \
  test_data/real/23andme/user1190_v5.txt \
  test_data/real/mhg/user1190.txt \
  test_data/real/livingdna/user1190.csv \
  test_data/real/myheritage/user1190.csv \
  test_data/transcoded/user1190_as_ancestrydna.txt \
  test_data/transcoded/user1190_as_ftdna.csv; do
  echo "=== $f ==="
  allelix analyze "$f" --exclude-snpedia --output /tmp/allelix-review/$(basename "$f").json 2>&1 | tail -3
done
```

Then compare annotation counts:

```bash
for f in /tmp/allelix-review/user1190_*.json; do
  echo "$(basename $f): $(python3 -c "import json; print(len(json.load(open('$f'))['annotations']))")"
done
```

**Expected:** All 6 files produce the same annotation count. Any
discrepancy is a parser or build-detection bug.

## 7. Wrong-allele safety — enrichment uses exact alt match (GH #18, #23, #42)

Verify that enrichment lookups use exact alt-allele matching, not
MAX-aggregated fallback and not complement-resolved coincidence. The
v2.0.1 ship closed three independent paths where a wrong-allele
number could attach to an annotation:

- **#18**: `resolve_strand` no longer falls back to complement
  matching at multi-allelic sites. CADD scores attached to alt-set
  annotations must come from a direct `(ref, alt)` match.
- **#23**: alt-less annotations (raw GWAS rows) no longer take the
  `MAX(af) GROUP BY rsid` fallback. They get enrichment only via the
  safe position-fallback path (rsID resolved on-the-fly via ClinVar,
  resolved tuple carries the user's specific alt). CADD enrichment
  now has the same position-fallback for symmetry.
- **#42**: ClinVar `CLNDN` list joined per record instead of
  index-paired with `CLNSIG`, eliminating Frankenstein pairings at
  multi-SCV variants.

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/enrichment_check.json
python3 -c "
import json
data = json.load(open('/tmp/allelix-review/enrichment_check.json'))
for a in data['annotations']:
    if a.get('am_pathogenicity') is not None and a.get('alt'):
        print(f\"{a['rsid']} alt={a['alt']} am={a['am_pathogenicity']:.3f} {a['am_class']}\")
" | head -20
```

**Expected:** AM scores correspond to the user's specific alt allele,
not the site-wide MAX. Spot-check a few rsIDs against the AlphaMissense
source data if available.

**#18 stronger invariant check.** For every alt-set annotation with a
stamped CADD score, the alt must appear directly in gnomAD's alts at
that rsID (no complement-resolved hits). Counted across the v2.0.1
HG002 gVCF battery: 578/578 direct, 0 via-complement.

**#42 ClinVar real-cache content gates.** After a `db update` against
the v2.2 per-SCV TSV loader. **These are ship-gates** in the same
sense §14 and §45's GWAS gates are: zero expected results mean
ship; any non-zero result means **do not tag**. The flagship
ClinVar loader gets the same gating posture as every other annotator.

**(a) Per-SCV row shape (Defect 4 — supersedes the pre-#42 join pin).**
The v2.2 loader stores one row per (variant, SCV submission), not one
row per variant with semicolons. Expect multiple rows per multi-SCV
rsID, each carrying a single SCV's condition exactly as the submitter
recorded it:

```bash
allelix-dev$ sqlite3 ~/.local/share/allelix/clinvar.GRCh38.sqlite \
  "SELECT rsid, COUNT(*) FROM clinvar_variants \
   WHERE rsid IN ('rs1800896', 'rs1063192') GROUP BY rsid;"
# Expected: ≥1 row per rsID. Multi-SCV variants will show >1 row
# (e.g. rs1063192 typically shows 2+). A single combined row with a
# semicolon-joined condition is a regression to the pre-#42 shape.
```

The pre-#42 Frankenstein-pairing safety question is preserved
structurally: each per-SCV row carries one submitter's
(significance, condition) pair, so there is no cross-product
pairing risk anymore. The semicolon-join in the old pin was the
loader's coping mechanism for the VCF's CLNSIG|CLNDN parallel
arrays; with per-SCV rows that mechanism is unnecessary.

**(b) Significance-sentinel guard (Defect 5 ship-gate).** ClinVar's
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
allelix-dev$ sqlite3 ~/.local/share/allelix/clinvar.GRCh38.sqlite \
  "SELECT clinical_significance, COUNT(*) FROM clinvar_variants \
   WHERE LOWER(clinical_significance) IN \
     ('-', '', 'not specified', 'no classification provided', \
      'other', 'association', 'association not found') \
   GROUP BY clinical_significance;"
# Expected: 0 rows. v2.2.1 (#116) added the three non-classification
# curatorial terms — "other", "association", "association not found".
# Keep this list in sync with _CLINVAR_PLACEHOLDER_CLNSIGS in
# databases/manager.py.

# Spot-check a known-affected rsID:
allelix-dev$ sqlite3 ~/.local/share/allelix/clinvar.GRCh38.sqlite \
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
tag." That's R-4 (CLNSIG vocabulary drift —
`tests/databases/test_clinvar_clnsig_drift.py` +
`allelix/data/clinvar_clnsig_snapshot.yaml`). When R-4 evolves
from snapshot-equality to "every value ∈ {known magnitude vocab,
placeholder set}", point it at the same caches this gate scans;
this §7b gate then graduates from "denylist confirmation" to a
spot-check of what R-4 enforces structurally.

### 7d. Strand-aware carrier matching (GH #14 strand-half, R-1 / #50, ADR-0035 PR 4)

v2.1.0 introduces strand-aware carrier matching: complement-strand
reads of the same biology now produce identical annotation sets to
forward-strand reads, while preserving the v2.0.1 #18 wrong-allele
multi-allelic safety invariant via the `Variant.ref` orientation
guard.

The CLAUDE.md R-1 mandatory test case (rs5742904 APOB forward C/C
vs coding G/G) isn't in the bundled test fixtures. The unit-test
suite pins the contract on rs1801133 (mock ClinVar: REF=G ALT=A
Pathogenic MTHFR) as the protocol-level substitute — same shape,
exercised against a known fixture. Re-run the dedicated test class
as a fast in-protocol sanity check:

```bash
.venv/bin/python -m pytest \
  tests/annotators/test_clinvar.py::TestStrandAwareCarrierMatch -v
```

**Expected:** three tests pass:

- `test_hom_ref_both_strands_no_annotation` — hom-ref forward G/G
  and coding C/C both produce zero ClinVar annotations (CLAUDE.md
  R-1 mandatory shape).
- `test_het_carrier_both_strands_same_annotation_set` — het carrier
  forward G/A and coding-strand C/T produce the **same single**
  ClinVar annotation (significance / condition / gene match across
  the two reads).
- `test_audit_18_multi_allelic_safety_preserved` — when
  `Variant.ref` disagrees with both the source ref and its
  complement, strand-flip abstains and direct-match-only applies.

For real-data validation against the Section 19 battery, see
`tests/utils/test_allele.py::TestStrandAwareCarrierMatch` (helper-
level, covers palindromic / orientation / multi-allelic branches)
and the wrong-allele safety invariant check at the end of
Section 19 — v2.1.0 produces a larger denominator (186/186 GRCh37,
173/173 GRCh38 vs v2.0.2's 107/107 + 113/113) because PR 2's per-
Annotation alt threading made more rows eligible for the exact-
`(rsid, alt)` CADD lookup. **All matches must remain 100% allele-
direct, 0 via-complement** — that's the v2.0.1 #18 invariant
holding against the broader denominator.

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
- Pop. Freq column shows gnomAD frequencies
- AM column shows AlphaMissense scores
- ClinPGx AM scores show dimmed caveat indicator
- Row borders are color-coded (red = pathogenic, green = benign)
- Zygosity column shows Heterozygous / Homozygous for each row
- CADD column scores are color-coded (red ≥30, orange ≥20, gray <20) with tooltips
- "Reading This Report" section is present
- Regulatory notice is present

### 8b. Terminal report (bare-min layout, GH #9)

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
"
```

**Expected:** Schema version 5 (v2.1.0; ADR-0035 Cluster B added
`variant.ref` internally and `annotation.trait` / `.p_value` /
`.phecode` to the JSON output). Multiple sources present. gnomAD and
AM enrichment counts > 0. On any file with GWAS hits, the structured
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

**Expected:** Per-chromosome concordance table. Coverage overlap stats.
High concordance expected (same biology, different format).

## 10. Config system

```bash
allelix config show
allelix config set license.commercial true
allelix config show
allelix analyze test_data/real/23andme/user1190_v5.txt 2>&1 | grep -i "snpedia\|skipping"
allelix config set license.commercial false
allelix config show
```

**Expected:** With `license.commercial = true`, SNPedia is excluded
from analysis automatically. After setting back to `false`, SNPedia
is included again.

## 11. CADD opt-in flow

CADD is disabled by default. Verify the opt-in path:

```bash
allelix config show | grep cadd
# Expected: sources.cadd = false

allelix db update --cadd
# Expected: CADD license confirmation prompt. Accept to download CADD cache.

allelix db status | grep -i cadd
# Expected: CADD shows "yes" in Ready column with version "v1.7"

allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/cadd_check.json
python3 -c "
import json
data = json.load(open('/tmp/allelix-review/cadd_check.json'))
has_cadd = sum(1 for a in data['annotations'] if a.get('cadd_phred') is not None)
print(f'With CADD score: {has_cadd}')
"
# Expected: Non-zero CADD enrichment count.
```

Verify commercial mode gates CADD:

```bash
allelix config set license.commercial true
allelix config show | grep -E "cadd|commercial"
# Expected: CADD is excluded when commercial mode is active (commercial_ok=False)
allelix config set license.commercial false
```

If full mode is available (pysam installed + local tabix file):

```bash
allelix config set options.cadd_full true
allelix analyze test_data/real/23andme/user1190_v5.txt 2>&1 | grep -i "cadd\|grch38"
# Expected: GRCh38-only guard — if input is not GRCh38, warning about skipping CADD full mode
allelix config set options.cadd_full false
```

## 12. Diff command

```bash
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/baseline.json
allelix analyze test_data/real/23andme/user1190_v5.txt --output /tmp/allelix-review/current.json --diff /tmp/allelix-review/baseline.json
```

**Expected:** Diff reports no changes (same input, same databases).

## 13. Database update signals

```bash
allelix db update
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
allelix db update --force
```

**Expected:** All annotators re-download and show green checkmarks.
Note: `--force` semantics differ by tier. Server-driven sources
override a "signal matches" skip; code-driven sources have no
signal-match path to override — `--force` is the only way to
re-trigger their download because pinned URLs are deterministic.
See ADR-0030.

## 14. GWAS Catalog real-data sanity (slow tests, auto-fetch GH #45)

These tests load the real GWAS Catalog (auto-fetched from EBI on
first run as of v2.0.2 — no manual download step) and verify that
the magnitude scoring formula produces bounded output.

```bash
python -m pytest tests/test_end_to_end.py -k "TestRealDataGwasSanity" -v
```

**Expected:** 2 tests **pass** (not skip). Default floor (9.0) keeps
output under 50 rows. Old floor (7.0) produces more output than new
floor. First run auto-fetches `test_data/gwas_catalog.zip` (~65 MB,
~30 s on a decent connection); subsequent runs are offline-fast.

A silent skip here is a ship-gate defect (GH #45 policy): either the
auto-fetch failed and was logged as `OSError`, or the test was
deselected. Investigate, do not tag.

## 15. Edge case files

```bash
# Build mismatch detection (analyze runs the build-detection pipeline; stats does not)
allelix analyze test_data/edge_cases/mhg_grch38_with_grch37_header.txt 2>&1 | grep -i "mismatch\|build"
# Expected: Build mismatch warning (header claims GRCh37, positions are GRCh38)

# P-A: canonical header tightening — this file should NOT be recognized as 23andMe
allelix stats test_data/edge_cases/23andme_lookalike_rejected_by_PA.txt 2>&1
# Expected: "No parser recognized" error

# Genes for Good — 23andMe-format export from a different service
allelix stats test_data/edge_cases/23andme_format_from_genes_for_good_service.txt
# Expected: Recognized as 23andMe format, stats displayed

# GRCh36 FTDNA file (analyze detects build from positions; stats shows parser default)
allelix analyze test_data/edge_cases/ftdna_grch36_positions.csv 2>&1 | grep -i "grch36\|build"
# Expected: GRCh36 detected. ClinVar skipped (no GRCh36 cache).

# Unsupported formats
allelix stats test_data/edge_cases/unsupported_decodeme.txt 2>&1
# Expected: "No parser recognized" (decodeme is genuinely unsupported)

# Previously-unsupported, now-supported: the 23andMe exome VCF is a valid
# VCFv4.1 file and is recognized by the VCF parser added after this fixture
# was named. Pin the regression — if "No parser recognized" comes back,
# VCF detection broke.
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

### 16b. gnomAD ref/alt resolution

```bash
allelix export plink test_data/real/23andme/user1190_v5.txt -o /tmp/allelix-review/user1190_gnomad --build grch37
python3 -c "
lines = open('/tmp/allelix-review/user1190_gnomad.bim').readlines()
with_alt = sum(1 for l in lines if l.strip().split('\t')[5] != '0')
mono = sum(1 for l in lines if l.strip().split('\t')[5] == '0')
print(f'With ref/alt: {with_alt}')
print(f'Monomorphic (A2=0): {mono}')
"
```

**Expected:** Majority of variants have ref/alt resolved (A2 != 0)
when gnomAD is available. Monomorphic count matches CLI output.

### 16c. Roundtrip with plink2 (optional)

If `plink2` is installed:

```bash
plink2 --bfile /tmp/allelix-review/user1190 --freq --out /tmp/allelix-review/freq_check
```

**Expected:** plink2 reads the files without error. Frequency report
produced.

## 17. Cleanup

```bash
rm -rf /tmp/allelix-review
```

Optionally remove downloaded databases to free ~16 GB (≈22 GB if you
also completed the CADD opt-in in Step 11):

```bash
rm -rf ~/.local/share/allelix/
```

(Default install is ~16 GB on disk: alphamissense ~7.8 GB +
gnomad ~6.1 GB + clinvar ~1.0 GB + gwas_catalog ~0.7 GB +
snpedia/gwas/pharmgkb ~0.6 GB. CADD adds ~5.8 GB on top.)

## 18. Upgrade-path verification (v2.0.1+ caches, GH #22 / #42)

When a user upgrades from an older minor (v2.0.0 → v2.0.2, or
v2.0.2 → v2.1.0), the next `allelix db update` should auto-invalidate
exactly the caches whose interpreter or schema version bumped at any
step in the upgrade window, and leave the rest signal-skipped.

```bash
# Pre-condition: a v2.0.0 cache already on disk
allelix db status
# Expected: ClinVar = "no" (CLINVAR_INTERPRETER_VERSION bumped 1 → 2 in v2.0.1)
#   if a v2.0.0 ClinVar cache predated the bump

allelix db update
# Expected:
#   clinvar: downloading…  ✓ clinvar ready
#   pharmgkb / gwas / snpedia / gnomad / alphamissense / cadd: already current
```

**Expected:** only ClinVar re-downloads (because v2.0.1 bumped
`CLINVAR_INTERPRETER_VERSION` to invalidate v2.0.0 caches so the
post-#42 loader runs). Other annotators are signal-matched and
skipped — no spurious re-downloads of the multi-GB gnomAD
(~2.7 GB compressed, ~6 GB on disk) or AlphaMissense (~1.8 GB
compressed, ~8 GB on disk) prebuilts.

**v2.0.2 specifics:** no new interpreter / schema bumps. The HF URL
move (#37, v2.0.1) is invisible — HF redirects from `dial481/...` to
`allelix/...`; existing pinned installs continue to fetch.

**v2.1.0 specifics:** no new interpreter / schema bumps in the cache
layer. `CLINVAR_INTERPRETER_VERSION` stays at 2 (v2.0.1 bump);
PharmGKB / gnomAD / AlphaMissense / CADD interpreter / schema
versions all unchanged. The ADR-0035 Cluster B work changed query
semantics, not loader output — `pharmgkb.sqlite` / `clinvar.*.sqlite`
content is byte-identical. JSON report `schema_version` bumped 4 → 5
at the cluster's first PR; the new fields (`annotation.trait`,
`annotation.p_value`, `annotation.phecode`) are additive and v4
consumers reading v5 output continue to succeed. Existing v2.0.2
caches stay current on upgrade — `allelix db update` will signal-skip
all annotators.

## 19. Gold-standard real-data VCF battery (~1.5 GB, optional but recommended)

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
# Build: GRCh38 (override). See "Floor invariants" and
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
# total annotation count dwarfs the GRCh38 run (an order of
# magnitude higher) because the upstream NIST GRCh37 benchmark
# VCF carries rsIDs in the ID column (~95% of rows) whereas the
# GRCh38 one has none — so rsID-keyed GWAS / SNPedia / ClinPGx
# lookups fire on nearly every row here.

# 3. HG00187 GATK-HC gVCF — real caller at WGS scale (GRCh37 build)
allelix analyze test_data/real/vcf/HG00187_gatkhc.g.vcf.gz \
  --build grch37 \
  --min-magnitude 0 \
  --output /tmp/allelix-review/hg00187_gatk.json \
  --report-format json
# This file is GATK-HC raw output (no rsIDs in ID column), so all
# hits come via position-based ClinVar + rsID-less resolution (#8)
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
  --sample HG00096 --build grch37 \
  --output /tmp/allelix-review/thousandG_HG00096.json
# Expected: chr22-only annotations for the chosen sample.

# 5. Synthetic edge-case fixtures (parser exercise, no DB lookups needed)
allelix analyze test_data/real/vcf/mock_vcf.vcf --build grch37 \
  --output /tmp/allelix-review/mock_vcf.json
allelix analyze test_data/real/vcf/mock_gvcf.g.vcf --build grch37 \
  --output /tmp/allelix-review/mock_gvcf.json
allelix analyze test_data/real/vcf/mock_vcf_rsidless.vcf --build grch37 \
  --output /tmp/allelix-review/mock_rsidless.json
```

### Ground-truth invariants — published HG002 spot checks + floors

Exact annotation counts are NOT pinned. Pinning them produced
stale-pin regressions three times (§7, §15, §19 historically) and the
last re-pin would have masked an over-filter regression of similar
shape to the one tagged this section's calibration. Instead the
battery's pass/fail is decided by **`test_data/check_ground_truth.py`**
against **`test_data/HG002_GROUND_TRUTH.yaml`** — the harness asserts:

1. **Floor invariants per file** (`floors:` in the YAML) — written /
   total / unique-key counts must clear floors sized at ~80% of
   v2.2.1 (de0de64, 2026-06-18) observed values. A drop below floor
   is an over-filter or vocab-leak regression: classify the removed
   annotations with a v2.2.0→v2.2.x baseline `--diff` run before
   ratcheting the floor down. **"Plausible drift" is not a license
   to re-pin.**

2. **Spot-check invariants per rsID** (`analyze:` in the YAML) —
   concrete published HG002 ground truth at MTHFR rs1801133, APOE
   rs7412 / rs429358, VEGFA rs2010963, Factor V rs6025: row counts,
   required sources (`{clinvar, pharmgkb}` for pharmacogenomic
   rsIDs), significance-pattern conformance, condition non-empty.
   The GIAB benchmark file IS the source of truth; the harness just
   verifies allelix reflects it faithfully. Unfalsifiable by drift.

3. **Vocabulary invariants** (`vocabulary:` in the YAML) — every
   emitted `significance` value must be in
   `clinvar_allowed ∪ pharmgkb_allowed`, with explicit forbidden
   placeholders the loader must filter. A novel value means the
   `_CLINVAR_PLACEHOLDER_CLNSIGS` filter regressed (a #116 / #101
   leak). The denylist form here is reactive; R-4
   (`tests/databases/test_clinvar_clnsig_drift.py`) is the durable
   structural form.

4. **Universal invariants** (`universal:` in the YAML) — every
   annotation has the required non-null fields, `source` ∈
   {clinvar, pharmgkb, gwas, snpedia}, `attribution` matches the
   `source`-to-display-name map, `magnitude` ∈ [0, 10].

Run the harness against each of the three battery outputs:

```bash
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/giab_grch38.json giab_grch38_benchmark
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/giab_grch37.json giab_grch37_benchmark
python3 test_data/check_ground_truth.py \
  /tmp/allelix-review/hg00187_gatk.json hg00187_gatkhc_gvcf
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
in place):** every annotation in the GIAB benchmark output must also
appear in the DeepVariant gVCF output at default filter. The gVCF
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
bench = json.load(open('/tmp/allelix-review/giab_grch38.json'))['annotations']
gvcf = json.load(open('/tmp/allelix-review/hg002_gvcf.json'))['annotations']
key = lambda a: (a['source'], a['rsid'], a.get('condition',''), a.get('description',''))
bk = {key(a) for a in bench}
gk = {key(a) for a in gvcf}
missing = bk - gk
print(f'benchmark rows: {len(bench)}  gvcf rows: {len(gvcf)}')
print(f'missing in gvcf at production filter: {len(missing)}')
assert not missing, f'SUPERSET BROKEN — {len(missing)} benchmark rows absent from gvcf'
print('SUPERSET HOLDS')
"
```

At `--min-magnitude 0` the superset relaxes slightly: a handful
(~5) of low-confidence variants (`clinvar_uncertain_significance`,
`clinvar_conflicting_classifications`) where DeepVariant disagrees
with the GIAB truth (calls `0/0` or `./.`) appear in the benchmark
but not the gVCF. Correct behavior — allelix is faithfully
reflecting each caller's view; the variants simply weren't called
by DeepVariant at those positions.

**Wrong-allele safety invariants (GH #18 / #23 / #42):** every CADD
score stamped on an alt-set annotation must come from a direct
`(ref, alt)` match in gnomAD (no complement-resolved hits). Per the
v2.0.1 verification on the HG002 gVCF: 578/578 (100%) of attached
CADD scores were allele-direct.

```bash
# Runs against whichever real VCF output is present.
python3 -c "
import json, os, sqlite3, contextlib, sys
for name in ['hg002_gvcf', 'hg00187_gatk', 'giab_grch37', 'giab_grch38']:
    p = f'/tmp/allelix-review/{name}.json'
    if not os.path.exists(p): continue
    rows_all = json.load(open(p))['annotations']
    direct = no_direct = 0
    with contextlib.closing(sqlite3.connect(os.path.expanduser('~/.local/share/allelix/gnomad.sqlite'))) as conn:
        for a in rows_all:
            if not a.get('alt') or a.get('cadd_phred') is None: continue
            rows = conn.execute('SELECT alt FROM gnomad_frequencies WHERE rsid=?', (a['rsid'],)).fetchall()
            if not rows: continue
            if a['alt'] in {r[0] for r in rows}: direct += 1
            else: no_direct += 1
    print(f'{name}: allele-direct CADD={direct}, via-complement={no_direct}')
    assert no_direct == 0, f'BUG in {name} — {no_direct} via-complement CADD scores'
"
```

**ClinVar per-SCV row shape (#42 v2.2):** multi-SCV variants surface
one row per submission, not a single semicolon-joined row. The
Frankenstein-pair safety question that the original audit asked
about is now structurally impossible — each row carries one
submitter's (significance, condition) pair from
`submission_summary.txt.gz`, so no cross-product pairing exists.
The legacy semicolon-joined pin from the VCF era is removed.
See §7 for the post-update real-cache ship-gate (per-SCV row shape
+ significance-sentinel scan).

### When to run this battery

- **Before every ship-gate tag.** v2.0.1's load-bearing chain-of-trust check.
- **After any change to** `_pipeline.py`, the VCF parser, `iter_clinvar_records`, or strand / enrichment logic.
- **Not for routine dev work** — unit tests cover the same logic against synthetic fixtures.

## Pass criteria

All of the following must be true:

- [ ] Unit test suite: **1,576 passed, 0 skipped** (v2.1.0 floor with
      `plink2` installed and GWAS auto-fetch succeeding). Skips are a
      ship-gate defect — investigate, don't ignore.
- [ ] Ruff lint + format: zero warnings (`ruff check .` and `ruff format --check .` both clean)
- [ ] `db update` downloads all enabled annotators without errors
- [ ] `db status` shows all annotators ready with version and record count; ClinPGx row labeled "ClinPGx" (not "PharmGKB")
- [ ] All 8 parser formats produce successful analysis (23andMe, AncestryDNA, FTDNA, MyHeritage, MHG / Tempus, Living DNA, FTDNA Illumina raw, VCF/gVCF)
- [ ] Cross-parser identity: same annotation count across all user1190 representations (6 array-based formats; VCF doesn't have a user1190 representation)
- [ ] VCF flagship feature: rsID-less VCF produces non-zero annotations via ClinPGx resolution (step 5i); multi-sample VCF requires `--sample` (step 5h)
- [ ] Build auto-detection: warning fires when no rsID + no header signal (step 5j); **chr-prefix contigs auto-infer GRCh38 without `--build` (GH #38)**
- [ ] **Wrong-allele safety invariants (GH #18, #23, #42)**: every alt-set CADD score has its alt directly in gnomAD's alts at that rsID; alt-less (raw GWAS) rows only get enrichment via the safe position-fallback path; ClinVar emits one row per (variant, SCV submission) — the per-SCV pairing structurally retires the Frankenstein-pair hazard the original audit raised (v2.2 per-SCV TSV loader)
- [ ] **ClinVar significance-sentinel ship-gate (Defect 5 / §7b)**: post-`db update` sentinel scan returns 0 rows on both clinvar.GRCh37.sqlite and clinvar.GRCh38.sqlite — `LOWER(clinical_significance) IN ('-', '', 'not specified', 'no classification provided', 'other', 'association', 'association not found')` (case-insensitive matches loader's `_normalize_clnsig`) — and rs137854557 carries a real classification. Set MUST match `_CLINVAR_PLACEHOLDER_CLNSIGS` in `databases/manager.py` — never broaden one without broadening the other (`not provided` is intentionally excluded; it maps to 2.0 in `_CLNSIG_MAGNITUDE`, safe below the floor). v2.2.1 (#116) added the three non-classification curatorial terms. Denylist form is reactive; durable allowlist form lives in R-4 (`test_clinvar_clnsig_drift.py`)
- [ ] **Terminal report (GH #9)**: bare-min columns only (`rsID | Gene? | Source | Significance | Mag | GT | Condition?`); Review Status / Zygosity / Freq / AM / CADD intentionally absent (still present in HTML/JSON)
- [ ] HTML report renders correctly in a browser; "Annotators:" subtitle uses display names (ClinPGx, not pharmgkb); enrichment columns (Review Status, Zygosity, Freq, AM, CADD) all present
- [ ] JSON report has schema version 6 with gnomAD + AM + CADD enrichment; structured GWAS fields populated (`annotation.trait` non-empty on GWAS rows, `p_value` parses as float when published, `phecode` populated when the upstream trait carried one); `panel_coverage` object appears when `--filter-file` is used (GH #75); `license_attributions[].source` shows "ClinPGx" with `source_url` `https://www.clinpgx.org`
- [ ] **Strand-aware carrier matching (GH #14 strand-half / R-1 / ADR-0035 PR 4)**: forward and reverse-strand reads of the same biological het carrier produce identical ClinVar annotation sets; hom-ref on either strand produces zero annotations; multi-allelic safety (variant.ref disagreement) abstains rather than guesses (sanity check at step 7d below)
- [ ] Config system correctly gates SNPedia on `license.commercial`
- [ ] CADD opt-in: `--cadd` downloads cache, license prompt shown, scores enriched
- [ ] CADD commercial gate: `license.commercial = true` excludes CADD
- [ ] Edge case files produce expected behavior
- [ ] `db update` (second run) skips already-current databases
- [ ] **Upgrade path (GH #22 / #42)**: upgrading from v2.0.0 to v2.0.2 against an existing cache re-downloads only ClinVar (interpreter-version invalidation from v2.0.1); upgrading from v2.0.2 to v2.1.0 against an existing cache signal-skips all annotators (no new interpreter / schema bumps in the ADR-0035 Cluster B work)
- [ ] GWAS Catalog slow tests pass (auto-fetch the fixture — silent skip is forbidden, GH #45)
- [ ] `methylation`, `pharmacogenomics`, `compare` subcommands produce output (pharmacogenomics --help says "ClinPGx-style sources")
- [ ] PLINK export produces valid .bed/.bim/.fam with correct magic and alignment
- [ ] PLINK export resolves ref/alt from gnomAD when available
- [ ] `allelix --version` reports the actual pyproject version even from a bare source checkout (GH #34)
- [ ] **Gold-standard VCF battery (step 19)**: GIAB benchmark + DeepVariant gVCF + GATK-HC GRCh37 all analyze cleanly; `check_ground_truth.py` exits 0 on all three battery files (`giab_grch38_benchmark`, `giab_grch37_benchmark`, `hg00187_gatkhc_gvcf`) — enforces floor counts, published HG002 spot-checks at rs1801133 / rs7412 / rs2010963 / rs6025, the significance vocabulary union, and per-row universal invariants; gVCF is a strict superset of the benchmark at production filter (0 missing); allele-direct CADD invariant holds (0 via-complement); ClinVar emits per-SCV rows (multiple rows per multi-SCV rsID, each with its submitter's exact condition) — semicolon-joined rows are a regression to the pre-#42 shape

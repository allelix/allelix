# Allelix Test Data Library

> **Note:** The `real/` and `transcoded/` directories are not committed
> to the repository. They are hosted as GitHub release assets. Run
> `scripts/fetch_testdata.sh` to download them.

Curated test data and frozen regression baselines for parser and annotator
verification.

## What's here

- **`real/`** — CC0 public-domain genotype files from openSNP. Real human
  genomes voluntarily released under [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/).
- **`transcoded/`** — Structural transcodes of `real/23andme/user1190_v5.txt`
  into formats that have no natural user1190 representation in the openSNP
  archive. Same biology, different format envelopes. Provenance documented
  in `transcoded/README.md`.
- **`baselines/`** — Frozen `allelix analyze` output keyed per release.
  Diff against the previous release to detect regressions. Generated with
  `--exclude-snpedia` so the JSON contains only CC0-compatible content.
- **`edge_cases/`** — Files that exercise specific parser quirks (P-A
  canonical-header tightening, P-B GRCh36 position detection, MHG/Tempus
  build-mismatch warning, unsupported formats). Documented in
  `edge_cases/README.md`.
- **`databases/`** — Pinned annotator databases. **Gitignored** (~15 GB).
  Restore via `allelix db update --data-dir test_data/databases`.

## Directory layout

```
test_data/
├── README.md                  # this file
├── real/                      # genotype files from CC0 openSNP donations
│   ├── mhg/                   # 1 file (user1190 transcoded — see transcoded/README.md)
│   ├── 23andme/               # 6 files (user1190 v5, user1 v1 build-36, +4 extras)
│   ├── ancestrydna/           # 7 files
│   ├── ftdna/                 # 7 files
│   ├── livingdna/             # 1 file (user1190 transcoded)
│   └── myheritage/            # 1 file (user1190 transcoded)
├── transcoded/                # user1190 → AncestryDNA / FTDNA format
├── databases/                 # pinned annotator databases — GITIGNORED
├── baselines/                 # frozen analyze outputs per release
│   ├── v1.0.0/
│   └── v1.1.0/
└── edge_cases/                # files that exercise specific parser quirks
```

The user1190 files in `real/mhg/`, `real/livingdna/`, and `real/myheritage/`
are all structurally transcoded from `real/23andme/user1190_v5.txt`. They
share the same rsIDs, positions, and genotypes — only the format envelope
differs. This is what makes the cross-parser identity test work across 6
formats.

## Licensing

| Source | License | Notes |
|---|---|---|
| openSNP genotype files (real/, transcoded/, edge_cases/, plus the synthetic source for baselines) | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) | Public domain. |
| Corpas family exome VCF (edge_cases/unsupported_23andme_exome_vcf.txt) | CC-BY / CC0 (see [Corpas et al. 2015](https://doi.org/10.1186/s12864-015-1973-7)) | Deliberately public-released by subjects. Attributed to satisfy CC-BY. |
| ClinVar (in databases/) | NCBI public domain | Restore via `allelix db update`. |
| GWAS Catalog (in databases/) | NCBI/EBI public domain | Same. |
| PharmGKB (in databases/) | [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) | Attribution required if reports are redistributed. |
| gnomAD (in databases/) | [ODbL](https://opendatacommons.org/licenses/odbl/1-0/) | Attribution required. |
| AlphaMissense (in databases/) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Attribution required. |
| SNPedia (in databases/) | [CC BY-NC-SA 3.0 US](https://creativecommons.org/licenses/by-nc-sa/3.0/us/) | Non-commercial only. Baselines here are generated with `--exclude-snpedia` so no SNPedia content reaches `baselines/`. |

## Pinned database versions (for regression testing)

```
clinvar.GRCh37.sqlite   2026-05-30
clinvar.GRCh38.sqlite   2026-05-30
pharmgkb.sqlite         2026-06-05
gwas.sqlite             2026-06-05
gnomad.sqlite           2026-06-05
alphamissense.sqlite    2026-06-05
snpedia.sqlite          pre-built cache from HuggingFace (source data scraped 2026-05-20)
```

Pinned databases are the constant for regression testing. When a new release
is verified, baselines are generated against THESE databases, not whatever
`allelix db update` pulls today. Database freshness is a separate concern
from code regression. Restore the pinned set via `allelix db update
--data-dir test_data/databases` then drop in the archived files.

## Usage

### Regenerate baselines for a new release

```bash
cd <repo-root>
TD=test_data
DB=$TD/databases
B=$TD/baselines/<new-version>
mkdir -p $B

# --exclude-snpedia keeps the baselines license-clean (SNPedia is non-commercial)
allelix analyze $TD/real/23andme/user1190_v5.txt  --data-dir $DB --exclude-snpedia --output $B/user1190_23andme_analyze.json
allelix analyze $TD/real/livingdna/user1190.csv   --data-dir $DB --exclude-snpedia --output $B/user1190_livingdna_analyze.json
allelix analyze $TD/real/myheritage/user1190.csv  --data-dir $DB --exclude-snpedia --output $B/user1190_myheritage_analyze.json
allelix analyze $TD/real/mhg/user1190.txt         --data-dir $DB --exclude-snpedia --output $B/user1190_mhg_analyze.json
allelix stats   $TD/real/mhg/user1190.txt                                                          > $B/user1190_mhg_stats.txt
```

Diff each JSON against the previous release. Any pathogenic ClinVar hit that
disappears or changes is a regression. New hits get spot-checked. **DO NOT
hardcode rsIDs as a regression set** — the baseline JSON files ARE the
regression set. Diff the full output, not a cherry-picked list.

### Cross-parser identity

user1190 exists in 6 representations across `real/` and `transcoded/`. All
six must produce identical annotation sets with the pinned databases
(currently 358 annotations each with `--exclude-snpedia`).

### Edge cases

- `edge_cases/mhg_grch38_with_grch37_header.txt` — synthetic MHG (Tempus)
  fixture with GRCh38 positions but a header claiming `build 37.1`. Exercises
  ADR-0021's auto-detect + mismatch warning path.
- `edge_cases/23andme_lookalike_rejected_by_PA.txt` — non-23andMe file that
  looks 23andMe-shaped. The v1.1.0 P-A canonical-header tightening correctly
  rejects it.
- `edge_cases/23andme_format_from_genes_for_good_service.txt` — Genes for
  Good service exports in 23andMe format. Auto-detect correctly identifies
  it as 23andMe.
- `edge_cases/ftdna_grch36_positions.csv` — FTDNA file with GRCh36
  coordinates. The FTDNA parser has no in-band build signal, so the file is
  silently labeled GRCh37 (tracked as roadmap item R-12).
- `edge_cases/unsupported_decodeme.txt` and `unsupported_23andme_exome_vcf.txt`
  — files in formats Allelix does not support. Auto-detect should return
  "no parser recognized" cleanly.

## Adding new files

- **New parser format** → add representative real files under
  `real/<format>/`. Document provenance in this README.
- **New regression edge case** → add to `edge_cases/` with an entry in
  `edge_cases/README.md` describing what it tests and the expected behavior.
- **New annotator** → add or refresh the `databases/` content (still
  gitignored), then regenerate baselines under a new release subdirectory.

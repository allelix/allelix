# ADR-0021: Genome build is detected from position data, not file headers

- **Date:** 2026-05-14
- **Status:** Accepted

## Context

MyHappyGenes/Tempus exports have been observed to declare "human reference build 37.1 coordinates" in their header while shipping positions on GRCh38. Cross-checking against NCBI's authoritative SNP positions on a small probe set spanning multiple chromosomes, and independently verifying via UCSC's GRCh38→GRCh37 chain, confirms the mislabel pattern is unambiguous: the header says GRCh37; the data is GRCh38.

The categorical risk: a false-positive pathogenic call on NIPA1 `rs104894490`. ClinVar's GRCh37 VCF records `REF=C ALT=G` at this position because chromosome 15's strand is inverted between GRCh37 and GRCh38. An MHG-format file at this position calls the genotype `G/G`. On GRCh38 forward strand (the actual build of the data) `G/G` is homozygous reference. On GRCh37 forward strand `G/G` looks like homozygous variant. The carrier rule (ADR-0007) works correctly — but it receives the wrong build's REF/ALT, and a healthy carrier of the reference allele looks pathogenic.

For ~99.6% of the genome the strand orientation is identical between GRCh37 and GRCh38 and a cross-build comparison silently agrees with the truth. For the ~0.4% of the genome where the strand was inverted between assemblies, REF and ALT swap, and any carrier-rule check produces an inverted result.

The deeper failure was Allelix trusted the file header. The first lesson is that providers mislabel headers. The second is that there's nothing intrinsic about the header that has to be trusted — coordinate differences between builds are typically tens of thousands to millions of bases. A handful of known-position SNPs is enough to identify the build unambiguously from the data alone.

## Decision

**Allelix detects the genome build of an input file from position data, not from the file header. When the header and the detected build disagree, Allelix uses the detected build and warns the user.**

### The detection algorithm

A hardcoded table in `allelix/utils/build_detect.py` maps a handful of well-known SNPs to their authoritative `(chromosome, position)` pairs in both GRCh37 and GRCh38. At parse time, the analyze pipeline scans variants until it finds a small number of rsIDs in the table and tallies which build's positions they match. The build with all (or all but a noise tolerance) of the matches wins.

Coordinates differ between GRCh37 and GRCh38 by tens of thousands to millions of bases for nearly every SNP on the autosomes — there is no overlap risk. A single matched rsID is unambiguous; multiple are confirmatory.

If the input contains no known SNPs from the table, detection returns `None` and Allelix falls back to the header (with a warning) or refuses to analyze (depending on a CLI flag). The default seed table includes ~10 SNPs spread across chromosomes 1, 10, 11, 12, 17, 19, 22 — common-array SNPs that virtually every consumer file carries.

### CLI surface

- `allelix analyze <file>` runs detection by default and prints a `[build]` line to stderr.
- `allelix analyze <file> --build grch37` forces the build. No detection, no warning.
- `allelix analyze <file> --build grch38` forces the build. No detection, no warning.
- `allelix db update` downloads BOTH builds' ClinVar VCFs by default. `--build grch37` or `--build grch38` flags exist to save bandwidth.

### Why not trust the header

Headers are wrong. MHG/Tempus exports confirmed it. Other providers may also be wrong. Even when right, headers conflate ambiguous reference identifiers — "build 37.1" could mean GRCh37 or NCBI 36, and there are old files in the wild that mismatch. Position data, by contrast, is a structured signal that's been authoritatively published in NCBI's reference for every common SNP.

This is the same principle as ADR-0016 (Data Classification Principle) applied to provenance: structured signals beat human-written labels.

### Why warn on mismatch

A header/data mismatch is itself diagnostic information. If a user's provider claims GRCh37 but ships GRCh38 (the observed MHG/Tempus pattern), that's a real-world data-quality issue the user may want to escalate to the provider. We surface it; we don't paper over it.

## Consequences

- **`Variant.build` is set per-file at analyze time, not per-parser.** Parsers may read the claimed build from the header (where one exists) but the analyze pipeline overrides with the detected build before annotators see the variant.
- **ClinVar annotator holds per-build SQLite caches** (`clinvar.GRCh37.sqlite` and `clinvar.GRCh38.sqlite`). `is_ready()` checks the cache for the build a variant carries. Missing build → a "build X cache not present; run db update" message.
- **`db update` defaults to downloading both builds.** Combined size ~200 MB compressed; user can opt out via `--build`. Once-per-month-ish refresh, acceptable.
- **The mock data generator supports both builds.** Three fixtures ship: clean GRCh37, clean GRCh38, GRCh38-with-GRCh37-header (replicating the MHG/Tempus mislabel pattern). Detection logic is tested against all three.
- **PharmGKB / CPIC are unaffected.** That pipeline is rsID + structured allele function — position-agnostic. No build dispatch needed.
- **Parser documentation tracks "header claims X" separately from "data is Y."** The observed MHG/Tempus mislabel pattern is noted wherever the format is referenced.
- **Existing v0.7.x ClinVar annotations are suspect** for any variant in a strand-inverted region between builds, when run against MHG/Tempus exports. After this ADR lands, users should regenerate with the auto-detected GRCh38 VCF and compare.
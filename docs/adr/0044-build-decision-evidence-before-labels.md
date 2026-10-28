# ADR-0044: The genome build is decided from the file's evidence before annotation

**Status:** Accepted
**Date:** 2026-10-12
**Supersedes:** ADR-0021's header fallback and its "scan until a few seeds
are found" detection; the chr-prefixed-contig inference added in v2.0.2.
**Keeps:** ADR-0037's confidence threshold, ADR-0043's position votes for
rsID-less records, and ADR-0025's GRCh36 fail-safe (as an explicit rule).

## Context

Up to 2.4.0, the build was detected from the known-SNP seeds in the first
100,000 rows only. Real exports are sorted by chromosome, so most seeds sit
far beyond that point. On full-size files detection therefore almost never
became confident, and the build came from the file's label, from the
`chr` naming of VCF contigs, or from the format default. When that label
was wrong or missing, every position-keyed lookup used the other build's
coordinates: real ClinVar findings were lost and false ones appeared,
including "pathogenic" rows the person does not carry. Three code paths
(analyze, the stats/extract/plink stream, and compare) each had their own
version of the logic.

Three further problems surfaced while fixing it:

- Contig naming is not evidence. UCSC hg19 uses `chr1`, and GRCh38 files with
  bare names are common.
- Most parsers could not say "this file states no build": they returned the
  GRCh37 default, so a missing label looked like a GRCh37 label.
- The seed table had wrong entries. rs80357906 was one base off in both
  builds, and two chr10 GRCh36 positions were 411 bases off, which left real
  GRCh36 files at a tentative 3/4.

## Decision

### Evidence

- **VCF contig lengths.** `##contig` lengths are decisive whatever the
  naming. Every nuclear contig (1-22, X, Y) differs in length between GRCh37
  and GRCh38, so any one of them decides. MT is the same length in both and
  is not evidence. Alt, decoy, random and unplaced contigs are ignored.
  Lengths that disagree with each other stop the run. A nuclear length that
  matches neither build stops the run with "unsupported assembly" (for
  example T2T-CHM13). When chrM is declared with the hg19 length (16,571 bp,
  not the rCRS 16,569 bp), MT records are not annotated and a note says so,
  whatever decides the build (contig lengths, `--build`, the assembly header
  or the seeds). *Amended (2.5):* in a GRCh37 VCF, mitochondrial records
  named `chrM` or `M` (UCSC hg19 names) are left out, counted and noted,
  unless the file declares that name with the 16,569 bp reference length.
  `MT` and `chrMT` are read. A declared hg19 length (16,571) leaves out
  every mitochondrial record, as above.
- **Known-SNP positions (seeds)** from the whole file. Arrays, and VCFs
  without contig lengths, are pre-scanned through the same parser before
  annotation; the scan stops early once every seed has been seen. Confidence
  is ADR-0037's: at least 3 matches and at least 80%. rsID-less records vote
  by position (ADR-0043).
- **Labels.** A VCF's `assembly=` or an array's build line. A parser reports
  `metadata["build"] = None` when the file states no build; only this rule
  applies a default.

### The rule

`--build grch37|grch38|grch36` overrides the label and the format default
and is recorded. `--build grch37` or `--build grch38` is still checked
against the file's own evidence: contig lengths of another build, or
confident known-SNP positions for another build, stop the run and the
message names both builds. Contrary positions that are not confident
(fewer than three, or under an 80% majority) only warn. `--build` with contig lengths that are mixed or unsupported runs as
forced. `grch36` is never checked: it runs the file with position-keyed
lookups off, as for a detected GRCh36 file (ADR-0025), so it cannot attach
another build's records, and a GRCh36 file that the rule stops has a
correct way through. `decide_build` accepts only GRCh36, GRCh37 or GRCh38 as an
override and raises `ValueError` for anything else, including callers of
the Python API. Otherwise:

| Evidence | Result |
|---|---|
| VCF contig lengths say X, label absent or X | X (`contig_lengths`) |
| VCF contig lengths say X, label says Y | stop |
| VCF contig lengths mixed, or unsupported | stop |
| No contig-length decision, seeds found but none at any build's coordinate (zero votes) | stop, label or not |
| Confident seeds say X | X (`seeds`); a different label gets the mismatch notice |
| Tentative seeds say GRCh36 | GRCh36 (`seeds`), the ADR-0025 fail-safe |
| No contig-length decision, seeds tied and GRCh36 among the tied builds | GRCh36 (`seeds`), the fail-safe, with the "inconclusive" warning, whatever the label |
| Label X, no seed inspected at all | stop: nothing checks the label (an mtDNA-only or X-only export, a small panel VCF) |
| Array label GRCh36, no seed inspected at all | GRCh36 (`file_label`): exempt from the stop above, like `--build grch36` (ADR-0025) |
| Label X, seeds tentatively X | X (`assembly_header` for a VCF, `file_label` for an array) |
| Label X, seeds tied and X is one of the tied builds | X, with the "inconclusive" warning |
| Label X, seeds tied and X is not one of them | stop |
| Label X, seeds tentatively Y | stop |
| VCF with no lengths, no label and no confident seeds | stop |
| Array with no label, seeds tentatively not the default | stop |
| Array with no label and no seed inspected at all | stop |
| Array with no label, at least one seed votes and the seeds' tentative winner is the default (a minority of seeds for another build, or at no build's coordinate, does not stop it) | the format default (`default`), with a visible note |
| Array with no label, seeds tied and the default is one of them | the default, with the note and the "inconclusive" warning |

GRCh36 is the one tentative winner that decides: it turns position lookups
off, so it cannot attach another build's records.

A tie means non-zero seed votes with no single winner; `detect_build`
reports the tied builds (`BuildDetectionResult.tied`). The build that would
otherwise decide is the label, the format default, or the contig-length
build. If it is one of the tied builds, the run proceeds on it, and the CLI
and the HTML report show a yellow "inconclusive" warning naming the tied
builds and the counts; JSON records them as `input.build_evidence.seeds_tied`.
If it is not, the run stops with the same message form as the other
contradictions. Under `--build`, a tie that leaves out the forced build
gets the `--build` warning instead; a tie that includes it gets nothing.

One exception comes first: when contig lengths did not decide (arrays, and
VCFs without lengths), a tie that includes GRCh36 takes the GRCh36
fail-safe, with the inconclusive warning, whether or not the label is among
the tied builds. A wrong GRCh36 costs only the position-keyed annotations
and says so; a wrong GRCh37 or GRCh38 on a real GRCh36 file attaches
another variant's records. When contig lengths decided, they already rule
GRCh36 out, and the rule above applies unchanged. With GRCh37 as the
format default, an unlabelled array's tie that excludes the default
always includes GRCh36, so it takes the fail-safe rather than stopping.

When contig lengths decide, there is no pre-scan. The seeds met during the
annotation pass are checked afterwards: any seed majority for another
build, or a tie that does not include the contig-length build, stops the
run (the realistic case is a lifted-over VCF that kept the
source assembly's contig lines). The stop message says the annotation pass
ran; nothing is written.

### One function for every command

`decide_file_build` (in `allelix/reports/_pipeline.py`, over `decide_build`
in `allelix/utils/build_detect.py`) is used by `analyze` and the focused
reports, by `stats`, `extract` and `export plink`, and by `compare`
(`--build1` / `--build2`). `compare` no longer takes a tentative winner.

### Stopping

"Wrong file" errors come first: a file no parser recognises, a
multi-sample VCF without `--sample`, a VCF of structural-variant /
copy-number records only, and a file with no genotype records each get
their own error on every command, whatever decides the build, never a
request for `--build`. When the header alone would stop the decision
(contig lengths that are mixed, unsupported or against `assembly=`), the
records are read up to the first usable one before the stop. A VCF whose
records were all dropped by the quality gate is not empty, and MT records
dropped for an hg19-length chrM are records.

A stop raises `BuildDecisionError` before any output exists. The message
names the file, says what is wrong and what evidence the file carries, and
asks for the build: "Tell Allelix which genome build the file uses (command
line: --build grch37 or --build grch38; for an old GRCh36 / hg18 file,
--build grch36, see docs/grch36-liftover.md)". The wording
works for a person on the command line and for a front end that offers the
same choice. No report, annotated VCF, diff, PLINK file set or `extract`
table is produced; `export plink` writes to temporary files and renames
them only on success.

### Notices on every command

Every command that decides a build (`analyze`, `stats`, `extract`, `export
plink`, `compare`) prints the same banner and notices for the same decision:
the label mismatch, the format-default note, the tie ("inconclusive")
warning, the hg19 chrM note and the GRCh36 notice. The wording lives in one
place (`allelix/reports/build_text.py`), which the HTML report also uses for
its mismatch and tie notices. A tie has its own source wording: it is never
described as decided by known-SNP positions.

### Provenance

`build_source` uses a fixed vocabulary: `override`, `contig_lengths`,
`seeds`, `assembly_header`, `file_label`, `default`.

- JSON report, schema 8 (ADR-0033): `input.build_source` and
  `input.build_evidence` (the file's label, the contig-length verdict, the
  seed counts, whether the label was overridden, MT records skipped). The
  diff loader accepts 8 and still loads 1 to 7.
- Annotated VCF, schema 0.3.0 (ADR-0036):
  `##ALLELIX_Build=<ID=GRCh38,Source=contig_lengths>`.
- HTML report: the source under the build ("forced with --build" when it
  was).

### Sources

- Seed positions: dbSNP b157 VCF releases for GRCh37.p13
  (GCF_000001405.25) and GRCh38.p14 (GCF_000001405.40); GRCh36 from dbSNP
  b130 on NCBI36/hg18 (UCSC `snp130` table). rs80357906 has no b130
  placement, so it has no GRCh36 entry. The corrected entries were
  cross-checked with the NCBI Variation API (RefSNP SPDI placements).
  `tests/fixtures/dbsnp_seed_placements.json` holds the snapshot; a test
  compares every entry with it.
- Contig lengths: the `##contig` headers of the GIAB v4.2.1 HG002 benchmark
  VCFs for GRCh37 (hs37d5) and GRCh38.

## Consequences

- A wrong or missing label can no longer silently select the other build;
  the run either resolves from evidence or stops.
- Arrays pay for one extra read of the file, also under `--build grch37` or
  `--build grch38` (the seeds check the forced build); `--build grch36`
  skips it. Measured on a one-CPU machine, before `--build` was checked:
  the pipeline with no annotators took a median of 7 s more over three runs
  on a 960,613-row 23andMe chip (14 s with `--build`, 21 s without); the
  decision alone took 1.4 to 3.5 s per real chip.
  A VCF with contig lengths pays nothing up front; the HG002 gVCF is
  decided from its header in milliseconds.
- Some inputs that used to complete now stop, by design: inconsistent VCF
  headers, rsID-less VCFs without contig lengths or an assembly,
  unlabelled chips whose few seeds point away from the default, and
  unlabelled chips that carry no seed at all, such as an X-chromosome-only
  export (all seeds are autosomal), and labelled files that carry no seed,
  such as an mtDNA-only export. The fix is `--build`, which then stops
  only if the file's contig lengths or confident seeds say another build.
- `stats` and `extract` gain `--build`; `compare` gains `--build1` and
  `--build2`. `stats` shows the decided build and its source.
- Parsers no longer fill in a default build, and the VCF parser exposes
  contig lengths instead of a chr-prefix flag.
- The committed fixtures are internally consistent (ADR-0015): contig
  lengths, `assembly=` and seed rows name one build.
  `tests/fixtures/mock_vcf_mixed_contigs.vcf` is the deliberate negative case.

## Not done

- A second evidence source from the rsID identity table: when it is
  installed, the GRCh37 and GRCh38 positions of thousands of rsIDs could
  confirm or rule out a build for files the seeds do not cover, such as
  chromosome-subset exports. This is a possible future change, not part of
  this decision.

## Amendment 2026-10-25 (2.5.2): positions that match no build

The zero-vote stop (table row "seeds found but none at any build's
coordinate") told the user to pass `--build`, and a run forced with `--build`
over such a file went ahead with no warning. The positions matching only in
the JSON's `build_evidence`, as `seeds_matched` 0 of N, was the only sign.

Forcing a build does not fix a file whose positions match none, so:

- **The stop** names the evidence and asks for the original file from the
  testing company. It does not suggest `--build`, and neither does
  `compare` (`--build1` / `--build2`) for this stop
  (`PositionsMatchNoBuildError`, `suggest_override = False`).
- **A run forced with `--build`** over such a file still runs, as forced,
  but warns. A **Build warning** notice names the evidence (0 of N
  known-SNP positions match any supported build) on the console of every
  command and in the HTML. The reports' `report_warnings` keep it as
  `build_override_unmatched`. This holds whether the seeds came from the
  pre-scan or from the annotation pass (contig lengths with `--build`), and
  in the indexed extract path.
- **Unchanged:** `--build grch36` stays exempt, and confident contrary
  evidence still stops. Matching or uninspected seeds give no warning.

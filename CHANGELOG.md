# Changelog

All notable changes are documented here. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.5.0] - 2026-10-19

After upgrading, run `allelix db update`. This release changes how
Allelix reads ClinVar (version 5 → 9, shown as `iv:9`), so existing
ClinVar caches are invalidated and the next `db update` rebuilds them
(about 860 MB to download). Until then `allelix analyze`
skips ClinVar, says so on the console ("Skipping ClinVar: out of date
(built for iv:5, needs iv:9). Run `allelix db update`."), and puts a
warning at the top of the report.

### Added

- `scripts/check_reference_sample.py`, an independent release gate on the
  reports of a healthy reference sample (the HG002 benchmark, the CC0
  arrays): no ClinVar Pathogenic or Likely pathogenic row is called
  homozygous, and none sits on a variant whose ClinVar aggregate is not
  Pathogenic or Likely pathogenic, except a row the report marks
  `outside_aggregate` (an outranked P/LP submission kept on purpose because
  the aggregate is not a disease classification, ADR-0045); such a row is
  counted and listed, and the homozygous rule still applies to it. It reads
  report JSON only and imports nothing from Allelix, so it does not depend
  on the tool's earlier output being right. Run it on reports made with the
  default filters. `--counts-only` prints no row. Not a CI job.
- Every GWAS Catalog row says which way the study's effect allele points
  (ADR-0024): JSON `effect_direction` (`increases`, `decreases`, `unknown`;
  empty on other sources), and plain words in the HTML detail panel and
  after the significance in the terminal ("C lowers odds"). The words say
  when the effect allele is the reference allele, so a person homozygous
  for a protective reference allele sees that it is protective. The odds
  ratio, or the beta's "unit increase/decrease", decides; anything else is
  `unknown`.
- Genes for Good's 23andMe-style text download is read (`genesforgood`
  parser). The imputed file is marked as imputed in every report (JSON
  `input.method` `array-imputed`), and the HTML's imputation note now
  says how each imputed input was made instead of always describing
  low-pass sequencing.
- deCODEme raw data (`deCODEme_scan.csv`) is read (`decodeme` parser).
  Calls on the minus strand are complemented to the plus strand; the
  file's Build 36 coordinates are detected from known-SNP positions and
  run on the GRCh36 path (position lookups off, ADR-0025).
- Every report names every source and its part in the run (ADR-0049).
  JSON `sources` lists each registered source as `used`,
  `ready_not_used`, or `skipped` with a reason: the cache's state
  (`missing`, `outdated`, `unreadable`, `corrupt`, `wal`,
  `readonly_wal`) or the user's choice (`disabled`, `license`,
  `excluded`, `flag`, `not_applicable_build`). The HTML lists the
  sources whose cache could not be used, with the fix, and the ones
  turned off in one line. `annotators` still lists only the sources that
  ran.
- When ClinVar did not run, the report says so at the top: a warning in
  the HTML, on the console and in JSON `report_warnings`
  (`{code, message, count}`, code `clinvar_skipped`). Before, a report
  without ClinVar read as complete.
- A report warning when three or more variants on chromosomes 1 to 22
  are called homozygous for an allele ClinVar classifies as Pathogenic or
  Likely pathogenic and are rare in gnomAD or have no gnomAD frequency
  (`multiple_homozygous_plp`, with the count). It is counted before the
  display filter. Several such calls in one file usually point to a
  genotyping or file problem (ADR-0050).
- A call whose plausibility was not assessed is marked: a muted ○ with the
  reason in the HTML and the terminal.
- The two warnings about a compressed input file as a whole, a bgzipped
  file that may be truncated and bytes ignored after the compressed data,
  are also report warnings: the HTML shows them at the top and JSON
  `report_warnings` lists them (codes `input_possibly_truncated` and
  `input_trailing_bytes`). Before, they were on the console only.
- The annotated VCF (`--vcf-out`, schema 0.4.0) has a
  `##ALLELIX_SkippedSources=` line naming each skipped source and its
  reason, or `none`.
- Zip downloads are read as they are, when they hold one genotype file:
  every format, the file streamed from the zip (never unpacked to disk),
  gzipped or not. 23andMe and AncestryDNA deliver their raw data as a zip;
  before, a zip stopped with "No parser recognized". A zip with several
  files, a password, a compression method Python cannot read, or another
  zip inside stops with a message saying to unzip it first. A file over
  4 GiB once unzipped (a gzipped file counts at its `.gz` size) stops,
  saying to unzip it or gzip it instead. A damaged zip, including a file
  whose data is longer or shorter than its zip entry declares, gets the
  "looks truncated or damaged" error. macOS archive
  entries (`__MACOSX/`, `.DS_Store`) and folders are not counted as files.
- Array reports say how many calls could not be checked against ClinVar
  (ADR-0051): vendor probes without an rsID (23andMe `i5000123`) at a
  locus where ClinVar's own classification is Pathogenic or Likely
  pathogenic, and calls written as `D` or `I` (a deletion or an
  insertion), with how many are at such an rsID. The console and the HTML
  input-quality notice print them; JSON has `input.probe_classes`
  (`i_probes_at_plp_loci`, `di_calls_at_plp_rsids`) and
  `input.array_quality` (`di_calls`, plus `unreadable_genotypes`,
  `not_genotyped` and `wrong_column_rows`, 0 until a parser counts them).
  Both are null for VCF input; `probe_classes` is also null when ClinVar
  did not run. Nothing new is matched.

### Changed

- The top of the HTML report is shorter. Still in view: a two-sentence
  disclaimer, and the warnings that change what to do (ClinVar not used or
  out of date, a file that may be cut short, a build warning, three or more
  rare homozygous Pathogenic calls, a forced `--format`, a source that could
  not be used, high-value no-calls, and the caution for imputed input).
  ClinVar's reason for not being used is said once, in its warning, not
  again under the sources not used. One line then names only what the
  section below holds: calls that could not be checked, and how many ClinVar
  rows were scored down or not shown; when it holds only general notes, the
  line says so. Every other notice (input quality, ClinVar rows scored down
  or not shown, strand, sources turned off, the full disclaimer) sits word
  for word in a collapsed "Technical details" section. Open it before
  printing or saving the page as PDF to include those notices. JSON and
  terminal output are unchanged.
- A file no reader recognises is rejected with a message that puts the
  original file first: "Allelix doesn't recognise the format of …. Use the
  original file exactly as your testing company provided it, not edited
  or re-saved." `--format <name>` still reads such a file, at the user's
  own risk. It now gets a warning before anything else in the terminal
  and on every report (the terminal summary, the top of the HTML report,
  JSON `report_warnings` code `forced_format`): Allelix cannot check a
  forced file's layout, strand or build label, so its results may be
  wrong. `--format` with a reader that recognises the file is not warned.
  Every command that takes `--format` behaves this way, and `compare`
  checks each of its two files on its own. Before, the rejection's only
  advice was to force a format, and a forced file was read without a
  word. `docs/formats.md` has a "Forcing a format" section.
- GWAS Catalog risk alleles and SNPedia genotypes are read on the forward
  strand against the site's alleles in the rsID identity table (ADR-0024).
  A risk allele that is among them is compared as written, else its
  complement is. A/T and C/G sites, rsIDs the table does not place at the
  call, and genotypes that fit neither reading abstain. Before, GWAS
  compared the letters as written, so at an A/T or C/G site a risk allele
  given on the other strand matched the wrong carriers. SNPedia pages with
  a stated orientation keep it, now checked against the site's alleles.
  One-row SNPedia sites of unknown orientation could not show they were
  palindromic and matched as plus. The JSON `strand_abstentions` block
  counts what abstained, per source and reason (GWAS association rows,
  SNPedia sites), and the console and HTML say so. Without the identity
  table, GWAS risk-allele rows and single-base SNPedia sites abstain, and
  so does GRCh36 input, which the identity table does not place: a GRCh36
  report loses those GWAS and SNPedia findings and says how many.
- Panel coverage (`--filter-file`) puts each requested rsID in one of
  four states (ADR-0051): `annotated`, `typed_no_annotation`, `no_call`
  or `not_on_chip` ("not in file" for VCF input). A no-call was counted as
  genotyped with no findings. A half call with a shown finding is
  annotated. JSON `panel_coverage` is `{requested, counts, variants}`:
  `counts` has the four states and sums to `requested`, and each entry of
  `variants` gives `rsid`, `gene`, `state`, `genotype` and
  `annotation_count`. `found`, `missing` and `no_findings` are gone. The
  HTML lists every requested rsID with its gene, genotype, state and
  number of annotations (a VCF half call's genotype is `T/.`, as in the
  annotations), and the terminal report prints the four counts, then the
  no-call and not-on-chip rsIDs.
- JSON values are machine codes (ADR-0050). `zygosity` is `homozygous`,
  `heterozygous`, `hemizygous`, `half_call` or `no_call` (it was
  "Homozygous", "Heterozygous", "Half call", "No Call"). `plausibility` is
  `not_assessed` where it was empty, and the new `plausibility_reason`
  says why: `no_frequency`, `gnomad_not_run`, `chromosome`, `half_call`
  or `no_allele`. HTML and the terminal show the display words; "No Call"
  is now "No call".
- A call the file declares as one copy, a one-letter array call (23andMe
  writes a male's X calls this way) or a one-index VCF genotype, is shown
  as Hemizygous with its one allele (`T`), not as Homozygous `TT`, and its
  plausibility counts one copy (ADR-0050). Findings do not change. A male
  X call written as two letters (AncestryDNA, FTDNA) is still Homozygous:
  Allelix does not infer sex.
- ClinVar rows are scored against ClinVar's own classification of the
  variant (ADR-0045). A Pathogenic or Likely pathogenic submission on a
  variant ClinVar calls "Conflicting classifications of pathogenicity" is
  capped at 4.0, so it leaves the default report; under a Benign or Likely
  benign classification the cap is 1.0 to 2.0. Drug response, risk factor
  and similar rows are not capped. JSON and the HTML detail panel show
  ClinVar's classification and how many submitters agree; the terminal
  shows it after a row that was scored down, for example
  `pathogenic (1 of 10; ClinVar: conflicting)`. Each report says how
  many Pathogenic or Likely pathogenic rows were scored down: 9 to 29 on
  three of the openSNP array files in `test_data/real` and the HG002
  benchmark VCF (TYR rs1042602, BBS2
  rs4784677 and HFE rs1799945, where 40 of 45 submitters say Pathogenic
  and 43 of 45 Pathogenic or Likely pathogenic, among others).
  These rows are still shown at `--min-magnitude 4` or lower.
- A ClinVar Pathogenic or Likely pathogenic row whose own submission has
  zero stars ("no assertion criteria provided" or "no classification
  provided") is scored 4.0, below the default floor, when gnomAD gives
  ClinVar's allele a frequency of 5% or more, the ACMG/AMP BA1 level
  (ADR-0052). The frequency is looked up for ClinVar's own REF and ALT at
  the rsID's GRCh38 locus, so it applies to GRCh37 files too. An
  allele that common does not cause a rare disease: CDKN2B-AS1 rs1063192,
  Likely pathogenic for a breast tumour from one zero-star submission, with
  a gnomAD frequency of about 0.65, was a homozygous Likely pathogenic
  finding in the HG002 benchmark report. The row says why, and each report
  counts these rows ("N ClinVar Pathogenic or Likely pathogenic rows not
  shown: zero-star classification on a common allele"). Rows with one star
  or more are unchanged, HFE C282Y among them. So is a row whenever the
  lookup is in doubt (no rsID, a locus or allele that does not agree
  between the identity table and both ClinVar builds, no gnomAD row) or
  gnomAD, the identity table or a ClinVar build is not installed.
- ClinVar classifications beside the five-tier scale have their own
  magnitude instead of 5.0: Pathogenic, low penetrance 7.5 (SERPINA1
  rs17580 and HFE rs1800562, for example); Likely pathogenic, low
  penetrance 6.5; established risk allele and confers sensitivity 6.0;
  likely risk allele, Affects and protective 5.0; the three VUS
  subdivisions 4.5, 4.0 and 3.5; uncertain risk allele 4.0. A
  classification with no magnitude is not shown, and the report counts
  it.
- JSON report `schema_version` 8 → 9. Each annotation gains
  `clinvar_aggregate`, `submitters_same`, `submitters_total`,
  `aggregate_capped`, `outside_aggregate`, `common_allele_capped` and
  `clinvar_allele_af`,
  plus `chromosome`, `position` and `ref`, so a finding without an rsID
  still names its variant. A top-level `clinvar` object gives
  `pathogenic_rows_scored_down`, `unrecognised_rows_hidden` and
  `common_allele_rows_scored_down`. For a VCF,
  `input.vcf_quality` gains `spanning_deletion_sites`, `half_calls`,
  `malformed_rows`, `no_genotype_rows`, `unusable_genotype_calls`,
  `symbolic_allele_calls` and `trimmed_matches`.

### Fixed

- `--filter-file`: a `#` after an entry now starts a comment, as a `#` at
  the start of a line already did, so `rs1801133  # MTHFR` keeps
  rs1801133. Before, the whole line was taken as a gene name that matched
  nothing, and the entry left the report with no warning. The help says
  so. A list saved with a byte order mark is still read in full.
- `analyze --category` takes only `clinical`, `pharma` or `trait` (in any
  case). A typo such as `clinicals` now stops with a message listing the
  three. Before, it gave an empty report with exit 0 and no warning. The
  help listed only two of the three.
- `config get --help`, `config set --help` and `db path --help` show their
  key lists and examples line by line. Before, they ran together into one
  paragraph with a stray `\b`, and the `db path` example did not run as
  shown.
- The GRCh36 warning said GWAS Catalog and SNPedia still run. It now says
  that ClinPGx still runs, and that GWAS Catalog risk alleles and
  single-base SNPedia genotypes are not assessed on GRCh36, as the counts
  printed with it show.
- A VCF whose header declares an imputation tool's field (Minimac's
  `IMPUTED`, `TYPED` or `R2`; Beagle's `DR2` or `IMP`; the per-call dosage
  `DS` that Minimac, Beagle and GLIMPSE write; GLIMPSE's `HS`) is now marked
  imputed, as Genes for Good's and MyHeritage's imputed files already are.
  Before, it was read as measured with no caution. Every report now says
  the genotypes are imputed and to confirm any finding by direct
  genotyping, and the HTML's "Reading This Report" explains imputation in
  place of the note on array probes. JSON `input.method` is `vcf-imputed`.
  Only the header's declarations count, never the values; genotype
  probabilities (`GP`) alone do not count. Measured VCFs (benchmark call
  sets, GATK and DeepVariant gVCFs) are unchanged.
- Four ways a file could be read wrongly with no warning are fixed:
  - A no-call written `00` or `NN` (or `0` or `N` in a one-letter or
    two-column file) is now a no-call in every array reader, as `--` is.
    Before, it was read as a call of "0" or "N": it matched nothing, but
    `stats` counted it as a homozygous call and the warning for
    clinically important SNPs with no call did not fire.
  - A `--filter-file` list saved with a byte order mark (as some
    Windows editors and spreadsheet "CSV UTF-8" exports do) lost its
    first entry. It is now read in full.
  - Genes for Good's imputed VCF download is marked imputed, as its text
    download already is, so every report gives the caution for imputed
    genotypes. It is recognised by its `##Genes_for_Good_version` header
    line together with phased calls.
  - In a GRCh37 VCF, a mitochondrial record named `chrM` or `M` (UCSC
    hg19's names for a 16,571 bp sequence that is not the reference ClinVar
    uses) is now left out, counted, with a note, unless the file declares
    that same name with the reference length (16,569). Before, it was read
    as the reference whenever no hg19 length was declared. Records named
    `MT` or `chrMT` are read as before, and a declared hg19 length still
    leaves out every mitochondrial record.
- The HTML report's "Understanding Magnitude Scores" legend now states the
  scores Allelix assigns. Its ClinVar, ClinPGx and GWAS sections dated from
  1.4.0 and disagreed with the scoring in many rows: for example, it said a
  single-submitter Pathogenic call scores 8 (every Pathogenic submission
  scores 9, whatever its review status), uncertain and conflicting 5 (4),
  risk factor and drug response 4 (6 and 6.5), and ClinPGx Level 2A 6 (7).
  The three sections are now generated from the scoring tables. They add the
  GWAS odds-ratio step and the 3.0 cap for an unnamed risk allele, ClinVar's
  own classification capping a submission, the zero-star rule for common
  alleles, and that a classification without a score is not shown. Scores
  and findings do not change.
- `--vcf-out` stamps, for each source, the row the report lists first at
  that variant (highest magnitude, then rsID). It used to stamp the last
  row, so a site the report shows as Pathogenic could carry a VUS in the
  annotated VCF (ADR-0036).
- `export plink` no longer reads a call on the other strand. A call whose
  alleles are not gnomAD's REF and ALT at the site, as written, is
  exported as missing and counted on the console ("N call(s) fit neither
  gnomAD allele at their site"). Before, such a call could be written as
  its complement: A/A at a C>T site came out as T/T.
- With VCF input, a person homozygous for a variant allele the reference
  genome itself carries (Factor V Leiden TT on GRCh37, and about 170 other
  ClinVar sites) now gets the finding. Their file has a `0/0` call or a
  reference block there, never a variant record. Allelix watches these
  sites and assesses a hom-ref call that passes the GQ gate as the
  homozygous variant genotype (JSON `input.vcf_quality.reference_carried_calls`).
  A plain VCF with no record at such a site cannot say which it is: the
  report counts those sites (`reference_carried_not_in_file`) and says a
  gVCF or an all-sites VCF would.
- A ClinVar normal-allele record now reads the same on both builds. A
  CASP12 rs497116 carrier of the normal allele G got the risk-factor row
  on GRCh38, where the record is written A>G, but not on GRCh37, where it
  is written G>G and was dropped at load. Such records are now kept on
  both builds (still without Pathogenic or Likely pathogenic
  classifications, as before). With a VCF on GRCh37, the GG genotype is a
  `0/0` call there, which the site watch above assesses. This changes the
  ClinVar cache (interpreter version 9). Such a record written G>G does
  not make a homozygous-reference call a ClinPGx or SNPedia non-finding:
  the reference allele there can be the one ClinVar classifies, so those
  sources show their rows for that call, as 2.4.2 did.
- A VCF record whose REF is the complement of ClinVar's REF no longer
  matches ClinVar. A VCF is forward strand, so such a record is
  mis-coordinated input. One-base records whose REF is none of ClinVar's
  at the position are counted (`input.vcf_quality.clinvar_ref_mismatches`),
  and the report asks you to check the file's build and coordinates.
- ClinPGx and SNPedia treat a call as homozygous reference by ClinVar's REF
  only where ClinVar places the rsID at the call's own position. Before, an
  rsID placed elsewhere could hide a finding.
- A ClinVar match from a multi-allelic or untrimmed VCF record carries the
  matched row's own rsID where the record's ID column names another rsID
  ClinVar holds at that site. A merged or alias rsID is kept, as before.
- A homozygous one-base deletion (`AC>A`, 1/1) is shown as `A/A`, no longer
  as `AA`, which looked like an SNV genotype.
- `--gwas-min-magnitude` replaced `--min-magnitude` for GWAS rows instead
  of taking the larger of the two, as ADR-0024 says. A GWAS floor
  set below the global one let GWAS rows in below the floor asked for. The
  must-include rsIDs keep their exemption from the GWAS floor and never go
  below `--min-magnitude`.
- A 23andMe export whose rows were re-saved with spaces instead of tabs
  (by a text editor or spreadsheet) is read. Before, the file was
  recognised as 23andMe by its banner and then every row was skipped, so
  the run stopped with "contains no genotype records". A row without a
  tab is now split on runs of spaces.
- With one ClinVar build installed (`allelix db update --build grch37`),
  `analyze` used no ClinVar at all and said it was not installed. It now
  uses the installed build; input on the other build skips ClinVar and
  names the missing build's cache (ADR-0049).
- With one ClinVar build installed, a person who does not carry Factor V
  Leiden (CC at F5 rs6025 on GRCh37) got a ClinPGx row for hormonal
  contraceptives and thrombosis (level 2B), and ClinPGx and SNPedia rows
  at other sites where the GRCh37 reference carries the variant allele
  read the wrong allele as the variant. The normal allele at those sites
  came only from the other build's ClinVar cache. Without that cache it
  now comes from the rsID identity table's reference base on the other
  build. With both ClinVar builds installed nothing changes.
- A source turned off in config, excluded, or turned off by a `--no-`
  flag is no longer reported as "not installed" when its cache is absent.
- The HTML note on source-specific magnitude thresholds named sources
  that did not run.
- `analyze --no-update` writes nothing to the data directory, so the
  directory can be read-only (ADR-0048). Before, it could add a column,
  stamp or rebuild a cache, build an index, create `config.toml`, or leave
  WAL side files. A gnomAD, AlphaMissense or CADD cache, installed in WAL
  mode, could not be read at all in a read-only directory. Every cache now
  opens read-only, and a WAL cache is never opened. The installer stores
  caches in rollback-journal mode, and the next `allelix db update`
  converts ones already installed, in place and without downloading, as
  does `analyze` without `--no-update`. Until then `analyze --no-update`
  skips a WAL cache with that reason, and `db status` shows it instead of
  a version. `db update` does not download a WAL cache again when it
  cannot convert it. Local migrations (a legacy version stamp, the GWAS
  re-ingest from its cached file, the SNPedia parse) run in the update
  step: `db update`, or `analyze` without `--no-update`.
- A source that cannot be used is named on the console with its reason
  and the fix: not installed or out of date (run `allelix db update`),
  cannot be read (check the file's permissions), damaged (run
  `allelix db update`, which replaces only that source), in WAL mode from
  an older Allelix (run `allelix db update` once to convert it; nothing is
  downloaded), or in WAL mode in a data directory you cannot write (run
  `allelix db update` once with write access). Before, every case said
  "run `allelix db update`". A cache whose main table has no rows is now
  damaged; before, ClinVar, ClinPGx, GWAS, gnomAD, AlphaMissense, CADD and
  rsID-identity caches with a current version tag and no rows were used as
  ready.
- A GRCh38 VCF without rsIDs no longer stalls on Ubuntu 24.04's SQLite
  (3.45). The lookup that finds rsIDs in gnomAD by position searched the
  cache by chromosome only. It now searches by chromosome and position,
  one chromosome per statement, and needs 4,063 statements instead of
  8,886. On the HG002 GRCh38 benchmark VCF with default settings, on a
  2-vCPU machine, the analysis did not finish in 30 minutes before this
  change and now takes about 7 to 8.5 minutes. Allelix also no longer
  builds an index on the gnomAD cache when it opens it, so a read-only
  cache built without that index no longer fails with "attempt to write
  a readonly database".
- VCF calls are no longer dropped without a word (ADR-0047). On the
  HG002 GRCh37 benchmark VCF (2026-10-10 ClinVar data) the unfiltered
  report gains 1,331 ClinVar rows on 978 variants, none removed; the
  default report gains one, OPRM1 rs71021021 (drug response, Tramadol).
  - A homozygous VCF deletion to one base (`AC>A` 1/1) now matches
    ClinVar. The indel-anchor rule meant for array calls also blocked it.
    HG002: 831 variants, among them NR2E3 rs11351249 (Likely pathogenic,
    shown at 4.0 under ClinVar's Conflicting classification) and OPRM1.
  - An allele written with bases it shares with REF (`GTT>GT,G`,
    `CA>TA`, or `ATG>ATGTG,A` 2/1) is trimmed, per carried allele, before
    ClinVar matching. HG002: 147 variants, among them a CFTR TG-tract
    Pathogenic submission (4.0, Conflicting). Indels are not left-aligned,
    so one the file places at a different position from ClinVar's still
    does not match.
  - A half call (`./1`) is assessed on its called allele by ClinVar and
    the GWAS Catalog, shown as `T/.` with zygosity "Half call"; ClinPGx and
    SNPedia, which need both alleles, skip it.
  - The calls the parser drops are counted with their reason on the
    console, in HTML and in JSON: malformed rows, rows without a genotype,
    genotypes naming an allele the record lacks or a symbolic allele such
    as `<NON_REF>` beside a real ALT, and spanning-deletion (`*`) calls,
    which are counted on one or both haplotypes and listed by position
    (516 in HG002). A no-call
    written `./.` is still reported as a no-call, not counted as a drop. A
    genotype naming the only ALT of a gVCF reference block (`<NON_REF>` or
    `<*>` alone) is not checked.
- VCF users now see ClinVar variants that have no rsID (ADR-0046).
  Allelix dropped them when building the ClinVar cache, so a VCF carrying
  one at the exact position got nothing. In the 2026-10-10 data they hold
  98,211 GRCh37 Pathogenic or Likely pathogenic rows on 91,391 variants,
  23.4% of the Pathogenic or Likely pathogenic variants ClinVar gives VCF
  alleles for. They are matched by chromosome, position, REF and
  ALT, for VCF input only: an array call's reference allele comes from a
  table keyed by its own rsID, which says nothing about a variant without
  one. Reports name such a finding by its position, for example
  `1:1000 A>G`, and `--vcf-out` writes `.` for its rsID.
- `--vcf-out` stamps ClinVar on multi-allelic records with no ID. It
  placed annotations by rsID, so a record with ID `.` whose ClinVar rows
  have an rsID got no stamp; it now finds the record by its position and
  the carried ALT, and stamps the right ALT (8 such records in HG002).
- Factor V Leiden (F5 rs6025) carriers see its Pathogenic submissions
  again. ClinVar's classification of the variant is an expert panel's
  "drug response", which outranks every disease submission, so Allelix
  showed only the drug response row. Where ClinVar's classification is
  not a disease classification, the outranked Pathogenic and Likely
  pathogenic submissions are kept, judged in their strongest review tier
  as ClinVar judges a variant: shown as submitted when that tier agrees
  (F5: 22 Pathogenic and 1 Pathogenic, low penetrance), capped at 4.0
  when it conflicts (MTHFR rs1801133), and not shown when it holds no
  Pathogenic or Likely pathogenic submission (VKORC1 rs9923231). In the
  2026-10-10 data this affects 46 GRCh37 variants with an rsID: 38 keep
  their Pathogenic and Likely pathogenic submissions and 8 show none.
- `allelix db update` no longer replaces a working ClinVar cache with an
  empty one when ClinVar changes the layout of its files. The update
  checks the column names in both ClinVar files and stops with an error
  that names the changed column. It builds the new cache beside the old
  one and replaces the old one only when the build finishes with records.
  Before, an inserted or renamed column made every record unreadable, and
  the update installed an empty cache that reported ready.
- A ClinVar placeholder such as "covers 10 genes, none of which curated
  to show dosage sensitivity" is no longer shown as the gene (rs4124874,
  for example).
- Submissions that ClinVar has flagged are no longer reported. They
  appeared only where every submission on a variant is flagged: 120 rows
  in the 2026-10-10 ClinVar data, all but one Pathogenic or Likely
  pathogenic.
- Array files (23andMe, AncestryDNA, Family Tree DNA, Living DNA,
  MyHeritage, MyHappyGenes) are read when gzipped, when they start with a
  UTF-8 byte order mark, or when their first line is blank. Before, these
  stopped with "No parser recognized" (some formats accepted one or two of
  the three). Family Tree DNA delivers its raw data as a `.csv.gz`, which
  can now be passed as it is, without unpacking it first.
- A compressed genotype file that was cut short or damaged is no longer
  read as complete or stopped with a bare "Aborted!". A file that ends in
  the middle of its compressed data, or holds corrupt data, stops with an
  error naming the file ("looks truncated or damaged"). A bgzipped VCF cut
  exactly between two blocks, which raised no error and was analysed with
  the records after the cut silently missing, now gets a warning that its
  end-of-file block is missing and the file may be truncated. Bytes after
  the end of the compressed data that are not compressed data (such as a
  web page saved onto the end of a download) are ignored with a warning
  instead of stopping the run after every record was read.
- A VCF whose chrM is declared with the hg19 length (16,571 bp) no longer
  gets mitochondrial findings at reference-sequence (rCRS) positions when
  its build comes from `--build`, the assembly header or the known-SNP
  positions. Its mitochondrial records are now left out with the same note
  as when the contig lengths decide the build.
- A file that states its build but carries none of the known-SNP positions
  (for example a 23andMe mtDNA-only export labelled build 37, or a small
  panel VCF with only an `assembly=` line) no longer runs on the label
  alone, which could attach another build's ClinVar records by position.
  It stops, names the label and the missing positions, and asks for
  `--build`. An array file labelled GRCh36 (Build 36) is exempt, as
  `--build grch36` is: it runs on the GRCh36 path, with position-keyed
  lookups off, so it cannot attach another build's records.
- `--build grch37` or `--build grch38` no longer skips the file's own
  evidence. Contig lengths of the other build, or confident known-SNP
  positions for it, stop the run and the message names both builds.
  Contrary positions that are not confident (fewer than three, or under
  an 80% majority) print a "Build warning" instead, as does a tie that
  leaves out the forced build. `--build grch36` is
  not checked, because it turns position-keyed lookups off. Arrays run
  under `--build` now take the same known-SNP pre-scan as without it.

### Documentation

- The README no longer says to pass `--cadd` to an analysis: `analyze` has
  no such flag. CADD runs once it is downloaded (`allelix db update --cadd`)
  and enabled in config (`allelix config set sources.cadd true`);
  `--no-cadd` turns it off for one run.
- The README roadmap no longer lists cross-source conflict surfacing,
  risk-allele display and a structured MTAG flag as shipped in v2.2. None
  of them has shipped.
- The README's IL10 rs1800896 example is removed. Current ClinVar has no
  Pathogenic classification for it, only "Uncertain risk allele" for
  leprosy susceptibility.
- `docs/formats.md` describes the 23andMe file, including the `file_id`,
  `signature` and `timestamp` lines that current v5 exports put before the
  "generated by 23andMe" banner. These files were already read correctly;
  a test now covers that layout.
- `docs/formats.md` says what happens when a MyHeritage file's preamble is
  deleted: the rest of the file has the same layout as a Family Tree DNA
  file, so it is read as one, and a low-pass export is then no longer
  marked as imputed.
- `docs/formats.md` says that every format reads Windows (CRLF), Linux
  (LF) and classic Mac OS (CR only) line endings, plain or gzipped. This
  was already the case; tests now cover each format with each ending.
- Test files: the AncestryDNA, FTDNA, Living DNA and MyHeritage test
  fixtures are rewritten from scratch with invented genotypes, and the
  exome VCF and the 23andMe lookalike in `test_data/edge_cases/` are
  replaced by synthetic files of the same shape. Reading is unchanged.
  The Sequencing.com structural-variant test VCF's `##fileDate` is a
  placeholder.

## [2.4.2] - 2026-10-13

After upgrading, run `allelix db update`. This release changes how
Allelix reads ClinVar: `CLINVAR_INTERPRETER_VERSION` bumps `4 → 5`, so
existing ClinVar caches are invalidated and the next `db update`
rebuilds them (about 860 MB to download). Until then `allelix analyze`
skips ClinVar and says so on the console ("Skipping unready annotators:
clinvar ... (run `allelix db update` to populate)").

### Fixed

- Array files (for example 23andMe v5 exports) no longer get a false
  homozygous pathogenic phenylketonuria (PAH) finding. Since 2.4.0,
  ClinVar matched array calls by position, so the probe ID rs62508729,
  which dbSNP has withdrawn and which reads TT on these chips, picked up
  the ClinVar record for rs62507269, a different variant at the same
  position. Two changes:
  - array calls whose rsID the rsID identity table does not hold no
    longer borrow a reference allele from another rsID at the same
    position;
  - an array call with no reference allele matches ClinVar by position
    only when ClinVar holds its rsID on that build.

  VCF records are unchanged: their reference allele comes from the file.
  Array calls carrying a merged rsID (for example rs10735810, merged into
  rs2228570) also lose their position match, as before 2.4.0. Their GWAS
  Catalog and SNPedia rows stay, but no longer show gnomAD, AlphaMissense,
  CADD or plausibility values. On three public array genomes with the
  2026-10-10 ClinVar snapshot, each file loses at most two false matches
  from withdrawn IDs and 200 to 290 ClinVar rows under merged IDs, of
  which at most three per file are not benign or likely benign.

- Factor V Leiden (rs6025) is no longer reported the wrong way round on
  GRCh37. The GRCh37 reference base at rs6025 is the Leiden allele, so
  ClinVar writes the Leiden variant there with the same reference and
  alternate allele (T and T), and writes its record for the normal allele
  (`c.1601=`) as T to C. Allelix dropped the first and matched the second:
  - people who do not carry Factor V Leiden (CC) were shown two ClinVar
    "Pathogenic" rows;
  - heterozygous carriers (CT) were shown the same two rows;
  - homozygous carriers (TT) were shown nothing from ClinVar, ClinPGx or
    SNPedia from 2.4.0, and before 2.4.0 the ClinVar row but nothing from
    ClinPGx or SNPedia.

  The same ClinVar data gives the false rows on 2.2.3, 2.3.1, 2.4.0 and
  2.4.1, and 2.4.0 and 2.4.1 also drop the Leiden row. GRCh38 was right.
  Now, on both builds, CC gets no ClinVar finding, and CT and TT get the
  Leiden row (ClinVar "drug response") with their zygosity. ClinPGx and
  SNPedia treat CC, not TT, as the normal genotype on GRCh37, as on
  GRCh38, and name T as the allele carried.

  The same applies where the GRCh37 reference carries another ClinVar
  variant allele (171 sites in the 2026-10-10 snapshot, most of them
  benign; for example DPYD rs1801265 and CYP3A5 rs776746), and at the 5
  sites where the GRCh38 reference does. ClinVar records for the normal
  allele no longer give a "Pathogenic" or "Likely pathogenic" row; their
  other rows are unchanged.

  This corrects the 2.4.0 release notes, which said such records were never
  reported. 2.4.0 hid only records written with the same reference and
  alternate allele, and on GRCh37 the rs6025 normal-allele record is
  written T to C.

  Known limitation, not new in this release: with VCF input on GRCh37, a
  person homozygous for a variant the GRCh37 reference carries (Factor V
  Leiden TT, and the other such loci) still gets no finding there, because
  the file has no variant record at that position. Array files are not
  affected. If a gVCF has no reference block at rs6025, the high-value
  list names rs6025 as not covered. A plain VCF with no record there, a
  homozygous-reference call, or a gVCF reference block covering it leaves
  rs6025 off the list, and the report reads as if it were genotyped with
  no findings.

## [2.4.1] - 2026-10-12

### Added

- Reports say where the build came from: forced with `--build`, the
  chromosome lengths, the known-SNP positions, the file's build line or
  VCF assembly, or the format default. The JSON report has
  `input.build_source` and `input.build_evidence` (schema version 8;
  `--diff` still reads older reports), the HTML report shows the source
  under the build, and the annotated VCF writes
  `##ALLELIX_Build=<ID=...,Source=...>` (VCF schema 0.3.0).
- `--build grch36` (and `--build1` / `--build2 grch36` for `compare`):
  runs an old GRCh36 / hg18 file that Allelix cannot identify by itself,
  with position-keyed lookups off, as for a detected GRCh36 file.
- `--build` for `stats` and `extract`, and `--build1` / `--build2` for
  `compare`.
- MyHeritage array files: the "human reference build 37" header sentence is
  now read as the file's build label.
- VCFs whose mitochondrial chromosome has the old hg19 length: mitochondrial
  records are skipped, with a note, because their positions do not match
  the reference sequence ClinVar uses.

### Removed

- `scripts/tag-release.sh`, a release-tagging helper, has been removed.
  Tag pushes are still checked against the `pyproject.toml` version by
  the pre-push hook in `.githooks/`.

### Fixed

- **Reports on the wrong genome build: files whose build label was wrong
  or missing.** Allelix looked for its build-detection SNPs only in the
  first 100,000 rows of a file. Real raw data files are sorted by
  chromosome, so on full-size files the build almost always came from the
  file's label, from the way a VCF names its chromosomes (`chr1` or `1`),
  or from the format's default. When that guess was wrong, the report used
  the other genome build's positions: real findings were lost, and
  findings the person does not carry appeared, including ClinVar
  "pathogenic" ones. This affected chip files whose header names the wrong
  build or none, UCSC hg19 VCFs with `chr`-named chromosomes, and GRCh38
  VCFs with plain chromosome names, in every release since 2.0.0.
  Allelix now decides the build before annotating, from the evidence in the
  file itself: a VCF's declared chromosome lengths, or well-known SNP
  positions across the whole file. The run stops instead of guessing when
  that evidence contradicts the file's label, or when the file carries no
  usable evidence: a VCF without chromosome lengths, without an assembly
  line and without enough of those SNPs, a chip with no build line and
  none of the SNPs (such as an X-chromosome-only export), or SNPs that sit
  at none of the supported builds' positions. The message says what was
  found and asks you to name the build with `--build grch37` or `--build
  grch38` (or `--build grch36` for an old GRCh36 / hg18 file). Nothing is
  written in that case.
  If you have reports from such files, run them again with this release.
- `stats`, `extract`, `export plink` and `compare` now decide the build the
  same way as `analyze`. `compare` no longer settles on a weak guess.
- Three entries in the build-detection SNP table were wrong (one SNP on
  both modern builds, two on GRCh36). They now match dbSNP. Old GRCh36
  FTDNA files that were detected only tentatively are now detected with
  confidence.
- The note about GRCh38 VCFs without rsIDs no longer claims GRCh37 gives
  more findings. It now appears only when the rsID identity table or
  gnomAD is not installed, says which, and tells you to run
  `allelix db update`.
- A failed `export plink` no longer leaves partial `.bed` / `.bim` / `.fam`
  files behind.
- `export plink` no longer writes a PLINK set with no variants. When every
  record is a no-call or an indel, or was dropped (for example mitochondrial
  records under an hg19-length chrM), it stops and says what was skipped.
- AlphaMissense and CADD caches without a schema version stamp are no longer
  treated as current. `allelix db update` replaces such a cache with the
  current one, as it already did for gnomAD.
- Very old GWAS Catalog and ClinPGx caches, made before Allelix recorded a
  version inside each cache, are accepted only when the version noted in
  them is the current one. Otherwise the GWAS cache is rebuilt from the
  download kept next to it, if that file is still there; if not, GWAS
  Catalog, like ClinPGx, is skipped until you run `allelix db update`.
- When an upgrade changes how Allelix reads ClinPGx data, `allelix analyze`
  no longer rebuilds the ClinPGx cache itself, which kept the old CPIC
  allele-function data. It skips ClinPGx and asks you to run
  `allelix db update`, which rebuilds the cache and fetches the CPIC data
  again. `analyze` makes no network request for it.
- The source distribution no longer ships a partial, unrunnable copy of the
  test suite (the top-level `tests/test_*.py` files without their fixtures).
  It contains the package, the licence, the README and the build files.
- Building Allelix from source now declares the setuptools version it needs
  (77 or later, for the SPDX licence field). With an older setuptools the
  build failed with a `project.license` configuration error.
- The GRCh36 warning no longer claims gnomAD results are complete: on
  GRCh36 files ClinVar, gnomAD, AlphaMissense and CADD give no values, and
  only ClinPGx, GWAS Catalog and SNPedia (matched by rsID) run. It also
  says how GRCh36 was chosen.

## [2.4.0] - 2026-10-11

### Fixed

- **Chip reports no longer show reference-allele records as findings.**
  ClinVar publishes records that describe the reference allele itself
  (written with `=`, for example `c.1601=`). On array input these were
  matched by rsID, so a person homozygous for the normal allele could be
  shown a "Pathogenic" ClinVar row whose allele was the reference base,
  for example at F5 or BRCA1. Such records are now never reported. To
  check an existing JSON report, re-run the same raw file with this
  release and compare: `allelix analyze RAW_FILE --diff OLD_REPORT.json`.
  ClinVar rows listed as removed whose ClinVar allele page is titled
  with `=` were affected. Other removals in that diff can come from
  database updates since the old report.

- **gVCF and other VCFs without rsIDs no longer drop ClinVar findings.**
  When a VCF record had no rsID, a ClinVar variant with several
  submissions for the same allele, or a call carrying two different
  non-reference alleles, could be skipped silently, while the same
  genotype on a genotyping chip was reported. This mainly affected runs
  without the gnomAD cache. ClinVar findings are now selected by
  chromosome, position, reference allele and alternate allele on every
  input type, never by rsID alone. The same genotype now produces the
  same ClinVar findings from a chip file and from a gVCF. Chip probes
  that carry a vendor ID instead of an rsID (such as 23andMe `i`
  probes) are never matched by position, because the strand their
  alleles are reported on is unknown. A chip rsID that ClinVar records at a different position is
  never given the record of another variant that shares its position.

- **ClinVar no longer matches by rsID at a different position.** When a
  file's positions disagreed with the genome build in use (for example
  a forced `--build`), a ClinVar row could still be attached by rsID.
  Rows must now match the variant's position in the build-specific
  ClinVar data.

- **VCF rows the variant caller rejected were reported as findings.**
  Allelix ignored the VCF FILTER column. A row the caller had marked as
  unreliable (for example `LowQual`, or a VQSR tranche in a filtered
  whole-genome VCF) was annotated and reported exactly like a passing
  call, and the annotated VCF output added Allelix fields to it. Such
  rows are now skipped by default and counted by filter name in the
  report, so you can see how many were left out. `--keep-filtered`
  includes them again; the annotated VCF then marks each one with
  `ALLELIX_FILTER`. A filtered row never receives Allelix fields
  otherwise. Genotype quality is now read as well: a call with GQ below
  20 (set with `--min-gq`), or with no reads, is reported as a no-call
  instead of as a genotype. JSON reports record these counts in
  `input.vcf_quality`.

- **MyHeritage and Family Tree DNA files re-saved through a spreadsheet
  are recognised again.** Opening an export in a spreadsheet and saving
  it pads every row with trailing commas. It can also leave an empty
  first row that pushes the "MyHeritage" line down, and remove the `#`
  from comment lines. Such files were not recognised, or their padded
  rows were skipped. They are now detected and parsed, including
  tab-delimited "save as text" copies. A genotype row with an empty
  genotype is still read as a no-call.

- **VCF files with only structural-variant or copy-number records now
  fail with a clear message.** Sequencing.com delivers a 30x genome as
  three VCF files. Two of them (`sv` and `cnv`) contain only records
  Allelix does not annotate, and analyzing either one used to finish
  with an empty report and no explanation. `allelix analyze` now stops
  and names the file to use instead (`snp-indel.genome.vcf.gz`).
  Structural-variant records in any VCF are counted in the report (JSON
  `input.vcf_quality.symbolic_rows`). Breakend records are no longer
  read as alleles.

- **`stats`, `extract` and `export plink` now recover rsIDs for VCF files
  without them.** Variant callers usually leave the VCF ID column empty.
  `analyze` recovered the rsIDs from position, but the other commands
  did not. As a result, `stats` missed high-value warnings, and `export
  plink` wrote every such variant without its reference and alternate
  allele. All four commands now use the same recovery step and the same
  installed databases. `stats` also reports high-value no-calls in a VCF
  the way `analyze` does, and gains a `--data-dir` option.

- **Genome build detection now works on VCF files without rsIDs.** Build
  detection looked for a set of well-known SNPs by rsID only, so it had
  no position evidence on variant-caller output with an empty ID column.
  A record at one of those SNPs' positions now counts as well. This
  matters most for GRCh37 files with `chr`-prefixed chromosome names,
  which were otherwise assumed to be GRCh38.

- **gVCF positions with no reads were treated as reference, and
  high-value no-calls were never listed for VCF input.** A gVCF reference
  block with no reads, or with low genotype quality, says nothing about
  the genotype. Allelix read it as a confident reference call. A block
  below the genotype-quality floor (`--min-gq`, default 20) or with no
  reads now makes the high-value positions it covers no-calls. So does a
  gVCF position that no record covers at all. The high-value no-call list
  is now filled in for VCF input, and each entry says why: for example a
  low-quality block, no coverage, a filtered row or a low-quality call. A
  plain VCF (without reference blocks) cannot tell reference from
  missing, so a high-value position without a record there is reported as
  unknown rather than listed. JSON reports record whether the input was a
  gVCF in `input.gvcf`.

- **VCF spanning-deletion alleles are no longer read as real alleles.**
  VCF 4.3 writes `*` for a base removed by a deletion on the other copy of
  the chromosome. Allelix treated it as an ordinary allele, so a report
  could show a genotype like `A/*` and `export plink` could write `*` into
  the `.bim` file. `*` is now handled like other symbolic alleles: a call
  that includes it is skipped and the character never appears in a
  genotype. Calls with `*` on one copy are counted, and the report says
  how many were not assessed. `export plink` also skips any allele that is not A, C, G or
  T, such as the `I`/`D` indel codes some genotyping chips use.

- **`allelix db update` now exits with an error when a database fails to
  install.** A failed download or refresh was reported on screen but the
  command still exited successfully, so scripts and scheduled jobs could
  not notice it. The remaining databases are still attempted, the failed
  ones are listed at the end, and the exit status is 1.

- **ClinPGx no longer reports false matches at palindromic SNP sites.**
  When matching a user's genotype to a pharmacogenomic source genotype,
  Allelix tries the complement strand if the direct alleles do not match.
  At palindromic sites — where the forward and reverse strands are
  indistinguishable (A↔T or C↔G pairs) — this step could produce a
  spurious match. The guard that was meant to prevent it did not account
  for multi-allelic sites or for the reference allele. Such sites are now
  detected and complement matching is suppressed at them; direct matches
  still work. A report affected by this issue shows a ClinPGx genotype
  the person does not carry — for example GG where the true call is CC.
  To check an existing report, re-run the same raw file with this release
  and compare: `allelix analyze RAW_FILE --diff OLD_REPORT.json`. A
  ClinPGx row listed as removed whose genotype is the complement of the
  retained row was affected. The same issue affected VCF and gVCF input
  in earlier releases wherever the reference allele was available.

### Added

- **Population frequency, missense and deleteriousness scores for GRCh37
  chip files.** A new pre-built table, derived from dbSNP, records where
  each rsID that Allelix can report sits on both GRCh37 and GRCh38 and how
  its alleles correspond. With it, a GRCh37 genotype (every major consumer
  chip format) is translated to the matching GRCh38 allele and receives
  gnomAD, AlphaMissense and CADD values again. Where the reference strand
  was flipped between assemblies the translation is explicit. Sites whose
  strand cannot be told from sequence (A/T, C/G, or sites where all four
  bases are known alleles) receive no values rather than a guess. The
  same table supplies the reference allele for chip data on both builds,
  so strand-aware matching now works on GRCh37 chips. It also recovers
  rsIDs for VCF files without rsIDs on GRCh37 as well as GRCh38. It is downloaded by `allelix db update`. Without it,
  GRCh37 analyses behave as before.

- **Plausibility flag for genotype calls.** Each finding is now compared
  with how common the allele is in gnomAD. A homozygous call for an allele
  rarer than 0.01% is marked *implausible*: two copies of an allele that
  rare are usually a genotyping error. A heterozygous call for an allele
  that gnomAD lists with no carriers is marked *suspicious*. Flagged rows
  are never hidden. The HTML report shows a warning icon next to the
  genotype with the reason, and the terminal table marks the genotype.
  Both thresholds can be changed with
  `allelix config set plausibility.homozygous_max_af` and
  `plausibility.heterozygous_max_af`. Calls without a gnomAD frequency are
  not assessed. Neither are homozygous X calls (a male's single X reads as
  homozygous) or Y and mitochondrial calls.

- **MyHeritage low-pass whole-genome exports are labelled as imputed.**
  MyHeritage kits processed from January 2026 are sequenced at about 2x
  coverage, and most genotypes are filled in by imputation, which is
  least reliable for rare variants. Allelix now reads the file's
  `##method` line and states the provenance in every report: "low-pass
  whole genome sequencing (about 2x coverage); genotypes are imputed",
  with advice to confirm any finding by direct genotyping. The HTML
  report's array caveat is replaced by an imputation caveat for these
  files. JSON reports carry it as `input.method`. The file's
  `##reference=build37` line is now used as the header build.

### Changed

- **JSON report schema version 7.** Every annotation carries a
  `plausibility` field: `plausible`, `suspicious`, `implausible`, or empty
  when not assessed. The `input` block gains `gvcf` and `vcf_quality`
  (rows excluded by the caller's filters, calls demoted for low genotype
  quality, and spanning-deletion calls not assessed; `null` for non-VCF
  input), and each `high_value_no_calls` entry gains a `reason`. The
  `input` block also gains `method`, the genotyping method a file
  declares (`null` when none). `--diff`
  accepts reports from earlier schema versions as before.

- **Annotated VCF output schema 0.2.0.** Rows the variant caller rejected
  are copied unchanged without Allelix fields. With `--keep-filtered`
  they are annotated and carry the new `ALLELIX_FILTER` field.

- **Every enrichment lookup is made by position and allele.** gnomAD,
  AlphaMissense and CADD values are looked up at the variant's GRCh38
  chromosome, position, reference and alternate allele, never by rsID.
  A GWAS association without a reported risk allele now shows the values
  for the single non-reference allele the person carries, and none when
  they carry two.

### Documentation

- Database sizes in the documentation have been updated to match
  current builds.

- New per-provider format guide: Sequencing.com's three deliverables
  and which one to use, MyHeritage's array and low-pass sequencing
  exports, and MyHappyGenes' build label.

- New architecture decision record: rsID recovery for rsID-less input
  is one shared step used by every command, and build detection matches
  well-known SNPs by position when the rsID is missing.

- New architecture decision record for the VCF quality gate: how
  caller filters, genotype quality, low-quality reference blocks and
  spanning deletions are handled, and how high-value positions are
  classified on gVCF and plain VCF input.

- New architecture decision record for the plausibility flag: which calls
  are assessed, the default thresholds, and why a missing gnomAD
  frequency is never treated as zero.

- New architecture decision record for the cross-build rsID identity
  table and position-keyed enrichment.

- New architecture decision record: ClinVar rows are identified by
  locus and alleles, and every annotation is checked at construction
  never to name an allele the person does not carry, and (for ClinVar,
  ClinPGx and SNPedia) never to name the reference allele. The GWAS
  Catalog allele is documented as the study's effect allele, which can
  be the reference allele.

### Infrastructure

- **Chip / gVCF parity harness.** A new test suite writes the same
  synthetic genotypes at the same positions as a consumer-array file
  and as an rsID-less GRCh38 gVCF, runs both through the analysis
  pipeline against the same ClinVar data (with and without gnomAD
  installed), and requires identical annotations: same alleles,
  significance, magnitude and frequency. Each case is also checked
  against the ClinVar rows the genotype actually carries, so two
  paths that agree on an empty result cannot pass. Disagreements on
  the current code are recorded as expected failures that the
  corresponding fix must remove.

## [2.3.1] - 2026-07-31

### Fixed

- **Cross-build enrichment is now fail-closed.** gnomAD v4.1 exome
  frequencies and coordinate resolution, AlphaMissense, and CADD are
  GRCh38-only in the current cache model. GRCh37 and GRCh36 analyses
  now skip all three sources and omit their provenance instead of
  allowing an rsID or coordinate collision to attach a valid GRCh38
  value to the wrong assembly or allele. ClinVar continues to dispatch
  to its build-specific GRCh37 or GRCh38 cache.

- **GWAS odds ratios are no longer confused with beta coefficients.**
  The effect-size magnitude modifier is applied only when the
  confidence-interval field positively identifies a bare odds ratio.
  Beta coefficients, unit-bearing values, missing or ambiguous effect
  types, and non-positive effects now receive no modifier. On the
  user1190 genome with the 2026-07-23 GWAS Catalog snapshot, this
  changes the default-filter GWAS output from 94 rows to 11; counts
  drift with future database snapshots, while the effect-type rule is
  the stable regression invariant. Beta-specific thresholds remain
  deliberately deferred until effect units can be modeled safely.

- **Position-based build detection requires meaningful evidence.** A
  detected build is authoritative only with at least three concordant
  seed-SNP matches and an 80% winning ratio. A 4/5 signal is confident;
  3/4 and 2/2 signals are tentative. The CLI and HTML report distinguish
  tentative evidence from a real detection even when the fallback build
  happens to be the same value. The existing GRCh36 fail-safe remains:
  a tentative GRCh36 winner stays on GRCh36 so no modern position cache
  is queried, while the diagnostics still label the evidence tentative.

- **SNPedia genotype matching honors source orientation.** Cache parser
  format 7 carries SNPedia's plus, minus, or unknown orientation into
  runtime matching. Known minus-strand rows are complemented to forward
  orientation. Unknown-orientation palindromic A/T and C/G genotypes
  abstain rather than risk a false match. Rebuild the SNPedia cache with
  `allelix db update` after upgrading. GWAS risk-allele strand
  normalization remains deferred pending a forward-strand data model.

- **Slow-suite ClinVar drift checks can no longer pass by skipping.** The
  workflow downloads the ClinVar cache explicitly, requires it in CI,
  selects the drift test by its exact node ID, emits JUnit XML, and fails
  unless the named drift test appears as executed and passed.

- **Annotated VCF output declares every emitted contig.** Before writing,
  Allelix scans the input records and synthesizes a missing
  `##contig=<ID=...>` declaration for each raw CHROM value. Existing
  declarations and record chromosome names remain unchanged, including
  mixed bare and `chr`-prefixed names, so `bcftools view` no longer warns
  about undeclared contigs created by an incomplete input header.

### Documentation

- Synchronized the README, security policy, example-report notes,
  GRCh36 guidance, test-data catalog, edge-case expectations, and full
  release protocol with v2.3.1 behavior. The real-data wrong-allele gate
  now requires a non-zero GRCh38 direct-match denominator and separately
  asserts that GRCh37 enrichment is zero. The protocol covers all nine
  parser formats and treats any release-gate skip as a failure.

## [2.3.0] - 2026-07-02

### Added

- **Annotated VCF output (`allelix analyze --vcf-out PATH`).** New
  output mode alongside the HTML and JSON reports: analyze a VCF/gVCF
  input as usual, and write the annotated form back out with
  annotations stamped as INFO fields and full provenance in the
  header. Positions Allelix as a lightweight, deterministic drop-in
  alternative to heavier annotation tools that pipeline users have
  been asking for. Every input header line is preserved verbatim;
  the injected block declares `##INFO=<>` fields for every source
  that produced an annotation, a provenance block (`##ALLELIX_Version`,
  `##ALLELIX_VCF_SchemaVersion`, `##ALLELIX_License`,
  `##ALLELIX_RunDate`, `##ALLELIX_Build`), per-database version pins
  (`##ALLELIX_Database=<Name=,Version=>`), and per-source attribution
  (`##ALLELIX_Attribution=<Name=,License=,URL=>`) mirroring the
  `license_attributions` block from the HTML/JSON renderers.
  Multi-allelic sites use per-allele `Number=A` semantics with
  comma-joined values in ALT-column order. Recovered rsIDs
  (`.` → `rs<N>` via ClinVar or gnomAD) land in the ID column with
  the pre-recovery value stashed under `ALLELIX_ORIGINAL_ID`. The
  output contract is versioned separately from the JSON report
  contract (`VCF_SCHEMA_VERSION="0.1.0"`, experimental — see
  ADR-0036); INFO-field churn during the 0.x window does not
  require the JSON `SCHEMA_VERSION` bump. Non-VCF input errors
  early; existing output paths are refused (no clobber, no `--force`
  in v2.3.0). Gzip-compressed input (`.vcf.gz`) is transparently
  handled, matching the parser. BCF binary output, `.vcf.gz`
  compressed output, consequence annotation coupling, and full
  VCF 4.3 §1.6 percent-encoding of INFO values are explicit non-goals
  for v2.3.0 and tracked as follow-ups.

- **README:** hosted-demo link (`[Demo](https://analyze.allelix.io)`)
  in the header navigation strip so a first-time reader can try the
  tool without a local install before deciding to `pip install`.

### Documentation

- **ADR-0036: Annotated VCF Output Contract.** Captures the two-pass
  writer architecture (pipeline builds the annotation dict; writer
  re-reads the input and stamps INFO — parser and `AnalysisResult`
  untouched), the per-allele `Number=A` decision (silent multi-allelic
  data loss the writer cannot afford), the schema-version namespace
  decision (independent from the JSON `SCHEMA_VERSION` governed by
  ADR-0033), and the explicit v2.3.0 non-goals.

- **`CONTRIBUTING.md` § Output-contract schema bumps.** New subsection
  under Release-content discipline cross-referencing both output
  contracts (JSON → ADR-0033, annotated VCF → ADR-0036) with the hard
  rule that any PR changing an output-contract field shape bumps the
  corresponding schema-version constant in the same PR. While
  `VCF_SCHEMA_VERSION` sits at `0.x`, INFO-shape churn is approved
  without ceremony; graduation to `1.0.0` waits for a follow-up ADR
  after real user feedback.

- **`test_data/FULL_TEST_PROTOCOL.md` § 20.** Manual verification
  battery for the annotated VCF output — smoke test, multi-allelic
  per-allele semantics, gVCF reference-block passthrough, rsID
  recovery visibility, multi-sample preservation, `bcftools view`
  round-trip, reproducibility from headers alone, non-VCF rejection,
  existing-output refusal, and gzip input handling.

### Infrastructure

- **CI: `bcftools view` gate.** The `lint-and-test` job now installs
  `bcftools` (~5 MB apt package) and runs three subprocess tests
  that validate `--vcf-out` output round-trips through a real VCF
  parser: (1) `bcftools view` accepts the annotated file with no
  writer-attributable warnings, (2) `bcftools view -h` lists every
  declared `##INFO=<ID=ALLELIX_*>` line plus the provenance block,
  (3) `bcftools query %INFO/ALLELIX_CLINVAR` extracts the value the
  writer stamped on a known-hit row (proves bcftools honors our
  `Number=A` declaration end-to-end). Any regression in the INFO
  syntax, provenance header format, or record-line structure now
  fails at PR time rather than at real-user integration.

## [2.2.3] - 2026-06-25

### Fixed

- **23andMe parser fails on newer v5 exports.** The parser anchored
  detection on the historical canonical first-line banner
  `"# This data file generated by 23andMe"`. A real user's June
  2026 v5 export uses an edited prose form
  `"# This data file is generated by 23andMe."` — with the extra
  word "is" and a trailing period — and the strict `startswith`
  check failed silently, falling the file through all parsers to
  "Unsupported format." The parser now matches against a tuple of
  both known canonical banner-line variants and is structured so
  that future vendor banner edits get a new entry without further
  detection-logic changes. Committed fixture
  `mock_23andme_v5_is_generated.txt` pins the v5 variant; both
  variants are tested in CI. The deliberate "no loose substring
  matching" policy (which keeps transcoded fixtures and
  third-party converter output from being misclassified as
  23andMe) is preserved. Banner-less 23andMe-shape files —
  WGSExtract output, Genetic Lifehacks converter output,
  privacy-conscious users who delete the comment block — remain
  out of scope for this parser by design and are tracked for a
  separate generic-TSV parser landing in v2.3.

- **gnomAD position resolver bound 1800 SQLite parameters per
  chunk.** The position-keyed resolver added with the v2.2.2
  rsID-recovery fix batched at the full `_BULK_BATCH_SIZE`
  constant (900) while binding two parameters per input pair
  (chrom + pos), so each chunk issued a 1800-param query. That
  was the only batched query in the codebase exceeding the
  ~999-param convention every other batched annotator
  deliberately respects. Modern SQLite (≥ 3.32, shipped with
  every Python version the project supports) raised the variable
  limit well past 1800, so the deviation was harmless at
  runtime — but it was a foot-gun for any future change that
  pushes `_BULK_BATCH_SIZE` higher or moves to a system
  `libsqlite3` of unknown vintage. Chunk size is now
  `_BULK_BATCH_SIZE // 2`, matching the established pattern in
  the annotator's own `bulk_lookup_by_alt` (identical 2-param
  shape, same divisor; `bulk_lookup_by_position` uses `// 4` for
  its 4-param shape, same convention scaled to the parameter
  count) and the other annotators' `_BATCH_CHUNK = 500`
  convention. Net runtime cost is negligible — twice as many
  smaller index-seek queries net out against half as many
  larger ones now that `idx_gnomad_position` covers both
  shapes — and three sampled measurements on a 5.2M-variant
  rsID-less DeepVariant gVCF showed user-CPU time identical
  across runs (within 1%). Convention compliance + removed
  foot-gun.

### Documentation

- **README "planned for v2.1" line on the data-lifecycle
  subcommand.** `allelix db clean` shipped in v2.2.0 and is
  documented as working elsewhere in the same README; the
  "planned for v2.1" sentence was a relic of the v2.0 README
  that never got refreshed at the v2.2.0 ship. Replaced with an
  accurate one-line note that the cache-side lifecycle ships
  and that managing report files written via `--output` remains
  the user's responsibility.

## [2.2.2] - 2026-06-20

### Fixed

- **rsID-less VCF coverage on `analyze --filter-file` and `extract`.**
  VCFs that don't stamp dbSNP identifiers in the ID column —
  DeepVariant gVCFs, Sequencing.com / DRAGEN output, GIAB GRCh38
  benchmarks — produced wellness-panel reports with most rsIDs
  classified as "not in file" even when the sample carried the
  underlying variants. The ClinVar position resolver already
  recovered an rsID for variants in ClinVar's curated subset, but
  most pharmacogenomic, GWAS-only, and SNPedia-only rsIDs aren't
  in ClinVar and stayed unresolved. The pipeline now runs a
  second-pass position resolver against gnomAD's full dbSNP
  coordinate map after ClinVar's pass, recovering rsIDs across
  the much broader gnomAD universe. `extract` was rewritten
  symmetrically — both the tabix fast-path and the sequential
  fallback now match by `(chrom, pos, ref, alt)` instead of an
  rsID equality check that always failed on rsID-less inputs.
  Real-world measurement on a DeepVariant HG002 gVCF with a
  22-rsID wellness panel: panel coverage went from 3/20 found
  to 8/22 found, with the rest reading honestly as "reference at
  this position or no variant call" instead of the misleading
  "not on chip" framing.

- **VCF parser tolerance for UTF-8 BOM and leading preamble lines.**
  `_open_vcf` now opens files with the `utf-8-sig` codec, which
  strips a leading UTF-8 BOM (`\xef\xbb\xbf`) transparently
  without affecting BOM-less inputs. `can_parse` and the header
  reader scan the first 10 non-blank lines for `##fileformat=VCF`
  rather than requiring it on line 1 — robust against tool
  banners, comment-line preambles, and other upstream
  pipeline quirks that previously caused silent format-detection
  failures. Four committed fixtures and 10 new tests cover the
  full matrix of {plain, gzipped} × {BOM only, preamble only,
  BOM + preamble}, plus a pathological-input rejection pin.

- **HTML panel-coverage wording on VCF / gVCF inputs.** State-3
  panel rsIDs (requested by the user, not present in the file)
  rendered as "not genotyped on this chip/platform" regardless
  of input format. For array inputs that's accurate; for VCF /
  gVCF inputs there is no chip, and the genuine ref/ref case
  (e.g., HG002 at rs429358 — no variant call because the sample
  is ref/ref) read as a platform-coverage failure. The label is
  now format-aware: array inputs keep the historical wording;
  VCF / gVCF inputs surface "reference at this position or no
  variant call." Both the per-row state-3 label and the
  summary-header count are updated.

- **Blank `Build` column in `allelix stats` on caller-output VCFs.**
  Caller outputs (DeepVariant, DRAGEN, GIAB GRCh38 benchmark)
  don't write `##reference=` in their headers, so the VCF
  parser returned an empty build string and the stats summary
  rendered a blank cell — reading as broken even though stats
  deliberately does not run build detection (that's an `analyze`
  job per the protocol). The cell now renders an em dash, plus
  a short hint pointing the user to `analyze` when chr-prefixed
  contigs suggest a likely build.

### Performance

- **gnomAD position-index lookup ~5x faster on WGS analyze.**
  The new rsID-less position resolver hit a missing covering
  index — `gnomad_frequencies` carried only an `rsid` index,
  so every position-keyed query was a full scan against 57
  million rows. Added `idx_gnomad_position` over `(chrom, pos)`
  to the canonical schema, with a lazy migration in
  `GnomadAnnotator._connection` so existing 2.7 GB caches don't
  need re-download: `CREATE INDEX IF NOT EXISTS` runs once on
  first use against a pre-v2.2.2 cache (a few minutes against
  the 57M rows), then stays a microsecond-cost no-op forever
  after. The per-batch SQL was also collapsed from one statement
  per chromosome to one combined row-value `IN` per batch — on
  WGS input each batch touches ~24 chromosomes, so the original
  shape was a 24× query multiplier the new index alone could
  not amortize. Real-world measurement on a 5.2 M-variant
  DeepVariant HG002 gVCF: ~165s → ~133s.

## [2.2.1] - 2026-06-18

### Fixed

- **VCF parser robustness on real-world inputs.** `can_parse` only
  caught `OSError`, so feeding the parser a binary file labeled
  `.vcf.gz` (e.g. a mislabeled tabix index that downloaded with the
  wrong filename) leaked a raw Python `UnicodeDecodeError` traceback
  instead of returning `False`/"No parser recognized." A real VCF
  carrying a non-UTF-8 byte in a tool-emitted header field
  (`##source=...`, `##CL=...` — extremely common in GATK / DeepVariant
  / bcftools output) hit the same fate. The parser now opens both
  gzip and plain paths with `errors="replace"` so real-world headers
  with stray latin-1 bytes parse cleanly (the replacement only
  touches header strings — VCF data lines are pure ASCII by spec so
  genotype fidelity is unaffected), and `can_parse` catches
  `UnicodeDecodeError` and `EOFError` (corrupt / truncated gzip)
  alongside `OSError` so format sniffing never leaks a traceback.
  Closed the underlying CI coverage gap that let this slip: four
  committed fixtures cover the full `{.vcf, .vcf.gz} × {UTF-8 only,
  non-UTF-8 header byte}` matrix, all under 400 bytes per fixture so
  every CI matrix leg exercises them. `.vcf.gz` parsing was already
  working in v2.2.0 — this is robustness on the same flagship path,
  not a missing feature.
- **ClinVar per-SCV duplicate rows on every report (regression from
  v2.2.0).** The per-SCV TSV migration correctly de-normalized the
  cache to one row per submission, but the annotator was mapping
  cache rows 1:1 to `Annotation` objects with no dedup. A variant
  with N identical SCV submissions rendered N identical rows in
  the terminal table, in the HTML report, and in the JSON
  `annotations` array. The annotator now collapses rows identical
  on (rsid, significance, condition, gene) into a single
  `Annotation`, preserving the strongest `review_status` (per
  ClinVar's documented star-rating order) and the union of
  `references` across the collapsed group. Distinct (significance,
  condition) pairs on the same variant survive — the per-SCV
  separation that the v2.2.0 cache redesign delivered to retire
  the Frankenstein-pair hazard is preserved, not undone.
- **Panel coverage accounting closes against the rendered report.**
  When `--filter-file` supplied a panel, `panel_coverage.found`
  reported the number of panel rsIDs in the input file but did
  not reconcile against the magnitude floor the renderer applied.
  Panel rsIDs whose every annotation scored below the analyze
  display floor were counted as "found" and then filtered out of
  every rendered surface — a flat lie. Now the panel annotated
  set is derived from the post-filter annotation list the
  renderer actually produces. Below-floor panel rsIDs drop into
  `no_findings` honestly; `requested = len(missing) +
  len(no_findings) + count(distinct displayed-annotation rsIDs)`
  becomes the pinned arithmetic invariant. `found` still means
  "in input file" (genotyped count) — same v6 contract; no
  schema change. The genotype-aware redesign for `panel_coverage`
  is reserved for v2.3.
- **FamFinder haploid MT/Y calls now parse as hemizygous.** The
  parser docstring claimed support for haploid rows with missing
  ALLELE2; the implementation required exactly 5 tab columns, so
  a real 4-column haploid line was warned and skipped, dropping
  every haploid MT/Y call. A 5-column line with an empty trailing
  ALLELE2 was normalized to `-`, surfacing as a half no-call that
  the analyze pipeline abstained on. Both shapes now parse as
  hemizygous `Variant` objects where `allele1 == allele2 == the
  called base` — consistent with the haploid convention used by
  the 23andMe parser. Empty ALLELE1 (lead allele absent) is the
  symmetry guard: stays as `-`, no haploid promotion.
- **ClinVar non-classification curatorial terms no longer surface
  as findings.** ClinVar's per-SCV TSVs expose three curatorial
  values — `other`, `association`, `association not found` —
  that the old summarized VCF data hid. They aren't
  classifications, aren't in the magnitude scoring table, and
  fell to the 5.0 default, surfacing as fake findings on real
  rsIDs. Added to the loader's ingest skip set alongside the
  empty / `-` / `not specified` / `no classification provided`
  placeholders that v2.2.0's earlier follow-up retired.
  Deliberately not added to the magnitude table (they aren't
  classifications) and not added to the benign set (they aren't
  benign).
- **`db clean` safety guard rejects generic project directories.**
  The pre-deletion guard accepted any directory containing a
  `config.toml` as an allelix data dir. `config.toml` is a
  generic filename used by Rust (Cargo), Hugo, Rocket, and
  others — a typo'd `--data-dir ~/my-rust-project` slipped past
  the guard simply because the project root had random content
  alongside a `config.toml`. Removed from the recognition prefix
  set; the remaining prefixes (`clinvar`, `gnomad`, `pharmgkb`,
  …) are all allelix-specific. The fresh-install case (an empty
  data dir with only `config.toml`) is unaffected — the empty
  cache list short-circuits to "Nothing to clean" before the
  guard runs.
- **CHANGELOG no longer contains dead issue links.** Bare
  issue-number references in the changelog auto-linked on GitHub
  to pages that do not exist in this repository, so every shipped
  reference rendered as a dead link. The previously documented
  `### Internal` section carried development notes users do not
  need. The released changelog content has been rewritten to stand
  on its own without issue numbers and without an Internal section.
  A structural guard (`tests/test_changelog_sanitize.py`) now runs
  on every fast-tier pytest pass and fails the build if either
  form reappears.

### Changed

- **Read-time placeholder skip set widened.** The ClinVar TSV
  loader's placeholder constant now lists seven values: empty,
  `-`, `not specified`, `no classification provided`, plus the
  three new non-classification curatorial terms. The protocol's
  significance-sentinel ship-gate scans the same seven values;
  loader and gate share a single source of truth.

## [2.2.0] - 2026-06-18

### Changed

- **ClinVar source switched from VCF to a per-SCV tab-delimited pair.**
  The loader now ingests `variant_summary.txt.gz` +
  `submission_summary.txt.gz` from NCBI's tab-delimited dump, emitting
  one cache record per (variant, SCV submission) instead of one combined
  row per variant. Structurally retires the CLNSIG|CLNDN
  Frankenstein-pair hazard of the VCF loader: each row carries one
  submitter's exact (significance, condition) pair, so there is no
  cross-product pairing risk. Record count roughly doubles (~5.79M →
  ~9.09M GRCh38) — expected from per-SCV granularity.
  `CLINVAR_INTERPRETER_VERSION` bumps `2 → 4` so existing caches
  auto-invalidate on the next `db update` (iv:3 was the in-flight bump
  for the initial loader; iv:4 forces the placeholder-skip fix described
  under Fixed). Old VCF loader + URL constants removed in the cluster's
  final stage.
- **JSON report `schema_version` bumped 5 → 6.** Adds
  the new top-level `panel_coverage` object when `--filter-file`
  supplies a rsID panel. Additive only — prior-version consumers still
  parse the new output. Diff loader's `_SUPPORTED_SCHEMA_VERSIONS` grows
  monotonically to include "6".
- **Focused subcommand `--min-magnitude` default lowered 5.0 → 3.0.**
  `allelix methylation` and `allelix pharmacogenomics`
  previously inherited the same magnitude floor as `allelix analyze`
  (5.0). That floor protects the general report from noise but filters
  out the exact signal a focused subcommand is supposed to surface —
  ClinPGx LoE 3 hits (mag 4.0) on COMT / MTHFR / MTR / CBS for
  `methylation`, similar pharmacogenomic hits for `pharmacogenomics`.
  New default 3.0 keeps LoE 3 in, keeps ClinVar benign (mag 1.0) out.
  `allelix analyze` default unchanged at 5.0.
- **Permission-resolution dedup.** `AllelixConfig.is_enabled`
  and the `config show` table-cell renderer used to each fetch the
  annotator class and call `check_permission()` inline. Both now
  delegate to a single `AllelixConfig.permission_for()` that returns
  the three-state `Permission | None`. No behavior change; eliminates
  drift surface.
- **Streaming ClinVar reference-lookup.**
  `ClinVarAnnotator._load_ref_lookup` switches from `fetchall()` to
  direct cursor iteration; peak memory drops from "full result list +
  output dict" to "one row + output dict". `iter_gwas_records`
  intentionally retained as buffering (dedup-on-min-p-value requires
  full input visibility); docstring now explicitly acknowledges this
  instead of claiming streaming.

### Added

- **Panel coverage in custom-panel reports.** When `--filter-file`
  supplies a rsID panel, the analyze command now distinguishes three
  states per requested rsID instead of silently dropping two of them:
  state 1 (genotyped + has annotations) renders normally in the main
  table; state 2 (genotyped + no annotations) and state 3 (not in input
  file) surface in a "Panel coverage" section. Terminal output emits a
  one-line warning when state 3 is non-empty
  (`⚠ N/M panel variants not found in input: rs1, rs2, rs3, …`); JSON
  payload carries a `panel_coverage` object
  (`{requested, found, missing, no_findings}`); HTML report shows the
  section between the high-value no-call banner and the education
  section with a per-rsID table. New `AnalysisResult.panel_coverage()`
  returns a typed dict — or `None` when no panel was supplied (zero
  behavior change on the default analyze run).
- **FTDNA FamFinder parser.** Third FTDNA file shape —
  tab-delimited with separate `ALLELE1` / `ALLELE2` columns instead of
  a concatenated `RESULT` (the FTDNA and FTDNA-Illumina-raw siblings).
  Detected by the `famfinder` substring in the file header plus the
  canonical 5-column header; requiring both signals avoids false
  positives against any future tab-delimited 5-column format that
  isn't a FamFinder file. Empty / `-` cells in either allele column
  map to the no-call convention (`-`). 25 new tests; synthetic
  fixture at `tests/fixtures/mock_ftdna_famfinder.txt`.
- **`allelix db clean` and `allelix db path` CLI subcommands.**
  `db clean` removes downloaded reference database caches (default
  `~/.local/share/allelix/`) with a size report and a confirmation
  prompt; `--keep-cadd` preserves the CADD opt-in cache (~5.8 GB,
  expensive to re-download); `--dry-run` previews; `--yes` skips the
  prompt for scripted use. `db clean` runs a pre-deletion safety
  guard: if the target directory contains content but no recognized
  allelix cache files, it refuses with a clear error even on the
  `--yes` path; `--force` bypasses for non-standard installs.
  `db path` prints the resolved data directory for shell capture
  (`$(allelix db path)`); `--check` verifies the directory is writable.
  Replaces the manual `rm -rf` recipe in the README "Managing your
  data" section.
- **Runtime degradation nudges.** Two once-per-run info
  lines surface previously-silent post-pipeline conditions: when array
  input runs without a usable gnomAD ref-resolution path (`--no-gnomad`
  or no gnomAD cache), every variant lands `ref=None` and the
  strand-flip path silently falls back to v2.0.x direct-match-only;
  analyze now emits a dim
  "strand-aware matching inactive — no reference context" notice.
  When the effective build is GRCh38 and a high fraction of
  input variants are rsID-less (the variant-caller `ID=.` case),
  analyze emits a yellow undercount warning recommending GRCh37 if
  available, pointing at the planned dbSNP-resolution work as the
  real long-term fix. Array inputs are silent on the GRCh38 nudge;
  VCF inputs are silent on the strand-aware nudge.
- **Methylation panel + magnitude floor retune.** Curated
  gene panel for the `methylation` subcommand with annotated
  citations; magnitude floor for the focused-subcommand path retuned
  to 3.0 (see Changed).
- **Documentation cross-links from README + CONTRIBUTING to
  `test_data/FULL_TEST_PROTOCOL.md`.** The full real-data
  release-validation battery is now linked from both the README's
  Development section and CONTRIBUTING's testing section, instead of
  being mentioned only at line ~461 of CONTRIBUTING.
- **`FULL_TEST_PROTOCOL.md §7b` ClinVar significance-sentinel
  ship-gate.** Post-`db update` SQL scan verifying zero rows match
  `_CLINVAR_PLACEHOLDER_CLNSIGS` in either build cache, with
  `LOWER()` for case-insensitive symmetry with the loader. Named
  ship-gate ("0 rows expected; any row → do not tag") and a
  forward-pointer to the CLNSIG vocabulary drift test as the durable
  allowlist backstop.

### Fixed

- **TSV loader filters ClinVar placeholder CLNSIGs.** ClinVar's
  `submission_summary.txt.gz` carries several placeholder values that
  mean "no clinical significance recorded." Values not in
  `_CLNSIG_MAGNITUDE` fell through to the dict's 5.0 default — equal to
  the analyze display floor — and surfaced as bogus annotations on real
  rsIDs (rs137854557 → Meningioma was the canary). The loader now
  filters the exact set that would otherwise default
  (`_CLINVAR_PLACEHOLDER_CLNSIGS = {"", "-", "not specified", "no
  classification provided"}`), with case-insensitive membership testing
  to mirror `_magnitude()`'s `_normalize_clnsig` normalization.
  `not provided` and `no_classification_for_the_single_variant` are
  intentionally kept (both map to 2.0 — safe below the floor; dropping
  them would lose real submitter records). The drift snapshot is updated
  to remove `no_classification_provided` so the live CLNSIG drift guard
  remains honest about what's filtered. `CLINVAR_INTERPRETER_VERSION`
  bumps `3 → 4` to force reingest of any iv:3 cache built by the pre-fix
  loader.
- **TSV loader stamps a real `version` value.** `db status` previously
  showed "ClinVar version: None" after the per-SCV rebuild;
  `stamp_remote_signal` inserts `version=NULL` by design (it's a
  freshness-only stamper meant to upsert remote_signal onto baked-in
  HuggingFace cache rows), and the TSV loader had no pre-existing row.
  Now stamps build date `YYYY-MM-DD` UTC — the practical equivalent of
  the old VCF loader's `##fileDate` header.
- **Terminal output emits the regulatory disclaimer.** HTML and JSON
  have carried the ADR-0003 "not medical advice / not a diagnosis"
  notice since v2.0.x; the terminal — the **default** analyze output —
  did not. Anyone running `analyze` without `--output` got a
  clinical-looking table with no disclaimer.
- **`FULL_TEST_PROTOCOL.md §7` per-SCV condition contract.** The
  flagship sanity block was still pinning the pre-TSV-loader VCF-era
  semicolon-joined condition shape (`rs1800896 → "Leprosy ... ;
  Hepatitis C ..."`). Doc rewritten to pin the per-SCV row shape
  (multi-row, never combined) and explain why the join-pin is
  structurally retired.
- **`FULL_TEST_PROTOCOL.md §15` stale exome-VCF expectation.**
  `unsupported_23andme_exome_vcf.txt` is a valid VCFv4.1 file and is
  parsed by the VCF parser added after the fixture was named. Doc split
  into "still unsupported" (decodeme) and "previously-unsupported /
  now-supported regression pin" (the 23andMe exome VCF, kept under its
  historical filename for git-mtime continuity).

## [2.1.2] - 2026-06-16

### Changed

- **`Annotator.__del__` safety net removed.** v2.1.1 finished the
  production-side audit; v2.1.2 finishes the test-side refactor and
  drops the destructor entirely. The pipeline's `contextlib.ExitStack`
  plus the CLI utility paths' explicit `try/finally close()` are now the
  only cleanup paths — deterministic, not GC-timed. `import contextlib`
  removed from `allelix/annotators/base.py` as a follow-on (no longer
  used).
- **Test-suite annotator instantiation routed through a `cm_stack`
  fixture.** 20 test methods in `tests/reports/test_pipeline.py`
  previously constructed annotators outside any context manager and
  relied on `__del__` for SQLite cleanup; that pattern surfaced as
  `PytestUnraisableExceptionWarning` once `__del__` was removed. New
  per-test `cm_stack` fixture (a `contextlib.ExitStack`) wraps every
  `enter_context(Annotator(...))` call so cleanup is deterministic at
  test-function teardown. Same shape as production's pipeline-ExitStack
  pattern.
- **PLINK gnomAD coordinate resolution moved out of the CLI.**
  `allelix/cli.py` had become a god-module with PLINK business
  logic — gnomAD coordinate lookup + multi-allelic disambiguation —
  sitting in the presentation layer. The CLI was already split into the
  `allelix/cli/` package before; the remaining business-logic leak
  was the PLINK ref/alt resolution loop in `allelix/cli/utility.py`. New
  `allelix.exporters.plink.resolve_ref_alt_via_gnomad()` owns the
  GnomadAnnotator lifecycle and the multi-allelic match-then-fallback
  disambiguation; the CLI handler now just calls it inside a small
  try/except for the user-facing warning. No behavior change — same
  ref/alt selection rules, same warning text on failure / unavailability.

### CI / Tooling

- **mypy added to the CI gate.** Backs the `py.typed` promise with
  a real type checker. Baseline is intentionally permissive (no strict
  mode, no `disallow_untyped_defs`) so the gate landed in CI immediately
  instead of negotiating dozens of remediation PRs; ratchet later by
  tightening individual flags. `mypy>=1.11` and `types-PyYAML>=6.0`
  added to dev extras; `[tool.mypy]` config in `pyproject.toml`; pysam
  and mwparserfromhell marked `ignore_missing_imports` (neither ships
  stubs). `.github/workflows/ci.yml` runs `mypy allelix/` after the
  ruff format check. Real-bug fixes that fell out of the baseline pass:
  `_pipeline.py` build-detection branches now assert non-None invariants
  before they `replace(build=…)`, and a `clinvar.py` loop reused the
  name `rows` for both a `dict.values()` iteration and a `dict.get()`
  result (renamed to disambiguate).

## [2.1.1] - 2026-06-16

Patch release. Focus: documentation drift —
SECURITY.md staleness, missing ADR index entry, sparse operational
docs, plus one structural refactor that eliminates the `is_must_include`
leak-by-omission risk at the model boundary. No user-facing behavior
changes; existing v2.1.0 caches and reports stay current.

### Security

- **SECURITY.md "Supported versions" table was stuck at 1.8.x while
  v1.9.0, v2.0.0, v2.0.1, v2.0.2, and v2.1.0 all shipped (caught
  during v2.1.1 prep).** Table updated to `2.1.x ✓` / `< 2.1 ✗`
  to match the "current minor release receives security fixes" policy
  line above it. Out-of-scope third-party data list also corrected:
  PharmGKB → `ClinPGx (formerly PharmGKB)` to match v2.0+ display-name
  conventions, plus AlphaMissense / CADD added (they were shipped after
  SECURITY.md was last touched). No reporting URL change.

### Changed

- **Centralized internal-field exclusion for public Annotation
  serialization.** Three production sites previously
  duplicated the inline `asdict(a) minus is_must_include` filter
  (JSON report's main annotations list, JSON diff "new" entries, the
  `diff_annotation_to_dict` change records). A fourth public-output
  path that forgot the filter would have leaked an internal boolean
  into reports without breaking anything obvious at PR time. New
  `allelix.models.annotation_to_public_dict()` helper + a
  module-level `_INTERNAL_ANNOTATION_FIELDS` frozenset are the single
  place to remember the rule. All three production serializers now
  route through the helper. Adding a new internal-only field is one
  line in the frozenset; every existing serializer picks up the
  exclusion automatically. No JSON schema change — same fields
  produced, same values, same shape. Test pin
  `TestAnnotation.test_annotation_to_public_dict_strips_internal_fields`
  asserts the contract.

### Documentation

- **Operational / data-lifecycle docs in the README.** New "Managing
  your data" section under the database-sizes block covers cache
  location (`$XDG_DATA_HOME/allelix`), per-file footprint (~16 GB
  default, ~22 GB with CADD), the disposability rule (caches are
  re-downloadable; safe to delete), the manual cleanup path
  (`rm -rf ~/.local/share/allelix/`), config.toml as the backup target
  (only non-reproducible state in the install), and the full uninstall
  procedure. New "Troubleshooting" section before "Architecture" covers
  six common failure states (3.10 ImportError, partial `db update`
  failure, "not ready" cache after a release upgrade, build mismatch
  warnings, the signal-driven staleness model, output-file sensitivity).
  Python 3.11+ requirement promoted from a buried sentence to the first
  line of Quickstart with the exact ImportError users see on 3.10
  (`cannot import name 'UTC' from 'datetime'` — grep-able). CADD table
  cell tightened from a full-paragraph cell to a one-line summary + a
  new "CADD modes" subsection that holds the cache-vs-full-mode
  explanation. Status banner refreshed from a stale `Latest: v2.0.0` to
  v2.1.0 with the ADR-0035 cluster summary. The `db clean` / `db path`
  CLI tooling for the same topic ships in v2.2.
- **`__del__` removal audit — partial closure.**
  v2.0.2 deferred `__del__` removal to v2.1 pending an audit of every
  callsite. v2.1.1 finishes the audit: zero production paths rely on
  `__del__` for SQLite-connection cleanup. `run_analysis` wires
  every primary + enrichment annotator into `contextlib.ExitStack`
  (v2.0.2 fix); `cli/utility.py:221, 476` use explicit `try /
  finally close()`; the annotators factory returns a list the
  pipeline immediately wraps. The blocker for actually removing
  `__del__` is ~27 test sites in `tests/reports/test_pipeline.py`
  that construct annotators outside a context manager; the
  pipeline closes them via its own ExitStack but the test still
  holds the closed-instance reference, and the sqlite3 Connection
  GC finalizer surfaces a `ResourceWarning` during the next test's
  startup which pytest's `error::PytestUnraisableExceptionWarning`
  filter converts to a hard failure. Verified empirically:
  removing `__del__` produced 3 such failures. `__del__` therefore
  stays as the safety net with an updated docstring documenting
  the partial-closure state; full removal is gated on the test
  refactor and lands in v2.1.2.
- **ADR-0035 added to the ADR index.** ADR-0035 itself (the v2.1.0
  `Variant.ref` and annotation-model design) shipped without updating
  `docs/adr/README.md`'s index. Index entry added with the
  supersedes note that ADR-0035's cluster carve-out partially
  supersedes ADR-0033's strict-per-field schema-version bump policy
  for documented coordinated clusters. Same docs-drift family as
  the SECURITY.md table fix and the stale README version banner —
  single-file additions that skip the index/inventory update step.
  The ADR index update is already step 4 of the "Writing a new ADR"
  checklist in `docs/adr/README.md`.

## [2.1.0] - 2026-06-16

### Added

- **Strand-aware carrier matching (ADR-0035).** ClinVar and ClinPGx
  carrier matching now accepts complement-strand reads of the same
  variant: a het carrier of an ALT scored on the forward strand and the
  same biology read on the coding strand produce identical annotation
  sets. The v2.0.1 multi-allelic safety invariant is preserved —
  strand-flip fires only when `Variant.ref` confirms reverse orientation
  (`variant.ref == complement(source_ref)`) AND the source `(ref, alt)`
  pair is non-palindromic. Without that orientation context
  (`variant.ref is None`) or against palindromic A/T and C/G sites, the
  match falls back to direct membership only — same behavior as v2.0.x.
  Two new helpers in `allelix.utils.allele`:
  `strand_aware_carrier_match` for per-row ClinVar-style `(ref, alt)`
  carrier checks; `strand_aware_genotype_match` for ClinPGx / SNPedia
  per-diploid matches. The reference test case (forward C/C and coding
  G/G both producing zero annotations on a het pathogenic site) is
  pinned via the fixture rsid, plus the het-carrier parity case (forward
  G/A and coding C/T producing the same annotation set) and a case where
  `variant.ref` disagrees with both the source ref and its complement,
  so strand-flip abstains.
- **Array-data `Variant.ref` population (ADR-0035).** The analysis
  pipeline now resolves the forward REF allele for array-format inputs via
  gnomAD's `bulk_resolve_coordinates(rsid)` lookup, mutating `Variant.ref`
  in place at batch-flush time. VCF inputs continue to carry `ref` from the
  REF column. The population step gates on `gnomad.is_ready()` —
  array runs without a current gnomAD cache continue to leave `ref=None`,
  and downstream consumers (SNPedia / ClinPGx alt threading; the new
  strand-aware carrier match) degrade gracefully to direct-match
  only. This is what unlocks allele-specific enrichment on 23andMe /
  AncestryDNA / MyHappyGenes / FTDNA / MyHeritage / LivingDNA inputs that
  previously got `alt=""` because `variant.ref` was None.
- **ClinPGx carrier rule extended to strand-aware (ADR-0035).**
  Previously `WHERE rsid = ? AND genotype = ?` SQL — direct exact-match
  only. Now fetches all rows for the rsid (typically ≤5) and applies
  `strand_aware_genotype_match` per row, so reverse-strand array reads
  match the source diploid when orientation context is available.


### Added

- **`Variant.ref` field (ADR-0035).** New optional field on the
  `Variant` dataclass carrying the reference allele on the forward
  strand. The VCF parser populates `ref` from the VCF REF column; array
  parsers leave it `None` for now (population via the existing
  `resolve_strand` path against gnomAD arrives together with the
  strand-aware carrier matching that consumes it). Required by
  downstream consumers (strand-aware comparison, plausibility flagging,
  future ACMG engine) to identify which side of the user's allele pair
  is the risk allele. `None` means "cannot resolve" — consumers degrade
  gracefully rather than guess.

### Changed

- **JSON `schema_version` bumped 4 → 5.** Enters the ADR-0035 cluster
  contract. The full ADR-0035 field set arriving within v5:
  `variant.ref`, `annotation.trait` / `annotation.p_value` /
  `annotation.phecode` (the structured GWAS fields). `annotation.alt`
  exists since v1.4.1; this release extends its population coverage from
  ClinVar to GWAS / SNPedia / ClinPGx without a schema change. v4
  consumers reading v5 output still succeed — every v4 field is
  unchanged. The diff loader accepts both v4 and v5 baselines.
- **Per-Annotation `alt` threading: GWAS / SNPedia / ClinPGx
  (ADR-0035).** The `Annotation.alt` field exists since v1.4.1 but until
  this release only ClinVar populated it. The v2.0.1 suppress-half
  removed the wrong-allele MAX-fallback for alt-less rows; this release
  closes the loop by routing per-row alt from upstream source data onto
  the annotation so rsID-bearing input VCFs once again get
  allele-specific gnomAD / AlphaMissense / CADD enrichment via the
  existing exact-`(rsid, alt)` lookup path. (a) GWAS Catalog: the
  catalog's per-association risk allele threads as `alt` when present;
  rows without a specified risk allele continue to emit `alt=""` and
  stay capped at magnitude 3.0 per ADR-0024. (b) SNPedia: alt is derived
  from `Variant.ref` + the matched diploid row via
  `derive_alt_from_diploid()` (new helper in `allelix.utils.allele`);
  empty when `Variant.ref` is None (array data without gnomAD-backed ref
  population). (c) ClinPGx: same derivation as SNPedia. The
  position-fallback enrichment path (used when rsID-less input VCFs
  resolve rsIDs on the fly via ClinVar) is unchanged; it was already
  allele-specific and stays that way. No schema-version change —
  `Annotation.alt` is the same field.
- **Structured GWAS fields on `Annotation`: `trait` / `p_value` /
  `phecode`.** v2.0.1's suppress-half merged two inconsistent
  `(PheCode ` delimiter strings into a shared constant; this release
  closes the structural half by promoting the values themselves out of
  the rendered `description` prose into typed `Annotation` fields. The
  GWAS annotator extracts the PheCode suffix (when the upstream trait
  carries one) at the ingestion boundary, populates the structured
  fields, and rebuilds the rendered `description` from the structured
  parts. Downstream code reads `a.trait` / `a.phecode` / `a.p_value`
  directly — the three regex-on-prose helpers (`_gwas_base_trait` /
  `_gwas_phecode_parent` / `_gwas_p_value`) and the shared
  `_PHECODE_DELIM` delimiter are deleted. Rollup behavior is unchanged;
  MTAG detection still uses `"(MTAG)" in a.description` (promotion to a
  structured flag is outside the ADR-0035 cluster manifest). **First
  change in the cluster to deliver observable v5 JSON fields** —
  `annotation.trait` / `.p_value` / `.phecode` now appear in every JSON
  report. v4 consumers continue to ignore the unknown fields; v5
  consumers can read them as the structured contract ADR-0035 advertised
  at the schema bump.

## [2.0.2] - 2026-06-16

### Added

- **High-value SNP no-call list: 16 additions.** The list grows
  from 12 to 28 entries, adding thrombophilia (Factor V Leiden,
  Prothrombin G20210A), hereditary hemochromatosis (HFE cluster:
  C282Y + H63D), cystic fibrosis carrier (CFTR F508del), sickle cell
  (HBB), BRCA2 Ashkenazi founder, completed warfarin / clopidogrel
  PGx profiles (CYP2C9*3, CYP2C19*17, VKORC1), CYP2D6*4 (the European
  LOF allele complementing the existing CYP2D6*10), CYP3A5*3
  (tacrolimus, CPIC Level A), TPMT cluster (\*3B / \*3C for
  thiopurines), NUDT15 (thiopurine toxicity, East-Asian-enriched),
  and LRRK2 G2019S (monogenic Parkinson's). All 16 verified on at
  least one consumer-array platform per SNPedia chip listings; all
  CPIC Level A pharmacogenomic entries are paired with their
  companion SNPs where the clinical action requires haplotype
  resolution. Two new clusters surface compound-het patterns the
  way APOE / MTHFR already do. Data-only change — no code
  modification; existing `Annotator.no_call_warnings` machinery
  consumes the expanded list automatically.
- **`__version__` reads `pyproject.toml` from a bare source checkout.**
  Previously, running Allelix from a `git clone` without
  `pip install -e .` made `__version__ == "0.0.0+local"`, which then
  flowed into the outbound HTTP `User-Agent` for downloads from NCBI /
  EBI / HuggingFace and the `allelix --version` output. The sentinel
  was never visible to users on a normal install but consistently leaked
  in dev runs. Now: the importlib-metadata read remains primary; on
  `PackageNotFoundError`, fall back to reading the version from
  `pyproject.toml` next to the package. Sentinel `0.0.0+local` remains
  only when both paths fail (genuinely broken state).
- **Build detection: infer GRCh38 from `chr`-prefixed contigs.**
  When rsID-based detection doesn't converge and the header carries no
  `##contig assembly=` tag, the VCF parser now records whether the
  contigs declared `chr1` / `chrX` style names and surfaces that as a
  tertiary build-detection signal. Modern variant callers (DeepVariant,
  DRAGEN, GATK HaplotypeCaller) overwhelmingly use `chr`-prefixed
  contigs on GRCh38; GRCh37 uses bare `1` / `X` / `MT`. Falls in
  priority after override → rsID detection → `##assembly` header,
  ahead of the previous bare-GRCh37 fallback. UX win for the
  increasingly common case of rsID-less WGS VCFs. The CLI build banner
  reports `Build: GRCh38 (inferred from chr-prefixed contig names; ...)`
  and replaces the blind-default warning with a positive informational
  message that still recommends `--build` for UCSC-hg19 edge cases.
  Detection covers every standard chromosome (1–22, X, Y, M, MT) via
  the contig regex — not just `chr1` / `chrX` as in the initial cut.

### Changed

- **Pipeline: enrichment annotators are stack-managed.**
  `run_analysis` previously context-managed only the `annotators=`
  list. The optional keyword params `gnomad=`, `alphamissense=`, and
  `cadd=` were used internally but never explicitly closed; their
  SQLite connections leaked at GC time and the leak was silenced by
  `Annotator.__del__`. All three are now entered into the same
  `ExitStack` as the primary annotators, so cleanup is deterministic.
  `__del__` is retained for now as a documented safety net — removing
  it surfaces residual leaks in code paths outside `run_analysis` that
  still need auditing. Full `__del__` removal tracked for v2.1.
- **Test suite: no more silent skips on missing fixtures.** The
  v2.0.1 ship-gate ran with "1525 passed, 2 skipped" — the 2 skips were
  `TestRealDataGwasSanity` silently bypassing when
  `test_data/gwas_catalog.zip` (~65 MB, gitignored) wasn't downloaded
  locally. That's a false-clean signal on a release. The fixture now
  auto-fetches from EBI on first use and caches under `test_data/`;
  subsequent runs are offline-fast. Two defensive `pytest.skip` guards
  in `test_pipeline.py` (a committed mock fixture that should always
  exist, and an `is_ready()` check immediately after a successful load)
  are converted to `assert` — silent skips would have hidden a real
  regression there. `CONTRIBUTING.md` documents the new "no silent
  skips on ship gates" policy.

## [2.0.1] - 2026-06-15

### Known limitations

- **ClinVar per-(classification, condition) pairing (planned for
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

- **HuggingFace asset URLs moved to `allelix/` org.** The gnomAD
  and AlphaMissense prebuilt-cache download URLs were previously
  hardcoded to `huggingface.co/datasets/dial481/...`; both datasets now
  live under the `allelix/` HF org to match the GitHub migration. HF
  serves redirects from the old paths, so existing pinned client
  installs continue to work, but new installs pull from the canonical
  org URL. Commit SHAs in the resolve paths are unchanged — the
  underlying asset bytes are identical and the SHA256 integrity checks
  pass without modification.
- **Terminal report: bare-min column set.** The terminal analyze
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
- `scripts/fetch_testdata.sh` now downloads from
  `allelix/allelix/releases/download/v2.0.0/test_data.tar.gz` (the new
  ~1.27 GB asset that includes the v2.0.0 VCF/gVCF additions: real
  GIAB HG002 GRCh38 + GRCh37 benchmarks, a real GATK-HC gVCF, a 1000
  Genomes chr22 multi-sample, plus the synthetic mocks). Previous URL
  pointed at `dial481/allelix/releases/download/v1.1.1/test_data.tar.gz`,
  which lacked all the VCF samples.

### Fixed

- **`Variant`: alleles normalized to uppercase at construction.**
  Reference databases (ClinVar, gnomAD, ClinPGx, ...) all ship uppercase
  alleles, and carrier matching is raw set membership against
  `{variant.allele1, variant.allele2}` — a lowercase user allele would
  silently miss every clinical annotator and a real carrier would get
  zero annotations with no warning. Production parsers all emit
  uppercase today, but a user-supplied filter file (custom panel) or
  future format variant could leak lowercase through. Normalization now
  happens in `Variant.__post_init__` so every emission point (parsers,
  filters, synthetic variants in tests, and any future code) is covered
  at the single chokepoint. No-call marker preserved as-is; multi-base
  indels uppercased in place. Strand-aware carrier matching (using
  `resolve_strand` at carrier-match time) is deferred to v2.1+ per
  ADR-0010.
- **SNPedia parser: allowlist guard on `raw_table` interpolation.**
  The `raw_table` name from `detect_raw_table` flows into three SQL
  queries via f-string interpolation (SQLite doesn't accept parameter
  binding for identifiers). Today the function only ever returns one
  of two literals or None — safe — but the signature is `str | None`
  with no constraint, so any future drift (config-driven name,
  scraped metadata) could open an injection path. Added a
  `_VALID_RAW_TABLES = frozenset({"_raw_pages", "pages"})` guard at
  the top of `_parse_raw_pages_inner`; an unexpected name raises
  `ValueError`. No behavior change on current paths.
- **GWAS rollup: consistent PheCode delimiter.**
  `_gwas_base_trait` was splitting on `" (PheCode "` (leading space)
  while `_gwas_phecode_parent` was searching for `"(PheCode "` (no
  leading space) — the two readers reaching into the same formatted
  description string with different delimiters. Promoted the delimiter
  to a module-level `_PHECODE_DELIM` constant used by both. The deeper
  structural fix (carrying trait / p-value / PheCode as structured
  fields on `Annotation` instead of re-parsing the formatted
  description, per ADR-0024 followup) is deferred to v2.1+.
- **`split_csv_line` docstring corrected.** The previous docstring
  claimed it handles single, double, and double-double-quoted fields
  per real CSV quoting rules. The implementation is
  ``line.split(",")`` + ``strip('"')``, which mishandles any quoted
  field containing a literal comma (silently dropped by callers'
  column-count guard). Docstring now states the limitation and the
  intended consumer set (FTDNA / MyHeritage / Living DNA, none of
  which embed commas in fields). No behavior change.
- **FTDNA `can_parse` rejects MyHeritage files.** FTDNA and
  MyHeritage exports share the same `RSID,CHROMOSOME,POSITION,RESULT`
  header and a byte-identical data section. Both `can_parse`
  implementations were accepting MyHeritage files, with routing masked
  only by the registry's listing order (MyHeritage first). A reorder
  or `--format ftdna` on a MyHeritage file would silently mislabel the
  source format and display name. FTDNA now scans the first 50 lines
  for a `MyHeritage` substring and returns `False` if found; detection
  is mutually exclusive regardless of registry order.
- **GWAS `_magnitude`: `<=` at p-value boundaries.** Strict `<`
  on the magnitude thresholds meant the canonical genome-wide-
  significance value `p = 5e-8` (which appears verbatim in the GWAS
  Catalog) landed in the suggestive bucket (4.0) rather than the
  significant bucket (6.0) — one full magnitude below a barely-
  significant 4.9e-8 hit. All boundary comparisons now use `<=`
  (`5e-100`, `5e-20`, `5e-8`, `5e-6`, `5e-4`). Practical surface
  area is small given the default `--gwas-min-magnitude 9.0`, but
  the scoring was simply wrong at the most meaningful threshold.
- **`diff`: validate baseline entries, magnitude tolerance, no implicit
  zero.** Two bugs in the report-diff path. (1)
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
- **ClinVar `condition` joins all CLNDN entries.**
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
- **Enrichment: drop MAX-fallback for alt-less annotations.**
  GWAS rows set `alt=""` because the catalog is rsID-keyed. The previous
  enrichment path took a `MAX(af) / MAX(am_pathogenicity) / MAX(phred)
  GROUP BY rsid` fallback for those rows, which at a multi-allelic site
  stamped the highest-frequency / highest-pathogenicity / highest-CADD
  alt's value next to the annotation as if it described the user's
  variant. Same wrong-allele hazard as the `resolve_strand` complement
  fallback (below). Symmetric fix: skip the enrichment number rather
  than show one that may belong to a different allele. The
  position-fallback path is allele-specific and stays — it pairs the
  alt-less annotation's rsID with the user's actual carried alt. CADD
  enrichment grew the same position-fallback path (it previously only had the MAX-fallback,
  which the suppress-half removed); for symmetry, alt-less rows whose
  rsID resolved on the fly now get a correct allele-specific CADD score
  rather than nothing. The full fix — carrying the user's alt onto
  every Annotation so GWAS rows can take the exact-alt path even when
  the input VCF already had rsIDs — is architectural and tracked for
  v2.1.
- **`resolve_strand`: drop complement-fallback.** At multi-allelic
  sites, the complement of the user's true forward allele can
  coincidentally equal a different alt at the same position. The
  previous fallback (when the user allele didn't match `{ref, alt}`
  directly, try `complement(user_allele) ∈ {ref, alt}`) would return
  that coincidental complement, and the CADD enrichment path then
  stamped a wrong-allele score onto an annotation that read as
  describing the user's variant. Example:
  `resolve_strand('A', 'C', 'T') -> 'T'` (the user's `A` is neither
  REF nor ALT but `complement('A') = 'T' = ALT`). The function now
  returns `None` whenever the user allele isn't directly in
  `{ref, alt}`. ADR-0010 already deferred proper strand handling; this
  makes the deferral explicit at the resolution layer. The CADD
  enrichment helper `_lookup_user_allele` is simplified — the
  complement loop is removed. Bounded blast radius (enrichment-only,
  not the carrier decision), but the previous behavior could read as
  describing the user's variant when it didn't.
- **`normalize_build_label`: anchored tokens, ambiguity rejection.**
  The previous implementation was first-substring-wins over a bare
  `"3X" in s` chain. Wrong results it gave: `2038-01-01` registered
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
- **LivingDNA parser: only inspect build-marker comment lines.**
  `get_metadata` previously fed every comment line into
  `normalize_build_label` and let the *last* match win, so a download-
  date comment could retag the file. It now skips comment lines that
  don't mention ``build`` / ``reference`` / ``genome``, and breaks on
  the first valid build-marker line — matching the format spec which
  puts ``# Human Genome Reference Build 37 (GRCh37.p13).`` near the
  top of the header.
- **gnomAD `is_ready()`: tagless caches fail.** The previous check
  was `tag == f"sv:{GNOMAD_SCHEMA_VERSION}" or not tag`, so any cache
  written before the schema-version mechanism existed was accepted as
  ready. That defeats invalidation: if `GNOMAD_SCHEMA_VERSION` is ever
  bumped, every legacy cache silently passes as the new version. The
  escape is dropped — `is_ready()` now requires the exact current tag.
  Production caches built via `install_prebuilt_cache` already stamp
  the current tag, so the only caches affected are pre-tag legacy ones,
  which now correctly require `db update`.
- **ClinVar `.md5` body validated as 32 hex chars.** When the
  ClinVar CDN returns a transient HTML error page, the previous parser
  treated the first whitespace-separated token (`<!DOCTYPE`) as the
  expected hash. That token then flowed into `verify_file_hash`, which
  on mismatch deletes the just-downloaded ~400 MB VCF and raises. So a
  five-second 503 upstream wiped the freshly fetched cache. The
  signal parser now requires the token to match `^[0-9a-fA-F]{32}$`
  before treating it as an md5 — non-hex, short hex, or HTML is logged
  and returns `None`, which `db update` handles as "freshness unknown,
  skip" rather than nuking the VCF.
- **`db update`: tagless caches re-download instead of being stamped.**
  A cache with no stored freshness signal almost always predates the signal
  mechanism, i.e. is old. The previous behavior was to stamp the live
  remote signal onto it on the next `db update` and call it current, which
  permanently marked stale data as fresh — every subsequent run saw
  `cached == remote` and skipped, so only `--force` could escape.
  `cached is None` now triggers `setup()` like a remote-signal mismatch
  does, with a distinct message (`cache predates the freshness signal;
  re-downloading…`) so the user knows why. The dead helper
  `_stamp_remote_signal` is removed.
- **GWAS / ClinPGx loaders: post-load row-count floor.** Closes the
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
- **Build detection: warn when inconclusive despite rsIDs.** Closes
  the silent-coords trap. When position-detection inspected known-rsID
  positions but couldn't pick a build (votes tied across builds or no row
  matched any), the pipeline silently routed through `header_build` — so a
  GRCh36 file with a GRCh37-mislabeled header used to silently get the
  GRCh37 ClinVar cache. Now fires a yellow warning naming the
  rsID-checks-attempted count, what the header claims, and how to force
  `--build`. The CLI-layer warning leaves `BuildDetectionResult`
  unchanged.

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
  variants) drove annotation count from 0 → 4,867.
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
  the pre-resolution empty rsid and missed).
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
- **ClinVar CLNSIG vocabulary drift CI test**
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
- **PharmGKB rebranded to ClinPGx**. Upstream
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

- **rsID resolver: indel-anchor disambiguation.** The position-keyed
  resolver matched ClinVar rows by ALT only. At positions where two co-located ClinVar entries share
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
- **VCFs from variant callers now annotated**. Pre-fix, real
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
  neon blue, unreadable against the dark background. Dark mode
  now uses `#93c5fd`; light mode uses `#1976d2` (matching the existing accent).
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
- **PLINK export.** `allelix export plink` converts any supported
  genotype format to PLINK1 binary (.bed/.bim/.fam) for downstream
  tools (plink2 PCA, ADMIXTURE, PRSice). Single-sample, SNP-major
  encoding. Uses gnomAD ref/alt for allele coding when available.
  No-call variants skipped.

### Fixed
- **Multi-allelic strand collision in PLINK coord selection.** At sites
  where one alt is the complement of another (e.g. ref=G, alts=A,T),
  single-pass coord selection picked the complement match over the
  forward match — same bug class as CADD. Fixed with a two-pass
  loop that prefers forward allele matches.

### Documentation
- **ADR-0034: Magnitude scoring scale and ceiling.** Formalizes the
  0-10 scale with practical ceiling of 9. Documents max-across-sources
  composite rule and reserves per-source scoring for v2.0.0.

## [1.6.1]

### Added
- **Zygosity column in all report outputs.** Every annotation row
  now shows `Heterozygous`, `Homozygous`, or `No Call` — derived from
  the genotype call. Appears in HTML, terminal, and JSON reports.
  Functional-medicine `+/−` notation deferred until the risk allele
  field lands (v2.1+).
- **CADD PHRED styling in HTML reports.** CADD scores are now
  color-coded by deleteriousness tier: red (≥30, top 0.1%), orange
  (≥20, top 1%), gray (10–20, top 10%), gray (<10, no tooltip).
  Tooltip shows the percentile tier on hover for scores ≥10.
- **Config file path in `config get` / `config set` output.**
  Both commands now print the resolved config file path, so users know
  which file is being read or written.

### Changed
- **JSON schema version bumped to 4.** New `zygosity` field on every
  annotation. Diff between v1–v4 reports still works.
- **Methylation panel expanded to 34 genes.** Added ACAT1, DHFR,
  GNMT, MAOA, NOS3, SUOX, VDR. Sorted alphabetically.

### Fixed
- **CADD multi-allelic scoring bug.** At multi-allelic positions,
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
- **Circular import.** `_versions.py` moved from
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
- **Multi-allelic enrichment accuracy.** gnomAD and AlphaMissense
  enrichment now uses exact alt-allele matching instead of `MAX()`
  aggregation at multi-allelic sites. Added `alt` field to the
  Annotation model; `bulk_lookup_by_alt()` on both enrichment
  annotators; pipeline splits exact-match and MAX-fallback paths.
- **Disk preflight multiplier.** Bumped from 5x to 6x for both
  gnomAD and AlphaMissense loaders. AlphaMissense compresses at 4.4x,
  so peak disk (gz + decompressed) is 5.4x gz — the old 5x check
  would greenlight a disk that ENOSPC'd at ~90% decompression.
- **Test suite disk usage.** db update tests were downloading real
  databases (7.8 GB AlphaMissense, 678 MB GWAS TSV) into pytest
  tmp_path because they lacked monkeypatches for all annotators. All
  db update tests now stub every annotator. Real-data GWAS tests
  delete the extracted TSV after SQLite load. Total pytest tmp_path
  reduced from ~3.4 GB to ~376 MB.
- **SQLite variable limit portability.** `bulk_lookup_by_alt()`
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
- **AlphaMissense gnomAD version stamping.** Build script stamps
  which gnomAD version provided the rsID mapping. Runtime warning on
  version mismatch: "AlphaMissense cache was built against gnomAD X
  but installed gnomAD is Y."
- **SNPedia HuggingFace download.** `db update` now downloads
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
- HTML report table overflows viewport, columns clipped on left.
  Added `overflow-x: auto` container, sticky rsID column,
  `max-width` on description cells, refs collapsed into `<details>`
  toggle, conditional Review Status column (hidden when all empty),
  stat card `flex-wrap`.
- AlphaMissense build script has zero unit-test coverage. Added
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
- **gnomAD population allele frequencies.** New `GnomadAnnotator`
  enriches report annotations with population frequency context from
  gnomAD v4.1 exomes (~16M rsIDs). Pre-built cache downloaded from
  HuggingFace via `db update`. Frequency column appears in terminal,
  HTML, and JSON reports when gnomAD data is available. `--no-gnomad`
  flag on `analyze`, `methylation`, `pharmacogenomics`, and `db update`
  to skip.
- **CPIC fallback for PharmGKB.** `db update` succeeds when the
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
  opt-in network access.
- `.gitignore` updated for GWAS Catalog test data.
- `scripts/fetch_testdata.sh` downloads GWAS Catalog associations
  from EBI FTP.

### Changed
- Pre-push hook reduced to version-tag check only; pytest removed
  (CI runs the full suite on every PR, pre-push pytest caused SSH timeouts).
- CI version-tag guard job added to `.github/workflows/ci.yml`.

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
- **Build mismatch warning in HTML report**. When the file header claims
  a different build than position data indicates, the HTML report renders a
  visible warning banner matching the CLI warning.
- **"Reading This Report" education section**. Static HTML block after
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
  builds.
- **No ClinVar GRCh36 cache.** `CLINVAR_SUPPORTED_BUILDS` is
  GRCh37/GRCh38 only. When a file is correctly detected as GRCh36 (e.g.
  23andMe with "build 36" in the header), `analyze` produces zero
  ClinVar annotations — PharmGKB, GWAS, and SNPedia still fire normally.
  The CLI prints a warning. Full GRCh36 annotation support requires
  either an NCBI GRCh36 VCF or liftover.

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
- Added roadmap feature: supplemental genotype file merging
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
  forward-strand assumption (deferred until strand-aware matching).
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

> **Genome build auto-detection (ADR-0021 + ADR-0022).** A real-world
> MyHappyGenes/Tempus export was confirmed to ship GRCh38 positions while its header claims "build
> 37.1." Cross-build REF/ALT comparison against ClinVar's GRCh37 VCF
> produced a false-positive pathogenic call on NIPA1 `rs104894490`.
> The carrier check (ADR-0007) was correct — it was matching against
> the wrong build's REF/ALT. The fix is structural: detect build from
> position data, hold per-build ClinVar caches, dispatch per variant.
> ADR-0022 documents the deliberate decision NOT to filter PharmGKB
> reference-genotype rows on non-CPIC genes.
>
> This release also fixes three fixture-layer defects and one
> `.gitignore` omission in the build-detection work.
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
- **Per-build ClinVar fixtures.** The earlier mock_clinvar.vcf
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
- **Test suite no longer hits the real CPIC API.**
  `test_db_update_with_file_url` monkey-patches
  `fetch_cpic_allele_functions` on the annotator module to return the
  same `MOCK_CPIC_LOOKUP` the rest of the suite uses. Run time on that
  single test dropped from ~4 s (network-bound) to 0.6 s (offline).
- **README refreshed to v0.7.1.** Status banner, supported-
  databases table, and architecture-decision summary reflect ADR-0020
  (CPIC API as the structured per-allele function source) and drop
  obsolete references to ADR-0013 / ADR-0014.
- **CPIC conflict-resolution policy flipped.** When the same
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
- **Retry-with-backoff on CPIC fetches.** `_http_get_json`
  retries up to 3 times with `(1, 2, 4) s` backoff on `URLError`,
  `TimeoutError`, and `JSONDecodeError`. A single transient TCP RST
  during `db update` no longer aborts the refresh.
- **Composite PharmGKB + CPIC freshness signal.**
  `PharmGKBAnnotator.fetch_remote_signal()` now returns
  `pgkb:<pgkb-signal>|cpic:<cpic-signal>` where the CPIC portion is
  the latest date from CPIC's `change_log` table. If CPIC publishes
  new allele functions while PharmGKB's zip is unchanged, the next
  `db update` detects it and refreshes. CPIC probe failure returns
  None (existing "can't verify, pass --force" UX).
- **Mutation-gap tests pinning critical policies.**
  `_classify_cpic_status` returning None for unknown statuses, and
  `fetch_cpic_allele_functions` skipping multi-base / non-ACGT values,
  both have direct pins. A future refactor can't silently break
  either without a test failure.
- **Zip cleanup wrapped in try/finally** in `PharmGKBAnnotator.setup()`.
  A failure between download and ingest no longer leaves the staged
  `clinicalAnnotations.zip` on disk.
- **`MOCK_CPIC_LOOKUP` shape invariants.** New
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
- **PharmGKB SNV per-allele function now structurally extracted
  (ADR-0018).** The DPYD cluster that v0.6.1 wrongly reported as
  findings (rs115232898 TT, rs1801266 GG, rs3918290 CC and ~17 others of
  the same shape) now classifies correctly. Reference-allele homozygotes
  do not emit.

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
  the reference-homozygote shapes v0.6.1 wrongly reported).

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
- Two prose-classifying regexes in `pharmgkb_loader.py` removed. Remaining
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
> fix; v0.5.0 adds the two PharmGKB fixes and is
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
- Minor fixes:
  - Removed dead `render_annotations` from `reports/terminal.py`;
    `render_terminal` is the single entry point. Test file rewritten to
    assert against `render_terminal` + `AnalysisResult` fixtures.
  - Dropped the unreachable `else` branch in
    `_run_analysis_command`; renderer dispatch is now a ternary.
  - `analyze --category` help text no longer lists "methylation"
    (no annotator emits that category — use `allelix methylation` instead).
  - HTML and JSON outputs are now atomic via a shared
    `atomic_write_text` helper in `allelix.reports` (`.tmp` then
    `os.replace`). A killed process mid-write leaves either the previous
    file or no file, never a half-written one.
  - `REGULATORY_NOTICE` moved from `json_report.py` to
    `allelix.reports`. Both `json_report` and `html` import from there.

### Added
- **Coverage**: PA-007 fixture row in the synthetic PharmGKB dump
  exercises the inner `_normalize_genotype` skip path
  (single-rsid annotation with only non-SNV per-allele rows).
- `tests/reports/test_init.py` covers `REGULATORY_NOTICE` presence and
  `atomic_write_text` (success path + simulated-rename-failure cleanup).

## [0.4.0] — 2026-05-11

### Added
- **Reports.** New `allelix/reports/_pipeline.py` builds an
  `AnalysisResult` once and hands it to format-specific renderers, so
  the CLI streams the file and queries each annotator exactly once per
  `analyze` invocation regardless of output format.
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
- **Strand-flip helpers** (partial): `allelix.utils.allele`
  ships `complement`, `flip_genotype`, `is_strand_ambiguous` for parser /
  annotator authors. Not wired into matching logic yet — opt-in.
- **ADR-0010**: documents the strand-flip decision and explicitly defers
  liftover to v1.0.0 with the rationale (chain-file weight, no real GRCh38
  consumer yet).

### Changed
- Annotators are entered into a `contextlib.ExitStack` inside the pipeline
  module, replacing the inline ExitStack in the CLI's `analyze` command.
  The same deterministic-cleanup contract applies across all report
  commands now.
- `allelix/cli.py` factored: `_resolve_parser`, `_ready_annotators`,
  `_run_analysis_command` are shared by `analyze`, `methylation`, and
  `pharmacogenomics`.

### Fixed
- `test_version_flag` no longer hardcodes a version string; it
  asserts the rendered output contains `__version__` and a new
  `test_pyproject_version_matches_metadata` test in `tests/test_version.py`
  pins that `pyproject.toml`'s version equals the installed package
  metadata. Catches the regression class where someone bumps
  pyproject.toml without reinstalling the editable package.

### Added
- PharmGKB loader gets a batch-flush spy test
  (`test_batched_insert_flushes`) — `INSERT_BATCH_SIZE = 3` against the
  9-record fixture asserts `[3, 3, 3]`. A mutation flipping the flush
  condition to `if False:` now fails the suite.
- Asymmetric no-call test (`Variant("A", "-")` and `("-", "A")`) for
  both ClinVar and PharmGKB annotators. Verified that the ClinVar
  test catches `if variant.is_no_call:` → `if variant.allele1 == "-":`.
- `_safe_float` direct unit tests cover the empty-string and ValueError
  branches in `pharmgkb_loader`.
- `test_leftover_tmp_file_is_cleared` covers the stale-`.tmp` cleanup
  branch in `load_pharmgkb_tsv`.

### Changed
- `pharmgkb_loader._open_directory` now imports `Path` at module top
  instead of inline-aliasing inside the function. `extractall`
  now has an explicit code comment about Python 3.11+'s built-in
  `..`/absolute-path sanitization.

## [0.3.0] — 2026-05-11

### Added
- **PharmGKB annotator** (first cut). `allelix db update` now also
  downloads PharmGKB's `clinicalAnnotations.zip` into `pharmgkb.sqlite`;
  `allelix analyze` layers in pharmacogenomic annotations alongside
  ClinVar's clinical-significance calls.
- `allelix.databases.pharmgkb_loader` parses the joined
  `clinical_ann.tsv` + `clinical_ann_alleles.tsv` view, normalizes 2-letter
  SNV genotypes (sorted, uppercased), and skips star alleles, multi-rsid
  composites, and indel genotypes (not supported in this release).
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
- Annotator test fixture now closes connections on teardown.
- `pytest.PytestUnraisableExceptionWarning` is escalated to error in
  `filterwarnings`, so future fixtures without teardown fail CI loudly.
- Tests covering the INSERT batch flush path, the download timeout path,
  the `parse_clinvar_version` end-of-file fall-through, and friendly
  `db update` failure output.
- The batch-flush test now spies on `executemany`
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


[2.5.0]: https://github.com/allelix/allelix/compare/v2.4.2...v2.5.0
[2.4.2]: https://github.com/allelix/allelix/compare/v2.4.1...v2.4.2
[2.4.1]: https://github.com/allelix/allelix/compare/v2.4.0...v2.4.1
[2.4.0]: https://github.com/allelix/allelix/compare/v2.3.1...v2.4.0
[2.3.1]: https://github.com/allelix/allelix/compare/v2.3.0...v2.3.1
[2.3.0]: https://github.com/allelix/allelix/compare/v2.2.3...v2.3.0
[2.2.3]: https://github.com/allelix/allelix/compare/v2.2.2...v2.2.3
[2.2.2]: https://github.com/allelix/allelix/compare/v2.2.1...v2.2.2
[2.2.1]: https://github.com/allelix/allelix/compare/v2.2.0...v2.2.1
[2.2.0]: https://github.com/allelix/allelix/compare/v2.1.2...v2.2.0
[2.1.2]: https://github.com/allelix/allelix/compare/v2.1.1...v2.1.2
[2.1.1]: https://github.com/allelix/allelix/compare/v2.1.0...v2.1.1
[2.1.0]: https://github.com/allelix/allelix/compare/v2.0.2...v2.1.0
[2.0.2]: https://github.com/allelix/allelix/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/allelix/allelix/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/allelix/allelix/releases/tag/v2.0.0
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

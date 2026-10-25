# Provider file formats

What each provider delivers, which file Allelix needs, and what Allelix
does with the provider-specific parts of the file. The general format
table is in the [README](../README.md#supported-formats).

## Compressed files and other wrapping

Every format below can be passed gzipped (`.gz`, or a gzipped file under
any name; Allelix checks the content, not the extension). A UTF-8 byte
order mark at the start of the file and blank lines before the first line
are ignored. Family Tree DNA, for example, delivers its raw data as a
`.csv.gz`, which Allelix reads without unpacking it first.

A zip download (23andMe and AncestryDNA deliver one) can be passed as it
is, when it holds one genotype file. The file is read from inside the zip,
never unpacked to disk, and may itself be gzipped. Folders and the hidden
entries macOS adds when it makes a zip (`__MACOSX/…`, `.DS_Store`) are not
counted as files. Allelix stops with a message saying what to do when the
zip:

- holds several files (it names them): unzip it and pass the genotype file;
- is password-protected, or uses a compression method Python cannot read
  (such as Deflate64, which Windows uses for very large zips): unzip it
  first;
- holds a file over 4 GiB once unzipped (for a gzipped file in the zip,
  that is the size of the `.gz` file itself): unzip it, or gzip it instead
  (a gzipped file has no size limit);
- holds another zip file: unzip it until you reach the genotype file;
- is damaged (its directory is missing, a file's checksum fails, or a
  file's data is longer or shorter than its zip entry declares): the same
  "looks truncated or damaged" error as a damaged gzip file.

Windows (CRLF), Linux (LF) and classic Mac OS (CR only) line endings are
all read, in every format, plain, gzipped or zipped.

In every array format a no-call may be `--`, `00` or `NN` (or `-`, `0`,
`N` in a one-letter or two-column file); all are read as no-calls.

A compressed file that was cut short (an interrupted download or copy) or
damaged is not read as if it were complete:

- If the file ends in the middle of its compressed data, or the data is
  corrupt, Allelix stops with an error naming the file ("looks truncated or
  damaged"). Download or copy the file again. This includes a gzip file cut
  within its first bytes, before any line can be read: it gets the same
  error, not "doesn't recognise the format".
- A bgzipped file (the usual `.vcf.gz`) that was cut exactly between two
  of its blocks still reads without error, so Allelix checks for the
  end-of-file block that every complete bgzipped file ends with. If it is
  missing, the records before the cut are used and a warning says the file
  may be truncated. The same check runs on a bgzipped file inside a zip.
- Bytes after the end of the compressed data that are not compressed data
  themselves (for example a web page saved onto the end of the file) are
  ignored with a warning; every complete record is read.

Both warnings are printed once on the console, and an `analyze` report
lists them at the top with its other warnings (JSON `report_warnings`,
codes `input_possibly_truncated` and `input_trailing_bytes`).

## Forcing a format

Allelix reads a file only when it recognises its source: each reader checks
the file's header and layout first. A file no reader recognises is
rejected: "Allelix doesn't recognise the format of …". Use the original
file exactly as your testing company provided it, not edited or re-saved
(opened and saved again in a spreadsheet, for example). Zipped or gzipped
is fine.

`--format <name>` reads a file with the reader you name even when that
reader does not recognise it (`analyze`, `stats`, `extract`, `export
plink`, `methylation`, `pharmacogenomics`; `--format1` and `--format2` for
`compare`, each file checked on its own). The reader names are
`myhappygenes`, `genesforgood`, `23andme`, `ancestrydna`, `livingdna`,
`decodeme`, `myheritage`, `ftdna`, `ftdna_famfinder`, `ftdna_illumina` and
`vcf`. This is an expert override, at your own risk. Allelix cannot check what it would normally check:

- **layout:** the reader takes the columns where its format puts them, so
  a file laid out differently is misread;
- **strand:** the calls are taken to be on the strand that format uses; a
  file on the other strand gives wrong ClinVar results, at A/T and C/G
  sites and elsewhere;
- **build label:** the build line the reader looks for may be missing or
  mean something else, so the format's default build can apply. Known-SNP
  positions still check the build, as for any file.

A forced file that its reader does not recognise gets a warning before
anything else in the terminal, and on every report the run writes: the
terminal summary, the top of the HTML report, and JSON `report_warnings`
(code `forced_format`). When the named reader does recognise the file,
`--format` only chooses between readers, and there is no warning.

Zip checks still apply under `--format`: a zip with several files, a
password-protected zip and a damaged zip stop as they do without it.

## Reference calls in VCFs

A homozygous-reference call (`0/0`, or a haploid `0`) is assessed where
the reference genome itself carries a ClinVar variant allele, such as
Factor V Leiden on GRCh37. Its input REF must agree with the known base
at that position; a disagreement is counted in
`input.vcf_quality.clinvar_ref_mismatches` and described on the console
and in the HTML input-quality notes. Reference sources are checked first:
single-base A/C/G/T REF evidence from the same identity-table build,
chromosome and position is used regardless of its rsID. If those rows
conflict with one another or with ClinVar's reference base, the call is
not reconstructed and `input.vcf_quality.reference_source_conflicts`
counts it. The note describes source ambiguity rather than blaming the
input REF. Duplicate alias rows do not multiply that count. Missing,
noncanonical, multi-base, or differently placed evidence does not create
a conflict. These rejection checks do not assign rsIDs or infer an array
probe's strand. Haploid calls keep one allele and
are shown as Hemizygous. Reconstructed calls follow the same mitochondrial
coordinate exclusions as other VCF calls.

Matching reference observations must also agree about ploidy. If usable
observations disagree between haploid and diploid, the site is left
unassessed in either row order. JSON schema 10 adds the count of
rejected matching observations in `input.vcf_quality.reference_ploidy_conflicts`,
with a note on the console and in HTML. REF-mismatching observations are
counted separately and do not establish ploidy. Reference-source conflicts
take precedence; alias rows never multiply observation counts. Repeated
and overlapping records retain only bounded base/ploidy counts per watched
position, including a separate state for unusable explicit REF bases.
With `--vcf-out`, the explicit reference records at watched positions are
also remembered, one entry each, so the annotated VCF can stamp them.

A reference block has ALT `<NON_REF>`, `<*>` or `.`, a homozygous-reference
genotype, and an `END` position. Its REF describes the block's start;
that base is checked only at positions explicitly written in REF, not
against a different base inside the block. Structural-variant alleles
such as `<DEL>`, `<DUP>` and breakends never provide reference coverage
over their `END` span. Their unassessed-record counts are still reported.

A VCF containing only reference calls is analysed and says
"This file has no non-reference calls." A VCF with no genotype records
still stops. `allelix stats` displays genotype-record, reference-call and
non-reference-call counts, before filtering. JSON schema 10 adds the same
counts as `genotype_records`, `hom_ref_calls` and `non_reference_calls`
under `input.vcf_quality`. Reference blocks count as records, not as one
record for every base they span. Non-reference counts include symbolic
ALT genotypes; they do not imply those structural calls were assessed.
Schema 10 also adds a nullable `reason` to panel coverage entries.
Reports in schemas 1–10 can be used as `--diff` baselines; versions 11
and later are rejected. `--diff` matches findings by source, variant, condition,
allele and genotype, whatever their order; a finding left on both sides with
a different genotype is listed as changed. An entry from an older report
without an allele is matched without allele and genotype. Each changed entry
in the JSON `diff` block carries `previous_genotype_match` (null when the
baseline entry had none), and the terminal shows `AG -> AA` when it differs.

## Annotated VCF: reference findings and re-annotation

`analyze --vcf-out` puts a ClinVar finding for an allele the reference
itself carries (see above) on the explicit `0/0` or haploid `0` record it
came from, in `ALLELIX_CLINVAR`. The field is `Number=A`, so the same
finding appears in each ALT slot, with `.` for symbolic ALTs; the field's
description says that on a reference record the finding is for the
reference allele. Only the record the analysis read the call from gets it,
so the selected sample (`--sample`), the quality rules (`--min-gq`, FILTER)
and the mitochondrial exclusions apply as in the report. A record whose REF
has more than one base gets it when that REF starts at or covers the site.
gnomAD, AlphaMissense and CADD values for the reference allele are not
written, since they would read as values for the listed ALT.

A gVCF reference block is never stamped, because one line covers many
positions, and a record with ALT `.` has no ALT slot. A finding from
either is in the JSON and HTML reports only.

When the input is itself an Allelix-annotated VCF, the old `ALLELIX_*` INFO
values, their `##INFO` declarations and the `##ALLELIX_` header lines are
removed first, on every row (filtered rows and reference blocks included),
and only the current result is written. Other INFO values keep their exact
text and order. Lookups match REF/ALT case-insensitively; the file keeps
its own case.

## Imputed VCFs

A VCF written by an imputation tool holds mostly genotypes filled in from
a reference panel, not measured. Allelix labels such a file imputed in
every report (JSON `input.method` is `vcf-imputed`), with the caution that
rare variants are the least reliable imputed calls. It decides from the
header alone, never from the values. A file is imputed when its header
declares one of these fields:

- `##INFO` `IMPUTED`, `TYPED` or `R2` (Minimac);
- `##INFO` `DR2` or `IMP` (Beagle);
- `##FORMAT` `DS`, the dosage (Minimac, Beagle, GLIMPSE), or `HS`
  (GLIMPSE).

A `##FORMAT` `GP` alone does not count, because a variant caller can add
it. Neither does GATK's `##INFO` `DS`, which means "downsampled". Genes for
Good's imputed VCF declares none of these fields; it is recognised as
described under [Genes for Good](#genes-for-good).

## Sequencing.com (30x whole genome)

A 30x whole-genome order delivers three VCF files. Use the first one:

| File | Contents | Use with Allelix |
|---|---|---|
| `…snp-indel.genome.vcf.gz` | The genotypes: SNVs and small indels, as a gVCF with reference blocks. Typically 150–300 MB compressed. | **Yes.** |
| `…sv.vcf.gz` | Structural variants, written as symbolic alleles (`<DEL>`, `<DUP>`, `<INV>`, `<INS>`) and breakends. | No |
| `…cnv.vcf.gz` | Copy-number calls, written as symbolic alleles. | No |

The genotype file is DRAGEN germline output with a header that
Sequencing.com rewrites. The start of the header looks like this:

```
##fileformat=VCFv4.2
##fileDate=January-01-2000 at 00:00 UTC
##source=Sequencing.com (30x WGS)
##dataSourceType=WGS
##dataAnalysisProvider=Sequencing.com
```

What Allelix does:

- **The `snp-indel.genome.vcf.gz` file** is read as a gVCF. The variant
  caller's FILTER and genotype quality are honoured. Reference blocks with
  low quality or no reads count as no-calls at high-value positions. See
  the VCF row of the README format table.
- **The SV and CNV files** contain only records Allelix does not annotate.
  `allelix analyze` stops with an error that names the
  `snp-indel.genome.vcf.gz` file, instead of producing an empty report.
- **File size.** The genotype file is too large for most upload forms. Run
  Allelix locally (`pip install allelix`, then `allelix analyze
  file.snp-indel.genome.vcf.gz`); it handles the full file in minutes on a
  laptop.

## MyHeritage

MyHeritage offers two kinds of raw data export. Both are the same CSV
with a header line (`RSID,CHROMOSOME,POSITION,RESULT`) and one genotype
per row. They differ in how the genotypes were produced and in their
preamble.

### Array export (kits processed before 2026)

Genotypes from a genotyping array. The preamble is a few comment lines
that name MyHeritage, for example:

```
# MyHeritage DNA raw data.
# For each SNP, we provide the identifier, chromosome number, base pair position and genotype.
# The genotype is reported on the forward (+) strand with respect to the human reference build 37.
RSID,CHROMOSOME,POSITION,RESULT
```

Allelix reads the "human reference build 37" sentence as the file's build
label (also when a spreadsheet stripped its `#`), and checks it against the
positions of well-known SNPs across the whole file (ADR-0044).

If the file was opened in a spreadsheet and saved again, it still works.
Allelix recognises all of these re-save changes:

- every row padded with trailing commas (`"rs4477212","1","82154","AA",,,`);
- an emptied first row (`,,,`) pushing the MyHeritage line down;
- comment lines that lost their `#`;
- a tab-delimited "save as text" copy.

### Low-pass whole genome sequencing (kits processed from January 2026)

MyHeritage now sequences kits at about 2x coverage and fills in the rest
by imputation. **The genotypes in this export are imputed**, not read
directly at each position. The file starts with a `##` preamble:

```
##fileformat=MyHeritage
##format=MHv1.0
##method=Low-pass Whole Genome Sequencing
##timestamp=2000-01-01 00:00:00 UTC
##reference=build37
RSID,CHROMOSOME,POSITION,RESULT
```

What Allelix does with it:

- **`##method`**: reports the input as "low-pass whole genome sequencing
  (about 2x coverage); genotypes are imputed" in the HTML report, the
  terminal output and the JSON `input.method` field (`low-pass-wgs-imputed`).
  The HTML report's array caveat is replaced by an imputation caveat.
  Imputation is accurate for common variants and least reliable for rare
  ones, so confirm any finding by direct genotyping before acting on it.
- **`##reference=build37`**: an explicit header build (GRCh37). It is the
  file's label in the build decision; confident position evidence can still
  override it, and contradicting tentative evidence stops the run (ADR-0044).
- **Not in the file:** Y-chromosome and mitochondrial genotypes. The X
  chromosome is included.

Only the CSV can be downloaded at present; a sequencing data file (CRAM)
has been announced but is not available.

### Keep the preamble

Keep the comment lines at the top of a MyHeritage file. Without them,
nothing in the file says MyHeritage: the header line and rows are laid out
exactly as in a Family Tree DNA file, so Allelix reads the file as Family
Tree DNA (the CSV layout, or the tab-delimited "Illumina raw" layout for a
tab-delimited copy). The genotypes are read the same way, but a low-pass export
loses its `##method` line, so the report no longer says the genotypes are
imputed, and the build the preamble states is lost with it.

## 23andMe

A tab-delimited text file: a comment block, then one genotype per row
(`rsid`, `chromosome`, `position`, `genotype`). Allelix recognises the file
by its banner line, `# This data file generated by 23andMe` or, in v5
exports from June 2026, `# This data file is generated by 23andMe.`

Current v5 exports start with three identification lines before the banner
(values replaced with placeholders here):

```
# file_id: 00000000-0000-0000-0000-000000000000
# signature: <64 hex characters>
# timestamp: 2026-01-01 00:00:00
#
# This data file is generated by 23andMe.
```

The banner does not have to be the first line; Allelix looks for it in the
first 50 comment lines. A file with the banner removed is not recognised
as 23andMe.

A row with no tab in it is split on spaces instead, so a file that a text
editor or spreadsheet re-saved with spaces between the columns still reads.

## Genes for Good

Genes for Good, a University of Michigan study now closed, gave each
participant a zip of their genotypes. It holds a measured and an imputed
set, each as a VCF and as a "23andMe-style" text file.

The text file has 23andMe's four tab-separated columns (`rsid`,
`chromosome`, `position`, `genotype`) under a comment block written by
PLINK. Allelix recognises it by its first line, which starts
`# Genes for Good v`. Insertions and deletions are written `D` and `I`, and
haploid calls (X and Y in males, and MT) as one letter. Uncalled sites are left
out of the file rather than written as `--`. The header names the
reference build, which the known-SNP positions check like any other label.

The imputed file says "imputed genotype calls" on its first line and has
about 7.7 million rows, against about 550,000 in the measured one. Its
genotypes are mostly not measured, so every report says so and the HTML
replaces the array caveat with a note on imputation.

The download's VCFs are read as VCF, and the imputed one is marked
imputed: its header has a `##Genes_for_Good_version` line and its calls
are phased. Allelix annotates a VCF's non-reference calls; a `0/0` record
counts only where the reference itself carries a ClinVar variant allele.
So most of the many `0/0` records of a measured VCF give no answer. The
VCFs' header states no build that Allelix reads, so they may need
`--build grch37`. The text file carries every call and states its build,
so it is the one to pass.

## deCODEme

deCODEme, deCODE genetics' consumer service (closed in 2013), delivered
`deCODEme_scan.csv`: no comment block, a column line
`Name,Variation,Chromosome,Position,Strand,YourCode`, then one SNP per
row. Allelix recognises the file by that column line.

`YourCode` is the genotype (`--` for a no-call) and `Strand` says which
strand it is written on. Every other format and annotation source uses the
plus strand, so a `-` row's call is complemented (`TG` becomes `AC`). This
was checked on real deCODEme files against real 23andMe exports: minus-strand
calls match only after complementing, plus-strand calls only as given.
`Variation` is written on the same strand as the call, so it is not used.
Chromosome `M` is read as MT.

The file states no build. Its coordinates are Build 36 (GRCh36), which the
known-SNP positions detect, so the run takes the GRCh36 path. ClinPGx
matches by rsID. GWAS Catalog risk alleles and single-base SNPedia
genotypes are not assessed: their alleles cannot be checked against the
strand without GRCh37 or GRCh38 positions, and the report says how many. ClinVar, gnomAD, AlphaMissense and
CADD give nothing (see [grch36-liftover.md](grch36-liftover.md) to lift
the file over first).

## Comparing and extracting VCF calls

`allelix compare first.vcf second.vcf` matches nonempty identifiers as
before. Unnamed calls use normalized chromosome, position, input REF and
the complete uppercase input ALT field when both effective file builds
are the same known build. The ALT list is literal: reordered lists or
different uncalled alternatives are distinct records. Quality-masked and
explicit no-calls retain their original allele identity. No liftover,
allele decomposition or genotype-based ALT inference is performed.

The coverage table separates parsed named/unnamed observations from
eligible unique named/unnamed keys. It labels shared unnamed variants
compared by position and counts omitted unnamed observations per file.
Unknown/different builds or missing input allele identity omit unnamed
calls. Nonempty vendor identifiers remain on the named path; named and
unnamed records do not match one another. Repeated identical keys retain
the existing last-observation convention, so conflicting duplicate
genotypes can depend on their order. Concordance percentages use shared
eligible keys; no shared keys is displayed without a percentage.

Indexed `allelix extract` preserves its normal build, sample, quality and
resolution rules. Any requested identifiers not actually found in the
index share one sequential pass, even if a coordinate mapping existed.
Already-found calls, including no-calls, are retained. Unavailable or
ineligible indexes use the normal sequential path. Actual index contigs
are matched through chromosome normalization, including M/MT/chrM/chrMT
and bare/prefixed autosomes; aliases are traversed in index order.

## MyHappyGenes (Tempus)

Tab-delimited, five columns. The header states the build as "37.1", but
real exports have shipped GRCh38 coordinates under that label. Allelix
therefore checks the label against the positions of well-known SNPs across
the whole file: confident positions win, with a warning when they disagree
with the header; weaker evidence against the header stops the run and asks
for `--build` (ADR-0021, ADR-0044).

## PLINK export

`allelix export plink` writes `<prefix>.bed`, `<prefix>.bim` and
`<prefix>.fam`, adding the endings to the whole prefix: `-o release.v1`
gives `release.v1.bed`. Without `-o`, the prefix is the input name minus one
genotype-file ending (`.vcf.gz`, `.vcf.bgz`, `.txt`, `.csv`, `.tsv`, `.vcf`,
`.zip` or `.gz`), so `my.genome.txt` gives `my.genome`.

A1/A2 come from the export's own build: on GRCh37, the rsID identity table's
GRCh37 REF/ALT for the record's rsID at the record's position, chosen for
each record on its own; on GRCh38, gnomAD. A call that fits neither allele
of that pair as written is exported as missing and counted. On GRCh37, a
site the identity table has no pair for uses gnomAD's pair when it is not
A/T or C/G and the call fits it as written: gnomAD holds GRCh38 alleles, and
a site GRCh38 writes on the other strand cannot fit, so it is never written
missing for that. Alleles are never complemented. Every other site keeps the
fallback coding (a homozygous call is written with A2 `0`); the console
counts the sites coded from the identity table, from gnomAD and with the
fallback. Other builds always use the fallback.

A variant without an rsID gets the ID `chr:pos:ref:alt` from the input, with
`.` for an unknown part and `:2`, `:3` and so on when that ID is already
taken by an ID the input supplies or an earlier generated one; the console
says how many. IDs the
input supplies are not changed.

## Output file safety

JSON/HTML reports, annotated VCF and PLINK exports reject destinations that
alias the input. Report and PLINK checks include the former fixed `.tmp`
names, although writers now allocate unique owned temporary files. Unrelated
old temp files survive both successful and failed writes.

Output content and format are unchanged. VCF remains streamed with UTF-8/LF;
PLINK retains its existing bytes and order. Publication is atomic per file,
with all PLINK content staged and closed before its first replacement. A
later PLINK replacement failure can leave a mix of complete new files and
unchanged prior files. See [ADR-0056](adr/0056-input-output-file-safety.md).

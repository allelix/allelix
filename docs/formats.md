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
`N` in a one-letter or two-column file); all are read as no-calls. An rsID
and genotype list is the exception: it takes only `--`, `-`, `00` and `0`
(see [rsID and genotype lists](#rsid-and-genotype-lists)).

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

A few characters typed before the `#` of a file's first line (for example
`X#AncestryDNA raw data download`, an edited download) are ignored, up to
8 of them, with a warning (`input_stray_characters`). The rest of the line
still has to be the vendor's banner, so this never makes Allelix read a
file it would not otherwise recognise.

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
`decodeme`, `myheritage`, `ftdna`, `ftdna_famfinder`, `ftdna_illumina`,
`vcf` and `rsid-genotype` (an rsID and genotype list, read only by
`analyze`, `methylation` and `pharmacogenomics`). This is an expert override, at your own risk. Allelix cannot check what it would normally check:

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

An explicit homozygous-reference call (`0/0`, `0|0`, or a haploid `0`;
ALT a base or `.`) that passes the quality gate is read by ClinPGx,
SNPedia and the GWAS Catalog as the same call from a chip file would be:
a GWAS association whose effect allele is the reference base, or a ClinPGx
or SNPedia row for that genotype. This applies at the placements of their
rsIDs, which come from the rsID identity table. A record is evidence only
for the variant it describes: a single-base change from any record with
that REF base at the position, an insertion, deletion or MNV only from a
record whose trimmed REF/ALT is the target's (no left-alignment), called
as its full REF, homozygous (one copy when the record is haploid). A record naming an rsID stands for that rsID
where the table places it; an unnamed record stands for every rsID placed
at its position that it is evidence for; a single-base record
naming an rsID stands for it at its own position even without the
identity table, and a placement elsewhere is counted as an identity
conflict. Where another record changes the position (a carried or unknown
ALT, as in a split multi-allelic site; a record the quality gate rejects;
a carried insertion at its anchor base; shared padding is no change), no
reference call is made. The `contradictory_calls` ploidy warning counts
records: every non-reference record, and an explicit reference record
only at a target site of an enabled source (a watched placement, a
single-base record naming a source rsID, or a ClinVar reference-carried
site), at the trimmed site of a padded single-base record and the written
position of any other. A reference record elsewhere is not counted. gVCF reference blocks are not expanded, and a
position the file does not mention is never read as reference (ADR-0064).

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
Reports in schemas 1–12 can be used as `--diff` baselines; versions 13
and later are rejected. `--diff` identifies a finding by its source, variant
(rsID, or locus without one), condition and allele, whatever the order of
either report. Its classification, magnitude, genotype and the text
generated from them are compared, not used to identify it: a ClinVar
reclassification or a reworded SNPedia, ClinPGx or GWAS text is one changed
finding, not a removed and a new one. An entry from an older report without
an allele is matched without the allele.

A baseline finding is removed only when the current run has no finding with
its identity at all, before the display floors. A finding still there but
now scored below `--min-magnitude` (a Pathogenic 9.0 now capped to 4.0, say)
is changed, shown with its old and new magnitude, and marked below the
display floor; it stays out of the annotation list. New findings are the
ones the report shows. Schema 11 adds `below_display_floor` to each changed
entry in the JSON `diff` block, which also carries `previous_genotype_match`
(null when the baseline entry had none); the terminal shows `AG -> AA` when
it differs. A SNPedia finding on a sequence allele from a report before
2.6.0 shows as changed here once (`GACGAC -> GAC/GAC`): its genotype now
keeps a slash between the alleles, as every source writes it. When the baseline was made from a different input file, the run
says so with a warning.

## Annotated VCF: reference findings and re-annotation

`analyze --vcf-out` writes a ClinVar finding for an allele the reference
itself carries (see above) in its own field, `ALLELIX_CLINVAR_REF`
(`Number=1`, `significance|rsID|gene`, VCF schema 0.6.0). It holds the
highest-ranked such finding for the record, and never sits in an ALT slot.
`ALLELIX_CLINVAR` stays `Number=A` and describes carried ALT alleles only;
a record can have both, as a `0/1` at a site with findings for REF and ALT.

The field goes on every record whose call carries the reference allele:
a het (`0/1`, `1|0`), a half call (`0/.`), a hom-ref (`0/0`, haploid `0`),
and a single-position record with ALT `.` (an all-sites VCF), whose site is
exact. A record whose REF has more than one base gets it when that REF
starts at or covers the site. gnomAD, AlphaMissense and CADD values for the
reference allele are not written, since they would read as values for an
ALT.

When the input's genotypes are imputed (a Genes for Good imputed VCF, or a
VCF whose header declares imputation fields), the header has one more
line, `##ALLELIX_InputMethod=<ID=array-imputed|vcf-imputed,Description="...">`
(VCF schema 0.7.0), as the JSON report's `input.method` does. A measured
input has no such line. `export plink` says so on the console for an
imputed input, since PLINK files carry no provenance.

A gVCF reference block (`END=`) is never stamped, because one line covers
many positions; its finding is in the JSON and HTML reports only.
ClinPGx, SNPedia and GWAS Catalog findings from an explicit reference call
(see Reference calls in VCFs) are not stamped either; they are in the
JSON, HTML and terminal reports only (ADR-0064).

Every Allelix field, ALT findings, enrichment and recovered rsIDs included,
goes only on the record the analysis read a call from that carries its
allele. The selected sample (`--sample`), the quality rules (`--min-gq`, DP
0, FILTER and the sample's FT) and the mitochondrial exclusions apply as in the report. A
second record at the same CHROM, POS, REF and ALT that gave no such call (a
low-quality, no-call or hom-ref duplicate) borrows nothing. Carrying the
allele is not enough either: a record gets only the findings its own call
produced, so a `1/1` duplicate of a `0/1` record does not take a finding for
the `0/1` genotype, in any source, enrichment or `ALLELIX_CLINVAR_REF` field,
and a record with another REF at the same position (`GA>A` beside `G>A`)
keeps its own finding.
A record whose call was rejected for contradicting its rsID's identity
placement (the `identity_locus_conflict` warning) gets none.

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
  caller's FILTER, the sample's FT and genotype quality are honoured. Reference blocks with
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

A file re-saved by a spreadsheet is read too: rows separated by spaces or
by commas instead of tabs, fields in quotes, comment lines padded with
commas or quoted. Keep the original download if you can; a spreadsheet
can change values it thinks are numbers or dates.

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

A row with no tab is split on commas as CSV; a row with neither tabs nor
commas is split on spaces, so a file that a text editor or spreadsheet
re-saved with spaces between the columns still reads.

### `D` and `I` calls

23andMe writes an indel call as letters, not bases: `D` for a deletion
and `I` for an insertion (`DD`, `DI`, `II`; AncestryDNA writes the two
alleles in separate columns). Allelix does not interpret these calls,
including those at known disease variants such as BRCA1 5382insC, BRCA2
6174delT and CFTR F508del. Allelix has found no published specification
of which allele `D` and `I` stand for at each probe, and has no sample
with a confirmed genotype to check one against. A `D`/`I` call at a
high-value rsID, or at one ClinVar classifies Pathogenic or Likely
pathogenic, is named in the report as not checked: gene, variant name,
rsID and call, in JSON `high_value_unchecked` (reason `di_call`), and
counted under the `report_warnings` code `di_calls_not_checked`. Other `D`/`I` calls are
only counted. Allelix neither confirms nor rules out carrier status for
these variants; that needs clinical testing.

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
are phased. A confident `0/0` record is read like the chip's
homozygous-reference call by ClinPGx, SNPedia and the GWAS Catalog, and by
ClinVar where the reference itself carries a ClinVar variant allele (see
Reference calls in VCFs). The
VCFs' header states no build that Allelix reads, so they may need
`--build grch37`. The text file carries every call and states its build,
so it is the one to pass.

## deCODEme

deCODEme, deCODE genetics' consumer service (closed in 2013), delivered
`deCODEme_scan.csv`: no comment block, a column line
`Name,Variation,Chromosome,Position,Strand,YourCode`, then one SNP per
row. Allelix recognises the file by that column line.

`YourCode` is the genotype (`--` for a no-call) and `Strand` says which
strand it is written on. Allelix compares calls on the plus strand, as
every other format writes them, so a `-` row's call is complemented (`TG`
becomes `AC`). Sources keep their own orientation: SNPedia pages say plus,
minus or unknown, and a minus page is complemented before matching. This
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

`allelix compare first.vcf second.vcf` matches each record by its allele
identity:
- **A named record** (one with an rsID) by its rsID, input REF and the
  complete input ALT field. The two records of a split multi-allelic site
  (`G>A` and `G>T` under one rsID) each meet their own partner, whatever
  the order of the records.
- **An unnamed record** by its chromosome, position, REF and ALT field,
  and only when both files have the same known build.

Each file's build is decided first, as for every command (ADR-0044), and
its records are read through that decision. Mitochondrial records the
decision leaves out (every MT record under an hg19-length chrM; in a
GRCh37 VCF, records named `chrM` or `M` unless that name is declared with
the 16,569 bp length) enter no match, count, concordance or chromosome
breakdown, and the file's notice counts them. A file whose calls were all
left out is not empty: it gets the notice and no eligible calls.

The alleles are trimmed first: REF and every ALT lose the bases they all
share at the end, then at the start, so `CG>CA` and `G>A` are one record.
There is no left-alignment. The ALT list is literal: reordered lists,
different uncalled alternatives, and a split record against a merged one
(`A` and `T` against `A,T`) are distinct records. Quality-masked and
explicit no-calls keep their original allele identity. No liftover,
allele decomposition or genotype-based ALT inference is performed.

A chip file has no REF or ALT. Its call is compared with the call under
the same rsID in another chip file. Against a VCF, it is compared with the
one record under that rsID whose REF and ALTs hold the call, directly or
on the other strand. When the VCF holds records under that rsID but none
or several of them hold the call, the call is not compared, and `compare`
says how many and why. A chip call whose rsID the VCF does not hold is
counted as only in its own file.

The coverage table separates parsed named/unnamed observations from
eligible unique named/unnamed keys. It labels shared unnamed variants
compared by position and counts omitted unnamed observations per file.
Unknown/different builds or missing input allele identity omit unnamed
calls. Nonempty vendor identifiers remain on the named path; named and
unnamed records do not match one another. Repeated records with the same
identity keep the last observation, so conflicting duplicates of the very
same record can depend on their order. Concordance percentages use shared
eligible keys; no shared keys is displayed without a percentage.

`allelix extract` shows confident reference calls too. An explicit `0/0`
(or `0|0`, haploid `0`, ALT `.`) record is shown as REF/REF. A gVCF reference
block covering the rsID's position is shown as that site's REF/REF. A
filtered, low-quality or zero-depth record is a no-call. In a VCF without
rsIDs, a no-call (`./.`, haploid `.`) or a half call carrying only REF
(`0/.`) at the rsID's position is shown too, as is a padded record whose
ALT trims to that position (ADR-0061).

When a file holds several records for one requested rsID, extract shows one,
by the same order with or without an index (ADR-0054):
1. a record at one of the rsID's placements (on GRCh38 gnomAD's, or the
   rsID identity table's when gnomAD does not place the rsID; on GRCh37 the
   identity table's) beats one elsewhere, even as a no-call; a padded
   record is at a placement its ALT trims to, and ties there with the same
   call written trimmed;
2. among those, a usable call beats a no-call or a rejected record;
3. then the header's contig order, POS and line order.

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

## rsID and genotype lists

Some services export only each marker's ID and its call, two tab-separated
columns (IYG files, for one):

```
# Comment lines are allowed
rs4477212	AG
rs3094315	C
rs3131972	--
```

The file has no chromosome, position or build. Allelix reads it as this
format only when every one of its first 50 data lines is an ID and a call
(one or two of A, C, G, T, D, I, or a no-call `--`, `-`, `00`, `0`; double
quotes around a field are removed), and at least three of those IDs, and
at least 80% of them, are rsIDs.
Two letters are a diploid call and one letter a haploid call. A column line,
a third column, or rows separated by spaces or commas make it another format.

Each row is placed on **GRCh38** from the rsID identity table (ADR-0069), so
the table must be installed (`allelix db update`); without it the run stops.
`--build grch38` is accepted and any other `--build` stops. A row is placed
at its rsID's locus when exactly one strand reading of its call fits the
site's alleles (REF plus every ALT): as written, or complemented to the
forward strand when only the complement fits. Each row is read on its own.
The report counts the rows placed, those complemented, and those not placed,
by reason:

| Reason | Meaning |
|---|---|
| `not_an_rsid` | the ID is not an rsID |
| `not_in_table` | the identity table holds no single GRCh38 locus for the rsID |
| `palindromic` | an A/T or C/G site: the strand cannot be told |
| `both_fit` | the call and its complement both fit a site with more alleles |
| `no_fit` | neither reading fits |
| `di_call` | a D/I call |

A row that is not placed has no chromosome or position (null in JSON). ClinVar,
gnomAD, AlphaMissense and CADD, which need a position, skip it. ClinPGx,
SNPedia and the GWAS Catalog still look it up by rsID. ClinPGx reads its
call as written, and so do SNPedia's indel and i-ID rows. SNPedia's
single-base genotypes and the GWAS Catalog's risk alleles keep their own
strand checks against the identity table's alleles (ADR-0066), and these
check the page or the risk allele, not the whole call that placement
checked. A SNPedia page with a stated orientation is used where that
orientation fits the table's alleles. A page without one, and a risk
allele, are read on the one strand that fits; where both readings fit or
neither does, or the table holds no single locus for the rsID, that source
abstains and counts it in `strand_abstentions`. A GWAS row without a risk
allele takes no strand check. So the reason a row was not placed does not
tell whether these sources report something for it, and an unplaced row is
not counted as an abstention by itself. Without a locus, a call that is
homozygous for the reference is not suppressed either, so findings for
that genotype can appear. A `di_call` row reaches no source: a D/I call is
counted and never interpreted (ADR-0070).

Only `analyze`, `methylation` and `pharmacogenomics` read these files.
`stats`, `extract`, `export` and `compare` stop with a message.

Two other layouts without positions (`SNP,Chr,Allele1,Allele2`, and
`SNP_ID`, `Allele1`, `Allele2` columns) are not read.

## PLINK export

`allelix export plink` writes `<prefix>.bed`, `<prefix>.bim` and
`<prefix>.fam`, adding the endings to the whole prefix: `-o release.v1`
gives `release.v1.bed`. Without `-o`, the prefix is the input name minus one
genotype-file ending (`.vcf.gz`, `.vcf.bgz`, `.txt`, `.csv`, `.tsv`, `.vcf`,
`.zip` or `.gz`), so `my.genome.txt` gives `my.genome`.

A1/A2 come from the export's own build where it has a pair (on GRCh37, the
gnomAD fallback below can supply one), chosen for each record on its own:
on GRCh37, the rsID identity table's GRCh37 REF/ALT for the record's rsID at
the record's position; on GRCh38, gnomAD's REF/ALT for the record's rsID at
the record's position, preferring a pair that is not A/T or C/G. Two records
with one rsID (split multi-allelic rows, or one rsID at two positions) never
share a pair. A call that carries two ALTs (`1/2`) fits no single pair; when
both its alleles are alleles of the site, it keeps its observed pair through
the fallback coding below. On GRCh37, a
site the identity table has no pair for uses gnomAD's pair when it is not
A/T or C/G and the call fits it as written: gnomAD holds GRCh38 alleles, and
a site GRCh38 writes on the other strand cannot fit, so it is never written
missing for that. Choosing the pair reads the call as written. Every other
site keeps the fallback coding (a homozygous call is written with A2 `0`);
the console counts the sites coded from the identity table, from gnomAD and
with the fallback. Other builds always use the fallback.

Against the chosen pair, a call is coded by the one strand reading that
fits it (ADR-0066):

- At a pair that is not A/T or C/G, the call is coded as written when its
  alleles are the pair's, and by its complement when only the complement's
  are (an array probe reported on the minus strand: `CC` at a G/A site is
  coded `GG`). The console counts the calls coded by their complement. A
  call that fits neither reading is exported as missing and counted.
- At an A/T or C/G pair the strand cannot be told: a call that fits as
  written keeps the fallback coding, and any other call is exported as
  missing and counted.

A third allele the pair does not list (for example G>C at a G/A site, called
on the minus strand) would be coded as its complement; that is the accepted
risk of reading minus-strand probes.

A variant without an rsID gets the ID `chr:pos:ref:alt` from the input, with
`.` for an unknown part and `:2`, `:3` and so on when that ID is already
taken by an ID the input supplies or an earlier generated one; the console
says how many. IDs the
input supplies are not changed.

## Output file safety

JSON/HTML reports, annotated VCF and PLINK exports reject destinations that
alias a file the run reads: the input, and for `analyze` the `--diff`
baseline and the `--filter-file`. Report and PLINK checks include the former fixed `.tmp`
names, although writers now allocate unique owned temporary files. Unrelated
old temp files survive both successful and failed writes.

Two outputs of one run never name the same file. `analyze` refuses an
`--output` and a `--vcf-out` that resolve to one file, and `export plink`
refuses a prefix whose `.bed`, `.bim` and `.fam` are not three files. That
covers the same path spelled two ways, a symlink (even a dangling one), a
linked parent folder and a hardlink. The check runs before any work starts,
and nothing is written.

A new output file gets the mode your umask gives any new file (`0644` under
the usual `022`); an existing report or PLINK file keeps its mode. `--vcf-out`
never replaces an existing file: the run stops instead.

Output content and format are unchanged. VCF remains streamed with UTF-8/LF;
PLINK retains its existing bytes and order. Publication is atomic per file,
with all PLINK content staged and closed before its first replacement. A
later PLINK replacement failure can leave a mix of complete new files and
unchanged prior files. See [ADR-0056](adr/0056-input-output-file-safety.md).

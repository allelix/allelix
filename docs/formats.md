# Provider file formats

What each provider delivers, which file Allelix needs, and what Allelix
does with the provider-specific parts of the file. The general format
table is in the [README](../README.md#supported-formats).

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
##fileDate=July-28-2026 at 02:21 UTC
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
##timestamp=2026-06-26 17:54:08 UTC
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

## MyHappyGenes (Tempus)

Tab-delimited, five columns. The header states the build as "37.1", but
real exports have shipped GRCh38 coordinates under that label. Allelix
therefore checks the label against the positions of well-known SNPs across
the whole file: confident positions win, with a warning when they disagree
with the header; weaker evidence against the header stops the run and asks
for `--build` (ADR-0021, ADR-0044).

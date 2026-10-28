# Converting GRCh36 Files to GRCh38

Allelix recognizes GRCh36 (hg18) position evidence. Three matching seed SNPs
and an 80% winning ratio make that evidence confident. If GRCh36 leads but
does not clear the confidence threshold, Allelix still fails closed to
GRCh36 rather than risk querying a GRCh37 or GRCh38 position cache; the
terminal and HTML diagnostics label the evidence tentative. ClinPGx
annotations, which match by rsID, still work. GWAS Catalog risk alleles
and single-base SNPedia genotypes are not assessed on GRCh36 (their strand
cannot be checked), and the report counts them. ClinVar is skipped because no GRCh36 cache is available. gnomAD,
AlphaMissense, and CADD are also skipped because their current Allelix
caches are GRCh38-only; using their values on another assembly would risk
wrong-allele enrichment.

If Allelix cannot tell the file's build by itself (for example an export
with no build line whose rows include none of the known SNPs it checks),
it stops and asks. For a GRCh36 file, pass `--build grch36`: the run then
uses GRCh36 with the same limits as above.

To get full ClinVar coverage, convert your file's coordinates to GRCh38 (or
GRCh37) using one of these tools, then re-run `allelix analyze` on the
converted file.

## UCSC liftOver (command-line)

Download the liftOver binary and chain file:

```bash
# Download liftOver (Linux)
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

# Download the GRCh36 → GRCh38 chain file
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
```

Convert a BED file of positions:

```bash
./liftOver input.bed hg18ToHg38.over.chain.gz output.bed unmapped.bed
```

For a genotype file, updating the positions alone is not enough:

1. Write one BED line per SNP. BED is zero-based and half-open, so a SNP
   at position 1,000 is `chr1  999  1000  rs123  0  +` (BED6: name, score,
   strand). Keep the rsID in the name column and write the strand as `+`,
   so the output shows where the mapping is inverted.
2. Lift over. Each output line gives the new chromosome, the new
   coordinates and the strand: `-` where the region is inverted in the new
   assembly. A line in `unmapped.bed` could not be lifted.
3. Rewrite each row's chromosome and position from the output (the
   position is the BED end, or start + 1). Drop the rows that did not map.
4. Where the output strand is `-`, complement that row's genotype (A↔T,
   C↔G). Allelix reads array calls on the forward strand of the build it
   decides, so an uncomplemented call there reads as a different allele.
   An A/T or C/G SNP looks the same either way, so check those rows most
   carefully.
5. Change the file's build line (or its `##reference` / `##contig` header
   for a VCF) to the new build.

For a VCF, CrossMap's `vcf` mode below does this allele-aware work itself
(it needs the target reference sequence). UCSC liftOver documentation:
https://genome.ucsc.edu/cgi-bin/hgLiftOver

## CrossMap (Python)

Install via pip:

```bash
pip install crossmap
```

Convert coordinates:

```bash
# Download chain file
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz

# Download and unpack the hg38 reference sequence (needed for a VCF;
# about 1 GB to download, 3 GB unpacked)
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Convert a VCF
CrossMap vcf hg18ToHg38.over.chain.gz input.vcf hg38.fa output.vcf

# Convert a BED
CrossMap bed hg18ToHg38.over.chain.gz input.bed output.bed
```

CrossMap documentation: https://crossmap.readthedocs.io/

## Notes

- A small fraction of positions (~0.1–1%) may fail to lift over due to
  structural rearrangements between assemblies. These are reported in the
  unmapped output file.
- After conversion, verify the build through the analysis pipeline:
  `allelix analyze converted_file.txt --output converted_report.html`
  should print a GRCh38 build banner. `allelix stats` runs the same build
  decision and shows the build and where it came from (ADR-0044).
- A converted array file usually still carries its old "build 36" header
  line. Confident GRCh38 positions override it with a mismatch notice; weaker
  evidence against it stops the run. A converted VCF must also get GRCh38
  `##contig` lines: if it keeps the GRCh36 or GRCh37 contig lengths,
  Allelix stops rather than mix builds.
- If fewer than three seed SNPs are present or the signal is mixed, the
  position evidence is tentative. After independently verifying a successful
  conversion, pass `--build grch38`.
- If your genotyping provider offers re-export on a newer build, that is
  simpler than liftover. Check your provider's download settings.

# Converting GRCh36 Files to GRCh38

Allelix recognizes GRCh36 (hg18) position evidence. Three matching seed SNPs
and an 80% winning ratio make that evidence confident. If GRCh36 leads but
does not clear the confidence threshold, Allelix still fails closed to
GRCh36 rather than risk querying a GRCh37 or GRCh38 position cache; the
terminal and HTML diagnostics label the evidence tentative. rsID-keyed
ClinPGx, GWAS Catalog, and SNPedia annotations can still work. ClinVar is
skipped because no GRCh36 cache is available. gnomAD,
AlphaMissense, and CADD are also skipped because their current Allelix
caches are GRCh38-only; using their values on another assembly would risk
wrong-allele enrichment.

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

For genotype files, extract positions into BED format first, lift over, then
update the positions in your genotype file. UCSC liftOver documentation:
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
  should print a GRCh38 build banner. `allelix stats` reports the parser's
  declared/default build and does not run position-based build detection.
- If fewer than three seed SNPs are present or the signal is mixed, detection
  remains tentative. After independently verifying a successful conversion,
  pass `--build grch38` rather than relying on a tentative signal.
- If your genotyping provider offers re-export on a newer build, that is
  simpler than liftover. Check your provider's download settings.

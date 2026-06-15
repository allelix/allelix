# Converting GRCh36 Files to GRCh38

Allelix detects GRCh36 (hg18) genotype files automatically. rsID-based
annotations (ClinPGx, GWAS Catalog, SNPedia, gnomAD) work normally, but
ClinVar annotations are skipped because ClinVar indexes by genomic position
and no GRCh36 cache is available.

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
- After conversion, verify the build is detected correctly:
  `allelix stats converted_file.txt` should show "Build: GRCh38".
- If your genotyping provider offers re-export on a newer build, that is
  simpler than liftover. Check your provider's download settings.

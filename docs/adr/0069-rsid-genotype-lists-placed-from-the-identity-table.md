# ADR-0069: rsID and genotype lists are placed on GRCh38 from the rsID identity table

**Status:** Accepted
**Date:** 2026-10-31
**Builds on:** ADR-0040 (the rsID identity table), ADR-0044 (the build decision), ADR-0066 (both strand readings together)

## Context

Some services export only each marker's ID and its call, two tab-separated
columns with no chromosome, position or build (IYG files, for one). Every
other reader needs a position: ClinVar, gnomAD, AlphaMissense and CADD are
keyed by locus, and the build decision (ADR-0044) reads known-SNP positions.
Before, such a file was refused as an unknown format.

The rsID identity table places each annotator rsID on GRCh37 and GRCh38 with
its REF and every ALT (ADR-0040). It can place these rows, provided the
call's strand can be told: the file does not say which strand it reports.

## Decision

**Detection is strict** (`rsid-genotype`, registered after every other
reader). Within the first 50 data lines, every line has exactly two
tab-separated fields (double quotes around a field are removed): an ID with
no spaces and a call of one or two of A, C, G, T, D, I, or a no-call (`--`,
`-`, `00`, `0`). At least 3 lines, and at least 80% of them, have an rsID.
Only `#` comments and blank lines may come between them; a column line, a
third column, a space- or comma-separated row or a call such as `NN` is not
this format. The two other position-free layouts seen (`SNP,Chr,Allele1,
Allele2` and `SNP_ID, Allele1, Allele2`) stay rejected: their source is not
known.

**GRCh38 only.** The file carries no build evidence, so its rows are placed
on GRCh38. `--build grch38` agrees; any other `--build` stops.

**The identity table is required.** Without it the run stops before
anything is read; an installed rsID-keyed source (ClinPGx, SNPedia, the GWAS
Catalog) does not make the run go ahead without placement.

**Each row is placed on its own, by ADR-0066's rule** (`forward_strand_reading`,
the one shared helper), against the site's forward alleles on GRCh38: REF
plus every ALT, at the rsID's single locus.

- Only the call as written fits: placed as given.
- Only its complement fits: placed, complemented to the forward strand, and
  counted (`complemented`).
- Both fit (`palindromic` at an A/T or C/G site, `both_fit` at a site with
  more alleles) or neither fits (`no_fit`): unplaced, counted.

The file's other rows never decide a row's strand: there is no file-strand
inference. A no-call is placed without a check, as it has no alleles to
orient. A row is also unplaced, counted, when its ID is not an rsID
(`not_an_rsid`), when the table holds no single GRCh38 locus for it
(`not_in_table`), or when the call is a D/I call (`di_call`).

**Unplaced rows still reach the rsID-keyed sources**, except a `di_call`
row: a D/I call reaches no source (ADR-0070). They keep chromosome
"" and position 0: position-keyed sources skip them, and ClinPGx, SNPedia
and the GWAS Catalog look them up by rsID under their own strand checks.
Those checks read a call with no position as written: SNPedia and the GWAS
Catalog check the source's alleles against the rsID's forward alleles in the
table (they abstain, counted as `not_in_identity_table`, when the table does
not hold the rsID), and ClinPGx matches the call exactly first. A finding
from an unplaced row publishes a null `chromosome` and `position`, not a
made-up locus.

**An unknown position is unknown, never a shared one.** Position 0 is no
locus, so two unplaced rows never stand at the same position:

- They raise no ploidy contradiction (ADR-0062).
- They do not split an annotation run; only a shared rsID or a real shared
  position does.
- A finding from an unplaced row finds its own call by rsID, for its record
  number and ploidy.

**Reports:** JSON `input.placement` holds `{build, rows, placed,
complemented, unplaced}`, with `unplaced` by reason; `input.build_source` is
`rsid_placement`. The terminal and HTML say how many rows were placed and why
the others were not.

**Only the analysis commands read these files** (`analyze`, `methylation`,
`pharmacogenomics`). `stats`, `extract`, `export` and `compare` stop on them
with a message that says so, rather than writing rows without a locus.

## Consequences

- An rsID and genotype list can be analyzed. Where a call fits both strands
  at its site, it is not placed. Its findings from the rsID-keyed sources
  may still appear, without a locus, but ClinVar and the frequency and score
  sources do not see it.
- An unplaced row's call is read as written by the rsID-keyed sources,
  although its strand was not told. A `both_fit` call reported on the minus
  strand can then match a source genotype on the wrong strand. This is the
  accepted cost of keeping unplaced rows visible; the console, HTML and JSON
  count every unplaced row by reason.
- The homozygous-reference suppression (ADR-0023) needs the call's locus, so
  it never applies to an unplaced row. An unplaced call that is homozygous
  for the reference can then show the rsID-keyed sources' findings for that
  genotype, which the same call with a position would not show.
- The identity table's coverage bounds what is placed: an rsID outside it
  stays unplaced.
- A third allele the table does not list, called on the minus strand, would
  be placed as its complement: the same accepted risk as ADR-0066.

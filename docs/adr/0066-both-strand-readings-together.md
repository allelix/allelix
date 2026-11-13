# ADR-0066: Both strand readings of a source's alleles are evaluated together

**Status:** Accepted
**Date:** 2026-10-31
**Amends:** ADR-0024 (the 2.5 abstain-first strand rule)

## Context

The GWAS Catalog's risk alleles and SNPedia's genotypes are read against the
site's forward alleles on the input's build: REF plus every ALT in the rsID
identity table (ADR-0040). ADR-0024's 2.5 rule read them in a fixed order:

1. The site is not A/T or C/G and the alleles are forward: read them as
   written, even where their complements are forward too.
2. Otherwise, if the complements are forward: complement them.
3. Abstain at A/T and C/G sites, where the table does not place the call,
   and where neither reading fits.

Rule 1 assumed that a source writes the forward strand whenever it can. At a
site with three or more alleles that holds a base and its complement (for
example forward `{A, G, T}`, which holds both A and T), a risk allele stated
on the reverse strand also fits as written. It was then read on the wrong
strand: the association went to the wrong genotype, and a carrier of the
real risk allele missed it.

Rule 3 also made a SNPedia page that states its orientation (plus or minus)
abstain at every A/T or C/G site. Yet such a page needs no reading: once its
genotype is turned to the forward strand, it either is among the site's
forward alleles or does not fit.

## Decision

**Both readings are evaluated together** (`forward_strand_reading` in
`allelix/utils/allele.py`), for the GWAS Catalog's risk allele and for
SNPedia rows without a stated orientation:

- Only the alleles as written are forward: read them as written (`direct`).
- Only their complements are forward: complement them (`complement`).
- Both are forward: the strand cannot be told from sequence, and the row
  abstains, counted as `palindromic` at an A/T or C/G site and as `both_fit`
  at a site with more alleles.
- Neither is: the row abstains, counted as `no_fit`.

**Each SNPedia page is read on its own.** A page's genotype is the only
evidence for its strand. Another page at the same rsID never decides it:
pooling the alleles of several pages would let an unrelated page settle,
or hide, a page's reading. Before this change, the alleles of every page
without a stated orientation were pooled into one reading. Abstentions are
counted once per site and reason.

**A stated orientation is independently known.** A SNPedia row whose page
states plus or minus is turned to the forward strand and kept wherever that
genotype is among the site's forward alleles, A/T, C/G and both-fit sites
included; otherwise it abstains as `no_fit` (`stated_orientation_fits`).

**ClinPGx is outside the rule.** It matches the user's genotype to the source
genotype exactly first, and tries the complement only when the call's REF
confirms reverse orientation at a site that is not A/T or C/G (ADR-0035).
Its genotypes are stated on the forward strand. Of the 8,304 two-base
genotypes in the ClinPGx cache, 4,117 show their strand against the identity
alleles: 4,114 read as written and 3 only as their complement. Applying the
both-readings rule there would make about half of its rows abstain (1,032 of
1,975 on the HG002 GRCh37 benchmark) where the exact match is correct.

**PLINK export uses the same reading** (`allelix/exporters/plink.py`,
`_pair_reading`). A call is coded against its site's REF/ALT pair, from the
identity table on GRCh37 or gnomAD on GRCh38:

- At a non-palindromic pair, the one reading that fits is used. As written,
  the call is coded directly. By its complement, the call is coded as the
  complemented alleles and counted on the console: "N call(s) coded by their
  complement". When neither reading fits, the call is written missing and
  counted, as before. Both cannot fit: a non-palindromic pair shares no base
  with its complement.
- At an A/T or C/G pair, the strand cannot be told: a call that fits as
  written keeps the fallback coding, and one that does not is written
  missing, as before.
- The pair selection is unchanged, and reads the call as written:
  - GRCh38 (gnomAD): among several records at the site, the first the call
    fits as written, preferring one that is not A/T or C/G. When none fits,
    the first record, unless every called allele is one of the site's (a
    `1/2` call), which gets no pair and the fallback coding.
  - GRCh37 (identity table): a single pair is used as is. With several, a
    pair is used only when exactly one fits the call as written; otherwise
    the call gets no pair and the fallback coding, so a minus-strand call
    at such a site is not complemented.
  - GRCh37, no identity pair: gnomAD's pair is used only when it is not A/T
    or C/G and the call fits it as written, because gnomAD's alleles are on
    the GRCh38 strand.

This reverses the 2.5 rule that export never reads a call on the other
strand, for non-palindromic pairs only. Arrays such as 23andMe report some
probes on the minus strand. A minus-strand `CC` at a G/A site was written
missing; it is now coded `GG`. The accepted risk: a third allele the pair
does not list (for example G>C) would be coded as its complement.

**Unchanged:** the forward alleles themselves; the `not_in_identity_table`
abstention where the table does not place the call (including GRCh36
input); rows without a risk allele (the 3.0 cap); SNPedia's 23andMe i-IDs
and indel or D/I genotypes, which keep the page's orientation with the
palindrome guard; ClinVar's own matching.

**Report.** `strand_abstentions` gains the `both_fit` reason for each
strand-checked source (`gwas`, `snpedia`). The console and HTML strand notice
name it: "alleles fit both strands at the site".

## Consequences

- Where both readings fit, a GWAS row abstains, counted under `both_fit` or
  `palindromic`, instead of being read as written. Where exactly one reading
  fits, the result is the same as under the 2.5 rule.
- SNPedia pages that state their orientation are read at A/T and C/G sites.
- Measured against 2.6 dev before this change, on HG002 GRCh38 chr22 and
  three array files:
  - At the default floors, 0 to 4 of the 0 to 11 displayed GWAS rows per
    file are lost, and 0 to 4 SNPedia rows gained.
  - With every floor at 0, 21.5% to 22.0% of the GWAS rows on the array
    files abstain as `both_fit` (5,123 to 6,646 rows per file), and 36% on
    chr22 (358 of 990). The identity table lists three or four alleles at
    many sites, so both readings fit there.
  - 2 to 5 GWAS rows per file appear. Each is a row the PheCode rollup
    (ADR-0024) used to hide: the best row of its rsID and parent phecode now
    abstains, so the next one is shown. No row on these files changes from
    one reading to the other.

## The trade-off

Abstaining at both-fit sites gives up rows that reading as written would
have got right, to avoid the ones it got wrong.

- **The Catalog is mostly forward.** Of the 343,557 distinct single-base risk
  alleles in the GWAS cache (GRCh38):
  - 190,675 fit only the direct reading and 2,459 only the complement;
  - 122,960 fit both at a site with more alleles, and 18,194 sit at A/T or
    C/G sites.

  Where a risk allele shows its strand, it is forward 98.7% of the time.
- **At both-fit sites, reading as written is wrong more often.** Where gnomAD
  can tell which ALTs are common, the risk allele resolves to the complement
  4.2% of the time. The test keeps the identity ALTs with gnomAD AF ≥ 0.001
  at the site. Of the 122,960 both-fit risk alleles, 20,573 have a gnomAD
  record; 15,478 of those then fit exactly one reading, 14,831 direct and
  647 complement. On one array file, 14 of the rows this change removes
  resolve to a complement the user does not carry: the old rows were on the
  wrong strand.
- **No frequency filter.** The local gnomAD cache is exomes only: 102,385 of
  the 122,960 (83%) have no record at their site. It would bring back none of
  the 9 default-floor rows lost on the files above, and it would make the
  result depend on whether a site is exonic. A common-allele site definition
  from a whole-genome frequency source would need its own ADR, superseding
  ADR-0024's "no allele-frequency rescue".

## Amendment 2.6.1 (2026-11-05): PLINK reads a VCF call as written

The complement reading in PLINK export is for chip probes, whose strand is
not known. A VCF or gVCF call is on the forward strand already, so
`export plink` reads it only as written (`forward_strand_calls` in
`allelix/exporters/plink.py`): at a pair that is not A/T or C/G, a VCF call
that does not fit the pair as written is exported as missing and counted,
not coded as its complement. In 2.6.0 a `T/T` call at an A/G site was
written as `A/A`. Chip calls keep the reading above, and the A/T and C/G
rule is unchanged.

## Amendment 2.7.0 (2026-11-09): ClinPGx's complement match reads the site

ClinPGx stays outside the rule for its exact match: a genotype that equals
the call as written matches, whatever the site reads. Its complement match
now also uses this reading, once per call: all the rsID's ClinPGx
genotypes together, read against the identity alleles at the call's locus
(for a call without a position, the rsID's one site).
Where only their complement fits, a single-base call (single-base carried
alleles, REF one base or not known) matches a row by its complement,
homozygotes included; an indel call's anchor letter is not evidence of a
sequence allele. Every other reading, and every other call, leaves the
complement match to the call's REF, as before (ADR-0035, amendment 2.7.0). No row abstains that
matched before, and no abstention is counted.

## Amendment 2.7.0 (2026-11-10): padded identity rows are read trimmed

The identity table stores a few single-base variants padded with a shared
base: `CG>CA` one position before `G>A`, or `GT>AT` at its own position
(20 such rows on GRCh37 and 70 on GRCh38 in the dbSNP b157 table, all
padded at the end). Read as stored, such a site had a multi-base allele
set at another position, so a GWAS risk allele or a SNPedia genotype
there found no forward alleles (`not_in_identity_table`) and an unnamed
call there recovered no rsID. Every identity lookup now reads such a row
as its trimmed single-base variant (`snv_form` in
`allelix/annotators/rsid_identity.py`, trimmed as the input is,
ADR-0047): the forward alleles above, the rows at a locus, each rsID's
placements and loci, the GRCh37-to-GRCh38 key translation, the sites
whose two builds' REFs differ and the normal allele on the other build
(so a padded reference-carried site gets the same normal-allele decision
as its trimmed form). A row padded before its variant is stored where a
lookup by the variant's position does not find it; such rows are found
once per run from the table's multi-base rows, whatever the padding
length (the b157 table has none). A genuine indel or MNV is read as
stored, and the strand readings and their abstentions are unchanged. The
reading happens when the table is read, so an installed cache needs no
rebuild.

## Amendment 2.7.0 RC2 (2026-11-10): ClinPGx reads a decided site one way

The ClinPGx amendment above says ClinPGx's exact match holds "whatever the
site reads", and that no row abstains that matched before. Both are
replaced for 2.7.0 (ADR-0035, amendment 2.7.0 RC2). Where this reading of
the rsID's ClinPGx genotypes at the call's site is `direct`, only the
exact match counts; where it is `complement`, only the complement does.
Elsewhere, a call that one genotype meets as written and another as its
complement meets neither, and is counted (`clinpgx_ambiguous_strand`).
SNPedia, the GWAS Catalog and placement keep the reading this ADR gives
them.

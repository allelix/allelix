# ADR-0046: ClinVar variants without an rsID are matched for VCF input only

**Status:** Accepted
**Date:** 2026-10-13
**Builds on:** ADR-0039 (locus-keyed ClinVar selection), ADR-0040 (rsID
identity table).

## Context

Until 2.5 the ClinVar loader skipped every variant whose `RS# (dbSNP)` is
`-1` or empty. ClinVar holds many such variants. In the 2026-10-10 data,
the GRCh37 cache gains 1,512,187 rows on 1,451,116 variants by keeping
them. 98,211 of those rows, on 91,391 variants, are P/LP: 23.4% of the
P/LP variants in the cache. A VCF that carries one of them
at the exact locus got nothing, with no warning.

Since ADR-0039, ClinVar rows are selected by locus, and a row annotates only
when its REF agrees with the call's REF and the user carries its ALT. For a
VCF record both alleles come from the file, so `(chrom, pos, ref, alt)`
identifies the variant without an rsID.

Array calls are different. Their REF is filled in from the rsID identity
table, from the rows of the call's own rsID (ADR-0040). A ClinVar row at the
same position without an rsID is not tied to that rsID: nothing confirms
that the probe measured the same variant, or on which strand.

## Decision

- The loader keeps a variant without an rsID, with `rsid` stored as `""`.
- `Variant.ref_from_input` is set by the VCF parser: REF was read from the
  input. A ClinVar row without an rsID annotates only such a call. Array
  calls never match it, with or without a REF from the identity table.
- The annotation keeps the input's own ID when ClinVar has none (a VCF ID
  such as `COSV...`), otherwise `""`.
- Lookups keyed by rsID skip `""`: the REF map behind `reference_for` and
  `normal_allele_for`, the locus resolver for rsID-less calls (it already
  ignored rows without an rsID), and the identity table's universe.
- A finding without an rsID is named by its locus:
  - Terminal, HTML and diff show `chrom:pos ref>alt` (for example
    `1:1000 A>G`) in place of the rsID.
  - HTML groups findings by that name, not by an empty rsID.
  - The diff keys on it, so two such findings with the same text do not
    collide.
  - Every JSON annotation carries `chromosome`, `position` and `ref`
    (schema 9), so a finding without an rsID still names its variant.
  - `--vcf-out` maps an annotation to its input record by locus and
    carried allele, not by rsID. The VCF INFO field writes `.` for no
    rsID.
  - `scripts/check_wrong_allele.py` places such a finding at its own JSON
    position instead of looking up an empty rsID.

## Consequences

- VCF users see ClinVar's rsID-less variants. On the HG002 GRCh37
  benchmark VCF (2026-10-10 data), a count from the cache gives the
  unfiltered report 917 more rows on 877 variants. 874 of them are benign
  or likely benign, and none is P/LP. The 20 rows that
  reach the default report are drug response (14 of them CYP2D6) and one
  established risk allele. Other genomes may carry P/LP variants without
  an rsID.
- The ClinVar cache grows by about 31%, from 4.83 to 6.34 million rows
  per build (about 1.2 GB per build). The ClinVar interpreter version is bumped
  (7 to 8), so an older cache without these rows is rebuilt.
- Array users see no change. A probe that measures such a variant still
  gets nothing from ClinVar until ClinVar or dbSNP links an rsID to it.
- JSON annotations are slightly larger: three fields more per row.

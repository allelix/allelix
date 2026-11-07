# ADR-0051: Count the array calls ClinVar cannot match, and panel coverage in four states

**Status:** Accepted
**Date:** 2026-10-16
**Supersedes:** the premise of ADR-0011 that arrays do not call indels (its
decision stands). **Builds on:** ADR-0039 (locus-keyed ClinVar selection),
ADR-0045 (ClinVar's aggregate), ADR-0050 (`is_pathogenic_aggregate`).

## Context

- **Two kinds of array call are never matched against ClinVar, and the
  report did not say so.**
  - A vendor probe without an rsID (23andMe `i5000123`) is not
    locus-matched (ADR-0039): nothing confirms that the probe tests the
    variant ClinVar records at that position.
  - An array writes an insertion or deletion as `D` or `I`. Those are not
    bases, so they never equal a ClinVar ALT.
  - A user could read "no pathogenic findings" as covering those sites. In
    one public 23andMe v5 file (openSNP), about 6,000 i-probes sit on loci
    with ClinVar Pathogenic / Likely pathogenic rows, and 970 calls are D/I.
- **ADR-0011's premise is false.** It says arrays do not call indels;
  23andMe does, as D/I codes. Its decision (a ClinVar indel row never
  matches a single-base array call) is still right: a single base cannot
  confirm an indel, and a D/I code cannot say which indel.
- **Array input had no parse counts.** VCF input has `input.vcf_quality`;
  an array file had nowhere to report what its parse saw.
- **Panel coverage folded no-calls into "found".** `panel_coverage` was
  `found` / `missing` / `no_findings`, so a no-call at a requested rsID
  read as "genotyped, no findings".

## Decision

- **Annotator run counts.** `Annotator.run_counts()` returns what the
  source counted during the run (empty by default). ClinVar returns
  `unrecognised_rows` (ADR-0045), `trimmed_matches` (ADR-0047),
  `i_probes_at_plp_loci` and `di_calls_at_plp_rsids`. The pipeline reads
  these in place of the annotator's attributes.
- **Probe classes.** P/LP is ClinVar's aggregate classification
  (`is_pathogenic_aggregate`, ADR-0050). Counted by call shape, whatever
  the vendor:
  - `i_probes_at_plp_loci`: called probes whose ID starts with `i`, at a
    locus where ClinVar has a row with a P/LP aggregate. A no-call is not
    counted.
  - `di_calls_at_plp_rsids`: calls with a `D` or `I` allele whose rsID has
    a ClinVar row with a P/LP aggregate, wherever ClinVar places it.
  - JSON `input.probe_classes`; null for VCF input and when ClinVar did
    not run on the input's build.
- **Array input stats.** `ArrayInputStats` in `parsers/base.py`, the array
  counterpart of `VcfQualityStats`, in JSON `input.array_quality` (null
  for VCF input):
  - `di_calls`: every call with a `D` or `I` allele, counted by the
    pipeline;
  - `unreadable_genotypes`, `not_genotyped`, `wrong_column_rows`: for
    parsers that can tell (a genotype field no rule reads, a reference
    call a converter filled in, a row with the wrong column count). A
    parser keeps them through `GenotypeParser.input_stats()`; the others
    leave them at 0.
- **Disclosure.** The console and the HTML input-quality notice say how
  many i-probes sit at P/LP loci and how many D/I calls there are, with
  the number at P/LP rsIDs.
- **Panel coverage.** Each requested rsID is in exactly one state:
  - `annotated`: at least one shown annotation (after the display filter).
    It wins over a no-call, so a half call with a finding is annotated.
  - `typed_no_annotation`: in the input and called, nothing shown.
  - `no_call`: in the input with either allele missing.
  - `not_on_chip`: not in the input. For VCF input it is shown as "not in
    file": the position may be reference, or have no variant call.

  JSON `panel_coverage` is `{requested, counts, variants}`. `counts` has
  the four states and sums to `requested`; each `variants` entry is
  `{rsid, gene, state, genotype, annotation_count}`, sorted by state, then
  rsID. A VCF half call's genotype is shown as the annotations show it
  (`T/.`). HTML lists every rsID with its state; the terminal report
  prints the four counts and flags the no-call and not-on-chip rsIDs.

### Amendment 2026-10-25 (2.5.2): D/I calls at high-value and P/LP sites are named

The aggregate D/I count stays. Every D/I call at a high-value rsID, or at
an rsID ClinVar classifies Pathogenic or Likely pathogenic, is also named:
its gene, the variant's name where the high-value list gives one
(`6174delT`, `5382insC`, `F508del`; new optional `name` field), its rsID,
and the call as written. It is marked as not checked: neither a finding
nor a normal result. The list is in JSON `report_warnings`
(`di_calls_not_checked`) and `high_value_unchecked`, in the HTML
high-value section, and in the console's high-value block. Carriers and
non-carriers get the same entry, because nothing is interpreted. That is
2.6. A D/I call at any other rsID is counted only.

## Consequences

- Schema 9 gains `input.array_quality` and `input.probe_classes`, and
  `panel_coverage` changes shape (`found`, `missing` and `no_findings` are
  gone). The CHANGELOG lists each.
- Nothing is matched that was not matched before; the counts only say what
  was not assessed.
- Matching i-probes or D/I calls would need a vendor probe map, which
  Allelix does not have; that is a separate decision.

## Amendment 2.6.0 (2026-11-03): D/I calls and ClinVar's run counts

- **D/I calls.** What this ADR left to 2.6 is settled by ADR-0070: a chip
  D/I call is interpreted by no source. It is still counted, and named as
  not checked at a high-value or P/LP rsID, as above.
- **Run counts.** ClinVar's `run_counts()` also returns `ref_mismatches`:
  one-base VCF REFs that are none of ClinVar's one-base REFs at the
  position (ADR-0035, `input.vcf_quality.clinvar_ref_mismatches`).

## Amendment 2.6.1 (2026-11-05): a chip call naming another rsID's allele

A chip call is matched to ClinVar only when every called base is one of
the alleles of its own rsID's SNV rows (REF and ALT one base each) in the
rsID identity table at the call's build and locus. An indel row is no
evidence (its one-base side says nothing about SNV alleles), and rows of
other rsIDs are never pooled in. At a locus that also holds another rsID's G>C row, a
`CC` call for a G/A rsID was taken as that rsID's C and could get its P/LP
finding (seen on 2.6.0). The guard establishes only that the call
disagrees with its rsID's alleles, not why: `CC` there may be `GG` read on
the other strand, but `AC` fits neither reading of G/A. Such a call is marked while its REF is filled
from the identity table (`foreign_alleles`, `_populate_ref_from_identity`)
and ClinVar does not match it. At a locus with an aggregate-P/LP row it is
counted in `input.probe_classes` as `foreign_allele_calls_at_plp_loci`
(JSON schema 13), and the console and HTML say so.

Nothing is complemented. A call with D/I alleles, a VCF record (its REF and
ALT come from the file), and a call whose own rsID has no SNV row in the
identity table at that build and locus (an indel-only rsID, unknown,
merged or alias IDs) are not judged. REF is still filled from the identity
table as before. Other
sources are unchanged: they read the call under its own rsID.

Every call the run could not check is now also named, not only counted,
where its site is high-value or ClinVar classifies it as Pathogenic or
Likely pathogenic: in JSON `high_value_unchecked`, on the console
and in HTML, with the reason. Besides `di_call` at an rsID and
`reference_carried_not_in_file`, the reasons are:

- `di_call` at a vendor probe without an rsID (23andMe `i…`), when an
  aggregate-P/LP row stands at the probe's position (named by probe ID and
  ClinVar's gene);
- `vendor_probe`: a base call at such a probe, the calls
  `i_probes_at_plp_loci` counts (named by probe ID and ClinVar's gene);
- `unplaced`: a row of an rsID and genotype list that placement could not
  place (ADR-0069), with why; P/LP is looked up by rsID on GRCh38;
- `identity_locus_conflict`: a call at a position the identity table does
  not give for its rsID (ADR-0063), with that position;
- `foreign_alleles`: a chip call that names a base its own rsID's SNV
  identity rows do not list, above, when it stands at a P/LP locus or its
  rsID is high-value.

Each entry says what ClinVar did not check and claims nothing about any
other source, which may be switched off or missing. Each site is named
once per reason. The counts (`probe_classes`,
`array_quality`, `placement`, the `identity_locus_conflict` warning) are
unchanged.

## Amendment 2.6.1 (2026-11-06): the call's own rsID's rows always match

The 2.6.1 guard above took a chip call's allele evidence from the identity
table alone and then matched the call to no ClinVar row at all. Two
defects followed: a call of its own rsID's Pathogenic allele that the
pinned identity snapshot did not list (rsO G/A, ClinVar G>T, call GT) lost
its finding; and where the identity table held no SNV row of the rsID
(absent, switched off, indel-only, or REF filled later by gnomAD) another
rsID's row could still match. Now:

- The evidence is the union of the identity table's SNV rows, gnomAD's rows
  of the rsID at the locus and ClinVar's own SNV rows of the rsID there
  (`_others_allowed`). The resolution step keeps what it learns while
  filling REF (`allele_evidence`); when that and ClinVar's rows do not hold
  every called base, an enabled, ready gnomAD is asked for the rsID directly
  (`set_allele_evidence_provider`), so the union does not depend on which
  source filled REF. A gnomAD switched off or missing supplies nothing,
  and gnomAD's positions are GRCh38, so for GRCh37 input it supplies
  nothing either: a call whose rsID only gnomAD knows can match another
  rsID's row on GRCh38 and is `unknown_alleles` on GRCh37. A gnomAD asked
  only for this evidence is recorded as used (ADR-0049).
- The call's own rsID's rows are always matched.
- Another rsID's row is matched only when every called base is in the
  evidence. A call naming a base outside it is `foreign_alleles`; a call
  whose rsID's alleles no source knows matches no other rsID's row
  (`unknown_alleles`). Both are counted at a P/LP locus
  (`foreign_allele_calls_at_plp_loci`, `unknown_allele_calls_at_plp_loci`)
  and named.
- An rsID call with no REF whose rsID neither ClinVar nor the identity
  table places at the locus, at a P/LP position, is counted and named as
  `unknown_alleles` too.

Naming: a D/I call ClinVar does not name by rsID (a vendor
probe's, or one under an rsID with no aggregate-P/LP row: absent, merged,
or held only as Benign or VUS) is named, as before,
by an aggregate-P/LP row at its own position, and now also by an
aggregate-P/LP indel row whose reference span (its first to last REF base)
covers the position or ends up to five bases before it (ClinVar left-aligns; a chip may give a later
position in a repeat). The span is the REF's length, so an insertion
reaches only its anchor and the slack, however long its ALT. Indels with a
REF up to 64 bases are found in a fixed window behind the call; longer
ones are read in one scan per build and looked up by chromosome, so a long
deletion covering the call is found and an unrelated long indel elsewhere
does not widen the search.
Such a call's line says which row named it and how it reaches the call,
naming that row's gene: a Pathogenic or Likely pathogenic variant at the
position, an indel whose REF span covers it, or an indel whose span ends a
stated number of bases (one to five) before it. The line says "covers"
only in the second case. Under an rsID ClinVar holds, the call is listed
under that rsID's own ClinVar gene, and the line adds that ClinVar does not
classify the rsID itself as Pathogenic or Likely pathogenic, so it does not
read as a P/LP rsID's D/I call. A vendor probe's call, or one under
an rsID ClinVar does not hold, is listed under the naming row's gene, and
its line says nothing of how ClinVar classifies the rsID. An
unplaced rsID-list D/I row is named once, as `di_call`. A
foreign-allele name mentions another rsID's row only where one stands at
the position, and the HTML banner says each line states what was not
checked.

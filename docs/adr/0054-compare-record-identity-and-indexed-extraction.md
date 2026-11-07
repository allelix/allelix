# ADR-0054: Comparison identity and indexed extraction completeness

Status: Accepted

## Comparison

Nonempty identifiers retain the named comparison path, including vendor
IDs. Unnamed records use a separate key namespace: normalized chromosome,
position, uppercase input REF and the complete literal uppercase input ALT
field. This path requires the same known effective build in both files.
Unknown/different builds or insufficient input allele identity omit those
observations and report counts per file. No named-to-unnamed matching,
liftover, allele decomposition or equivalence of reordered ALT lists is
introduced.

The parser retains the input ALT field independently of genotype and
quality masking in optional internal model metadata. The default is null,
old positional constructors remain valid, and this metadata is excluded
from model equality and public report serialization. Report schemas and
writers do not change.

Input counts describe parsed named and unnamed observations, including
duplicates. Eligible counts describe unique keys. Existing last-observation
map semantics remain for duplicate keys; contradictory duplicates are not
adjudicated and can depend on order. Distinct keys are order-independent.
Concordance categories and their shared-eligible-key denominator retain
the existing genotype/strand/no-call rules.

## Extraction

After eligible indexed lookup, pending requests are `wanted - found`.
One sequential pass resolves that whole pending set, including unplaced
requests, empty placement lists and mapped fetch misses. Indexed hits,
including no-calls, are excluded from the pending set and are not replaced.
The existing shared parser, resolver, sample and build decisions remain
authoritative; an index changes the lookup path, not those policies.

Actual index contig names are grouped through existing chromosome
normalization and traversed in index order. This handles prefixes in
either direction and M/MT/chrM/chrMT. Seed checks use the same aliases.
Repeated exact intervals share parsed-row results within one extraction
and sample binding. Position filtering still excludes overlapping rows
whose start is elsewhere. Output retains one row per requested identifier.

## Verification

Default offline tests exercise indexed behavior with a deterministic
index double. Tiny real BGZF/tabix cases require the optional `vcf-index`
extra; environments without it skip those cases. The full release runner
runs all tests and rejects skips, and its prepared environment includes
that extra. Ordinary default CI does not install the extra and does not
establish real-index coverage. No mandatory runtime dependency or workflow
change is introduced here.

**Superseded for CI (since 2.5.1; noted in 2.5.2):** CI installs
`.[dev,vcf-index]` and fails the gate on any skipped test, so ordinary CI runs the real BGZF/tabix cases
too. The paragraph above records the decision as it was made.

## Amendment 2026-10-25 (2.5.2): reference calls and one record order for extract

Two gaps in `extract`:
- **Paths disagreed.** The indexed path fetched rows at the placements gnomAD
  maps a requested rsID to; the sequential path kept the first record in the
  file. With records at two positions, the paths could show different calls.
- **Confident reference calls were missing.** A `0/0` record was outside the
  extract stream and printed "not in file", while panel coverage showed the
  same record as typed.

Both paths now share one stream rule and one order (supervisor rulings):

- **The stream** includes, besides the calls `analyze` reads:
  - explicit reference records (`0/0`, `0|0`, haploid `0`, ALT `.`), as
    REF/REF;
  - gVCF reference blocks covering a placement, END inclusive, as that
    site's own REF/REF, never the block's anchor base;
  - filtered, low-quality and zero-depth records, as no-calls.

  Such a record counts for a requested rsID when its ID is that rsID, or when
  it starts at (explicit) or covers (block) one of the rsID's placements and
  carries no other rsID. Only `extract` reads these records; analysis and
  export streams are unchanged.
- **Placements:** on GRCh38, gnomAD's (the rsID identity table's when gnomAD
  is not ready). On GRCh37, the identity table's GRCh37 placements, never
  GRCh38 coordinates.
- **The order** among an rsID's records:
  1. A record at one of its placements beats one elsewhere. A no-call at a
     placement still beats a usable call elsewhere.
  2. Among records at a placement, a usable call beats a no-call or a
     rejected record.
  3. Then file order: the header's contig order, then POS, then line order.
     A file without header contigs (arrays) uses line order alone.
- **The indexed path** fetches every mapped placement, covering blocks
  included, and applies the same order. It never takes a placement by
  gnomAD's row order. A request with no record at any placement goes to the
  sequential pass, which applies the same order.
- **Early stop:** the sequential pass stops once every request is settled,
  except when contig lengths decided the build. An rsID without placements
  is settled by its first call. One with placements is settled only by a
  usable call at its best-ranked placement (the header's contig order, then
  POS), since no later record can rank before it in any record order. A
  usable call at another placement does not settle it, and in a file
  without header contigs nothing does, so the pass reads to the end.
  When MT is excluded (chrM has the hg19 length), MT records and placements
  are left out of this decision as they are of the choice: an excluded call
  settles nothing, and an rsID placed only on MT is not settled early.

Build, sample and quality rules are unchanged; they decide which records are
usable.

## Amendment 2026-10-27 (2.5.2): identity placements beside gnomAD, and unplaced calls

Two more gaps in `extract` on a VCF without rsIDs:
- **gnomAD alone placed GRCh38 requests.** With gnomAD ready, an rsID that
  gnomAD does not hold had no placement, even when the rsID identity table
  places it, so a call there printed "not in file".
- **Unnamed calls were dropped.** A no-call (`./.`, haploid `.`) or a half
  call carrying only REF (`0/.`) is a call the rsID recovery cannot name, and
  it was never matched to the placement it sits at, even one gnomAD maps.

Now:
- **Placements** on GRCh38 are gnomAD's, and the identity table's GRCh38
  placements for every rsID gnomAD does not place (all of them when gnomAD is
  not ready). An rsID gnomAD places keeps gnomAD's placements alone. GRCh37 is
  unchanged. The indexed path fetches the same placements.
- **A call that carries no ALT** (a no-call, or a half call whose one allele
  is REF) and so has no rsID after recovery counts for every requested rsID
  placed at its position, on both paths. A call carrying an ALT counts only
  when the recovery names it as the rsID: left unnamed, it is another
  variant. The order above then chooses among the candidates.

## Amendment 2026-10-29 (2.6.0): padded calls at a placement

A padded record stands at its POS and at the site each ALT trims to
(ADR-0061). In the choice above, a call at a placement ranks by that
placement's position, not its record's POS. So a padded call and the same
call written trimmed tie, and the one earlier in the file wins. The early
stop uses the same rule. The call shown keeps its record's own columns.

## Amendment 2026-10-30 (2.6.0): comparison records keep their alleles; chip calls take an explaining record

Supervisor ruling, 2026-10-30.

**Named records.** A named record with a literal REF and ALT is keyed by:
- its rsID;
- its trimmed REF;
- its trimmed ALTs, in their literal order.

Trimming drops the bases REF and every ALT share at the end, then at the
start, while each keeps one base (as in ADR-0047 and ADR-0061; no
left-alignment). The position is not part of the named key, so a named
record still matches across builds. As a result:
- the records of a split multi-allelic site each have their own partner,
  whatever the order of the records in either file;
- a record without a partner counts as file-only;
- a split record and a merged record (`A` and `T` against `A,T`) are
  different records.

**Unnamed records** are keyed by their trimmed position, REF and ALTs, so a
padded and a trimmed spelling of the same record match. Calls are compared
in the trimmed form of their record's alleles.

**Chip calls.** A chip call has no REF or ALT, so it is keyed by its rsID.
- **Against another chip,** it matches the other file's chip call with
  that rsID, as before.
- **Against a VCF,** it is matched with the one record under its rsID
  whose declared alleles (trimmed REF and ALTs) hold the call, directly or
  complemented. The partner is chosen by declared alleles, never by whether
  the genotypes agree.
- **When the VCF has records under its rsID but none or more than one of
  them does,** the call is not compared, is not counted as file-only, and
  `compare` prints the count and the reason. A chip call whose rsID the VCF
  does not hold at all is counted as file-only.

This replaces, for distinct records, the rule in Comparison that keeps the
last observation per rsID. Duplicates of the exact same key still keep the
last observation.

## Amendment 2.6.0 (2026-11-02): excluded calls are not compared

A call that the file's build decision leaves out (ADR-0044: MT under an
hg19-length chrM; a UCSC-named mitochondrion in a GRCh37 VCF) is in no
comparison key, coverage count, concordance or chromosome breakdown. A
declared rCRS mitochondrion, and MT on GRCh38, are compared as before.

## Amendment 2.6.1 (2026-11-05): a gVCF's reference-block ALT is not part of a record

A gVCF lists `<NON_REF>` (or `<*>`) as an extra ALT on every variant
record. It names no allele the record carries, so `compare` leaves it out
of a record's alleles before trimming (`_trimmed` in `allelix/compare.py`):
`rs1 A G,<NON_REF>` in a gVCF and `rs1 A G` in a plain VCF are one record
and are compared, and so are the same unnamed records in two files of one
build. In 2.6.0 the named records got different keys and were each counted
as in one file only; 2.5.2, which keyed a named record by its rsID alone,
paired them. Every other ALT stays literal, as above.

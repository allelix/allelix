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

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

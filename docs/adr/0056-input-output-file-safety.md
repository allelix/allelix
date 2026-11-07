# ADR-0056: Protect input files and own staged output resources

## Status

Accepted

## Context

A report destination equal to its input replaces that input. Fixed report and
PLINK temp names can overwrite or delete unrelated files, including an input
named like a former temp. Direct VCF output can leave a partial destination on
a streaming failure. CLI safety checks must precede cache and output mutation.

## Decision

Validate every requested report/VCF destination and all three actual PLINK
paths before cache-directory setup, refresh or output allocation. Use resolved
identity for lexical/new paths and samefile for existing hard links; real
stat/resolve errors propagate. Refuse input aliases of former fixed report and
PLINK temps as well, without opening or deleting those old temp paths. Direct
report/VCF APIs validate their known input path. The iterator-only PLINK API
cannot infer an input origin; its CLI supplies that boundary check. The CLI's
existing VCF no-overwrite policy stays intact.

Allocate unique same-directory temps with mkstemp, prefix `.` plus destination
name plus `.`, suffix `.tmp`. Consume the returned descriptors, finish writing
and close successfully before os.replace. Stage all three PLINK files before
first publication. Keep new staging files private and apply only existing
ordinary permission bits (0777) to complete content before publication, never
copying setuid/setgid/sticky bits. Supported replacement retains those bits;
permission failures preserve the previous destination.

Clean only allocated descriptors and still-owned unpublished temp paths on
ordinary failures and interrupts. Successful publication ends ownership of the
old staging pathname; cleanup must leave any later unrelated file there alone.
On Windows, a still-owned temp made readonly for publication may need its own
readonly attribute cleared after replacement fails. Never change the final
file's permissions to make cleanup succeed, and preserve the primary error.

## Consequences

Prepublication failure preserves all prior outputs and creates no partial
newly named output. PLINK publication is atomic per file: later rename failure
can leave complete already-published files alongside unchanged prior files.
Cleanup never deletes those published outputs as purported rollback.

JSON/HTML serialization, streaming VCF UTF-8/LF, PLINK content/order/counts and
schema versions are unchanged. Filename/build/allele changes and output/output
collisions remain separate work. There is no claim of preserving ownership,
ACLs, xattrs or link topology, crash durability, a three-file transaction, or
protection from hostile concurrent path swaps. An uncatchable kill may leave
an owned temp. Cross-platform behavior requires platform-specific verification.

## Amendment 2026-10-25 (2.5.2): outputs that alias each other; new-file mode

- **Two destinations of one run must be two files.** `analyze --output` and
  `--vcf-out`, and the three PLINK outputs of a prefix, are compared before
  any work (`same_destination`): resolved paths, which follow symlinks even
  when they dangle, and `samefile` for existing hardlinks. A collision stops
  the run, and nothing is written. The direct PLINK writer checks too, before
  it reads its input.
- **New-file mode.** A published output that did not exist before gets
  `0666 & ~umask`, the mode the user would get from any new file, instead of
  the staging temp's `0600`. An existing destination keeps its mode, as
  before.

## Amendment 2.6.0 (2026-11-02): cache writers own their staging

The caches did not follow this decision. Every writer of a cache staged
under a name made from the target alone: `<cache>.tmp` (GWAS Catalog,
ClinPGx, and every prebuilt install: gnomAD, AlphaMissense, CADD, SNPedia,
the rsID identity table), `<file>.part` (every download), `<cache>.partial`
(`load_clinvar_tsv` called directly), and fixed download names in the data
directory (`gnomad.sqlite.gz` and the other prebuilt archives,
`gwas_catalog_associations.zip` and its extracted table,
`clinicalAnnotations.zip`). The offline identity builder wrote its target in
place, and the offline gnomAD, AlphaMissense and CADD builders
(`scripts/build_*_cache.py`) staged at `<output>.tmp`. Two refreshes of one cache (two `db update`s, or an `analyze` that
refreshes a stale source beside another) could then publish or delete each
other's bytes: in the outside review of RC3, a successful GWAS refresh
installed the other refresh's empty file, losing the previous cache.

Now each operation owns its staging, as outputs do above:
- a database is staged in a unique hidden file beside it
  (`databases/staging.py`, `owned_staging`), moved onto it only when
  complete, and removed on any failure or interruption;
- a setup downloads into a directory of its own in the data directory
  (`owned_directory`) and verifies, extracts and loads from there, so it
  reads the bytes it fetched; the directory is removed when the setup ends.
  The GWAS Catalog table a refresh loaded is then moved into the data
  directory for a later re-ingest (ADR-0024), and ClinPGx's archive is no
  longer left there;
- the offline builders (the identity table's, and gnomAD's, AlphaMissense's
  and CADD's) stage their output the same way.

No run removes a staging file it did not create, so a file left by a killed
process stays until the user deletes it (`db clean` removes the data
directory's contents). A cache staging file is created with the mode a new
file gets under the user's umask, as before, not private first: caches hold
no personal data, and the mode is that of every cache installed before.
Production ClinVar setup already staged everything in its own temporary
directory (ADR-0053) and is unchanged. Two refreshes that both succeed
publish their own complete caches in turn; which one is installed last is
not ordered.

A cache and everything that names its freshness belong to one refresh.
Unique staging alone did not give that where a refresh publishes two things,
or updates a cache it did not just build:
- GWAS Catalog: the database and the table kept for re-ingest were moved
  into place one after the other, so a refresh finishing between the two
  left one refresh's database beside another's table, and a later re-ingest
  built that table's rows under the database's freshness signal. A GWAS
  cache now records the SHA-256 of the table it was built from (an optional
  `content_hash` field, as ClinPGx's, ADR-0068). A re-ingest takes its own
  link or copy of the kept table before it reads the cache's signal and
  hash, rebuilds only when the table matches a recorded hash, and keeps the
  signal and hash read with it; otherwise it rebuilds nothing and
  `db update` downloads again. A cache with no recorded hash re-ingests its
  own table only, because the table is published after the database
  (ADR-0024, amendment 2026-11-03).
- ClinPGx: the unchanged-content shortcut stored its new signal on whatever
  cache the path named after its hash check. It now stores it only while
  that cache still holds the content it checked, in the same statement; a
  cache another refresh installed meanwhile keeps its own signal (ADR-0068,
  amendment 2026-11-03).

## Amendment 2.6.1 (2026-11-05): the default PLINK prefix never replaces a file

Without `--output-prefix`, `export plink` derives its prefix from the
input's stem, so its three destinations are files the user did not name.
When any of them already exists, a link (even a dangling one) included,
the command stops before any work, names it and writes nothing
(`export_plink_cmd` in `allelix/cli/utility.py`). In 2.6.0 an unrelated
`genome.bed` beside `genome.vcf` was replaced. An explicit
`--output-prefix` keeps the rule above: its files are replaced, atomically
per file.

## Amendment 2.6.1 (2026-11-05): `db clean` removes only allelix's names

`db clean` removed every entry in the data directory except `config.toml`,
behind a guard that only asked for one name starting with a source prefix.
CADD full-mode files, which the user downloads and places there and
`db update` cannot fetch, were deleted with the caches, `--keep-cadd`
included, and so was any other file kept there.

It now removes only names allelix creates in the data directory: each
cache and its SQLite side files (`-wal`, `-shm`, `-journal`); GWAS's kept
table and the ClinVar build log; staging entries, both this ADR's
per-operation names (`.<name>.<16 hex digits><suffix>` and the
`.<source>-<8 characters>` directories) and the shared `<name>.tmp`,
`.part` and `.partial` names of earlier releases; and the downloads
earlier releases left in the data directory (the prebuilt `.gz` files, the
GWAS and ClinPGx zips, ClinVar's VCF). A name given to a file counts only
for a file or a symbolic link, and a link is removed, never followed; a
directory is removed only when it has the staging-directory shape. Every
other entry is kept, and listed apart from `config.toml` and what
`--keep-cadd` keeps.
`--keep-cadd` keeps the CADD cache's names. The guard is gone, since
nothing allelix did not create is removed; `--force` is accepted and has
no effect. A staging entry left by a killed process is therefore removed by
`db clean` without `--force`.

## Amendment 2.6.1 (2026-11-05): signals end a command through its cleanup

A `db update` stopped by SIGTERM or SIGHUP left its staging entries (which
can be several GB), since Python's default action for those signals ends
the process without running `finally` blocks. SIGTERM and SIGHUP now end
every command as a normal exit, with status 128 + the signal number, so each
owner removes its own entries. A signal the caller ignores (SIGHUP under
`nohup`) or handles is left alone; during the cleanup only the handlers
installed here are set to ignore a second signal, and they are restored
afterwards.

A process killed outright (SIGKILL, a crash) still removes nothing, and no
other run removes what it left: it cannot be told from a live operation's.
`db clean` removes it, since it removes only names allelix creates (above).
`db clean` cannot tell that an update is running, so it is not to be run
during one. Telling a dead operation's entries from a live one's (a lock the
live operation holds, and a sweep at the start of the next) is left for a
later release.

## Amendment 2.6.1 (2026-11-06): staging-shaped directories, legacy names, Windows

`db clean` removes a directory of the staging shape (`.<source>-<8
characters>`) only when every entry in it is an exact name allelix
downloads there (the ClinVar TSV pair, the prebuilt `.sqlite.gz` files, the
GWAS zip and its table, the ClinPGx zip) or a partial file of one. The
shape alone is not enough: a user's dated backup such as
`.clinvar-20261101/` has it too. GWAS's extracted table is named by EBI's
zip, so it cannot be proven allelix's either. A staging directory holding
anything else, a cache built there by a killed update included, is kept and
listed as possibly an interrupted update's leftovers, to delete by hand.

v0.7.x's single `clinvar.sqlite`, before the per-build files, is one of
allelix's names, with its side files. No other name of a release before 2.0
is written into the data directory: the CHANGELOG's other file names are
read inside a download, fetched without being kept, or were never written.

On Windows a forced end (TerminateProcess, as `taskkill /F` and Task
Manager use, and SIGTERM there), closing the console and Ctrl-Break do not
run the SIGTERM/SIGHUP cleanup of the amendment above. A normal exit, an
error and Ctrl-C do. After a forced end, as after SIGKILL, `db clean`
removes what is left.

## Amendment 2.6.1 RC3 (2026-11-06): what a killed update's staging directory holds

`db clean` removes a staging directory when it holds only names allelix
writes there: its downloads and their partial files (amendment above), a
cache's partial file while it is built there (`.<cache>.<16
hex>.partial`, with its side files), and an empty `extracted/`, which a
GWAS refresh leaves once its table is moved out. So an update killed while
it downloads, loads the GWAS table or builds the first ClinVar build leaves
only allelix's names, and `db clean` removes the directory. A staging
directory holding anything else is kept and listed, to delete by hand: a
finished candidate under its own name (the first build of a two-build
ClinVar update, SNPedia's staged pages while they parse), a table still in
`extracted/`, a user's backup. Marking a finished candidate
as allelix's own, so that it too can be removed, is left for a later
release.

# Contributing to Allelix

Allelix is an open-source genotype analysis toolkit licensed under AGPL-3.0-or-later.
Contributions are welcome.

## Development Setup

```bash
git clone https://github.com/allelix/allelix.git
cd allelix
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,vcf-index]"
git config core.hooksPath .githooks
pre-commit install --hook-type pre-commit
```

The `vcf-index` extra installs pysam, which the indexed-VCF tests need; it
is the same install as CI's, so a local run is the CI run.

System tools outside pip: ``bcftools`` (fast suite; its tests skip
locally when it is missing, and CI installs it) and ``plink2`` (one
integration test in the slow suite; skips when missing):

- ``bcftools`` — annotated-VCF round-trip validation for ``--vcf-out``
  (ADR-0036). ``apt install bcftools`` on Debian /
  Ubuntu, ``brew install bcftools`` on macOS. Tests skip gracefully
  when absent locally; CI installs it unconditionally.
- ``plink2`` — PLINK1 binary export round-trip integration test
  (marked ``@pytest.mark.integration``, excluded from fast CI, runs
  in the manually dispatched slow suite).

Run the test suite:

```bash
pytest
```

Tests marked `@slow` exercise real-data invariants against fixtures
that are gitignored due to size. These fixtures are fetched on first
use and cached; no manual download step is required for the GWAS
catalog.

For the full real-data release-validation battery — all parser formats,
VCF/gVCF, PLINK round-trip, edge cases, and upgrade-path checks — see
[`test_data/FULL_TEST_PROTOCOL.md`](test_data/FULL_TEST_PROTOCOL.md).
This is the project's primary answer to "do you validate against real
data" and is what every release run executes.

Two fixture sources, cached under `test_data/`:

- **GWAS Catalog** (`test_data/gwas_catalog.zip`, ~65 MB) — auto-fetched
  by the `TestRealDataGwasSanity` fixture in
  `tests/test_end_to_end.py` on first run. Cached forever after.
- **Real genotype files** (`test_data/real/`, `test_data/transcoded/`) —
  used by cross-parser identity checks. Fetched via:

      scripts/fetch_testdata.sh

  One-time; subsequent runs detect the existing data dir.

### Silent skips are forbidden on ship gates

A test must not `pytest.skip()` because an optional fixture is missing.
Either:

1. **Auto-fetch the fixture** when absent (preferred — see
   `TestRealDataGwasSanity` for the pattern).
2. **Mark `@pytest.mark.integration`** and document the external
   precondition (e.g., `plink2` binary, live `~/.local/share/allelix/`
   ClinVar cache). The marker makes the dependency explicit; the
   ship-gate procedure must ensure those preconditions are met before
   tagging.

`pytest.skip` for a committed mock fixture that "should always be
present" is also forbidden — convert to `assert` so a real regression
is caught instead of hidden.

### Run the full suite locally

    pytest                      # runs everything: fast + slow
                                # (auto-fetches the GWAS zip if missing)

`pytest` runs everything locally, fast and slow; it fetches the GWAS zip
if missing. Hosted CI runs the fast tier
(`-m "not slow and not integration"`) on pull requests targeting `main`,
manual dispatch and publication refs. Development pushes and pull requests
targeting `dev` do not start hosted CI. The
slow and integration workflow is manually dispatched for release validation;
it has no scheduled run.

### Development checks and release gates

During development, run the affected regression tests, lint, formatting
and type checks locally. Regression tests must demonstrate the defect
before the fix and pass afterwards. Use an adequately provisioned local
machine for the full suite.

Before merging to `main` or publishing, require a successful hosted CI run
for the finished candidate, through the `main` pull request or manual
dispatch, and record the resolved commit SHA. All matrix jobs, warning checks and
coverage requirements must pass. Run the full release-validation
battery, including slow and integration tests and their required tools,
on that same candidate. Inspect and resolve failures and missing
prerequisites before publishing.

A changed candidate needs validation of the new commit. Batch related
development work before dispatching another hosted matrix. Each hosted
workflow cancels a superseded run of the same workflow and ref; a
canceled run never supplies release evidence. Publication and version
checks remain enabled.

Lint and format:

```bash
ruff check .
ruff format .
```

## Coding Standards

- Python 3.11+. Use `from __future__ import annotations`.
- Type hints on all signatures. No bare `Any`.
- Google-style docstrings on public classes and functions.
- Ruff enforces linting and formatting. Zero warnings before commit.
- Every file starts with the AGPL-3.0-or-later license header and copyright.

## How to Add a Parser

Parsers live in `allelix/parsers/`. Each parser is a single file that implements
the `GenotypeParser` abstract base class. The parser's job is to read a vendor's
genotype file format and yield normalized `Variant` objects.

### Step 1: Create the parser file

Create `allelix/parsers/vendorname.py`:

```python
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parser for VendorName genotype files."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from allelix.models import Variant
from allelix.parsers._helpers import no_call_allele, open_genotype_text
from allelix.parsers.base import GenotypeMetadata, GenotypeParser

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)


class VendorNameParser(GenotypeParser):
    name: ClassVar[str] = "vendorname"
    display_name: ClassVar[str] = "VendorName"
    file_extensions: ClassVar[list[str]] = [".txt"]
    url: ClassVar[str] = "https://vendorname.com"

    def can_parse(self, file_path: Path) -> bool:
        """Check for the vendor's signature in the first few lines."""
        with open_genotype_text(file_path) as fh:
            for line in fh:
                if line.startswith("# VendorName"):
                    return True
                if not line.startswith("#"):
                    break
        return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Yield Variants from the file. Stream, don't load into memory."""
        with open_genotype_text(file_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                # Skip the header row
                if line.startswith("rsid"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    logger.warning("Skipping malformed line: %s", line.strip())
                    continue
                yield Variant(
                    rsid=parts[0],
                    chromosome=parts[1],
                    position=int(parts[2]),
                    allele1=no_call_allele(parts[3]),
                    allele2=no_call_allele(parts[4]),
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Extract metadata from comment headers."""
        sample_id = ""
        with open_genotype_text(file_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if "Sample ID" in line:
                    sample_id = line.split("\t")[-1].strip()
        return GenotypeMetadata(
            format=self.name,
            sample_id=sample_id,
            build="GRCh37",
        )
```

Key rules:

- `can_parse()` must be fast. Only look at header/comment lines.
- `parse()` yields `Variant` objects one at a time (streaming).
- Malformed lines log a warning and skip. Never crash the whole parse.
- Open files with `open_genotype_text`. Test CRLF, LF and CR-only endings,
  plain and compressed.
- No-calls use `"-"` as the allele value (matches `allelix.models.NO_CALL_MARKER`).
  Pass each allele through `no_call_allele()` so a vendor's `0` or `N`
  becomes a no-call.

### Step 2: Register the parser

Add your parser to `allelix/parsers/__init__.py`:

```python
from allelix.parsers.vendorname import VendorNameParser

PARSERS: list[GenotypeParser] = [
    # ... existing parsers ...
    VendorNameParser(),
]
```

Order matters: auto-detection tries each parser's `can_parse()` in order.
Put a more specific parser before a broader one that would also accept its
files (for example Genes for Good before 23andMe).

### Step 3: Add a test fixture

Create `tests/fixtures/mock_vendorname.txt` with synthetic data. Include:

- Comment lines matching the vendor's format
- At least one known rsID with a specific genotype (for annotation tests)
- A no-call line
- An edge case (blank line, extra whitespace, etc.)

**The synthetic fixture IS the format spec.** It must exercise every
documented field-shape variant — concatenated vs. separate alleles,
no-calls (full and partial where applicable), haploid calls on MT/Y if
the format permits them, comment-line handling, header-detection
signals, chromosome normalization cases. Test coverage of this fixture
is what gates the PR; if a future ambiguity is discovered against a
real file, that's a bug, file an issue and extend the synthetic fixture
to cover it.

**Real-data validation is OPPORTUNISTIC, not a merge gate.** Public
CC0 files for the new format may not exist in our corpus (openSNP /
Personal Genome Project don't host every vendor). Don't block your PR
waiting for one. The pattern when a public file later surfaces:

1. Add the public file under `test_data/real/<vendor>/` (gitignored —
   joins the GitHub release asset on the next bump).
2. Add a cross-parser identity test if biology overlaps an existing
   donor (see the user1190 transcoded set for the pattern).
3. If the file exposes a real-vs-synthetic ambiguity, extend the
   synthetic fixture FIRST (the spec gets richer), then ship the
   real-file test as a regression pin.

Annotators follow the same rule — synthetic mock loader fixtures cover
the parse / interpreter contract; real-data validation runs against
the pinned databases via `test_data/databases/` and the `@slow`
integration set documented above.

All `tests/fixtures/` files are synthetic (produced by mock data
generators or hand-written). Real-data integration tests use CC0
public-domain openSNP genotype files fetched via
`scripts/fetch_testdata.sh`.

### Step 4: Write tests

Create `tests/parsers/test_vendorname.py`:

```python
from allelix.parsers.vendorname import VendorNameParser


class TestCanParse:
    def test_recognizes_vendor_format(self, tmp_path):
        f = tmp_path / "sample.txt"
        f.write_text("# VendorName\nrsid\tchr\tpos\ta1\ta2\nrs1\t1\t100\tA\tG\n")
        assert VendorNameParser().can_parse(f)

    def test_rejects_other_format(self, tmp_path):
        f = tmp_path / "other.txt"
        f.write_text("# OtherVendor\ndata\n")
        assert not VendorNameParser().can_parse(f)


class TestParse:
    def test_yields_variants(self, tmp_path):
        f = tmp_path / "sample.txt"
        f.write_text("# VendorName\nrsid\tchr\tpos\ta1\ta2\nrs1\t1\t100\tA\tG\n")
        variants = list(VendorNameParser().parse(f))
        assert len(variants) == 1
        assert variants[0].rsid == "rs1"

    def test_handles_no_call(self, tmp_path):
        f = tmp_path / "sample.txt"
        f.write_text("# VendorName\nrsid\tchr\tpos\ta1\ta2\nrs1\t1\t100\t-\t-\n")
        variants = list(VendorNameParser().parse(f))
        assert variants[0].is_no_call
```

### Step 5: Run tests

```bash
pytest tests/parsers/test_vendorname.py -v --no-cov
ruff check allelix/parsers/vendorname.py tests/parsers/test_vendorname.py
```

## How to Add an Annotator

Annotators live in `allelix/annotators/`. Each annotator queries a reference
database and returns `Annotation` objects for variants the user carries.

### Step 1: Create the annotator file

Create `allelix/annotators/mydb.py`:

```python
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Annotator for MyDB reference database."""

from __future__ import annotations

import sqlite3
from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.databases.cache_state import CacheStatus, cache_state, connect_readonly
from allelix.models import Annotation

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Variant


class MyDBAnnotator(Annotator):
    name: ClassVar[str] = "mydb"
    display_name: ClassVar[str] = "MyDB"
    attribution: ClassVar[str] = "MyDB"
    requires_download: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
        spdx="CC-BY-4.0",
        license_url="https://example.com/mydb/license",
        attribution_text="MyDB variant data.",
        source_url="https://example.com/mydb",
        commercial_ok=True,
    )

    def __init__(self, data_dir: Path) -> None:
        super().__init__(data_dir)
        self._db_path = data_dir / "mydb.sqlite"
        self._conn: sqlite3.Connection | None = None

    def setup(self) -> None:
        """Download and ingest the database. Idempotent."""
        # Download from source, parse into SQLite cache
        ...

    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return annotations for variants the user carries.

        MUST check both rsID AND genotype. Presence in the database
        is not enough -- verify the user carries the flagged allele.
        """
        if variant.is_no_call:
            return []
        conn = self._connection()
        rows = conn.execute(
            "SELECT alt, significance, condition, gene FROM mydb_variants WHERE rsid = ?",
            (variant.rsid,),
        ).fetchall()

        annotations: list[Annotation] = []
        carrier_alleles = {variant.allele1, variant.allele2}
        for alt, significance, condition, gene in rows:
            if alt not in carrier_alleles:
                continue
            annotations.append(
                Annotation(
                    source=self.name,
                    rsid=variant.rsid,
                    significance=f"mydb_{significance}",
                    category="clinical",
                    magnitude=5.0,
                    description=f"MyDB: {significance}",
                    attribution=self.attribution,
                    genotype_match=f"{variant.allele1}{variant.allele2}",
                    condition=condition or "",
                    gene=gene or "",
                )
            )
        return annotations

    def is_ready(self) -> bool:
        return self.cache_status().ready

    def cache_status(self) -> CacheStatus:
        # Read-only check: the file, its main table with a row, and a
        # database_versions row named "mydb" (ADR-0048).
        return cache_state(self._db_path, "mydb", tables=("mydb_variants",))

    def version(self) -> str | None: ...

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def fetch_remote_signal(self) -> str | None:
        return None

    def cached_remote_signal(self) -> str | None:
        return None

    def record_count(self) -> int | None:
        return None

    def _connection(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = connect_readonly(self._db_path)
        return self._conn
```

Key rules:

- Every annotation must be **source-attributed**. Set `attribution` to your
  database name. Never omit it.
- **Check genotype, not just rsID.** The user must carry the flagged allele.
  A variant existing in the database means nothing if the user has the
  reference (normal) allele.
- **Readiness comes from `cache_status()`, never from the file existing**
  (ADR-0048). `is_ready()` returns `self.cache_status().ready`. A file
  that exists can still be missing a table, empty, or built by older code.
- **Analysis never writes a cache.** Open it with `connect_readonly()`.
  `setup()` downloads and builds it during `db update`; work that brings an
  older cache up to date from local data belongs in `prepare_cache()`,
  which also runs only in the update step.
- Implement `close()` to release SQLite connections.
- No-calls return empty (no annotation possible without a genotype).

### Step 2: Register the annotator

Add to `allelix/annotators/__init__.py` in `get_annotators()`:

```python
from allelix.annotators.mydb import MyDBAnnotator

def get_annotators(data_dir, ...):
    # ... existing annotators ...
    mydb = MyDBAnnotator(data_dir)
    return [clinvar, pharmgkb, gwas, snpedia, gnomad, alphamissense, cadd, rsid_identity, mydb]
```

The real list has eight sources: ClinVar, ClinPGx (`pharmgkb`), GWAS
Catalog, SNPedia, gnomAD, AlphaMissense, CADD and the rsID identity table.
Also add your class to `_ANNOTATOR_CLASS_LIST` in the same file, which
builds the `_ANNOTATOR_CLASSES` name registry.

### Step 3: Write tests

Create `tests/annotators/test_mydb.py` with a fixture that builds a small
SQLite database in a `tmp_path`. Test:

- Carrier of the flagged allele triggers annotation
- Homozygous reference does not trigger
- No-call does not trigger
- Unknown rsID returns empty
- `attribution` field is set correctly on all results
- `close()` releases the connection

### Step 4: Run tests

```bash
pytest tests/annotators/test_mydb.py -v --no-cov
ruff check allelix/annotators/mydb.py tests/annotators/test_mydb.py
```

## Architecture Notes

- Parsers are stateless. Annotators hold database connections.
- All annotators run on every variant (unlike parsers, which are exclusive).
- The `Annotation.significance` field is always source-prefixed
  (`clinvar_pathogenic`, not `pathogenic`).
- Reports never assert significance directly. They attribute: "ClinVar
  classifies this as pathogenic", not "this is pathogenic."
- The resolved user data directory (default
  `~/.local/share/allelix/`) holds the local database caches and
  `config.toml`. It is populated by `allelix db update`; tests use
  isolated temporary data directories.

## Hooks and CI

Two hooks run locally:

- **pre-commit** (managed by pre-commit framework): `ruff check` + `ruff format --check`
- **pre-push** (raw hook in `.githooks/`): blocks tag pushes where the tag doesn't match `pyproject.toml`

CI (`.github/workflows/ci.yml`) fires on:

- Manual dispatch via the Actions UI (`workflow_dispatch`) for a finished candidate
- Pull requests targeting `main` (before merge)
- Push to `main`
- Push of a `v*` tag

Routine development pushes and pull requests targeting `dev` use local
checks and review. They do not trigger a hosted matrix.

CI runs ruff (lint + format), mypy, and the fast pytest suite (synthetic
fixtures only) on a Python 3.11 / 3.12 / 3.13 matrix. A skipped test fails
the run, so every fast test runs on every Python. `@slow` and
`@integration`-marked tests are deselected in CI
(`-m "not slow and not integration"`) because the runner doesn't
have the ~1.5 GB real-genotype fixture set or the ~19 GB default annotator
database cache — fetching either would push CI runtime past 45 minutes
per development update and consume GitHub Actions minutes already covered
by deliberate local release validation.

### Where the slow / real-data battery actually runs

The complete battery runs locally at ship time against the exact RC
commit, following `test_data/FULL_TEST_PROTOCOL.md`. It includes the
full pytest suite, fresh caches, both genome builds, all parser formats,
and the GIAB / 1000 Genomes corpus. A release does not ship without it.

`.github/workflows/slow.yml` runs only by manual `workflow_dispatch`. It downloads a current ClinVar cache, executes the
CLNSIG drift sentinel by exact node ID, emits JUnit XML, and fails unless
that named test is present and passed; it then runs the slow + integration
marker set. Every RC requires a manual dispatch against the same commit
the local battery verified. An absent or skipped drift testcase is a
release failure, not a green result.

### Why the complete real-data battery is not in per-push CI

- The large public-genome corpus and all production-size caches exceed
  a practical per-push Actions budget.
- Manually dispatched `slow.yml` checks live-source drift and integration
  behavior for the release candidate.
- The local RC protocol supplies the larger, build-paired end-to-end
  signal and runs against the exact release commit before tagging.

To fire the release slow gate, run
`gh workflow run slow.yml --ref <release-branch>` and confirm the
resulting run's `headSha` equals the RC commit. The ordinary fast gate
can still be dispatched with `gh workflow run ci.yml --ref <branch>`.

### Release-content discipline: no drift-vulnerable counts without a snapshot caveat

Raw annotation counts, timing measurements, and `N/M found` panel
ratios drift across ClinVar / GWAS Catalog / ClinPGx refreshes.
A figure measured on ship day will not match a figure measured
the next week — even with byte-identical code — because the
reference databases change underneath. A user (or reviewer)
reading a release note and re-running the analysis on current
DBs will see the numbers diverge and assume a regression.

**Hard rule: no raw drift-vulnerable counts in any of the
following surfaces unless wrapped with an explicit
snapshot caveat:**

- `CHANGELOG.md` entries (the `[Unreleased]` block and every
  shipped release section)
- Commit messages and annotated tag bodies
- GitHub Release notes on the public repo

Patterns that count as drift-vulnerable include (non-exhaustive):
`N/M found`, `M,NNN variants`, `~NNNs → ~NNNs`, `N,NNN annotations`,
`went from X to Y` measurements anchored to a specific data set.

**Two acceptable shapes when a measurement is illustrative:**

1. **Wrap with snapshot tag** — *"On the GIAB GRCh38 benchmark
   (DB snapshot 2026-06-20): X annotations → Y after the fix.
   Counts drift with each ClinVar / GWAS refresh; the §19
   ground-truth harness is the canonical regression signal."*
2. **Quote the harness output** — *"`✓ giab_grch38_benchmark:
   N annotations, M unique keys, all invariants hold`"* — the
   harness floor invariants in `test_data/HG002_GROUND_TRUTH.yaml`
   are drift-tolerant by design.

The deterministic signals that don't drift and can be quoted
freely:

- **Cross-parser identity**: same individual transcoded to
  N vendor formats produces N-way identical annotation sets
  (verified by symmetric set diff of `(source, rsid, sig)`
  tuples — see `test_data/transcoded/`).
- **Same-sample superset**: gVCF analyze output is a strict
  superset of the GIAB benchmark output at default magnitude
  filter — `0 missing` is the pass condition.
- **Wrong-allele safety**: a non-zero allele-direct enrichment
  denominator on both builds, 0 non-direct (complement-resolved)
  scores, and every GRCh37 enrichment value traceable through
  the rsID identity table.
- **Harness invariants**: vocabulary union holds, floor counts
  cleared, spot-checks on published HG002 carriers (rs1801133 /
  rs7412 / rs2010963 / rs6025) pass.

This applies to `CHANGELOG.md`, commit messages and release notes
alike: check release-content drafts for
drift-vulnerable numeric tokens in unsnapshotted release sections
before publishing.

### Output-contract schema bumps

Allelix ships two versioned output contracts and each one has its own
independent schema version. A change to one contract does not
implicitly bump the others — the schemas evolve on independent axes
so consumers of one format aren't forced to re-validate when only
another format churned.

| Contract | Constant | Governing ADR |
|---|---|---|
| JSON report | `allelix.reports.json_report.SCHEMA_VERSION` | ADR-0033 |
| Annotated VCF output | `allelix.reports.vcf.VCF_SCHEMA_VERSION` | ADR-0036 |

**Hard rule.** Any PR that adds, removes, renames, or changes the
value grammar of a field in an output contract must bump the
corresponding schema-version constant in the same PR and note the
change in the PR body.

Field shapes that trigger a bump:

- JSON report — any field emitted by `render_json` (annotations,
  panel_coverage, license_attributions, filters, top-level metadata)
- Annotated VCF — any Allelix-declared `##INFO=<>` field, any
  `##ALLELIX_*` provenance line format, the `##ALLELIX_Attribution=<>`
  grammar, or `ALLELIX_ORIGINAL_ID` semantics

`VCF_SCHEMA_VERSION` is currently at `0.x` (experimental). While it
sits at `0.x`, INFO-shape churn is approved without ceremony because
the version conveys instability to downstream readers. The `1.0.0`
graduation is a follow-up ADR after real user feedback; it locks the contract and revokes the
churn freedom. `SCHEMA_VERSION` is stable and every field addition
bumps it.

### Coverage gate: CI is the source of truth

Local `pytest` and the CI matrix can report different `--cov-branch`
totals (~0.1-0.2pp drift) because newer Python versions execute a
handful of lines that 3.11 / 3.12 don't (CPython internals, `__future__`
import paths, error-message construction). The coverage gate
(`--cov-fail-under` in `pyproject.toml`) is calibrated against the **CI
matrix**, not against local measurement.

When raising the floor:

1. Confirm CI is green on `dev` at the candidate threshold.
2. Pull the CI run's reported coverage:
   `gh run view <run_id> --log-failed | grep TOTAL`
   (or `--log` on a passing run; `TOTAL` line is the project number).
3. Set the new floor slightly below the lowest CI number across the
   3.11 / 3.12 / 3.13 matrix. Local numbers are not the source of truth.

The floor-history comment in `pyproject.toml` records the CI-measured
baseline at each ratchet so the next bump calibrates correctly.

## Pull Request Checklist

- [ ] Tests pass: `pytest`
- [ ] Lint clean: `ruff check .`
- [ ] Format clean: `ruff format --check .`
- [ ] License header on new files
- [ ] No private or identifying genetic data in fixtures
- [ ] Source attribution on all annotations

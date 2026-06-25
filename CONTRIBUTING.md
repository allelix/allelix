# Contributing to Allelix

Allelix is an open-source genotype analysis toolkit licensed under AGPL-3.0-or-later.
Contributions are welcome.

## Development Setup

```bash
git clone https://github.com/allelix/allelix.git
cd allelix
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
git config core.hooksPath .githooks
pre-commit install --hook-type pre-commit
```

Run the test suite:

```bash
pytest
```

Tests marked `@slow` exercise real-data invariants against fixtures
that are gitignored due to size. As of v2.0.2 (GH #45), these fixtures
are **auto-fetched on first use** and cached locally — no manual
download step is required for the GWAS catalog.

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

CI runs the same `pytest` invocation. The auto-fetch handles fixture
availability so CI and local runs are identical.

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
        with open(file_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("# VendorName"):
                    return True
                if not line.startswith("#"):
                    break
        return False

    def parse(self, file_path: Path) -> Iterator[Variant]:
        """Yield Variants from the file. Stream, don't load into memory."""
        with open(file_path, encoding="utf-8") as fh:
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
                    allele1=parts[3],
                    allele2=parts[4],
                )

    def get_metadata(self, file_path: Path) -> GenotypeMetadata:
        """Extract metadata from comment headers."""
        sample_id = ""
        with open(file_path, encoding="utf-8") as fh:
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
- No-calls use `"-"` as the allele value (matches `allelix.models.NO_CALL_MARKER`).

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
Put more specific parsers before generic ones.

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
pytest tests/parsers/test_vendorname.py -v
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
            "SELECT alt, significance, condition, gene "
            "FROM mydb_variants WHERE rsid = ?",
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
        db_path = self.data_dir / "mydb.sqlite"
        return db_path.exists()

    def version(self) -> str | None:
        ...

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
            self._conn = sqlite3.connect(self.data_dir / "mydb.sqlite")
        return self._conn
```

Key rules:

- Every annotation must be **source-attributed**. Set `attribution` to your
  database name. Never omit it.
- **Check genotype, not just rsID.** The user must carry the flagged allele.
  A variant existing in the database means nothing if the user has the
  reference (normal) allele.
- Implement `close()` to release SQLite connections.
- No-calls return empty (no annotation possible without a genotype).

### Step 2: Register the annotator

Add to `allelix/annotators/__init__.py` in `get_annotators()`:

```python
from allelix.annotators.mydb import MyDBAnnotator

def get_annotators(data_dir, ...):
    # ... existing annotators ...
    mydb = MyDBAnnotator(data_dir)
    return [clinvar, pharmgkb, gwas, snpedia, mydb]
```

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
pytest tests/annotators/test_mydb.py -v
ruff check allelix/annotators/mydb.py tests/annotators/test_mydb.py
```

## Architecture Notes

- Parsers are stateless. Annotators hold database connections.
- All annotators run on every variant (unlike parsers, which are exclusive).
- The `Annotation.significance` field is always source-prefixed
  (`clinvar_pathogenic`, not `pathogenic`).
- Reports never assert significance directly. They attribute: "ClinVar
  classifies this as pathogenic", not "this is pathogenic."
- The `data/` directory at project root is the local database cache.
  It is gitignored and populated by `allelix db update`.

## Hooks and CI

Two hooks run locally:

- **pre-commit** (managed by pre-commit framework): `ruff check` + `ruff format --check`
- **pre-push** (raw hook in `.githooks/`): blocks tag pushes where the tag doesn't match `pyproject.toml`

CI (`.github/workflows/ci.yml`) fires on:

- Push to `main` or `dev`
- Push of a `v*` tag
- Pull requests targeting `main` or `dev` (so feature PRs to `dev` are
  coverage / mypy / matrix-gated at PR time, not only post-merge)
- Manual dispatch via the Actions UI (`workflow_dispatch`)

CI runs ruff (lint + format), mypy, and the fast pytest suite (synthetic
fixtures only) on a Python 3.11 / 3.12 matrix. `@slow` and
`@integration`-marked tests auto-skip in CI because the runner doesn't
have the ~1.5 GB real-genotype fixture set or the ~15 GB annotator
database cache — fetching either would push CI runtime past 45 minutes
per push and burn GitHub Actions minutes that the local-run discipline
already covers.

### Where the slow / real-data battery actually runs

**Locally, at ship time, by the developer.** `GITHUB_WORKFLOW.md` Phase
0c (the pre-Phase-1 gate) requires the full pytest suite + the real-data
battery from `test_data/FULL_TEST_PROTOCOL.md` to pass against the dev
tree before any squash. The reviewer asks for that report at ship time
and blocks the release without it. That's the gate that catches real-data
regressions; CI runs the fast checks, the ship procedure runs the
expensive ones.

### Why we don't put the slow battery in CI

- The slow tests don't change on every commit. They only matter when
  parser / annotator / exporter code changes or before a ship — both
  cases are already covered by the local-run discipline.
- A 45-minute CI run per push-to-main would burn ~225 minutes/month for
  zero new signal beyond what the ship gate already produces.
- The two failure modes CI-on-slow-tests would defend against —
  developer forgets to run locally, or a mid-cycle regression goes
  unnoticed until ship — are both caught at Phase 0c, which is enforced
  by the reviewer.

When in doubt, you can manually fire the existing CI workflow against
any branch via `gh workflow run ci.yml --ref <branch>` or the "Run
workflow" button in the Actions UI. The fast suite catches everything
CI is supposed to catch.

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
- Commit messages (release commit on `dev`, private squash on
  `main`, public squash on the public mirror, annotated tag
  body)
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
- **Wrong-allele safety**: 0 via-complement CADD scores.
- **Harness invariants**: vocabulary union holds, floor counts
  cleared, spot-checks on published HG002 carriers (rs1801133 /
  rs7412 / rs2010963 / rs6025) pass.

The `tests/test_changelog_sanitize.py` CI guard enforces this on
`CHANGELOG.md` — drift-vulnerable numeric tokens in unsnapshotted
release sections fail the build. The commit-message and release-
notes paths are policy-enforced (reviewer responsibility): scan
release-content drafts before publishing.

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
   3.11 / 3.12 matrix. Local numbers are not the source of truth.

The floor-history comment in `pyproject.toml` records the CI-measured
baseline at each ratchet so the next bump calibrates correctly.

## Pull Request Checklist

- [ ] Tests pass: `pytest`
- [ ] Lint clean: `ruff check .`
- [ ] Format clean: `ruff format --check .`
- [ ] License header on new files
- [ ] No private or identifying genetic data in fixtures
- [ ] Source attribution on all annotations

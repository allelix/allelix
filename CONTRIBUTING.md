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

Tests marked `@slow` require `test_data/gwas_catalog.zip` (~150 MB
compressed, ~400 MB unzipped). Fetch it with:

    scripts/fetch_testdata.sh

Slow tests are skipped automatically when the fixture is absent. CI does
not fetch this fixture and skips slow tests — the CI suite uses small
synthetic fixtures only.

### Run the full suite locally

Before pushing, run the complete test suite with slow tests included:

    scripts/fetch_testdata.sh   # one-time download
    pytest                      # runs everything: fast + slow

This is the only place slow tests run. CI uses small synthetic fixtures
and skips slow tests to keep runs fast and disk-friendly. If you add or
change anything that touches real-data parsing paths, verify it locally
with the full suite before pushing.

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

All `tests/fixtures/` files are synthetic (produced by mock data
generators). Real-data integration tests use CC0 public-domain openSNP
genotype files available via `scripts/fetch_testdata.sh`.

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

CI runs the fast test suite (synthetic fixtures only) on every push to
main and on pull requests. Slow tests that require real-data fixtures run
locally only. There is no pre-push pytest gate — run the full suite
yourself before pushing.

## Pull Request Checklist

- [ ] Tests pass: `pytest`
- [ ] Lint clean: `ruff check .`
- [ ] Format clean: `ruff format --check .`
- [ ] License header on new files
- [ ] No private or identifying genetic data in fixtures
- [ ] Source attribution on all annotations

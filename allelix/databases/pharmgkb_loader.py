# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinPGx clinical-annotation download, parse, and load into SQLite.

ClinPGx publishes a `clinicalAnnotations.zip` containing two TSVs:

- `clinical_annotations.tsv`: one row per clinical annotation
  (id, variant/haplotypes, gene, drug(s), phenotype(s), level of evidence,
  score, phenotype category, …)
- `clinical_ann_alleles.tsv`: per-genotype rows for each annotation
  (annotation id, genotype/allele, annotation text, allele function)

This loader joins the two on annotation id and emits one record per
(rsid, genotype) pair. Star alleles, multi-rsid composites, and indel
genotypes are skipped — they require haplotype reconstruction.

See ADR-0009 for the genotype-matching rationale.
"""

from __future__ import annotations

import contextlib
import csv
import logging
import os
import re
import sqlite3
import tempfile
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

from allelix.databases.schema import PHARMGKB_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

# PharmGKB rebranded to ClinPGx in 2026 (https://www.clinpgx.org).
# `api.pharmgkb.org` still answers via 301 → `api.clinpgx.org` and both
# endpoints 303 to the same `s3.pgkb.org` archive, but the new domain is
# the right one to call directly.
PHARMGKB_CLINICAL_URL = "https://api.clinpgx.org/v1/download/file/data/clinicalAnnotations.zip"
PHARMGKB_DB_FILENAME = "pharmgkb.sqlite"

INSERT_BATCH_SIZE = 5_000

CLINICAL_ANN_FILENAME = "clinical_annotations.tsv"
CLINICAL_ANN_ALLELES_FILENAME = "clinical_ann_alleles.tsv"

# Structural format validation only — NOT prose classification.
# Per ADR-0016, regex is permitted for ID format checking and shape
# validation; it's forbidden as input to any classification decision.
_RSID_RE = re.compile(r"^rs\d+$")
_TWO_LETTER_GENOTYPE_RE = re.compile(r"^[ACGT]{2}$")

# ADR-0020 (v0.9.0): per-allele function lives in the structured CPIC API,
# fetched into `pharmgkb_allele_function` at db-update time and queried as
# a join. The filter is: for the user's `(rsid, genotype)`, look up each
# base in the lookup; if every base maps to Normal function, the row is a
# non-finding. No regex, no prose parsing, no description classification.
#
# Function class enumeration mirrors CPIC's structured field. Values
# outside this set are treated as not-Normal (variant) and emit the row.
FUNCTION_CLASS_NORMAL = "normal"
FUNCTION_CLASS_DECREASED = "decreased"
FUNCTION_CLASS_NO_FUNCTION = "no_function"
FUNCTION_CLASS_INCREASED = "increased"
FUNCTION_CLASS_UNKNOWN = "unknown"

# Schema migration. v0.5.x lacks `function_class`; v0.6.x lacks the
# `pharmgkb_allele_function` table. `schema_is_current()` returns False on
# either, so `db update` automatically refreshes into the v0.9.0 schema.
_REQUIRED_PHARMGKB_COLUMNS = frozenset(
    {
        "rsid",
        "genotype",
        "gene",
        "drugs",
        "phenotype",
        "phenotype_category",
        "annotation_text",
        "level_of_evidence",
        "score",
        "pgkb_annotation_id",
        "allele_function",
        "function_class",
        "is_nonfinding",
    }
)
_REQUIRED_PHARMGKB_TABLES = frozenset({"pharmgkb_annotations", "pharmgkb_allele_function"})


def classify_function(allele_function: str | None) -> str:
    """Map ClinPGx's `Allele Function` field to a stable enum string.

    The structured field is authoritative (ADR-0016). When it's empty, we
    record `unknown` rather than guess from prose — the user sees the row
    and decides what to do with it.
    """
    if not allele_function:
        return FUNCTION_CLASS_UNKNOWN
    value = allele_function.strip().lower()
    if "no function" in value:
        return FUNCTION_CLASS_NO_FUNCTION
    if "decreased" in value:
        return FUNCTION_CLASS_DECREASED
    if "increased" in value:
        return FUNCTION_CLASS_INCREASED
    if "normal" in value:
        return FUNCTION_CLASS_NORMAL
    return FUNCTION_CLASS_UNKNOWN


def is_nonfinding_for_row(
    allele_function: str | None,
    annotation_text: str | None = None,  # kept for back-compat; unused.
    *,
    rsid: str | None = None,
    genotype: str | None = None,
    allele_function_lookup: dict[tuple[str, str], str] | None = None,
) -> bool:
    """Decide whether a row is a non-finding (ADR-0020, v0.9.0).

    The filter is a join, not a text classifier:

    1. **ClinPGx's structured `Allele Function` column** (ADR-0016).
       Authoritative on the rare row where ClinPGx populates it for an
       SNV genotype (most SNV rows have it empty).

    2. **CPIC per-allele function lookup** (ADR-0020). For the row's
       `(rsid, genotype)`, look up each user-carried base in the
       `pharmgkb_allele_function` table. If every base maps to
       `Normal function`, the row is a non-finding. If any base is
       non-Normal — or absent from the lookup for an rsid that HAS
       entries — the row emits.

    If neither tier has data for an rsid, the row emits (rows are never
    silently suppressed without structured evidence).
    """
    function_class = classify_function(allele_function)
    if function_class != FUNCTION_CLASS_UNKNOWN:
        return function_class == FUNCTION_CLASS_NORMAL

    if rsid and genotype and allele_function_lookup is not None:
        lookup_result = is_nonfinding_by_allele_lookup(rsid, genotype, allele_function_lookup)
        if lookup_result is not None:
            return lookup_result

    return False


def is_nonfinding(function_class: str) -> bool:
    """Structured-only non-finding check (back-compat shim for tests).

    Returns True iff function_class == 'normal'. Production code should
    use `is_nonfinding_for_row()` which also handles the empty-field
    prose fallback per ADR-0017.
    """
    return function_class == FUNCTION_CLASS_NORMAL


def schema_is_current(db_path: Path) -> bool:
    """True iff the cache has the v0.7.0 ClinPGx schema.

    v0.7.0 requires both:
      - all v0.6.0 columns on pharmgkb_annotations
      - the new pharmgkb_allele_function table (ADR-0018)
    """
    if not db_path.exists():
        return False
    try:
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            try:
                tables = {
                    row[0]
                    for row in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
                }
                if not _REQUIRED_PHARMGKB_TABLES.issubset(tables):
                    return False
                cols = {row[1] for row in conn.execute("PRAGMA table_info(pharmgkb_annotations)")}
            except sqlite3.DatabaseError:
                return False
            return _REQUIRED_PHARMGKB_COLUMNS.issubset(cols)
    except sqlite3.DatabaseError:
        return False


def is_nonfinding_by_allele_lookup(
    rsid: str,
    genotype: str,
    allele_function_lookup: dict[tuple[str, str], str],
) -> bool | None:
    """Per-allele structured classifier (ADR-0018).

    Returns True if every allele in the user's genotype is either absent
    from the CPIC lookup (Normal-by-absence) or explicitly classified as
    Normal function. Returns False if ANY allele has a flagged non-Normal
    function. Returns None if the lookup has no entries for this rsid at
    all (callers fall back to prose).
    """
    if len(genotype) != 2:
        return None
    rsid_has_entries = any(k[0] == rsid for k in allele_function_lookup)
    if not rsid_has_entries:
        return None
    for allele in set(genotype.upper()):
        function = allele_function_lookup.get((rsid, allele))
        # Under ADR-0020, the CPIC source classifies every allele ClinPGx
        # cares about — Normal for reference, non-Normal for variant. An
        # allele MISSING from the lookup at an rsid that otherwise has
        # entries is an uncharacterized base; never silently suppressed.
        if function != FUNCTION_CLASS_NORMAL:
            return False
    return True


def _normalize_genotype(raw: str) -> str | None:
    """Return a sorted 2-letter SNV genotype, or None if not an SNV diploid call."""
    cleaned = raw.replace(":", "").replace(";", "").replace("/", "").strip().upper()
    if not _TWO_LETTER_GENOTYPE_RE.match(cleaned):
        return None
    return "".join(sorted(cleaned))


def _is_single_rsid(variant_field: str) -> bool:
    """True if the Variant/Haplotypes field is a single rsid."""
    return bool(_RSID_RE.match(variant_field.strip()))


def _open_directory(zip_or_dir: Path) -> tuple[Path, tempfile.TemporaryDirectory | None]:
    """Return a directory path containing the TSVs.

    If `zip_or_dir` is a directory, return it (no cleanup needed).
    If it's a ZIP, extract to a temp dir and return (path, tempdir to clean).
    """
    if zip_or_dir.is_dir():
        return zip_or_dir, None
    tmp = tempfile.TemporaryDirectory(prefix="allelix-pharmgkb-")
    # Python 3.11+ zipfile.extractall sanitizes "../" and absolute paths in
    # member names. The project pins requires-python >= 3.11 (pyproject.toml).
    with zipfile.ZipFile(zip_or_dir) as zf:
        zf.extractall(tmp.name)
    return Path(tmp.name), tmp


def iter_pharmgkb_records(
    zip_or_dir: Path,
    allele_function_lookup: dict[tuple[str, str], str] | None = None,
) -> Iterator[dict[str, object]]:
    """Yield one record per (rsid, genotype) pair from a clinical annotations dump.

    Skips:
      - rows whose Variant/Haplotypes is not a single rsid (star alleles,
        multi-variant composites)
      - per-allele rows whose Genotype/Allele is not a 2-letter SNV genotype
        (indels, star alleles)
    """
    dir_path, tmp = _open_directory(zip_or_dir)
    try:
        annotations: dict[str, dict[str, str]] = {}
        ann_tsv = dir_path / CLINICAL_ANN_FILENAME
        alleles_tsv = dir_path / CLINICAL_ANN_ALLELES_FILENAME
        if not ann_tsv.exists() or not alleles_tsv.exists():
            raise FileNotFoundError(
                f"ClinPGx dump missing required TSVs in {dir_path}: "
                f"need {CLINICAL_ANN_FILENAME} + {CLINICAL_ANN_ALLELES_FILENAME}"
            )

        with ann_tsv.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                ann_id = row.get("Clinical Annotation ID", "").strip()
                variant = row.get("Variant/Haplotypes", "").strip()
                if not ann_id or not _is_single_rsid(variant):
                    continue
                annotations[ann_id] = {
                    "rsid": variant,
                    "gene": row.get("Gene", "").strip(),
                    "drugs": row.get("Drug(s)", "").strip(),
                    "phenotype": row.get("Phenotype(s)", "").strip(),
                    "phenotype_category": row.get("Phenotype Category", "").strip(),
                    "level_of_evidence": row.get("Level of Evidence", "").strip(),
                    "score": row.get("Score", "").strip(),
                }

        with alleles_tsv.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                ann_id = row.get("Clinical Annotation ID", "").strip()
                if ann_id not in annotations:
                    continue
                normalized = _normalize_genotype(row.get("Genotype/Allele", ""))
                if normalized is None:
                    continue
                meta = annotations[ann_id]
                allele_function = row.get("Allele Function", "").strip()
                function_class = classify_function(allele_function)
                annotation_text = row.get("Annotation Text", "").strip()
                yield {
                    "rsid": meta["rsid"],
                    "genotype": normalized,
                    "gene": meta["gene"],
                    "drugs": meta["drugs"],
                    "phenotype": meta["phenotype"],
                    "phenotype_category": meta["phenotype_category"],
                    "annotation_text": annotation_text,
                    "level_of_evidence": meta["level_of_evidence"],
                    "score": _safe_float(meta["score"]),
                    "pgkb_annotation_id": ann_id,
                    "allele_function": allele_function,
                    "function_class": function_class,
                    "is_nonfinding": is_nonfinding_for_row(
                        allele_function,
                        annotation_text,
                        rsid=meta["rsid"],
                        genotype=normalized,
                        allele_function_lookup=allele_function_lookup,
                    ),
                }
    finally:
        if tmp is not None:
            tmp.cleanup()


def _safe_float(value: str) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_pharmgkb_tsv(
    zip_or_dir: Path,
    db_path: Path,
    source_url: str = "",
    version: str = "",
    remote_signal: str | None = None,
    allele_function_lookup: dict[tuple[str, str], str] | None = None,
) -> int:
    """Load a ClinPGx clinical-annotations dump into a fresh SQLite cache atomically.

    Writes to a `.tmp` SQLite then `os.replace`s onto `db_path`. A failed
    mid-parse leaves the previous cache (if any) intact.

    `allele_function_lookup` is the structured `(rsid, base) → function_class`
    table that drives the non-finding filter (ADR-0020). Production fetches
    it from CPIC's API; tests inject a synthetic dict directly. When None
    the loader falls back to an empty lookup — every row emits.

    `remote_signal` is the value `fetch_remote_signal` returned at the time
    of this download; stored so the next `db update` can detect remote
    changes without re-downloading.
    """
    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    resolved_version = version or datetime.now(UTC).strftime("%Y-%m-%d")
    lookup = allele_function_lookup or {}

    try:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.executescript(PHARMGKB_SCHEMA)

            # Populate the per-allele function table (ADR-0020) first.
            # The lookup arrives pre-built from cpic_loader.fetch_cpic_allele_functions
            # (production) or a test fixture (unit tests).
            af_insert_sql = (
                "INSERT INTO pharmgkb_allele_function "
                "(rsid, allele, function_class, source) "
                "VALUES (?, ?, ?, 'cpic_api')"
            )
            for (rsid, allele), function_class in lookup.items():
                conn.execute(af_insert_sql, (rsid, allele, function_class))

            insert_sql = (
                "INSERT INTO pharmgkb_annotations "
                "(rsid, genotype, gene, drugs, phenotype, phenotype_category, "
                "annotation_text, level_of_evidence, score, pgkb_annotation_id, "
                "allele_function, function_class, is_nonfinding) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            )
            batch: list[tuple] = []
            count = 0
            for record in iter_pharmgkb_records(zip_or_dir, lookup):
                batch.append(
                    (
                        record["rsid"],
                        record["genotype"],
                        record["gene"],
                        record["drugs"],
                        record["phenotype"],
                        record["phenotype_category"],
                        record["annotation_text"],
                        record["level_of_evidence"],
                        record["score"],
                        record["pgkb_annotation_id"],
                        record["allele_function"],
                        record["function_class"],
                        int(bool(record["is_nonfinding"])),
                    )
                )
                if len(batch) >= INSERT_BATCH_SIZE:
                    conn.executemany(insert_sql, batch)
                    count += len(batch)
                    batch.clear()
            if batch:
                conn.executemany(insert_sql, batch)
                count += len(batch)
            from allelix.databases._versions import PHARMGKB_INTERPRETER_VERSION

            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    "pharmgkb",
                    source_url,
                    resolved_version,
                    datetime.now(UTC).isoformat(),
                    count,
                    remote_signal or "",
                    f"iv:{PHARMGKB_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()
        os.replace(tmp_path, db_path)
        return count
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                logger.warning("Could not remove failed temp DB %s", tmp_path)
        raise

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Parse raw SNPedia wiki markup into structured genotype rows.

Called automatically by the SNPedia annotator when raw pages exist but
the structured ``snpedia_genotypes`` table does not. Can also be invoked
standalone via ``scripts/parse_snpedia.py``.
"""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from datetime import UTC, datetime

import mwparserfromhell

logger = logging.getLogger(__name__)

_PARSER_VERSION = 6


def _parse_title_prefix(title: str) -> tuple[str, str] | None:
    """Extract (prefix, number) from a title like 'Rs12345' or 'I4000178'.

    Returns None if the title doesn't start with Rs or I followed by digits.
    """
    if title.startswith("Rs") or title.startswith("I"):
        prefix = "Rs" if title.startswith("Rs") else "I"
        rest = title[len(prefix) :]
        digits = []
        for ch in rest:
            if ch.isdigit():
                digits.append(ch)
            else:
                break
        if digits:
            return prefix, "".join(digits)
    return None


def _parse_title_alleles(title: str) -> tuple[str, str] | None:
    """Extract alleles from a title like 'Rs12345(A;G)' or 'I4000178(C;T)'.

    Returns None if the title doesn't contain a valid (allele;allele) suffix.
    """
    paren_start = title.find("(")
    if paren_start == -1 or not title.endswith(")"):
        return None
    inner = title[paren_start + 1 : -1]
    semi = inner.find(";")
    if semi == -1:
        return None
    a1 = inner[:semi].strip()
    a2 = inner[semi + 1 :].strip()
    if a1 and a2:
        return a1, a2
    return None


def _tmpl_param(tmpl: object, name: str) -> str:
    """Extract a named parameter from a mwparserfromhell template."""
    # mwparserfromhell ships no type stubs (see [[tool.mypy.overrides]]); the
    # parameter is typed as `object` at the boundary, then unwrapped here.
    if tmpl.has(name):  # type: ignore[attr-defined]
        return str(tmpl.get(name).value).strip()  # type: ignore[attr-defined]
    return ""


_STRUCTURED_SCHEMA = """
CREATE TABLE IF NOT EXISTS snpedia_genotypes (
    rsid TEXT NOT NULL,
    allele1 TEXT NOT NULL,
    allele2 TEXT NOT NULL,
    magnitude REAL,
    repute TEXT,
    summary TEXT,
    gene TEXT,
    scraped_at TEXT
);

CREATE INDEX IF NOT EXISTS idx_snpedia_rsid_alleles
    ON snpedia_genotypes(rsid, allele1, allele2);

CREATE UNIQUE INDEX IF NOT EXISTS idx_snpedia_genotype_dedup
    ON snpedia_genotypes(rsid, allele1, allele2, COALESCE(summary, ''));

CREATE TABLE IF NOT EXISTS database_versions (
    name TEXT PRIMARY KEY,
    source_url TEXT NOT NULL,
    version TEXT,
    downloaded_at TEXT NOT NULL,
    record_count INTEGER NOT NULL,
    remote_signal TEXT,
    local_version_tag TEXT
);
"""


def parser_is_current(conn: sqlite3.Connection) -> bool:
    """Return True if the cache was built by the current parser version.

    Checks ``local_version_tag`` first.  If absent, falls back to the
    legacy ``|pv:N`` suffix in ``remote_signal`` and migrates the tag
    in-place — avoiding a full re-parse just to populate the column.
    """
    tag = f"pv:{_PARSER_VERSION}"
    try:
        row = conn.execute(
            "SELECT local_version_tag FROM database_versions WHERE name='snpedia'"
        ).fetchone()
        if row and row[0] == tag:
            return True
    except sqlite3.OperationalError:
        pass
    try:
        row = conn.execute(
            "SELECT remote_signal FROM database_versions WHERE name='snpedia'"
        ).fetchone()
        if row and row[0] and f"|pv:{_PARSER_VERSION}" in row[0]:
            from allelix.databases.manager import _ensure_local_version_tag_column

            _ensure_local_version_tag_column(conn)
            clean_signal = row[0].split("|pv:")[0]
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = ? WHERE name = 'snpedia'",
                (clean_signal, tag),
            )
            conn.commit()
            return True
    except sqlite3.OperationalError:
        pass
    return False


def _dedupe_existing(conn: sqlite3.Connection) -> int:
    """Collapse pre-existing duplicate rows in old caches. Returns rows removed."""
    before = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
    conn.execute("""
        DELETE FROM snpedia_genotypes
        WHERE rowid NOT IN (
            SELECT MIN(rowid) FROM snpedia_genotypes
            GROUP BY rsid, UPPER(allele1), UPPER(allele2), COALESCE(summary, '')
        )
    """)
    after = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
    return before - after


# GH #12: identifier allowlist for the raw_table f-string interpolation
# below. The interpolation cannot be parameterized (SQLite doesn't support
# bind variables for identifiers); the allowlist gives a programmatic
# guarantee that only these two literals can reach the SQL.
_VALID_RAW_TABLES: frozenset[str] = frozenset({"_raw_pages", "pages"})


def detect_raw_table(conn: sqlite3.Connection) -> str | None:
    """Return the name of the raw pages table, or None if absent."""
    tables = {
        row[0]
        for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    }
    if "_raw_pages" in tables:
        return "_raw_pages"
    if "pages" in tables:
        return "pages"
    return None


def has_structured_table(conn: sqlite3.Connection) -> bool:
    """Return True if snpedia_genotypes exists and has rows."""
    try:
        count = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
        return count > 0
    except sqlite3.OperationalError:
        return False


def parse_raw_pages(db_path: str, *, verbose: bool = False) -> int:
    """Parse raw wiki markup into structured genotype rows.

    Returns the number of structured rows created.
    """
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        return _parse_raw_pages_inner(conn, verbose=verbose)


def _parse_raw_pages_inner(conn: sqlite3.Connection, *, verbose: bool = False) -> int:
    """Inner parser logic. Caller owns the connection lifecycle."""
    raw_table = detect_raw_table(conn)
    if raw_table is None:
        return 0
    # GH #12: `raw_table` flows into three SQL queries via f-string
    # interpolation because SQLite doesn't support parameterized
    # identifiers. Today it's safe — `detect_raw_table` only ever returns
    # one of two literals or None — but the function's `str | None` return
    # type doesn't pin that. A future edit (config-driven table name,
    # scraped metadata) could drift it into an injection path. Allowlist
    # explicitly so the guarantee outlives memory of the original design.
    if raw_table not in _VALID_RAW_TABLES:
        msg = f"unexpected raw table name: {raw_table!r}"
        raise ValueError(msg)

    if verbose:
        logger.info("Parsing SNPedia raw pages from '%s' table", raw_table)

    conn.execute("DROP INDEX IF EXISTS idx_snpedia_genotype_dedup")

    has_table = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='snpedia_genotypes'"
    ).fetchone()
    if has_table:
        deduped = _dedupe_existing(conn)
        if deduped:
            logger.info("Backfill dedupe: removed %d duplicate row(s)", deduped)

    conn.executescript(_STRUCTURED_SCHEMA)
    conn.execute("DELETE FROM snpedia_genotypes")

    # Build gene map from SNP pages
    gene_map: dict[str, str] = {}
    snp_rows = conn.execute(
        f"SELECT title, content FROM {raw_table} WHERE category = 'snp'"
    ).fetchall()
    print(f"  Building gene map from {len(snp_rows)} SNP pages...", flush=True)
    for title, content in snp_rows:
        parsed_prefix = _parse_title_prefix(title)
        if not parsed_prefix or not content:
            continue
        prefix, num = parsed_prefix
        snp_key = f"{prefix.lower()}{num}"
        try:
            wikicode = mwparserfromhell.parse(content)
            for template in wikicode.filter_templates():
                tname = template.name.strip().lower()
                if tname in ("rsnum", "snp"):
                    gene = _tmpl_param(template, "Gene")
                    if gene:
                        gene_map[snp_key] = gene
                    break
                if tname == "23andme snp":
                    gene = _tmpl_param(template, "Gene_s")
                    if gene:
                        gene_map[snp_key] = gene
                    break
        except Exception:
            logger.debug("Failed to parse SNP page %s", title, exc_info=True)
            continue

    print(f"  Gene map: {len(gene_map)} mappings built.", flush=True)

    # Parse genotype pages
    genotype_rows = conn.execute(
        f"SELECT title, content, scraped_at FROM {raw_table} WHERE category = 'genotype'"
    ).fetchall()
    print(f"  Parsing {len(genotype_rows)} genotype pages...", flush=True)

    batch: list[tuple[str, str, str, float | None, str | None, str | None, str | None, str]] = []

    for title, content, scraped_at in genotype_rows:
        parsed_prefix = _parse_title_prefix(title)
        if not parsed_prefix or not content:
            continue

        prefix, num = parsed_prefix
        snp_id = f"{prefix.lower()}{num}"

        try:
            wikicode = mwparserfromhell.parse(content)
        except Exception:
            logger.debug("Failed to parse genotype page %s", title, exc_info=True)
            continue

        templates = [
            t for t in wikicode.filter_templates() if t.name.strip().lower() == "genotype"
        ]
        if not templates:
            continue

        tmpl = templates[0]

        allele1 = _tmpl_param(tmpl, "allele1").upper()
        allele2 = _tmpl_param(tmpl, "allele2").upper()
        if not allele1 or not allele2:
            title_alleles = _parse_title_alleles(title)
            if not title_alleles:
                continue
            allele1, allele2 = title_alleles[0].upper(), title_alleles[1].upper()
            if not allele1 or not allele2:
                continue

        if allele1 > allele2:
            allele1, allele2 = allele2, allele1

        mag_str = _tmpl_param(tmpl, "magnitude")
        magnitude: float | None = None
        if mag_str:
            try:
                magnitude = float(mag_str)
            except ValueError:
                magnitude = None

        repute = _tmpl_param(tmpl, "repute") or None
        summary = _tmpl_param(tmpl, "summary") or None
        snp_gene = gene_map.get(snp_id) or None

        batch.append((snp_id, allele1, allele2, magnitude, repute, summary, snp_gene, scraped_at))

        if len(batch) >= 1000:
            conn.executemany(
                "INSERT OR IGNORE INTO snpedia_genotypes "
                "(rsid, allele1, allele2, magnitude, repute, summary, gene, scraped_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                batch,
            )
            batch.clear()

    if batch:
        conn.executemany(
            "INSERT OR IGNORE INTO snpedia_genotypes "
            "(rsid, allele1, allele2, magnitude, repute, summary, gene, scraped_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            batch,
        )

    row_count = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]

    date_row = conn.execute(f"SELECT MIN(scraped_at) FROM {raw_table}").fetchone()
    scrape_date = date_row[0][:10] if date_row and date_row[0] else "unknown"

    existing_signal = ""
    try:
        sig_row = conn.execute(
            "SELECT remote_signal FROM database_versions WHERE name = 'snpedia'"
        ).fetchone()
        if sig_row and sig_row[0]:
            existing_signal = sig_row[0].split("|pv:")[0]
    except sqlite3.OperationalError:
        pass

    conn.execute("DELETE FROM database_versions WHERE name = 'snpedia'")
    conn.execute(
        "INSERT INTO database_versions "
        "(name, source_url, version, downloaded_at, record_count, "
        "remote_signal, local_version_tag) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            "snpedia",
            "https://bots.snpedia.com/api.php",
            f"scraped {scrape_date} ({row_count} genotypes)",
            datetime.now(UTC).isoformat(),
            row_count,
            existing_signal,
            f"pv:{_PARSER_VERSION}",
        ),
    )

    conn.commit()
    return row_count

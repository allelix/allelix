# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
#!/usr/bin/env python3
"""One-time SNPedia archival scraper.

Downloads every page from Category:Is_a_snp and Category:Is_a_genotype
into a single SQLite database. Stores raw wiki markup with no parsing
or filtering. Resumable: skips titles already in the database.

Output defaults to the allelix data directory (same location used by
``allelix db update``). Override with --output.

Usage:
    python scripts/scrape_snpedia.py
    python scripts/scrape_snpedia.py --output /path/to/snpedia.sqlite
"""

from __future__ import annotations

import argparse
import functools
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterator

print = functools.partial(print, flush=True)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from allelix.databases import resolve_data_dir  # noqa: E402

API_URL = "https://bots.snpedia.com/api.php"
USER_AGENT = "allelix-scraper/1.0 (archival; contact@redneck.travel)"
CM_BATCH_SIZE = 500
CONTENT_BATCH_SIZE = 50
COMMIT_EVERY = 100
REQUEST_DELAY = 1
MAX_RETRIES = 8
RETRY_BACKOFF = (2, 5, 10, 20, 30, 60, 60, 60)

SCHEMA = """
CREATE TABLE IF NOT EXISTS pages (
    title TEXT PRIMARY KEY,
    category TEXT,
    content TEXT,
    scraped_at TEXT
);

CREATE TABLE IF NOT EXISTS enumerated_titles (
    title TEXT NOT NULL,
    category TEXT NOT NULL,
    PRIMARY KEY (title, category)
);
"""

CATEGORIES = [
    ("Category:Is_a_snp", "snp"),
    ("Category:Is_a_genotype", "genotype"),
]


def api_get(params: dict[str, str]) -> dict[str, Any]:
    """GET request to the SNPedia MediaWiki API with retry."""
    qs = "&".join(f"{k}={urllib.parse.quote(str(v))}" for k, v in params.items())
    url = f"{API_URL}?{qs}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_exc = None
    for attempt in range(MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 502, 503) or attempt == MAX_RETRIES:
                raise
            last_exc = exc
        except (TimeoutError, OSError) as exc:
            if attempt == MAX_RETRIES:
                raise
            last_exc = exc
        delay = RETRY_BACKOFF[min(attempt, len(RETRY_BACKOFF) - 1)]
        print(f"  retry {attempt + 1}/{MAX_RETRIES} in {delay}s ({last_exc})")
        time.sleep(delay)
    raise last_exc


def enumerate_titles(category_title: str) -> Iterator[str]:
    """Yield all page titles from a MediaWiki category."""
    cmcontinue = None
    while True:
        params = {
            "action": "query",
            "list": "categorymembers",
            "cmtitle": category_title,
            "cmlimit": str(CM_BATCH_SIZE),
            "format": "json",
        }
        if cmcontinue:
            params["cmcontinue"] = cmcontinue
        data = api_get(params)
        for member in data.get("query", {}).get("categorymembers", []):
            yield member["title"]
        if "continue" not in data:
            break
        cmcontinue = data["continue"]["cmcontinue"]
        time.sleep(REQUEST_DELAY)


def fetch_content_batch(titles: list[str]) -> dict[str, str]:
    """Fetch page content for up to 50 titles."""
    params = {
        "action": "query",
        "prop": "revisions",
        "rvprop": "content",
        "format": "json",
        "titles": "|".join(titles),
    }
    data = api_get(params)
    results = {}
    for page in data.get("query", {}).get("pages", {}).values():
        title = page.get("title", "")
        revisions = page.get("revisions", [])
        content = revisions[0].get("*", "") if revisions else ""
        results[title] = content
    return results


def main() -> None:
    """Scrape all of SNPedia into a local SQLite archive."""
    default_output = str(resolve_data_dir() / "snpedia.sqlite")
    parser = argparse.ArgumentParser(description="Scrape all of SNPedia into SQLite.")
    parser.add_argument(
        "--output",
        default=default_output,
        help=f"Output SQLite file (default: {default_output})",
    )
    args = parser.parse_args()

    conn = sqlite3.connect(args.output)
    conn.executescript(SCHEMA)

    existing = set()
    for row in conn.execute("SELECT title FROM pages"):
        existing.add(row[0])
    if existing:
        print(f"Resuming: {len(existing)} pages already in database")

    # Phase 1: enumerate titles from both categories (saved to DB for resume)
    saved_titles = {}
    for row in conn.execute("SELECT title, category FROM enumerated_titles"):
        saved_titles[row[0]] = row[1]

    if saved_titles:
        print(f"Resuming with {len(saved_titles)} previously enumerated titles")
        all_titles = saved_titles
    else:
        all_titles = {}
        for cat_title, cat_label in CATEGORIES:
            print(f"Enumerating {cat_title}...")
            count = 0
            batch = []
            for title in enumerate_titles(cat_title):
                if title not in all_titles:
                    all_titles[title] = cat_label
                    batch.append((title, cat_label))
                count += 1
                if len(batch) >= 500:
                    conn.executemany(
                        "INSERT OR IGNORE INTO enumerated_titles (title, category) VALUES (?, ?)",
                        batch,
                    )
                    conn.commit()
                    batch = []
                if count % 5000 == 0:
                    print(f"  {count} titles...")
            if batch:
                conn.executemany(
                    "INSERT OR IGNORE INTO enumerated_titles (title, category) VALUES (?, ?)",
                    batch,
                )
                conn.commit()
            print(f"  {count} titles from {cat_title}")

    print(f"Total unique titles: {len(all_titles)}")

    # Phase 2: fetch content for titles not already in database
    to_fetch = [(t, c) for t, c in all_titles.items() if t not in existing]
    print(f"Need to fetch: {len(to_fetch)} pages")

    total = len(all_titles)
    fetched = len(existing)
    failed_titles = []
    inserts_since_commit = 0

    for batch_start in range(0, len(to_fetch), CONTENT_BATCH_SIZE):
        batch = to_fetch[batch_start : batch_start + CONTENT_BATCH_SIZE]
        batch_titles = [t for t, _ in batch]
        batch_categories = {t: c for t, c in batch}

        try:
            contents = fetch_content_batch(batch_titles)
        except Exception as exc:
            print(f"  FAILED batch at {batch_start}: {exc}")
            failed_titles.extend(batch_titles)
            time.sleep(REQUEST_DELAY)
            continue

        now = datetime.now(UTC).isoformat()
        for title in batch_titles:
            content = contents.get(title, "")
            category = batch_categories[title]
            conn.execute(
                "INSERT OR IGNORE INTO pages (title, category, content, scraped_at) "
                "VALUES (?, ?, ?, ?)",
                (title, category, content, now),
            )
            inserts_since_commit += 1

        if inserts_since_commit >= COMMIT_EVERY:
            conn.commit()
            inserts_since_commit = 0

        fetched += len(batch)
        if fetched % 1000 < CONTENT_BATCH_SIZE:
            pct = fetched / total * 100
            print(f"Fetched {fetched}/{total} pages ({pct:.1f}%)")

        time.sleep(REQUEST_DELAY)

    conn.commit()

    # Summary
    print("\n=== COMPLETE ===")
    rows = conn.execute("SELECT category, COUNT(*) FROM pages GROUP BY category").fetchall()
    total_pages = 0
    for cat, count in rows:
        print(f"  {cat}: {count}")
        total_pages += count
    print(f"  TOTAL: {total_pages}")

    empty = conn.execute(
        "SELECT COUNT(*) FROM pages WHERE content IS NULL OR content = ''"
    ).fetchone()[0]
    print(f"  Empty/stub pages: {empty}")

    i_count = conn.execute("SELECT COUNT(*) FROM pages WHERE title LIKE 'I%'").fetchone()[0]
    print(f"  i-prefixed entries: {i_count}")

    if failed_titles:
        print(f"  FAILED fetches: {len(failed_titles)}")
        for t in failed_titles[:20]:
            print(f"    {t}")
        print("\n  Run again to retry failed pages (resume is automatic).")
    else:
        print("  Failed fetches: 0")

    conn.close()
    print(f"\nDatabase: {args.output}")


if __name__ == "__main__":
    main()

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""CPIC API → per-allele function lookup.

CPIC (Clinical Pharmacogenetics Implementation Consortium) publishes
structured per-allele functional status assignments via its public
PostgREST API. Each allele (haplotype or single variant) carries a
`clinicalfunctionalstatus` chosen from a small enumeration:

    Normal function | Decreased function | No function |
    Increased function | Uncertain function | (gene-specific tags)

The ClinPGx filter is a join: for an annotation row matching the
user's `(rsid, genotype)`, look up each base of the user's genotype
in the per-allele function table. If both alleles map to
`Normal function`, the row is a non-finding (the user does not carry
the studied variant). Otherwise the row emits.

ADR-0020 documents this as the canonical structured source. Three
CPIC tables are joined client-side:

    sequence_location.dbsnpid (rsid)
        ↔ allele_location_value.locationid
            ↔ allele_location_value.alleledefinitionid
                ↔ allele.definitionid
                    → allele.clinicalfunctionalstatus

Result: `(rsid, base) → function_class` for every CPIC-curated
variant. Genes outside CPIC's scope have no entries; the filter
treats absence as "no opinion" and emits the row.
"""

from __future__ import annotations

import json
import logging
import re
import time
import urllib.error
import urllib.request

from allelix import __version__
from allelix.databases.manager import PROBE_ERRORS

logger = logging.getLogger(__name__)

CPIC_API_BASE = "https://api.cpicpgx.org/v1"
CPIC_TIMEOUT_SECONDS = 60
CPIC_MAX_ROWS = 99_999  # PostgREST defaults to limit=1000 without explicit Range.

# Retry on transient failures. CPIC's PostgREST API is generally
# reliable but TCP RSTs and brief 5xx blips do happen; one retry burst
# saves the user from a manual `db update --force` rerun. Backoff is
# capped low because the loader runs interactively at db-update time.
CPIC_RETRY_ATTEMPTS = 3
CPIC_RETRY_BACKOFF_SECONDS: tuple[float, ...] = (1.0, 2.0, 4.0)

USER_AGENT = f"allelix/{__version__} (+https://github.com/allelix/allelix)"

# CPIC's clinicalfunctionalstatus enumeration → Allelix's function_class enum.
# Anything not mapped here is treated as "not Normal" (i.e., a variant) and
# the row containing that allele emits. We never coerce an unknown status
# into Normal - silent suppression is the failure mode v0.5-v0.8 kept
# producing.
FUNCTION_CLASS_NORMAL = "normal"
FUNCTION_CLASS_DECREASED = "decreased"
FUNCTION_CLASS_NO_FUNCTION = "no_function"
FUNCTION_CLASS_INCREASED = "increased"
FUNCTION_CLASS_UNCERTAIN = "uncertain"

_CPIC_TO_FUNCTION_CLASS: dict[str, str] = {
    "normal function": FUNCTION_CLASS_NORMAL,
    "decreased function": FUNCTION_CLASS_DECREASED,
    "no function": FUNCTION_CLASS_NO_FUNCTION,
    "increased function": FUNCTION_CLASS_INCREASED,
    "possibly increased function": FUNCTION_CLASS_INCREASED,
    "uncertain function": FUNCTION_CLASS_UNCERTAIN,
}


def _classify_cpic_status(status: str | None) -> str | None:
    """Map a CPIC `clinicalfunctionalstatus` string to a function_class enum.

    Returns None for empty/unrecognized strings — the caller skips the
    row rather than guessing.
    """
    if not status:
        return None
    return _CPIC_TO_FUNCTION_CLASS.get(status.strip().lower())


def fetch_cpic_remote_signal(api_base: str = CPIC_API_BASE) -> str | None:
    """Return a freshness signal for CPIC's data, or None on failure.

    ClinPGx's bulk-download Last-Modified header tells us nothing
    about CPIC's curation database. CPIC publishes a `change_log` table
    with one row per curated change; the most recent date is a stable,
    cheap freshness proxy. The signal format is `lastchange:{date}`.

    Never retries — this is the lightweight probe used at `db update`
    freshness-check time. Persistent CPIC outages should NOT block the
    user; returning None signals "can't verify" and the CLI prints
    "pass --force to refresh anyway" rather than aborting.
    """
    url = f"{api_base}/change_log?select=date&order=date.desc&limit=1"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=CPIC_TIMEOUT_SECONDS) as response:
            rows = json.loads(response.read().decode("utf-8"))
    except PROBE_ERRORS as exc:
        if hasattr(exc, "close"):
            exc.close()
        logger.warning("CPIC freshness probe failed: %s", exc)
        return None
    # Valid JSON of another shape (an object, a scalar) is an unusable reply too.
    if not isinstance(rows, list) or not rows or not isinstance(rows[0], dict):
        return None
    date = rows[0].get("date")
    if not date:
        return None
    return f"lastchange:{date}"


class CpicIncompleteError(ValueError):
    """A CPIC response that does not hold the whole table it was asked for."""


_CONTENT_RANGE = re.compile(r"(?:(\d+)-(\d+)|\*)/(\d+|\*)")

# The type of each column read from the three tables. A value of another type
# (a number where a string is read, a list where a key is) is no CPIC table
# either; null is allowed, and the row is skipped (2.6.1).
_COLUMN_TYPES: dict[str, dict[str, type]] = {
    "sequence_location": {"id": int, "dbsnpid": str},
    "allele_location_value": {"alleledefinitionid": int, "locationid": int, "variantallele": str},
    "allele": {"definitionid": int, "clinicalfunctionalstatus": str},
}


def _check_whole_table(url: str, rows: object, content_range: str | None) -> None:
    """Raise CpicIncompleteError unless ``rows`` is the whole table (2.6.1).

    The response must be a list of objects whose columns have their types
    (``_COLUMN_TYPES``), its range must start at row 0 and cover
    every row of the body, and the body must hold the stated total. A reply
    with no total ("*", or no Content-Range at all) is not proof of a whole
    table: a server can cap a reply below the window asked for. Each request
    asks for the exact total (``Prefer: count=exact``), and CPIC's API gives
    it for all three tables (``0-991/992``, ``0-3159/3160``, ``0-1377/1378``
    on 2026-11-05), so a reply without one is incomplete.
    """
    endpoint = url.split("?", 1)[0].rsplit("/", 1)[-1]
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        # Each row is read with .get(): a list of anything but objects is no
        # table either (2.6.1).
        raise CpicIncompleteError(f"CPIC {endpoint}: the response is not a list of row objects")
    for row in rows:
        for column, kind in _COLUMN_TYPES.get(endpoint, {}).items():
            value = row.get(column)
            if value is not None and (not isinstance(value, kind) or isinstance(value, bool)):
                raise CpicIncompleteError(
                    f"CPIC {endpoint}: {column} is a {type(value).__name__}, not a {kind.__name__}"
                )
    total: int | None = None
    if content_range is not None:
        match = _CONTENT_RANGE.fullmatch(content_range.strip())
        if match is None:
            raise CpicIncompleteError(
                f"CPIC {endpoint}: unreadable Content-Range {content_range!r}"
            )
        first, last, stated = match.groups()
        if first is not None and (int(first) != 0 or int(last) - int(first) + 1 != len(rows)):
            raise CpicIncompleteError(
                f"CPIC {endpoint}: range {content_range} does not match "
                f"the {len(rows)} rows returned"
            )
        if first is None and rows:
            raise CpicIncompleteError(
                f"CPIC {endpoint}: range {content_range} with {len(rows)} rows returned"
            )
        if stated != "*":
            total = int(stated)
    if total is None:
        raise CpicIncompleteError(f"CPIC {endpoint}: the reply gives no total row count")
    if len(rows) != total:
        raise CpicIncompleteError(f"CPIC {endpoint}: {len(rows)} of {total} rows returned")


def _http_get_json(url: str, timeout: float = CPIC_TIMEOUT_SECONDS) -> list[dict]:
    """Fetch a CPIC PostgREST endpoint and return the JSON body.

    Sends `Range: 0-N` to bypass PostgREST's default 1000-row cap.
    Retries up to `CPIC_RETRY_ATTEMPTS` times with exponential backoff
    on transient transport failures. Raises the last exception
    after the final attempt so the caller surfaces a clear failure
    instead of silently producing an empty lookup.
    """
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
            "Range-Unit": "items",
            "Range": f"0-{CPIC_MAX_ROWS}",
            # Without it PostgREST answers a slice with 200 and "0-N/*", so a
            # partial table cannot be told from a whole one (2.6.1).
            "Prefer": "count=exact",
        },
    )
    last_error: Exception | None = None
    for attempt in range(CPIC_RETRY_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                rows = json.loads(response.read().decode("utf-8"))
                _check_whole_table(url, rows, response.headers.get("Content-Range"))
                return rows
        except PROBE_ERRORS as exc:
            if hasattr(exc, "close"):
                exc.close()
            last_error = exc
            if attempt + 1 < CPIC_RETRY_ATTEMPTS:
                backoff = CPIC_RETRY_BACKOFF_SECONDS[attempt]
                logger.warning(
                    "CPIC fetch failed (attempt %d/%d): %s — retrying in %.1fs",
                    attempt + 1,
                    CPIC_RETRY_ATTEMPTS,
                    exc,
                    backoff,
                )
                time.sleep(backoff)
    assert last_error is not None  # loop runs at least once
    raise last_error


def fetch_cpic_allele_functions(
    api_base: str = CPIC_API_BASE,
) -> dict[tuple[str, str], str]:
    """Build the `(rsid, base) → function_class` lookup from CPIC's API.

    ADR-0020: this IS the ClinPGx non-finding filter's data source.
    Three CPIC tables are fetched and joined client-side:

      - `sequence_location` (id, dbsnpid)
      - `allele_location_value` (alleledefinitionid, locationid, variantallele)
      - `allele` (definitionid, clinicalfunctionalstatus)

    Only single-base alleles (A/C/G/T) are emitted; CPIC's tables also
    contain multi-base haplotype components which don't apply to the
    SNV genotype-matching path (ADR-0009).

    On network failure, raises `urllib.error.URLError` (or similar).
    The caller decides whether to abort `db update` or fall back to a
    cached lookup.
    """
    seq_url = f"{api_base}/sequence_location?dbsnpid=not.is.null&select=id,dbsnpid"
    loc_url = (
        f"{api_base}/allele_location_value?select=alleledefinitionid,locationid,variantallele"
    )
    allele_url = f"{api_base}/allele?select=definitionid,clinicalfunctionalstatus"

    sequence_locations = _http_get_json(seq_url)
    location_values = _http_get_json(loc_url)
    alleles = _http_get_json(allele_url)

    location_to_rsid: dict[int, str] = {}
    for row in sequence_locations:
        loc_id = row.get("id")
        rsid = row.get("dbsnpid")
        if loc_id is not None and rsid:
            location_to_rsid[loc_id] = rsid

    allele_to_function: dict[int, str] = {}
    for row in alleles:
        definition_id = row.get("definitionid")
        function_class = _classify_cpic_status(row.get("clinicalfunctionalstatus"))
        if definition_id is not None and function_class is not None:
            allele_to_function[definition_id] = function_class

    out: dict[tuple[str, str], str] = {}
    for row in location_values:
        rsid = location_to_rsid.get(row.get("locationid"))  # type: ignore[arg-type]
        function_class = allele_to_function.get(row.get("alleledefinitionid"))  # type: ignore[arg-type]
        base = (row.get("variantallele") or "").strip().upper()
        if not rsid or function_class is None:
            continue
        if len(base) != 1 or base not in "ACGT":
            continue
        # Conflict policy: when the same (rsid, base) appears under multiple
        # allele definitions with different function classes, prefer the
        # non-Normal classification. Suppressing happens only when EVERY
        # base maps to Normal, so when CPIC's own data has a Normal-vs-non-
        # Normal conflict the safe choice is "treat as variant and emit"
        # — never silently suppress a real variant. In practice CPIC's
        # tables are internally consistent; this is defense in depth.
        prev = out.get((rsid, base))
        if prev is None:
            out[(rsid, base)] = function_class
        elif prev != function_class and FUNCTION_CLASS_NORMAL in (prev, function_class):
            out[(rsid, base)] = function_class if prev == FUNCTION_CLASS_NORMAL else prev
        # else: both classifications agree, or both non-Normal — keep first.
    return out

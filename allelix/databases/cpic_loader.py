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
import time
import urllib.error
import urllib.request

from allelix import __version__

logger = logging.getLogger(__name__)

CPIC_API_BASE = "https://api.cpicpgx.org/v1"
CPIC_TIMEOUT_SECONDS = 60
CPIC_MAX_ROWS = 99_999  # PostgREST defaults to limit=1000 without explicit Range.

# M-1: retry on transient failures. CPIC's PostgREST API is generally
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

    M-2: ClinPGx's bulk-download Last-Modified header tells us nothing
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
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        if hasattr(exc, "close"):
            exc.close()
        logger.warning("CPIC freshness probe failed: %s", exc)
        return None
    if not rows or not isinstance(rows[0], dict):
        return None
    date = rows[0].get("date")
    if not date:
        return None
    return f"lastchange:{date}"


def _http_get_json(url: str, timeout: float = CPIC_TIMEOUT_SECONDS) -> list[dict]:
    """Fetch a CPIC PostgREST endpoint and return the JSON body.

    Sends `Range: 0-N` to bypass PostgREST's default 1000-row cap.
    Retries up to `CPIC_RETRY_ATTEMPTS` times with exponential backoff
    on transient transport failures (M-1). Raises the last exception
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
        },
    )
    last_error: Exception | None = None
    for attempt in range(CPIC_RETRY_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
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

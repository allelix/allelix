# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A CPIC table is used only when the response holds all of it.

CPIC's PostgREST API answers a range request with a slice. Without
``Prefer: count=exact`` a slice comes back as HTTP 200 with
``Content-Range: 0-0/*``, so the status cannot tell a partial table from a
whole one; with it, a slice is 206 with ``0-0/1378`` and the whole table 200
with ``0-1377/1378`` (observed on api.cpicpgx.org, 2026-11-05). A partial
``allele_location_value`` could drop a Decreased-function row and hide a
ClinPGx finding. Every endpoint now asks for the total and rejects a body
that is not the whole table; the rejection takes the existing CPIC-unavailable
path. Synthetic replies only; nothing reaches the network.
"""

from __future__ import annotations

import io
import json
from typing import TYPE_CHECKING

import pytest

from allelix.databases import cpic_loader

if TYPE_CHECKING:
    import urllib.request

_RS = "rs900000002"
_SEQUENCE = [{"id": 1, "dbsnpid": _RS}]
_LOCATIONS = [
    {"alleledefinitionid": 1, "locationid": 1, "variantallele": "A"},
    {"alleledefinitionid": 2, "locationid": 1, "variantallele": "A"},
]
_ALLELES = [
    {"definitionid": 1, "clinicalfunctionalstatus": "Normal function"},
    {"definitionid": 2, "clinicalfunctionalstatus": "Decreased function"},
]
_TABLES = {
    "/sequence_location?": _SEQUENCE,
    "/allele_location_value?": _LOCATIONS,
    "/allele?": _ALLELES,
}


# The healthy replies have the shape CPIC's API gives for the three requests
# (status 200, ``Content-Range: 0-<n-1>/<n>``), observed live on 2026-11-05:
# sequence_location 0-991/992, allele_location_value 0-3159/3160,
# allele 0-1377/1378.
class _Reply(io.BytesIO):
    def __init__(self, rows, *, status=200, content_range=None):
        super().__init__(json.dumps(rows).encode())
        self.status = status
        self.headers = {} if content_range is None else {"Content-Range": content_range}


def _whole(rows):
    return _Reply(rows, content_range=f"0-{len(rows) - 1}/{len(rows)}" if rows else "*/0")


@pytest.fixture
def cpic(monkeypatch):
    """Serve the three tables; ``replies[endpoint]`` overrides one endpoint's reply."""
    replies: dict[str, object] = {}
    requests: list[urllib.request.Request] = []

    def urlopen(request, **kwargs):
        requests.append(request)
        for endpoint, rows in _TABLES.items():
            if endpoint in request.full_url:
                override = replies.get(endpoint)
                return override(rows) if callable(override) else _whole(rows)
        raise AssertionError(request.full_url)

    monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(cpic_loader, "CPIC_RETRY_BACKOFF_SECONDS", (0.0, 0.0, 0.0))
    return replies, requests


def test_the_whole_tables_resolve_the_conflict_to_decreased(cpic):
    lookup = cpic_loader.fetch_cpic_allele_functions()

    assert lookup[(_RS, "A")] == "decreased"


def test_every_request_asks_for_the_exact_total(cpic):
    _, requests = cpic
    cpic_loader.fetch_cpic_allele_functions()

    assert len(requests) == 3
    for request in requests:
        assert request.get_header("Prefer") == "count=exact"


@pytest.mark.parametrize("endpoint", list(_TABLES))
def test_a_partial_slice_of_any_table_is_rejected(cpic, endpoint):
    replies, _ = cpic
    # One row short of the table: an empty slice of a one-row table.
    replies[endpoint] = lambda rows: _Reply(
        rows[:-1],
        status=206,
        content_range=f"0-{len(rows) - 2}/{len(rows)}" if len(rows) > 1 else f"*/{len(rows)}",
    )

    with pytest.raises(cpic_loader.CpicIncompleteError):
        cpic_loader.fetch_cpic_allele_functions()


@pytest.mark.parametrize(
    ("description", "reply"),
    [
        ("200 with a short body", lambda rows: _Reply(rows[:1], content_range="0-1/2")),
        ("range not from 0", lambda rows: _Reply(rows, status=206, content_range="1-2/3")),
        ("range shorter than the body", lambda rows: _Reply(rows, content_range="0-0/2")),
        ("total larger than the range", lambda rows: _Reply(rows, content_range="0-1/3")),
        ("unparsable range", lambda rows: _Reply(rows, content_range="items 0-1/2")),
        ("not a list", lambda rows: _Reply({"rows": rows}, content_range="0-1/2")),
    ],
)
def test_a_response_that_is_not_the_whole_table_is_rejected(cpic, description, reply):
    replies, _ = cpic
    replies["/allele_location_value?"] = reply

    with pytest.raises(cpic_loader.CpicIncompleteError):
        cpic_loader.fetch_cpic_allele_functions()


@pytest.mark.parametrize(
    "reply",
    [
        lambda rows: _Reply(rows, content_range=f"0-{len(rows) - 1}/*"),
        lambda rows: _Reply(rows),
    ],
    ids=["unknown-total", "no-content-range"],
)
def test_a_reply_without_a_total_is_incomplete_even_when_whole(cpic, reply):
    """The body here is in fact the whole table, but nothing in the reply says
    so: a server can cap a reply below the window asked for."""
    replies, _ = cpic
    replies["/allele_location_value?"] = reply

    with pytest.raises(cpic_loader.CpicIncompleteError, match="no total"):
        cpic_loader.fetch_cpic_allele_functions()


def test_a_table_larger_than_the_window_is_rejected(cpic, monkeypatch):
    replies, _ = cpic
    monkeypatch.setattr(cpic_loader, "CPIC_MAX_ROWS", 1)
    replies["/allele_location_value?"] = lambda rows: _Reply(
        rows, status=206, content_range="0-1/5"
    )

    with pytest.raises(cpic_loader.CpicIncompleteError):
        cpic_loader.fetch_cpic_allele_functions()


def test_an_empty_table_is_whole(cpic):
    replies, _ = cpic
    replies["/sequence_location?"] = lambda rows: _Reply([], content_range="*/0")

    assert cpic_loader.fetch_cpic_allele_functions() == {}


def test_an_incomplete_response_is_retried_before_it_is_rejected(cpic):
    replies, requests = cpic
    replies["/allele?"] = lambda rows: _Reply(rows[:1], status=206, content_range="0-0/2")

    with pytest.raises(cpic_loader.CpicIncompleteError):
        cpic_loader.fetch_cpic_allele_functions()

    tries = [r for r in requests if "/allele?" in r.full_url]
    assert len(tries) == cpic_loader.CPIC_RETRY_ATTEMPTS

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for CPIC API → per-allele function lookup.

ADR-0020: CPIC's structured per-allele clinicalfunctionalstatus is the
canonical source for the ClinPGx non-finding filter. These tests
exercise the loader without touching the real CPIC API — `urllib`'s
urlopen is monkey-patched to return canned JSON.
"""

from __future__ import annotations

import io
import json
from contextlib import contextmanager

import pytest

from allelix.databases import cpic_loader
from allelix.databases.cpic_loader import (
    FUNCTION_CLASS_DECREASED,
    FUNCTION_CLASS_INCREASED,
    FUNCTION_CLASS_NO_FUNCTION,
    FUNCTION_CLASS_NORMAL,
    FUNCTION_CLASS_UNCERTAIN,
    _classify_cpic_status,
    fetch_cpic_allele_functions,
    fetch_cpic_remote_signal,
)


def _stub_responses(seq_locs: list[dict], loc_values: list[dict], alleles: list[dict]):
    """Build a urlopen replacement that returns canned JSON per CPIC endpoint."""

    payloads = {
        "sequence_location": seq_locs,
        "allele_location_value": loc_values,
        "allele": alleles,
    }

    @contextmanager
    def fake_urlopen(request, timeout=None):
        url = request.full_url if hasattr(request, "full_url") else str(request)
        for key, payload in payloads.items():
            if f"/{key}?" in url or url.endswith(f"/{key}"):
                yield io.BytesIO(json.dumps(payload).encode("utf-8"))
                return
        raise AssertionError(f"unexpected URL: {url}")

    return fake_urlopen


class TestClassifyCpicStatus:
    """M-4: mutation-gap pin. Unknown statuses must NEVER coerce to Normal."""

    def test_known_statuses_map(self):
        assert _classify_cpic_status("Normal function") == FUNCTION_CLASS_NORMAL
        assert _classify_cpic_status("Decreased function") == FUNCTION_CLASS_DECREASED
        assert _classify_cpic_status("No function") == FUNCTION_CLASS_NO_FUNCTION
        assert _classify_cpic_status("Increased function") == FUNCTION_CLASS_INCREASED
        assert _classify_cpic_status("Possibly Increased function") == FUNCTION_CLASS_INCREASED
        assert _classify_cpic_status("Uncertain function") == FUNCTION_CLASS_UNCERTAIN

    def test_status_lookup_is_case_insensitive(self):
        assert _classify_cpic_status("normal function") == FUNCTION_CLASS_NORMAL
        assert _classify_cpic_status("NO FUNCTION") == FUNCTION_CLASS_NO_FUNCTION

    def test_unknown_status_returns_none(self):
        """Critical contract: unknown statuses MUST NOT map to Normal. The
        failure mode v0.5-v0.8 kept producing was silent suppression via
        unknown-treated-as-default. Pin it here.
        """
        assert _classify_cpic_status("Malignant Hyperthermia associated") is None
        assert _classify_cpic_status("Some Future Status") is None
        assert _classify_cpic_status("typo function") is None

    def test_empty_status_returns_none(self):
        assert _classify_cpic_status(None) is None
        assert _classify_cpic_status("") is None
        assert _classify_cpic_status("   ") is None


class TestFetchCpicAlleleFunctions:
    """ADR-0020 end-to-end with mocked urlopen."""

    def test_three_way_join(self, monkeypatch):
        seq_locs = [
            {"id": 1, "dbsnpid": "rs116855232"},
            {"id": 2, "dbsnpid": "rs1800559"},
        ]
        loc_values = [
            {"alleledefinitionid": 100, "locationid": 1, "variantallele": "C"},
            {"alleledefinitionid": 101, "locationid": 1, "variantallele": "T"},
            {"alleledefinitionid": 200, "locationid": 2, "variantallele": "C"},
        ]
        alleles = [
            {"definitionid": 100, "clinicalfunctionalstatus": "Normal function"},
            {"definitionid": 101, "clinicalfunctionalstatus": "No function"},
            {"definitionid": 200, "clinicalfunctionalstatus": "Normal function"},
        ]
        monkeypatch.setattr(
            cpic_loader.urllib.request, "urlopen", _stub_responses(seq_locs, loc_values, alleles)
        )

        result = fetch_cpic_allele_functions()
        assert result == {
            ("rs116855232", "C"): FUNCTION_CLASS_NORMAL,
            ("rs116855232", "T"): FUNCTION_CLASS_NO_FUNCTION,
            ("rs1800559", "C"): FUNCTION_CLASS_NORMAL,
        }

    def test_multi_base_alleles_skipped(self, monkeypatch):
        """M-4: mutation-gap pin. CPIC's tables also contain haplotype
        components like 'CT' or 'CGTACG' that aren't valid SNV genotype
        bases. The loader must skip them — otherwise the lookup would
        contain garbage keys that mis-classify real user genotypes.
        """
        seq_locs = [{"id": 1, "dbsnpid": "rs1"}]
        loc_values = [
            {"alleledefinitionid": 100, "locationid": 1, "variantallele": "C"},
            {"alleledefinitionid": 101, "locationid": 1, "variantallele": "CT"},
            {"alleledefinitionid": 102, "locationid": 1, "variantallele": ""},
            {"alleledefinitionid": 103, "locationid": 1, "variantallele": "ACGT"},
            {"alleledefinitionid": 104, "locationid": 1, "variantallele": "N"},
        ]
        alleles = [
            {"definitionid": d, "clinicalfunctionalstatus": "Normal function"}
            for d in (100, 101, 102, 103, 104)
        ]
        monkeypatch.setattr(
            cpic_loader.urllib.request, "urlopen", _stub_responses(seq_locs, loc_values, alleles)
        )

        result = fetch_cpic_allele_functions()
        assert result == {("rs1", "C"): FUNCTION_CLASS_NORMAL}

    def test_unknown_status_alleles_skipped(self, monkeypatch):
        """An allele whose clinicalfunctionalstatus is not in the enum
        (e.g., gene-specific labels like 'Malignant Hyperthermia
        associated') doesn't get a function_class — so it doesn't end up
        in the lookup either. Downstream filter then emits the row.
        """
        seq_locs = [{"id": 1, "dbsnpid": "rs1800559"}]
        loc_values = [
            {"alleledefinitionid": 100, "locationid": 1, "variantallele": "C"},
            {"alleledefinitionid": 101, "locationid": 1, "variantallele": "T"},
        ]
        alleles = [
            {"definitionid": 100, "clinicalfunctionalstatus": "Normal function"},
            {"definitionid": 101, "clinicalfunctionalstatus": "Malignant Hyperthermia associated"},
        ]
        monkeypatch.setattr(
            cpic_loader.urllib.request, "urlopen", _stub_responses(seq_locs, loc_values, alleles)
        )

        result = fetch_cpic_allele_functions()
        # C is Normal; T is not in the lookup (unknown status). User CC →
        # both Normal → non-finding. User CT or TT → T missing → finding.
        assert result == {("rs1800559", "C"): FUNCTION_CLASS_NORMAL}

    def test_dbsnp_null_locations_filtered(self, monkeypatch):
        """sequence_location rows with null dbsnpid don't carry an rsid;
        the resulting allele entries would have no key. Skip them.
        """
        seq_locs = [
            {"id": 1, "dbsnpid": None},
            {"id": 2, "dbsnpid": "rs1"},
        ]
        loc_values = [
            {"alleledefinitionid": 100, "locationid": 1, "variantallele": "C"},
            {"alleledefinitionid": 101, "locationid": 2, "variantallele": "G"},
        ]
        alleles = [
            {"definitionid": 100, "clinicalfunctionalstatus": "Normal function"},
            {"definitionid": 101, "clinicalfunctionalstatus": "Normal function"},
        ]
        monkeypatch.setattr(
            cpic_loader.urllib.request, "urlopen", _stub_responses(seq_locs, loc_values, alleles)
        )

        result = fetch_cpic_allele_functions()
        assert result == {("rs1", "G"): FUNCTION_CLASS_NORMAL}

    def test_conflict_prefers_non_normal(self, monkeypatch):
        """M-3 pin: when CPIC's data has a Normal-vs-non-Normal conflict for
        the same (rsid, base), prefer the non-Normal classification. The
        filter only suppresses when EVERY base is Normal — biasing
        conflicts toward non-Normal ensures we never silently suppress a
        real variant just because one CPIC row happened to flag it Normal.
        """
        seq_locs = [{"id": 1, "dbsnpid": "rs_conflict"}]
        loc_values = [
            {"alleledefinitionid": 100, "locationid": 1, "variantallele": "G"},
            {"alleledefinitionid": 101, "locationid": 1, "variantallele": "G"},
        ]
        alleles = [
            {"definitionid": 100, "clinicalfunctionalstatus": "Normal function"},
            {"definitionid": 101, "clinicalfunctionalstatus": "Decreased function"},
        ]
        monkeypatch.setattr(
            cpic_loader.urllib.request, "urlopen", _stub_responses(seq_locs, loc_values, alleles)
        )

        result = fetch_cpic_allele_functions()
        assert result == {("rs_conflict", "G"): FUNCTION_CLASS_DECREASED}

    def test_network_error_propagates_after_retries(self, monkeypatch):
        """M-1: persistent network failure must raise after exhausting
        retries — never silently return an empty lookup (which would
        treat every row as a finding and flood the user with noise).
        """
        import urllib.error

        attempts = []

        @contextmanager
        def boom(*args, **kwargs):
            attempts.append(1)
            raise urllib.error.URLError("simulated transport failure")
            yield  # pragma: no cover - generator contract

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", boom)
        # Speed up the test — no real sleeping.
        monkeypatch.setattr(cpic_loader.time, "sleep", lambda _s: None)

        with pytest.raises(urllib.error.URLError, match="simulated"):
            fetch_cpic_allele_functions()
        assert len(attempts) == cpic_loader.CPIC_RETRY_ATTEMPTS


class TestFetchCpicRemoteSignal:
    """M-2: lightweight CPIC freshness probe via change_log table."""

    def test_returns_lastchange_date(self, monkeypatch):
        @contextmanager
        def fake(request, timeout=None):
            yield io.BytesIO(b'[{"date": "2026-05-11"}]')

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", fake)
        assert fetch_cpic_remote_signal() == "lastchange:2026-05-11"

    def test_network_failure_returns_none(self, monkeypatch):
        import urllib.error

        @contextmanager
        def boom(*args, **kwargs):
            raise urllib.error.URLError("simulated")
            yield  # pragma: no cover

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", boom)
        assert fetch_cpic_remote_signal() is None

    def test_timeout_returns_none(self, monkeypatch):
        @contextmanager
        def slow(*args, **kwargs):
            raise TimeoutError("timed out")
            yield  # pragma: no cover

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", slow)
        assert fetch_cpic_remote_signal() is None

    def test_malformed_json_returns_none(self, monkeypatch):
        @contextmanager
        def garbage(request, timeout=None):
            yield io.BytesIO(b"not json {")

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", garbage)
        assert fetch_cpic_remote_signal() is None

    def test_empty_result_returns_none(self, monkeypatch):
        @contextmanager
        def empty(request, timeout=None):
            yield io.BytesIO(b"[]")

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", empty)
        assert fetch_cpic_remote_signal() is None

    def test_missing_date_field_returns_none(self, monkeypatch):
        @contextmanager
        def no_date(request, timeout=None):
            yield io.BytesIO(b'[{"version": 1}]')

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", no_date)
        assert fetch_cpic_remote_signal() is None


class TestHttpRetryBehavior:
    """M-1: transient failures recover via retry; persistent ones surface."""

    def test_retries_then_succeeds(self, monkeypatch):
        import urllib.error

        attempts = {"n": 0}

        @contextmanager
        def flaky(request, timeout=None):
            attempts["n"] += 1
            if attempts["n"] < 3:
                raise urllib.error.URLError(f"transient error #{attempts['n']}")
            yield io.BytesIO(b'[{"id": 1, "dbsnpid": "rs1"}]')

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", flaky)
        monkeypatch.setattr(cpic_loader.time, "sleep", lambda _s: None)

        result = cpic_loader._http_get_json("https://example/test")
        assert result == [{"id": 1, "dbsnpid": "rs1"}]
        assert attempts["n"] == 3

    def test_timeout_retried(self, monkeypatch):
        attempts = {"n": 0}

        @contextmanager
        def slow(request, timeout=None):
            attempts["n"] += 1
            if attempts["n"] == 1:
                raise TimeoutError("first call timed out")
            yield io.BytesIO(b"[]")

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", slow)
        monkeypatch.setattr(cpic_loader.time, "sleep", lambda _s: None)

        result = cpic_loader._http_get_json("https://example/test")
        assert result == []
        assert attempts["n"] == 2

    def test_malformed_json_retried(self, monkeypatch):
        attempts = {"n": 0}

        @contextmanager
        def garbage_then_good(request, timeout=None):
            attempts["n"] += 1
            if attempts["n"] == 1:
                yield io.BytesIO(b"not json {")
            else:
                yield io.BytesIO(b'[{"ok": true}]')

        monkeypatch.setattr(cpic_loader.urllib.request, "urlopen", garbage_then_good)
        monkeypatch.setattr(cpic_loader.time, "sleep", lambda _s: None)

        result = cpic_loader._http_get_json("https://example/test")
        assert result == [{"ok": True}]
        assert attempts["n"] == 2

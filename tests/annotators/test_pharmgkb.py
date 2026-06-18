# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the ClinPGx annotator."""

from __future__ import annotations

import contextlib
import urllib.error
from typing import TYPE_CHECKING
from unittest.mock import patch

import pytest

from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def annotator(pharmgkb_data_dir: Path):
    ann = PharmGKBAnnotator(pharmgkb_data_dir)
    try:
        yield ann
    finally:
        ann.close()


class TestSetupAndStatus:
    def test_unconfigured_is_not_ready(self, tmp_path: Path):
        ann = PharmGKBAnnotator(tmp_path)
        assert ann.is_ready() is False
        assert ann.version() is None
        assert ann.record_count() is None

    def test_configured_is_ready(self, annotator: PharmGKBAnnotator):
        assert annotator.is_ready() is True
        assert annotator.version() is not None
        assert annotator.record_count() == 16


class TestSignalGuard:
    def test_setup_aborts_when_signal_fetch_fails(self, tmp_path: Path, monkeypatch):
        """setup() raises RuntimeError when remote signal is None."""
        ann = PharmGKBAnnotator(tmp_path)
        monkeypatch.setattr(ann, "fetch_remote_signal", lambda: None)
        with pytest.raises(RuntimeError, match="cannot verify remote freshness signal"):
            ann.setup()


class TestInterpreterVersionStamp:
    """Verify the cache self-heals on interpreter version bumps."""

    def test_stale_stamp_reingests_from_cached_zip(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        """Simulates PHARMGKB_INTERPRETER_VERSION bumping: stale iv:0 → auto-reingest."""
        import sqlite3

        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db_path = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db_path, remote_signal="test-signal")

        zip_path = tmp_path / "clinicalAnnotations.zip"
        import zipfile

        with zipfile.ZipFile(zip_path, "w") as zf:
            for f in mock_pharmgkb_dir.iterdir():
                zf.write(f, f.name)

        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'iv:0' WHERE name = 'pharmgkb'"
            )
            conn.commit()

        ann = PharmGKBAnnotator(tmp_path)
        try:
            assert ann.is_ready() is True
            assert ann.record_count() == 16
        finally:
            ann.close()

    def test_stale_stamp_no_zip_returns_false(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        """Stale stamp without a retained ZIP → not ready."""
        import sqlite3

        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db_path = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db_path, remote_signal="test-signal")

        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'iv:0' WHERE name = 'pharmgkb'"
            )
            conn.commit()

        ann = PharmGKBAnnotator(tmp_path)
        try:
            assert ann.is_ready() is False
        finally:
            ann.close()

    def test_legacy_no_stamp_self_heals(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        """Pre-mechanism cache (no local_version_tag) gets stamped without reingest."""
        import sqlite3

        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        db_path = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, db_path, remote_signal="test-signal")

        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = NULL WHERE name = 'pharmgkb'"
            )
            conn.commit()

        ann = PharmGKBAnnotator(tmp_path)
        try:
            assert ann.is_ready() is True
        finally:
            ann.close()


class TestGenotypeMatching:
    """ADR-0009: ClinPGx matches the user's exact normalized diploid call."""

    def test_heterozygous_match_triggers(self, annotator: PharmGKBAnnotator):
        # rs1801133 stored AG, MHG fixture has G/A (normalizes to AG).
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        a = results[0]
        assert a.significance == "pharmgkb_loe_2a"
        assert a.attribution == "ClinPGx"
        assert a.source == "pharmgkb"
        assert a.category == "pharma"
        assert a.gene == "MTHFR"
        assert a.genotype_match == "AG"
        assert a.magnitude == 7.0
        assert "methotrexate" in a.description.lower()

    def test_homozygous_alt_match_triggers(self, annotator: PharmGKBAnnotator):
        # rs4680 stored AA, MHG fixture has A/A.
        v = Variant("rs4680", "22", 19963748, "A", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].significance == "pharmgkb_loe_3"
        assert "tramadol" in results[0].description.lower()

    def test_genotype_mismatch_does_not_trigger(self, annotator: PharmGKBAnnotator):
        # rs1801133 has rows for AG/AA/GG. User with TT (impossible biologically
        # but tests the lookup) must produce no annotation.
        v = Variant("rs1801133", "1", 11796321, "T", "T")
        assert annotator.annotate(v) == []

    def test_no_call_does_not_trigger(self, annotator: PharmGKBAnnotator):
        v = Variant("rs1801133", "1", 11796321, "-", "-")
        assert annotator.annotate(v) == []

    def test_asymmetric_no_call_does_not_trigger(self, annotator: PharmGKBAnnotator):
        """r-2: one good allele + one no-call must short-circuit before lookup."""
        v_left = Variant("rs1801133", "1", 11796321, "-", "A")
        v_right = Variant("rs1801133", "1", 11796321, "A", "-")
        assert annotator.annotate(v_left) == []
        assert annotator.annotate(v_right) == []

    def test_unknown_rsid_does_not_trigger(self, annotator: PharmGKBAnnotator):
        v = Variant("rs999000111", "1", 1000, "A", "G")
        assert annotator.annotate(v) == []

    def test_indel_does_not_trigger(self, annotator: PharmGKBAnnotator):
        # ClinPGx v0.3.0 doesn't model indels — _normalize_genotype rejects them.
        v = Variant("rs113993960", "7", 117199644, "CTT", "C")
        assert annotator.annotate(v) == []

    def test_alt_threading_with_ref_heterozygous(self, annotator: PharmGKBAnnotator):
        """ADR-0035 PR 2: REF=A, user is A/G → annotation.alt = "G"."""
        v = Variant("rs1801133", "1", 11796321, "G", "A", ref="A")
        results = annotator.annotate(v)
        assert results
        assert results[0].alt == "G"

    def test_alt_threading_without_ref_stays_empty(self, annotator: PharmGKBAnnotator):
        """Array data without Variant.ref leaves alt empty (no guess)."""
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert results
        assert results[0].alt == ""


class TestAttribution:
    """ADR-0003: significance source-prefixed; attribution always present."""

    def test_attribution_consistent(self, annotator: PharmGKBAnnotator):
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        for a in annotator.annotate(v):
            assert a.attribution == "ClinPGx"
            assert a.significance.startswith("pharmgkb_")
            assert a.category == "pharma"
            assert a.description.startswith("ClinPGx:")


class TestRegistryMetadata:
    def test_class_attributes(self):
        assert PharmGKBAnnotator.name == "pharmgkb"
        assert PharmGKBAnnotator.display_name == "ClinPGx"
        assert PharmGKBAnnotator.attribution == "ClinPGx"
        assert PharmGKBAnnotator.requires_download is True


class TestNonFindingSuppression:
    """ADR-0020 (fallback tier): the cache's CPIC-derived `is_nonfinding`
    flag still suppresses reference homozygotes when ClinVar has no REF
    data for the rsid. The bare `annotator` fixture has no ClinVar
    provider wired, so this exercises exactly that fallback path.

    The ClinVar-REF primary tier is exercised separately in
    `TestClinvarRefPrimaryFilter`.
    """

    def test_reference_homozygote_suppressed(self, annotator: PharmGKBAnnotator):
        # PA-008 GG: rs900000010 G → Normal function in the mock CPIC lookup.
        v = Variant("rs900000010", "1", 1, "G", "G")
        assert annotator.annotate(v) == []

    def test_carrier_emits(self, annotator: PharmGKBAnnotator):
        # PA-008 AA: rs900000010 A → Decreased function. Both alleles A → finding.
        v = Variant("rs900000010", "1", 1, "A", "A")
        assert annotator.annotate(v)

    def test_real_carrier_still_fires(self, annotator: PharmGKBAnnotator):
        # rs1801133 AG: A is Decreased, G is Normal → at least one non-Normal → finding.
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        assert annotator.annotate(v)


class TestClinvarRefPrimaryFilter:
    """ADR-0023: the primary non-finding filter is the ClinVar REF check.

    Universal across all genes — CFTR (where CPIC's vocabulary is drug-
    specific), MTHFR/F2/F5 (where CPIC has no data at all), and DPYD
    (where CPIC has proper function classes) all flow through the same
    code path: if ClinVar says REF=X and the user is X/X, suppress.
    """

    def _annotator_with_ref(
        self, pharmgkb_data_dir: Path, refs: dict[str, str]
    ) -> PharmGKBAnnotator:
        ann = PharmGKBAnnotator(
            pharmgkb_data_dir,
            clinvar_ref_provider=lambda rsid, _build: refs.get(rsid),
        )
        return ann

    def test_homozygous_reference_suppressed_via_clinvar_ref(self, pharmgkb_data_dir: Path):
        """rs1801133 GG with ClinVar REF=G → user is hom-ref → suppress.
        The cache row has `is_nonfinding=0` (because CPIC marks A as
        Decreased, G as Normal → GG is non-finding via CPIC too), but
        the test pins that even without CPIC agreeing the REF check is
        decisive.
        """
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs1801133": "G"})
        try:
            v = Variant("rs1801133", "1", 11796321, "G", "G")
            assert ann.annotate(v) == []
        finally:
            ann.close()

    def test_heterozygous_carrier_emits_via_clinvar_ref(self, pharmgkb_data_dir: Path):
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs1801133": "G"})
        try:
            v = Variant("rs1801133", "1", 11796321, "G", "A")
            results = ann.annotate(v)
            assert results, "heterozygous carrier should emit"
            assert results[0].genotype_match == "AG"
        finally:
            ann.close()

    def test_cftr_class_leak_is_now_filtered(self, pharmgkb_data_dir: Path):
        """The user's real-world CFTR leak: CPIC's CFTR table uses
        'ivacaftor responsive' (no Normal entries) so the v0.7.1 join
        had zero CFTR coverage. With ClinVar REF as primary, a CFTR
        rsid with REF=G and user genotype GG is suppressed cleanly.

        Synthesizes the failure shape using the fixture's MTHFR
        rsid + the mock CPIC lookup's existing entries — the structural
        invariant is the same: when CPIC's lookup would say "no usable
        Normal entry," ClinVar REF still suppresses correctly.
        """
        # rs900000010 GG: CPIC fallback would suppress (mock has G=Normal).
        # rs1801133 GG: ClinVar REF=G is decisive — even if CPIC said
        # otherwise, the REF check fires first.
        refs = {"rs1801133": "G"}
        ann = self._annotator_with_ref(pharmgkb_data_dir, refs)
        try:
            v = Variant("rs1801133", "1", 11796321, "G", "G")
            assert ann.annotate(v) == []
        finally:
            ann.close()

    def test_clinvar_missing_falls_through_to_cpic(self, pharmgkb_data_dir: Path):
        """ADR-0023 fallback: rsids ClinVar doesn't know about still
        use the pre-computed `is_nonfinding` flag from the cache.
        """
        # No REF provided for rs900000010 → fall through to CPIC.
        # rs900000010 GG is a non-finding per the mock CPIC lookup.
        ann = self._annotator_with_ref(pharmgkb_data_dir, {})
        try:
            v = Variant("rs900000010", "1", 1, "G", "G")
            assert ann.annotate(v) == []
        finally:
            ann.close()

    def test_multi_base_clinvar_ref_falls_through(self, pharmgkb_data_dir: Path):
        """Indel REFs (e.g. CTT) can't suppress single-base SNV
        genotypes — fall through to CPIC instead.
        """
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs900000010": "CTT"})
        try:
            v = Variant("rs900000010", "1", 1, "G", "G")
            # Falls through to CPIC, which says GG is Normal → still suppressed.
            assert ann.annotate(v) == []
        finally:
            ann.close()

    def test_genotype_match_uses_user_diploid(self, pharmgkb_data_dir: Path):
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs1801133": "G"})
        try:
            v = Variant("rs1801133", "1", 11796321, "A", "A")
            results = ann.annotate(v)
            assert results
            # Sorted diploid representation regardless of input order.
            assert results[0].genotype_match == "AA"
        finally:
            ann.close()

    def test_cpic_normal_suppresses_even_when_clinvar_says_carrier(self, pharmgkb_data_dir: Path):
        """Regression: rs1801265 GG — CPIC assigns Normal function to both
        alleles → is_nonfinding = 1 → must be suppressed regardless of
        ClinVar REF. ClinVar REF=A means the user is NOT hom-ref, but
        the CPIC filter is additive and independently sufficient.
        """
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs1801265": "A"})
        try:
            v = Variant("rs1801265", "1", 97515839, "G", "G")
            assert ann.annotate(v) == []
        finally:
            ann.close()

    def test_non_cpic_gene_hom_ref_suppressed_by_clinvar(self, pharmgkb_data_dir: Path):
        """Regression: rs900000020 GG on a gene with zero CPIC coverage.
        ClinVar REF=G → user is hom-ref → suppressed by the ClinVar
        check even though CPIC can't help.
        """
        ann = self._annotator_with_ref(pharmgkb_data_dir, {"rs900000020": "G"})
        try:
            v = Variant("rs900000020", "1", 1, "G", "G")
            assert ann.annotate(v) == []
        finally:
            ann.close()


class TestSchemaMigration:
    """Pre-v0.6.0 caches lack `function_class` → is_ready returns False so
    db update refreshes into the v0.6.0 schema with structured classification
    (ADR-0016).
    """

    def test_v05x_cache_reports_not_ready(self, tmp_path: Path):
        import sqlite3

        # v0.5.x shape: has is_nonfinding + is_somatic but no function_class.
        db = tmp_path / "pharmgkb.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(
                """
                CREATE TABLE pharmgkb_annotations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    rsid TEXT NOT NULL, genotype TEXT NOT NULL, gene TEXT,
                    drugs TEXT, phenotype TEXT, phenotype_category TEXT,
                    annotation_text TEXT, level_of_evidence TEXT, score REAL,
                    pgkb_annotation_id TEXT, allele_function TEXT,
                    is_nonfinding INTEGER NOT NULL, is_somatic INTEGER NOT NULL
                );
                CREATE TABLE database_versions (
                    name TEXT PRIMARY KEY, source_url TEXT NOT NULL,
                    version TEXT, downloaded_at TEXT NOT NULL,
                    record_count INTEGER NOT NULL, remote_signal TEXT
                );
                """
            )
            conn.execute(
                "INSERT INTO database_versions VALUES (?, ?, ?, ?, ?, ?)",
                ("pharmgkb", "old", "v0.5.x", "2025-01-01T00:00:00", 1, "lm:old"),
            )
            conn.commit()

        ann = PharmGKBAnnotator(tmp_path)
        try:
            assert ann.is_ready() is False
        finally:
            ann.close()


class TestRemoteSignal:
    """Composite freshness signal: ClinPGx (ETag/Last-Modified) + CPIC (M-2)."""

    def test_etag_preferred_when_present(self, annotator: PharmGKBAnnotator, monkeypatch):
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(
            pharmgkb_module,
            "head_request_headers",
            lambda url: {"ETag": '"abc123"', "Last-Modified": "old"},
        )
        monkeypatch.setattr(
            pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:2026-05-11"
        )
        assert annotator.fetch_remote_signal() == 'pgkb:etag:"abc123"|cpic:lastchange:2026-05-11'

    def test_falls_back_to_last_modified(self, annotator: PharmGKBAnnotator, monkeypatch):
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(
            pharmgkb_module,
            "head_request_headers",
            lambda url: {"Last-Modified": "Wed, 21 Oct 2025 07:28:00 GMT"},
        )
        monkeypatch.setattr(
            pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:2026-05-11"
        )
        assert annotator.fetch_remote_signal() == (
            "pgkb:lm:Wed, 21 Oct 2025 07:28:00 GMT|cpic:lastchange:2026-05-11"
        )

    def test_returns_none_when_no_signal_headers(self, annotator: PharmGKBAnnotator, monkeypatch):
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(
            pharmgkb_module, "head_request_headers", lambda url: {"Server": "nginx"}
        )
        monkeypatch.setattr(
            pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:2026-05-11"
        )
        assert annotator.fetch_remote_signal() is None

    def test_returns_none_on_network_error(self, annotator: PharmGKBAnnotator, monkeypatch):
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(pharmgkb_module, "head_request_headers", lambda url: None)
        monkeypatch.setattr(
            pharmgkb_module, "fetch_cpic_remote_signal", lambda: "lastchange:2026-05-11"
        )
        assert annotator.fetch_remote_signal() is None

    def test_returns_unavailable_when_cpic_probe_fails(
        self, annotator: PharmGKBAnnotator, monkeypatch
    ):
        """R-5: CPIC probe failure is non-fatal — signal carries cpic:unavailable."""
        from allelix.annotators import pharmgkb as pharmgkb_module

        monkeypatch.setattr(
            pharmgkb_module,
            "head_request_headers",
            lambda url: {"ETag": '"abc123"'},
        )
        monkeypatch.setattr(pharmgkb_module, "fetch_cpic_remote_signal", lambda: None)
        signal = annotator.fetch_remote_signal()
        assert signal is not None
        assert "cpic:unavailable" in signal
        assert "pgkb:etag:" in signal

    def test_cached_round_trip(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        load_pharmgkb_tsv(
            mock_pharmgkb_dir,
            tmp_path / "pharmgkb.sqlite",
            source_url="test",
            remote_signal='etag:"v1"',
        )
        ann = PharmGKBAnnotator(tmp_path)
        try:
            assert ann.cached_remote_signal() == 'etag:"v1"'
        finally:
            ann.close()


class TestCloseable:
    def test_close_releases_connection(self, annotator: PharmGKBAnnotator):
        annotator.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
        assert annotator._conn is not None
        annotator.close()
        assert annotator._conn is None

    def test_close_is_idempotent(self, annotator: PharmGKBAnnotator):
        annotator.close()
        annotator.close()

    def test_context_manager(self, pharmgkb_data_dir: Path):
        ann = PharmGKBAnnotator(pharmgkb_data_dir)
        with ann as bound:
            assert bound is ann
            ann.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
        assert ann._conn is None


class TestCpicFallback:
    """R-5: ClinPGx setup succeeds when CPIC API is down."""

    def test_setup_succeeds_when_cpic_down(self, tmp_path: Path, mock_pharmgkb_dir: Path):
        ann = PharmGKBAnnotator(tmp_path)
        signal = "pgkb:etag:test123|cpic:unavailable"
        with (
            patch.object(ann, "fetch_remote_signal", return_value=signal),
            patch(
                "allelix.annotators.pharmgkb.download",
                side_effect=lambda url, dest: _fake_download(mock_pharmgkb_dir, url, dest),
            ),
            patch(
                "allelix.annotators.pharmgkb.fetch_cpic_allele_functions",
                side_effect=urllib.error.URLError("CPIC is down"),
            ),
        ):
            ann.setup()
        assert ann.is_ready()
        assert ann.version() is not None

    def test_signal_cpic_down_returns_unavailable(self):
        ann = PharmGKBAnnotator.__new__(PharmGKBAnnotator)
        with (
            patch(
                "allelix.annotators.pharmgkb.head_request_headers",
                return_value={"ETag": '"test-etag"'},
            ),
            patch(
                "allelix.annotators.pharmgkb.fetch_cpic_remote_signal",
                return_value=None,
            ),
        ):
            signal = ann.fetch_remote_signal()
        assert signal is not None
        assert "cpic:unavailable" in signal
        assert "pgkb:etag:" in signal

    def test_nonfinding_filter_degraded_without_cpic(
        self, tmp_path: Path, mock_pharmgkb_dir: Path
    ):
        ann = PharmGKBAnnotator(tmp_path)
        signal = "pgkb:etag:test123|cpic:unavailable"
        with (
            patch.object(ann, "fetch_remote_signal", return_value=signal),
            patch(
                "allelix.annotators.pharmgkb.download",
                side_effect=lambda url, dest: _fake_download(mock_pharmgkb_dir, url, dest),
            ),
            patch(
                "allelix.annotators.pharmgkb.fetch_cpic_allele_functions",
                side_effect=TimeoutError("CPIC timeout"),
            ),
        ):
            ann.setup()
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        with ann:
            results = ann.annotate(v)
        assert len(results) > 0


def _fake_download(mock_dir: Path, _url: str, dest: Path) -> None:
    """Create a zip from the mock fixture dir, simulating a real download."""
    import zipfile

    with zipfile.ZipFile(dest, "w") as zf:
        for tsv in mock_dir.glob("*.tsv"):
            zf.write(tsv, arcname=tsv.name)


class TestBatchAnnotateParity:
    """batch_annotate(vs) must return identical results to flatmap(annotate, vs).

    Pins the base.py contract: a chunked SQL override must produce the
    same annotations in the same order as the per-variant path.
    """

    def test_parity_mixed_carriers_and_non_carriers(self, annotator: PharmGKBAnnotator):
        variants = [
            Variant("rs1801133", "1", 11796321, "G", "A"),  # het carrier
            Variant("rs4680", "22", 19963748, "A", "A"),  # hom alt carrier
            Variant("rs1801133", "1", 11796321, "T", "T"),  # genotype mismatch
            Variant("rs999000111", "1", 1000, "A", "G"),  # unknown rsid
            Variant("rs1801133", "1", 11796321, "-", "-"),  # no-call
            Variant("rs113993960", "7", 117199644, "CTT", "C"),  # indel (rejected)
        ]
        per_variant = [a for v in variants for a in annotator.annotate(v)]
        batched = list(annotator.batch_annotate(variants))
        assert per_variant == batched

    def test_parity_empty_input(self, annotator: PharmGKBAnnotator):
        assert list(annotator.batch_annotate([])) == []

    def test_parity_clinvar_homref_filter_applied(self, pharmgkb_data_dir: Path):
        """ClinVar REF homref check (ADR-0023) must apply identically in batch."""
        refs = {"rs1801133": "G"}
        ann = PharmGKBAnnotator(
            pharmgkb_data_dir,
            clinvar_ref_provider=lambda rsid, _build: refs.get(rsid),
        )
        try:
            variants = [
                Variant("rs1801133", "1", 11796321, "G", "G"),  # homref, suppressed
                Variant("rs4680", "22", 19963748, "A", "A"),  # carrier, fires
            ]
            per_variant = [a for v in variants for a in ann.annotate(v)]
            batched = list(ann.batch_annotate(variants))
            assert per_variant == batched
        finally:
            ann.close()

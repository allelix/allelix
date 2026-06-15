# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the GWAS Catalog annotator."""

from __future__ import annotations

import contextlib
import shutil
import sqlite3
from pathlib import Path

import pytest

from allelix.annotators.gwas import (
    _EXCLUDED_TRAIT_CATEGORIES,
    _MUST_INCLUDE_RSIDS,
    _UNKNOWN_RISK_ALLELE_MAG_CAP,
    GWASCatalogAnnotator,
    _magnitude,
)
from allelix.databases.gwas_loader import (
    GWAS_DB_FILENAME,
    _is_metabolite_ratio,
    _is_uncharacterized_analyte,
    _parse_risk_allele,
    classify_gwas_trait,
    iter_gwas_records,
    load_gwas_tsv,
    schema_is_current,
)
from allelix.models import Annotation, Variant
from allelix.reports._pipeline import AnalysisResult

FIXTURES_DIR = Path(__file__).parents[1] / "fixtures"


class TestLoaderHelpers:
    """Unit tests for loader helper functions."""

    def test_parse_risk_allele_standard(self) -> None:
        assert _parse_risk_allele("rs1801133-A") == "A"

    def test_parse_risk_allele_question_mark(self) -> None:
        assert _parse_risk_allele("rs1801133-?") is None

    def test_parse_risk_allele_empty(self) -> None:
        assert _parse_risk_allele("") is None

    def test_parse_risk_allele_no_dash(self) -> None:
        assert _parse_risk_allele("rs1801133") is None

    def test_parse_risk_allele_multi_base(self) -> None:
        assert _parse_risk_allele("rs123-AT") is None


class TestLoaderDeduplication:
    """The loader keeps only the best p-value per (rsid, trait) pair."""

    def test_dedup_keeps_lowest_pvalue(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        mthfr_homocysteine = [
            r for r in records if r["rsid"] == "rs1801133" and r["trait"] == "Homocysteine levels"
        ]
        assert len(mthfr_homocysteine) == 1
        assert mthfr_homocysteine[0]["p_value"] == 2e-15


class TestLoaderIngestion:
    """SQLite ingestion from the mock TSV."""

    def test_load_creates_db(self, gwas_data_dir: Path) -> None:
        assert (gwas_data_dir / GWAS_DB_FILENAME).exists()

    def test_schema_is_current(self, gwas_data_dir: Path) -> None:
        assert schema_is_current(gwas_data_dir / GWAS_DB_FILENAME)

    def test_record_count(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir)
        try:
            assert annotator.record_count() == 8
        finally:
            annotator.close()


class TestMagnitudeScoring:
    """P-value and effect size magnitude mapping."""

    def test_genome_wide_significant(self) -> None:
        assert _magnitude(1e-10, None) == 6.0

    def test_strong_gwas_signal(self) -> None:
        assert _magnitude(1e-25, None) == 7.0

    def test_hyper_significant(self) -> None:
        assert _magnitude(1e-150, None) == 8.0

    def test_suggestive(self) -> None:
        assert _magnitude(1e-6, None) == 4.0

    def test_nominal(self) -> None:
        assert _magnitude(1e-4, None) == 3.0

    def test_weak(self) -> None:
        assert _magnitude(0.01, None) == 2.0

    def test_none_pvalue(self) -> None:
        assert _magnitude(None, None) == 2.0

    def test_high_or_adds_one(self) -> None:
        assert _magnitude(1e-10, 3.5) == 7.0

    def test_moderate_or_adds_half(self) -> None:
        assert _magnitude(1e-10, 2.5) == 6.5

    def test_or_capped_at_nine(self) -> None:
        assert _magnitude(1e-150, 5.0) == 9.0

    def test_protective_or_adds_one(self) -> None:
        assert _magnitude(1e-10, 0.2) == 7.0

    def test_protective_moderate_adds_half(self) -> None:
        assert _magnitude(1e-10, 0.4) == 6.5


class TestSetupAndStatus:
    """Annotator lifecycle: ready, version, close."""

    def test_unconfigured_is_not_ready(self, tmp_path: Path) -> None:
        annotator = GWASCatalogAnnotator(tmp_path)
        assert not annotator.is_ready()

    def test_configured_is_ready(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir)
        assert annotator.is_ready()

    def test_version_returns_string(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir)
        assert annotator.version() is not None


class TestSignalGuard:
    def test_setup_aborts_when_signal_fetch_fails(self, tmp_path: Path) -> None:
        """setup() raises RuntimeError when remote signal is None."""
        ann = GWASCatalogAnnotator(tmp_path)
        ann.fetch_remote_signal = lambda: None  # type: ignore[assignment]
        with pytest.raises(RuntimeError, match="cannot verify remote freshness signal"):
            ann.setup()


class TestAutoReingest:
    """Auto-reingest from cached TSV when categorizer version bumps."""

    def test_is_ready_auto_reingests_when_categorizer_bumped(
        self, tmp_path: Path, mock_gwas_tsv: Path
    ) -> None:
        """Stale cv: tag + TSV present → auto-reingest succeeds."""
        db_path = tmp_path / GWAS_DB_FILENAME
        tsv_dest = tmp_path / "gwas_catalog_associations.tsv"
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock")
        shutil.copy2(mock_gwas_tsv, tsv_dest)
        assert schema_is_current(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'cv:0' WHERE name = 'gwas'"
            )
            conn.commit()
        assert not schema_is_current(db_path)
        annotator = GWASCatalogAnnotator(tmp_path)
        assert annotator.is_ready() is True
        assert schema_is_current(db_path)

    def test_is_ready_returns_false_when_tsv_missing(
        self, tmp_path: Path, mock_gwas_tsv: Path
    ) -> None:
        """Stale tag + no TSV → can't auto-reingest, returns False."""
        db_path = tmp_path / GWAS_DB_FILENAME
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'cv:0' WHERE name = 'gwas'"
            )
            conn.commit()
        annotator = GWASCatalogAnnotator(tmp_path)
        assert annotator.is_ready() is False


class TestLegacyMigration:
    """One-shot migration from |cv: in remote_signal to local_version_tag."""

    def test_legacy_cv_in_remote_signal_migrates(
        self, tmp_path: Path, mock_gwas_tsv: Path
    ) -> None:
        """Pre-v1.5.0 cache with |cv:N in remote_signal self-heals."""
        from allelix.databases.gwas_loader import _CATEGORIZER_VERSION

        db_path = tmp_path / GWAS_DB_FILENAME
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock", remote_signal="etag:abc")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = NULL "
                "WHERE name = 'gwas'",
                (f"etag:abc|cv:{_CATEGORIZER_VERSION}",),
            )
            conn.commit()
        annotator = GWASCatalogAnnotator(tmp_path)
        assert annotator.is_ready() is True
        from allelix.databases.manager import get_database_info

        info = get_database_info(db_path, "gwas")
        assert info is not None
        assert info["remote_signal"] == "etag:abc"
        assert info["local_version_tag"] == f"cv:{_CATEGORIZER_VERSION}"

    def test_cached_remote_signal_returns_clean(self, tmp_path: Path, mock_gwas_tsv: Path) -> None:
        """After migration, cached_remote_signal returns clean signal."""
        db_path = tmp_path / GWAS_DB_FILENAME
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock", remote_signal="etag:xyz")
        annotator = GWASCatalogAnnotator(tmp_path)
        assert annotator.cached_remote_signal() == "etag:xyz"

    def test_stamp_missing_file_returns_false(self, tmp_path: Path) -> None:
        from allelix.annotators.gwas import _stamp_existing_gwas_cache

        assert _stamp_existing_gwas_cache(tmp_path / "nonexistent.sqlite") is False

    def test_stamp_no_row_returns_false(self, tmp_path: Path) -> None:
        from allelix.annotators.gwas import _stamp_existing_gwas_cache
        from allelix.databases.schema import GWAS_SCHEMA

        db = tmp_path / "empty.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.executescript(GWAS_SCHEMA)
        assert _stamp_existing_gwas_cache(db) is False

    def test_stamp_already_current_returns_true(self, tmp_path: Path, mock_gwas_tsv: Path) -> None:
        from allelix.annotators.gwas import _stamp_existing_gwas_cache

        db_path = tmp_path / GWAS_DB_FILENAME
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock")
        assert _stamp_existing_gwas_cache(db_path) is True

    def test_stamp_stale_tag_returns_false(self, tmp_path: Path, mock_gwas_tsv: Path) -> None:
        from allelix.annotators.gwas import _stamp_existing_gwas_cache

        db_path = tmp_path / GWAS_DB_FILENAME
        load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'cv:0' WHERE name = 'gwas'"
            )
            conn.commit()
        assert _stamp_existing_gwas_cache(db_path) is False

    def test_has_current_columns_on_bad_db(self, tmp_path: Path) -> None:
        from allelix.annotators.gwas import _has_current_gwas_columns

        bad = tmp_path / "bad.sqlite"
        bad.write_text("not a database", encoding="utf-8")
        assert _has_current_gwas_columns(bad) is False


class TestGenotypeMatching:
    """Carrier rule: only fire when the user carries the risk allele."""

    def test_het_carrier_fires(self, gwas_data_dir: Path) -> None:
        """rs1801133 G/A — risk allele A, should fire for Homocysteine levels."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Homocysteine levels" in traits
        finally:
            annotator.close()

    def test_hom_alt_fires(self, gwas_data_dir: Path) -> None:
        """rs4680 A/A — risk allele A, should fire."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs4680",
                chromosome="22",
                position=19963748,
                allele1="A",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Pain sensitivity" in traits
        finally:
            annotator.close()

    def test_hom_ref_does_not_fire(self, gwas_data_dir: Path) -> None:
        """rs121918506 G/G — risk allele T, user doesn't carry it."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs121918506",
                chromosome="17",
                position=7674222,
                allele1="G",
                allele2="G",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert len(results) == 0
        finally:
            annotator.close()

    def test_no_call_does_not_fire(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="-",
                allele2="-",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert len(results) == 0
        finally:
            annotator.close()

    def test_unknown_rsid_returns_empty(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs999999999",
                chromosome="1",
                position=1,
                allele1="A",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert len(results) == 0
        finally:
            annotator.close()

    def test_unknown_risk_allele_fires_with_cap(self, gwas_data_dir: Path) -> None:
        """rs1801133 has a second GWAS entry with risk allele '?' (Height).

        ADR-0024: unknown-risk-allele magnitude capped at 3.0 so it
        doesn't pass typical --min-magnitude thresholds.
        """
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            height_hits = [a for a in results if a.condition == "Height"]
            assert len(height_hits) == 1
            assert "risk allele not specified" in height_hits[0].description
            assert height_hits[0].magnitude <= _UNKNOWN_RISK_ALLELE_MAG_CAP
            homocysteine_hits = [a for a in results if a.condition == "Homocysteine levels"]
            assert len(homocysteine_hits) == 1
            assert height_hits[0].magnitude < homocysteine_hits[0].magnitude
        finally:
            annotator.close()

    def test_unknown_risk_allele_cap_value_is_3_0(self) -> None:
        """ADR-0024 pin: the magnitude cap for unknown risk allele is 3.0."""
        assert _UNKNOWN_RISK_ALLELE_MAG_CAP == 3.0

    def test_unknown_risk_allele_caps_high_evidence(self) -> None:
        """ADR-0024: even genome-wide significant p-value is capped at 3.0
        when the risk allele is unknown.
        """
        base_mag = _magnitude(1e-15, None)
        assert base_mag > _UNKNOWN_RISK_ALLELE_MAG_CAP
        capped = min(base_mag, _UNKNOWN_RISK_ALLELE_MAG_CAP)
        assert capped == _UNKNOWN_RISK_ALLELE_MAG_CAP


class TestAttribution:
    """Source, category, significance labels."""

    def test_source_is_gwas(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert all(a.source == "gwas" for a in results)
        finally:
            annotator.close()

    def test_category_is_trait(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert all(a.category == "trait" for a in results)
        finally:
            annotator.close()

    def test_significance_label(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert all(a.significance == "gwas_association" for a in results)
        finally:
            annotator.close()

    def test_attribution_label(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert all(a.attribution == "GWAS Catalog" for a in results)
        finally:
            annotator.close()

    def test_references_include_pubmed(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            homocysteine = next(a for a in results if a.condition == "Homocysteine levels")
            assert any(r.startswith("pubmed:") for r in homocysteine.references)
            assert any(r.startswith("gwas:") for r in homocysteine.references)
        finally:
            annotator.close()


class TestCloseable:
    """Context manager and close semantics."""

    def test_close_is_idempotent(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir)
        annotator.close()
        annotator.close()

    def test_context_manager(self, gwas_data_dir: Path) -> None:
        with GWASCatalogAnnotator(gwas_data_dir) as annotator:
            assert annotator.is_ready()


class TestTraitClassifier:
    """Unit tests for classify_gwas_trait()."""

    def test_mondo_uri_is_disease(self) -> None:
        uri = "http://purl.obolibrary.org/obo/MONDO_0005319"
        assert classify_gwas_trait("myopathy", uri) == "disease"

    def test_oba_uri_is_measurement(self) -> None:
        uri = "http://purl.obolibrary.org/obo/OBA_0000060"
        assert classify_gwas_trait("body height", uri) == "other_measurement"

    def test_body_height_is_body_measurement(self) -> None:
        assert classify_gwas_trait("body height", "") == "body_measurement"

    def test_body_mass_index_is_body_measurement(self) -> None:
        assert classify_gwas_trait("body mass index", "") == "body_measurement"

    def test_cholesterol_is_lipid(self) -> None:
        assert classify_gwas_trait("HDL cholesterol measurement", "") == "lipid_measurement"

    def test_triglyceride_is_lipid(self) -> None:
        assert classify_gwas_trait("triglyceride measurement", "") == "lipid_measurement"

    def test_platelet_count_is_hematological(self) -> None:
        assert classify_gwas_trait("platelet count", "") == "hematological_measurement"

    def test_red_blood_cell_is_hematological(self) -> None:
        assert classify_gwas_trait("red blood cell count", "") == "hematological_measurement"

    def test_blood_pressure_is_other_measurement(self) -> None:
        assert classify_gwas_trait("blood pressure", "") == "other_measurement"

    def test_homocysteine_is_other_measurement(self) -> None:
        assert classify_gwas_trait("homocysteine measurement", "") == "other_measurement"

    def test_educational_attainment_is_behavioral(self) -> None:
        assert classify_gwas_trait("educational attainment", "") == "behavioral"

    def test_breast_cancer_is_cancer(self) -> None:
        assert classify_gwas_trait("breast carcinoma", "") == "cancer"

    def test_neoplasm_is_cancer(self) -> None:
        assert classify_gwas_trait("neoplasm", "") == "cancer"

    def test_warfarin_is_drug_response(self) -> None:
        assert classify_gwas_trait("warfarin dose", "") == "drug_response"

    def test_rheumatoid_arthritis_is_immune(self) -> None:
        assert classify_gwas_trait("rheumatoid arthritis", "") == "immune"

    def test_coronary_artery_disease_is_cardiovascular(self) -> None:
        assert classify_gwas_trait("coronary artery disease", "") == "cardiovascular"

    def test_type_2_diabetes_is_metabolic(self) -> None:
        assert classify_gwas_trait("type 2 diabetes mellitus", "") == "metabolic"

    def test_alzheimer_is_neurological(self) -> None:
        assert classify_gwas_trait("Alzheimer's disease", "") == "neurological"

    def test_generic_disease_catchall(self) -> None:
        assert classify_gwas_trait("myopathy", "") == "disease"

    def test_measurement_suffix_catchall(self) -> None:
        assert classify_gwas_trait("some unknown measurement", "") == "other_measurement"

    def test_unknown_trait_is_other(self) -> None:
        assert classify_gwas_trait("something totally unknown", "") == "other"

    def test_disease_priority_over_measurement(self) -> None:
        """Multi-trait with both disease and measurement keywords."""
        assert classify_gwas_trait("breast cancer, body height", "") == "cancer"

    def test_mondo_priority_over_keywords(self) -> None:
        """MONDO URI trumps keyword-based classification."""
        result = classify_gwas_trait("body height", "http://purl.obolibrary.org/obo/MONDO_0000001")
        assert result == "disease"

    def test_empty_trait_is_other(self) -> None:
        assert classify_gwas_trait("", "") == "other"


class TestTraitFiltering:
    """ADR-0024: default trait-category filtering excludes noise categories."""

    def test_excluded_categories_frozen(self) -> None:
        assert isinstance(_EXCLUDED_TRAIT_CATEGORIES, frozenset)
        assert "body_measurement" in _EXCLUDED_TRAIT_CATEGORIES
        assert "lipid_measurement" in _EXCLUDED_TRAIT_CATEGORIES
        assert "hematological_measurement" in _EXCLUDED_TRAIT_CATEGORIES
        assert "other_measurement" in _EXCLUDED_TRAIT_CATEGORIES
        assert "behavioral" in _EXCLUDED_TRAIT_CATEGORIES

    def test_disease_categories_not_excluded(self) -> None:
        for cat in (
            "disease",
            "cancer",
            "drug_response",
            "immune",
            "cardiovascular",
            "metabolic",
            "neurological",
            "other",
        ):
            assert cat not in _EXCLUDED_TRAIT_CATEGORIES

    def test_filter_excludes_measurement_traits(self, gwas_data_dir: Path) -> None:
        """Default filtering drops homocysteine measurement and body height."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Homocysteine levels" not in traits
            assert "Height" not in traits
        finally:
            annotator.close()

    def test_filter_keeps_disease_traits(self, gwas_data_dir: Path) -> None:
        """Default filtering keeps cancer and drug response traits."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            v = Variant(
                rsid="rs80357906",
                chromosome="17",
                position=43057063,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Breast cancer" in traits
        finally:
            annotator.close()

    def test_filter_keeps_drug_response(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            v = Variant(
                rsid="rs1799853",
                chromosome="10",
                position=94942290,
                allele1="C",
                allele2="T",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Warfarin dose requirement" in traits
        finally:
            annotator.close()

    def test_filter_keeps_mondo_disease(self, gwas_data_dir: Path) -> None:
        """Statin-induced myopathy has MONDO_ URI -> classified as disease."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            v = Variant(
                rsid="rs4149056",
                chromosome="12",
                position=21178615,
                allele1="T",
                allele2="C",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Statin-induced myopathy" in traits
        finally:
            annotator.close()

    def test_unfiltered_returns_all(self, gwas_data_dir: Path) -> None:
        """filter_traits=False disables trait-category filtering."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Homocysteine levels" in traits
            assert "Height" in traits
        finally:
            annotator.close()


class TestLoaderTraitCategory:
    """Loader populates trait_category from MAPPED_TRAIT and MAPPED_TRAIT_URI."""

    def test_records_have_trait_category(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        assert all("trait_category" in r for r in records)

    def test_records_have_mapped_trait_uri(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        assert all("mapped_trait_uri" in r for r in records)

    def test_breast_carcinoma_classified_as_cancer(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        breast = [r for r in records if r["trait"] == "Breast cancer"]
        assert len(breast) == 1
        assert breast[0]["trait_category"] == "cancer"

    def test_homocysteine_classified_as_measurement(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        hc = [r for r in records if r["trait"] == "Homocysteine levels"]
        assert len(hc) == 1
        assert hc[0]["trait_category"] == "other_measurement"

    def test_body_height_classified_as_body_measurement(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        height = [r for r in records if r["trait"] == "Height"]
        assert len(height) == 1
        assert height[0]["trait_category"] == "body_measurement"

    def test_myopathy_mondo_classified_as_disease(self, mock_gwas_tsv: Path) -> None:
        records = list(iter_gwas_records(mock_gwas_tsv))
        myop = [r for r in records if r["trait"] == "Statin-induced myopathy"]
        assert len(myop) == 1
        assert myop[0]["trait_category"] == "disease"


class TestStructuralNoiseDetection:
    """ADR-0024 step 1.5: metabolite ratios, uncharacterized analytes, UKB body-comp."""

    @pytest.mark.parametrize(
        ("trait", "expected"),
        [
            ("whole body water mass", "body_measurement"),
            ("body water mass", "body_measurement"),
            ("impedance of whole body", "body_measurement"),
            ("impedance of arm (left)", "body_measurement"),
            ("impedance of leg (right)", "body_measurement"),
            ("impedance of trunk", "body_measurement"),
            ("UKB data field 23100 (whole body water mass)", "body_measurement"),
        ],
        ids=[
            "whole-body-water",
            "body-water-mass",
            "impedance-whole-body",
            "impedance-arm",
            "impedance-leg",
            "impedance-trunk",
            "ukb-data-field-231",
        ],
    )
    def test_ukb_body_composition(self, trait: str, expected: str) -> None:
        assert classify_gwas_trait(trait, "") == expected

    @pytest.mark.parametrize(
        "trait",
        [
            "cholesterol-to-phospholipid ratio in small VLDL",
            "triglyceride-to-phospholipid ratio in large HDL",
            "free cholesterol-to-esterified cholesterol ratio",
        ],
        ids=["cholesterol-phospholipid", "triglyceride-phospholipid", "free-esterified"],
    )
    def test_metabolite_ratio_classified_as_measurement(self, trait: str) -> None:
        assert classify_gwas_trait(trait, "") == "other_measurement"

    @pytest.mark.parametrize(
        "trait",
        [
            "X-12345 level",
            "X-09876 level in plasma",
        ],
        ids=["x-analyte-basic", "x-analyte-in-plasma"],
    )
    def test_uncharacterized_analyte_classified_as_measurement(self, trait: str) -> None:
        assert classify_gwas_trait(trait, "") == "other_measurement"

    def test_metabolite_ratio_helper_positive(self) -> None:
        assert _is_metabolite_ratio("cholesterol-to-phospholipid ratio in small vldl")

    def test_metabolite_ratio_helper_negative(self) -> None:
        assert not _is_metabolite_ratio("cholesterol measurement")

    def test_uncharacterized_analyte_helper_positive(self) -> None:
        assert _is_uncharacterized_analyte("x-12345 level")

    def test_uncharacterized_analyte_helper_negative(self) -> None:
        assert not _is_uncharacterized_analyte("albumin level")

    def test_disease_not_misrouted_by_ratio(self) -> None:
        """A disease trait containing '-to-' must not be caught by the ratio filter."""
        assert classify_gwas_trait("response to metformin", "") == "drug_response"

    def test_disease_not_misrouted_by_body_kw(self) -> None:
        """'body' in a disease context must not trigger body_measurement."""
        uri = "http://purl.obolibrary.org/obo/MONDO_0005015"
        assert classify_gwas_trait("whole body irradiation response", uri) == "disease"

    def test_alzheimer_not_misrouted(self) -> None:
        assert classify_gwas_trait("Alzheimer's disease", "") == "neurological"


class TestMustInclude:
    """ADR-0024: must-include rsID allowlist bypasses per-source floor."""

    def test_must_include_constant(self) -> None:
        assert isinstance(_MUST_INCLUDE_RSIDS, frozenset)
        assert "rs10737680" in _MUST_INCLUDE_RSIDS
        assert "rs11209026" in _MUST_INCLUDE_RSIDS
        assert "rs9271366" in _MUST_INCLUDE_RSIDS

    def test_must_include_flag_set_on_carrier(self, gwas_data_dir: Path) -> None:
        """Must-include rsIDs get is_must_include=True when carried."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            assert len(results) > 0
            assert all(not a.is_must_include for a in results)
        finally:
            annotator.close()

    def test_source_floor_bypass(self) -> None:
        """Must-include annotations bypass per-source magnitude floor."""
        must_include_ann = Annotation(
            source="gwas",
            rsid="rs10737680",
            significance="gwas_association",
            category="trait",
            magnitude=6.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="AC",
            is_must_include=True,
        )
        normal_ann = Annotation(
            source="gwas",
            rsid="rs9999999",
            significance="gwas_association",
            category="trait",
            magnitude=6.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="AC",
            is_must_include=False,
        )
        result = AnalysisResult(
            file_path=Path("/dev/null"),
            parser_name="test",
            parser_display_name="Test",
            sample_id="TEST",
            build="GRCh38",
            total_variants=1,
            skipped_count=0,
            annotators_used=[],
            annotations=[must_include_ann, normal_ann],
        )
        filtered = result.filter(source_min_magnitudes={"gwas": 9.0})
        rsids = {a.rsid for a in filtered}
        assert "rs10737680" in rsids, "must-include should bypass source floor"
        assert "rs9999999" not in rsids, "normal annotation should be blocked by source floor"

    def test_global_min_magnitude_still_applies(self) -> None:
        """Global --min-magnitude still filters must-include annotations."""
        must_include_ann = Annotation(
            source="gwas",
            rsid="rs10737680",
            significance="gwas_association",
            category="trait",
            magnitude=4.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="AC",
            is_must_include=True,
        )
        result = AnalysisResult(
            file_path=Path("/dev/null"),
            parser_name="test",
            parser_display_name="Test",
            sample_id="TEST",
            build="GRCh38",
            total_variants=1,
            skipped_count=0,
            annotators_used=[],
            annotations=[must_include_ann],
        )
        filtered = result.filter(min_magnitude=5.0, source_min_magnitudes={"gwas": 9.0})
        assert len(filtered) == 0, "global min_magnitude should still apply to must-include"

    def test_trait_filter_still_applies(self, gwas_data_dir: Path) -> None:
        """Trait-category filter still applies to must-include rsIDs."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="G",
                allele2="A",
                build="GRCh38",
            )
            results = annotator.annotate(v)
            traits = {a.condition for a in results}
            assert "Homocysteine levels" not in traits
        finally:
            annotator.close()


class TestRegistryMetadata:
    """Class attributes match the annotator contract."""

    def test_name(self) -> None:
        assert GWASCatalogAnnotator.name == "gwas"

    def test_display_name(self) -> None:
        assert GWASCatalogAnnotator.display_name == "GWAS Catalog"

    def test_requires_download(self) -> None:
        assert GWASCatalogAnnotator.requires_download is True


class TestBatchAnnotateParity:
    """batch_annotate(vs) must return identical results to flatmap(annotate, vs).

    Pins the base.py contract: a chunked SQL override must produce the
    same annotations in the same order as the per-variant path.
    """

    def test_parity_mixed_carriers_and_non_carriers(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            variants = [
                Variant("rs1801133", "1", 11796321, "G", "A", build="GRCh38"),  # carrier
                Variant("rs4680", "22", 19963748, "A", "A", build="GRCh38"),  # carrier
                Variant("rs121918506", "17", 7674222, "G", "G", build="GRCh38"),  # non-carrier
                Variant("rs999000111", "1", 1000, "A", "T", build="GRCh38"),  # unknown rsid
                Variant("rs1801133", "1", 11796321, "-", "-", build="GRCh38"),  # no-call
            ]
            per_variant = [a for v in variants for a in annotator.annotate(v)]
            batched = list(annotator.batch_annotate(variants))
            assert per_variant == batched
        finally:
            annotator.close()

    def test_parity_empty_input(self, gwas_data_dir: Path) -> None:
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=False)
        try:
            assert list(annotator.batch_annotate([])) == []
        finally:
            annotator.close()

    def test_parity_with_trait_filter_enabled(self, gwas_data_dir: Path) -> None:
        """Trait-category filter must apply identically in batch path."""
        annotator = GWASCatalogAnnotator(gwas_data_dir, filter_traits=True)
        try:
            variants = [
                Variant("rs1801133", "1", 11796321, "G", "A", build="GRCh38"),
                Variant("rs4680", "22", 19963748, "A", "A", build="GRCh38"),
            ]
            per_variant = [a for v in variants for a in annotator.annotate(v)]
            batched = list(annotator.batch_annotate(variants))
            assert per_variant == batched
        finally:
            annotator.close()

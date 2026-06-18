# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared pytest fixtures."""

from __future__ import annotations

import contextlib
import sqlite3
from pathlib import Path

import pytest

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.gwas_loader import load_gwas_tsv
from allelix.databases.pharmgkb_loader import (
    FUNCTION_CLASS_DECREASED,
    FUNCTION_CLASS_NO_FUNCTION,
    FUNCTION_CLASS_NORMAL,
    load_pharmgkb_tsv,
)
from allelix.databases.schema import CLINVAR_SCHEMA

# Synthetic ClinVar cache rows used by clinvar_data_dir and
# all_annotators_data_dir. Pre-stage-C this content lived in
# tests/fixtures/mock_clinvar_grch{37,38}.vcf and was parsed by
# iter_clinvar_records. Stage C removes the VCF loader, so the
# data is encoded directly in Python and inserted into the cache
# via SQLite. The biology is identical to the prior VCFs — every
# test that depends on these records keeps its semantics.
#
# Format: (rsid, chromosome, position, ref, alt, clinical_significance,
#          condition, gene, review_status, allele_id)
_MOCK_CLINVAR_ROWS_GRCH37: tuple[tuple, ...] = (
    (
        "rs1801133",
        "1",
        11856378,
        "G",
        "A",
        "Pathogenic",
        "MTHFR_deficiency",
        "MTHFR",
        "criteria_provided,_single_submitter",
        100001,
    ),
    (
        "rs1801131",
        "1",
        11854476,
        "T",
        "G",
        "Likely_pathogenic",
        "Hyperhomocysteinemia",
        "MTHFR",
        "criteria_provided,_single_submitter",
        100002,
    ),
    (
        "rs1801394",
        "5",
        7870860,
        "A",
        "G",
        "Likely_benign",
        "Folate_metabolism_disorder",
        "MTRR",
        "criteria_provided,_single_submitter",
        100003,
    ),
    (
        "rs4680",
        "22",
        19951271,
        "G",
        "A",
        "Drug_response",
        "Methylphenidate_response",
        "COMT",
        "criteria_provided,_single_submitter",
        100004,
    ),
    (
        "rs1799853",
        "10",
        96702047,
        "C",
        "T",
        "Drug_response",
        "Warfarin_response",
        "CYP2C9",
        "criteria_provided,_single_submitter",
        100005,
    ),
    (
        "rs4149056",
        "12",
        21331549,
        "T",
        "C",
        "Drug_response",
        "Statin-induced_myopathy",
        "SLCO1B1",
        "criteria_provided,_single_submitter",
        100006,
    ),
    (
        "rs80357906",
        "17",
        41209080,
        "G",
        "A",
        "Pathogenic",
        "Hereditary_breast_and_ovarian_cancer_syndrome",
        "BRCA1",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100007,
    ),
    (
        "rs121918506",
        "17",
        7577538,
        "G",
        "T",
        "Pathogenic",
        "Li-Fraumeni_syndrome",
        "TP53",
        "criteria_provided,_single_submitter",
        100008,
    ),
    (
        "rs999999999",
        "1",
        100,
        "A",
        "T",
        "Benign",
        "Synthetic_test_only",
        "TESTGENE",
        "no_assertion_provided",
        100009,
    ),
    (
        "rs113993960",
        "7",
        117199644,
        "CTT",
        "C",
        "Pathogenic",
        "Cystic_fibrosis",
        "CFTR",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100010,
    ),
    # ADR-0021 strand-inversion pin: GRCh37 NIPA1 REF=C ALT=G.
    (
        "rs104894490",
        "15",
        23060816,
        "C",
        "G",
        "Pathogenic",
        "Hereditary_spastic_paraplegia_6",
        "NIPA1",
        "criteria_provided,_single_submitter",
        100011,
    ),
    # rs1065852 multi-allelic split (G->A drug_response, G->C benign).
    (
        "rs1065852",
        "22",
        42526694,
        "G",
        "A",
        "Drug_response",
        "Codeine_response",
        "CYP2D6",
        "criteria_provided,_single_submitter",
        100020,
    ),
    (
        "rs1065852",
        "22",
        42526694,
        "G",
        "C",
        "Benign",
        "Synthetic_benign_pair",
        "CYP2D6",
        "criteria_provided,_single_submitter",
        100021,
    ),
    # GH #111 fixtures: multi-SCV variants needed to pin #42's read-side
    # behavior. Pre-#42 each variant was one row; post-#42 (per-SCV TSV
    # loader) each SCV submission gets its own row, so these mirror the
    # real cache shape and exercise the per-SCV semantics fast-tier
    # tests previously couldn't reach (root-cause coverage gap behind
    # #106 and #109).
    #
    # rs999000111 — three identical SCV rows for the same (variant,
    # sig, condition, gene) quadruple. On today's dev tip the
    # annotator emits three identical Annotation objects (the #109
    # regression — visible in terminal/HTML/JSON). The #109 fix
    # collapses to one Annotation while preserving the strongest
    # review_status and the union of references.
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100030,
    ),
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_single_submitter",
        100031,
    ),
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_single_submitter",
        100032,
    ),
    # rs999000222 — two SCV rows for the same variant with distinct
    # (significance, condition) pairs. #109's dedup MUST preserve
    # both rows (this is exactly what #42 was built to surface;
    # collapsing them would be a Frankenstein-pair regression).
    (
        "rs999000222",
        "1",
        200100,
        "C",
        "G",
        "Pathogenic",
        "Disease_A_test",
        "MULTI_SCV_CONFLICT",
        "criteria_provided,_single_submitter",
        100040,
    ),
    (
        "rs999000222",
        "1",
        200100,
        "C",
        "G",
        "Likely_benign",
        "Disease_B_test",
        "MULTI_SCV_CONFLICT",
        "criteria_provided,_single_submitter",
        100041,
    ),
    # rs999000333 — sub-floor: the only annotation maps to magnitude
    # 2.0 (Likely_benign), below the 5.0 analyze floor. On today's
    # dev tip this rsid lands in panel_coverage.found and vanishes
    # from every rendered surface — the #106 accounting lie. The
    # #106 patch reroutes it to no_findings.
    (
        "rs999000333",
        "1",
        200200,
        "T",
        "C",
        "Likely_benign",
        "Sub_floor_test",
        "SUB_FLOOR",
        "criteria_provided,_single_submitter",
        100050,
    ),
    # rs999000444 — unmapped CLNSIG long-tail. "protective" is a
    # real ClinVar term but is NOT in _CLNSIG_MAGNITUDE → falls to
    # the 5.0 default. Pinned here for #108 (v2.3) to retire the
    # silent 5.0-default surface; not actionable in v2.2.1.
    (
        "rs999000444",
        "1",
        200300,
        "G",
        "A",
        "protective",
        "Unmapped_clnsig_test",
        "UNMAPPED_CLNSIG",
        "criteria_provided,_single_submitter",
        100060,
    ),
    # rs999000555 — multi-allelic site with two SCV rows that share
    # (significance, condition, gene) but differ ONLY on alt. Pinned
    # for #109 Finding 1: the dedup key MUST include alt or these
    # rows collapse and the surviving annotation mis-attaches its
    # `alt` field downstream — breaking the exact-(rsid, alt) lookup
    # gnomAD / AlphaMissense / CADD use for allele-specific
    # enrichment. Carrier biology: a het user (A/C) at this site
    # carries both alts and should receive two distinct annotations,
    # one per alt — the #18 wrong-allele safety case applied to dedup.
    (
        "rs999000555",
        "1",
        200400,
        "T",
        "A",
        "Pathogenic",
        "Multi_alt_test",
        "MULTI_ALT_DEDUP",
        "criteria_provided,_single_submitter",
        100070,
    ),
    (
        "rs999000555",
        "1",
        200400,
        "T",
        "C",
        "Pathogenic",
        "Multi_alt_test",
        "MULTI_ALT_DEDUP",
        "criteria_provided,_single_submitter",
        100071,
    ),
)

_MOCK_CLINVAR_ROWS_GRCH38: tuple[tuple, ...] = (
    (
        "rs1801133",
        "1",
        11796321,
        "G",
        "A",
        "Pathogenic",
        "MTHFR_deficiency",
        "MTHFR",
        "criteria_provided,_single_submitter",
        100001,
    ),
    (
        "rs1801131",
        "1",
        11794419,
        "T",
        "G",
        "Likely_pathogenic",
        "Hyperhomocysteinemia",
        "MTHFR",
        "criteria_provided,_single_submitter",
        100002,
    ),
    (
        "rs1801394",
        "5",
        7870973,
        "A",
        "G",
        "Likely_benign",
        "Folate_metabolism_disorder",
        "MTRR",
        "criteria_provided,_single_submitter",
        100003,
    ),
    (
        "rs4680",
        "22",
        19963748,
        "G",
        "A",
        "Drug_response",
        "Methylphenidate_response",
        "COMT",
        "criteria_provided,_single_submitter",
        100004,
    ),
    (
        "rs1799853",
        "10",
        94942290,
        "C",
        "T",
        "Drug_response",
        "Warfarin_response",
        "CYP2C9",
        "criteria_provided,_single_submitter",
        100005,
    ),
    (
        "rs4149056",
        "12",
        21178615,
        "T",
        "C",
        "Drug_response",
        "Statin-induced_myopathy",
        "SLCO1B1",
        "criteria_provided,_single_submitter",
        100006,
    ),
    (
        "rs80357906",
        "17",
        43057063,
        "G",
        "A",
        "Pathogenic",
        "Hereditary_breast_and_ovarian_cancer_syndrome",
        "BRCA1",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100007,
    ),
    (
        "rs121918506",
        "17",
        7674222,
        "G",
        "T",
        "Pathogenic",
        "Li-Fraumeni_syndrome",
        "TP53",
        "criteria_provided,_single_submitter",
        100008,
    ),
    (
        "rs999999999",
        "1",
        100,
        "A",
        "T",
        "Benign",
        "Synthetic_test_only",
        "TESTGENE",
        "no_assertion_provided",
        100009,
    ),
    (
        "rs113993960",
        "7",
        117559590,
        "CTT",
        "C",
        "Pathogenic",
        "Cystic_fibrosis",
        "CFTR",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100010,
    ),
    # ADR-0021 strand-inversion pin: GRCh38 NIPA1 REF=G ALT=A.
    (
        "rs104894490",
        "15",
        22812251,
        "G",
        "A",
        "Pathogenic",
        "Hereditary_spastic_paraplegia_6",
        "NIPA1",
        "criteria_provided,_single_submitter",
        100011,
    ),
    # rs1065852 multi-allelic split.
    (
        "rs1065852",
        "22",
        42130692,
        "G",
        "A",
        "Drug_response",
        "Codeine_response",
        "CYP2D6",
        "criteria_provided,_single_submitter",
        100020,
    ),
    (
        "rs1065852",
        "22",
        42130692,
        "G",
        "C",
        "Benign",
        "Synthetic_benign_pair",
        "CYP2D6",
        "criteria_provided,_single_submitter",
        100021,
    ),
    # GH #111 fixtures (mirror of GRCh37 set above — see those rows
    # for the rationale per variant). Same synthetic chr1 coordinates
    # for cross-build symmetry; the per-SCV semantics being pinned
    # are build-independent.
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_multiple_submitters,_no_conflicts",
        100030,
    ),
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_single_submitter",
        100031,
    ),
    (
        "rs999000111",
        "1",
        200000,
        "A",
        "T",
        "Pathogenic",
        "Cystic_fibrosis_test",
        "MULTI_SCV_AGREE",
        "criteria_provided,_single_submitter",
        100032,
    ),
    (
        "rs999000222",
        "1",
        200100,
        "C",
        "G",
        "Pathogenic",
        "Disease_A_test",
        "MULTI_SCV_CONFLICT",
        "criteria_provided,_single_submitter",
        100040,
    ),
    (
        "rs999000222",
        "1",
        200100,
        "C",
        "G",
        "Likely_benign",
        "Disease_B_test",
        "MULTI_SCV_CONFLICT",
        "criteria_provided,_single_submitter",
        100041,
    ),
    (
        "rs999000333",
        "1",
        200200,
        "T",
        "C",
        "Likely_benign",
        "Sub_floor_test",
        "SUB_FLOOR",
        "criteria_provided,_single_submitter",
        100050,
    ),
    (
        "rs999000444",
        "1",
        200300,
        "G",
        "A",
        "protective",
        "Unmapped_clnsig_test",
        "UNMAPPED_CLNSIG",
        "criteria_provided,_single_submitter",
        100060,
    ),
    (
        "rs999000555",
        "1",
        200400,
        "T",
        "A",
        "Pathogenic",
        "Multi_alt_test",
        "MULTI_ALT_DEDUP",
        "criteria_provided,_single_submitter",
        100070,
    ),
    (
        "rs999000555",
        "1",
        200400,
        "T",
        "C",
        "Pathogenic",
        "Multi_alt_test",
        "MULTI_ALT_DEDUP",
        "criteria_provided,_single_submitter",
        100071,
    ),
)


def _build_synthetic_clinvar_cache(
    db_path: Path,
    build: str,
    *,
    source_url: str = "test://mock",
    remote_signal: str | None = None,
) -> None:
    """Build a per-build ClinVar SQLite cache from the synthetic rows above.

    Stage C (#42) removed load_clinvar_vcf; this helper replaces it for
    tests that just need a populated cache. The rows are identical to
    what the prior mock_clinvar_grch{37,38}.vcf fixtures parsed into.
    Stamps the current CLINVAR_INTERPRETER_VERSION so the annotator's
    is_ready() check passes without invoking stamp_existing_clinvar_cache.
    """
    if db_path.exists():
        db_path.unlink()
    rows = _MOCK_CLINVAR_ROWS_GRCH37 if build == "GRCh37" else _MOCK_CLINVAR_ROWS_GRCH38
    record_name = clinvar_record_name(build)
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, "
            "alt, clinical_significance, condition, gene, review_status, "
            "allele_id) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, "
            "downloaded_at, record_count, remote_signal, local_version_tag) "
            "VALUES (?, ?, '20260101', '2026-06-08', ?, ?, ?)",
            (
                record_name,
                source_url,
                len(rows),
                remote_signal,
                f"iv:{CLINVAR_INTERPRETER_VERSION}",
            ),
        )
        conn.commit()


@pytest.fixture
def build_synthetic_clinvar_cache():
    """Provide the synthetic ClinVar cache builder to tests outside conftest.

    Stage C (#42) replaced the VCF loader with a TSV loader; production
    tests that previously called `load_clinvar_vcf(...)` against a VCF
    to pre-populate a cache for unit tests now request this fixture and
    call it as `build_synthetic_clinvar_cache(db_path, build)`. Same
    biology (13 ClinVar rows per build), same per-build dispatch,
    no VCF parsing path involved.
    """
    return _build_synthetic_clinvar_cache


@pytest.fixture(autouse=True)
def _bypass_loader_row_floors(monkeypatch):
    """GH #19: production loaders enforce a row-count floor against truncated
    downloads (``GWAS_MIN_ROWS``, ``PHARMGKB_MIN_ROWS``). Mock fixtures only
    have a handful of rows, so every test that runs an annotator ``setup()``
    or auto-reingest path against mock data would trip the floor. Patch both
    constants to 0 across the test suite. Loader-level tests that need to
    exercise the floor explicitly pass ``min_rows=`` themselves (see
    ``tests/databases/test_*_loader.py::TestMinRowsFloor``).
    """
    from allelix.annotators import gwas as gwas_module
    from allelix.annotators import pharmgkb as pharmgkb_module

    monkeypatch.setattr(gwas_module, "GWAS_MIN_ROWS", 0)
    monkeypatch.setattr(pharmgkb_module, "PHARMGKB_MIN_ROWS", 0)


FIXTURES_DIR = Path(__file__).parent / "fixtures"

# ADR-0020: the structured per-allele function lookup the ClinPGx filter
# joins against. In production it's fetched from CPIC's API at db-update
# time; in tests we inject a deterministic dict so the filter is
# self-contained and offline. Entries here mirror real CPIC classifications
# for the rsids the MHG fixture carries, plus a few synthetic rsids that
# exercise the non-finding suppression path end-to-end.
MOCK_CPIC_LOOKUP: dict[tuple[str, str], str] = {
    # MTHFR C677T — both bases classified so the GG row stores as non-finding
    # and AG/AA emit. (CPIC doesn't actually publish MTHFR; this is a test
    # fixture choice — not a claim about real CPIC coverage.)
    ("rs1801133", "G"): FUNCTION_CLASS_NORMAL,
    ("rs1801133", "A"): FUNCTION_CLASS_DECREASED,
    # COMT (synthetic — see above note).
    ("rs4680", "G"): FUNCTION_CLASS_NORMAL,
    ("rs4680", "A"): FUNCTION_CLASS_DECREASED,
    # CYP2C9*2 (rs1799853).
    ("rs1799853", "C"): FUNCTION_CLASS_NORMAL,
    ("rs1799853", "T"): FUNCTION_CLASS_DECREASED,
    # SLCO1B1*5 (rs4149056).
    ("rs4149056", "T"): FUNCTION_CLASS_NORMAL,
    ("rs4149056", "C"): FUNCTION_CLASS_DECREASED,
    # PA-008 synthetic: G reference, A decreased. GG → non-finding.
    ("rs900000010", "G"): FUNCTION_CLASS_NORMAL,
    ("rs900000010", "A"): FUNCTION_CLASS_DECREASED,
    # Pins for the three v0.7.0/v0.8.0 production leakers (matches real CPIC).
    ("rs1800559", "C"): FUNCTION_CLASS_NORMAL,
    ("rs1800559", "T"): FUNCTION_CLASS_DECREASED,
    ("rs116855232", "C"): FUNCTION_CLASS_NORMAL,
    ("rs116855232", "T"): FUNCTION_CLASS_NO_FUNCTION,
    # DPYD rs1801265: CPIC assigns Normal function to BOTH alleles.
    # Regression pin: GG must be suppressed by the CPIC is_nonfinding flag
    # even when the user is not homozygous-reference per ClinVar.
    ("rs1801265", "G"): FUNCTION_CLASS_NORMAL,
    ("rs1801265", "A"): FUNCTION_CLASS_NORMAL,
}


@pytest.fixture
def mock_mhg_path() -> Path:
    """Path to the committed synthetic MyHappyGenes fixture (clean GRCh38).

    Generate it with `python tests/generate_mock_data.py` if missing.
    """
    path = FIXTURES_DIR / "mock_myhappygenes.txt"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}. Run: python tests/generate_mock_data.py")
    return path


@pytest.fixture
def mock_23andme_path() -> Path:
    """Path to the committed synthetic 23andMe fixture."""
    path = FIXTURES_DIR / "mock_23andme.txt"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_ancestrydna_path() -> Path:
    """Path to the committed synthetic AncestryDNA fixture."""
    path = FIXTURES_DIR / "mock_ancestrydna.txt"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_ftdna_path() -> Path:
    """Path to the committed synthetic FTDNA fixture."""
    path = FIXTURES_DIR / "mock_ftdna.csv"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_ftdna_illumina_path() -> Path:
    """Path to the committed synthetic FTDNA Illumina raw fixture."""
    path = FIXTURES_DIR / "mock_ftdna_illumina.txt"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_ftdna_famfinder_path() -> Path:
    """Path to the committed synthetic FTDNA FamFinder fixture."""
    path = FIXTURES_DIR / "mock_ftdna_famfinder.txt"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_myheritage_path() -> Path:
    """Path to the committed synthetic MyHeritage fixture."""
    path = FIXTURES_DIR / "mock_myheritage.csv"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_livingdna_path() -> Path:
    """Path to the committed synthetic Living DNA fixture."""
    path = FIXTURES_DIR / "mock_livingdna.csv"
    if not path.exists():
        pytest.fail(f"Mock fixture missing: {path}")
    return path


@pytest.fixture
def mock_mhg_grch37_path() -> Path:
    """ADR-0021 fixture: clean GRCh37 positions, GRCh37 header."""
    path = FIXTURES_DIR / "mock_myhappygenes_grch37.txt"
    if not path.exists():
        pytest.fail(
            f"GRCh37 mock fixture missing: {path}. Run: "
            "`python tests/generate_mock_data.py --build grch37 "
            "--output tests/fixtures/mock_myhappygenes_grch37.txt`"
        )
    return path


@pytest.fixture
def mock_mhg_mislabeled_path() -> Path:
    """ADR-0021 fixture: GRCh38 positions, header claims GRCh37 (real MHG bug)."""
    path = FIXTURES_DIR / "mock_myhappygenes_mislabeled.txt"
    if not path.exists():
        pytest.fail(
            f"Mislabeled mock fixture missing: {path}. Run: "
            "`python tests/generate_mock_data.py --build grch38 --header-build grch37 "
            "--output tests/fixtures/mock_myhappygenes_mislabeled.txt`"
        )
    return path


@pytest.fixture
def clinvar_data_dir(tmp_path: Path) -> Path:
    """Build a fresh data dir with populated per-build ClinVar caches.

    ADR-0021 + ADR-0015: GRCh37 and GRCh38 caches built from the
    synthetic row tables at the top of this module. NIPA1 (rs104894490)
    has REF=C ALT=G on GRCh37 and REF=G ALT=A on GRCh38 — tests that
    exercise the strand-inverted regression case observe DIFFERENT
    results across caches.
    """
    for build in ("GRCh37", "GRCh38"):
        _build_synthetic_clinvar_cache(
            tmp_path / clinvar_db_filename(build),
            build,
            source_url=f"test://mock-{build}",
        )
    return tmp_path


@pytest.fixture
def mock_pharmgkb_dir() -> Path:
    """Path to the synthetic ClinPGx clinical-annotations directory."""
    path = FIXTURES_DIR / "mock_pharmgkb"
    if not path.exists():
        pytest.fail(
            f"Mock ClinPGx fixture missing: {path}. Run: python tests/generate_pharmgkb_fixture.py"
        )
    return path


@pytest.fixture
def mock_cpic_lookup() -> dict[tuple[str, str], str]:
    """Synthetic CPIC per-allele function lookup for tests (ADR-0020)."""
    return dict(MOCK_CPIC_LOOKUP)


@pytest.fixture
def pharmgkb_data_dir(tmp_path: Path, mock_pharmgkb_dir: Path) -> Path:
    """Build a fresh data dir with a populated ClinPGx SQLite cache."""
    db_path = tmp_path / "pharmgkb.sqlite"
    load_pharmgkb_tsv(
        mock_pharmgkb_dir,
        db_path,
        source_url="test://mock-pharmgkb",
        allele_function_lookup=dict(MOCK_CPIC_LOOKUP),
    )
    return tmp_path


@pytest.fixture
def mock_gnomad_gz() -> Path:
    """Path to the gzipped mock gnomAD SQLite fixture."""
    path = FIXTURES_DIR / "mock_gnomad.sqlite.gz"
    if not path.exists():
        pytest.fail(
            f"Mock gnomAD fixture missing: {path}. Run: python tests/generate_mock_data.py"
        )
    return path


@pytest.fixture
def mock_gwas_tsv() -> Path:
    """Path to the synthetic GWAS Catalog associations TSV."""
    path = FIXTURES_DIR / "mock_gwas_catalog.tsv"
    if not path.exists():
        pytest.fail(f"Mock GWAS Catalog fixture missing: {path}.")
    return path


@pytest.fixture
def gwas_data_dir(tmp_path: Path, mock_gwas_tsv: Path) -> Path:
    """Build a fresh data dir with a populated GWAS Catalog SQLite cache."""
    db_path = tmp_path / "gwas.sqlite"
    load_gwas_tsv(mock_gwas_tsv, db_path, source_url="test://mock-gwas")
    return tmp_path


@pytest.fixture
def all_annotators_data_dir(
    tmp_path: Path,
    mock_pharmgkb_dir: Path,
    mock_gwas_tsv: Path,
) -> Path:
    """Build a fresh data dir with all annotators ready (ClinVar + ClinPGx + GWAS).

    Per-build ClinVar caches built from the synthetic row tables (ADR-0021).
    """
    for build in ("GRCh37", "GRCh38"):
        _build_synthetic_clinvar_cache(
            tmp_path / clinvar_db_filename(build),
            build,
            source_url=f"test://mock-{build}",
        )
    load_pharmgkb_tsv(
        mock_pharmgkb_dir,
        tmp_path / "pharmgkb.sqlite",
        source_url="test://mock-pharmgkb",
        allele_function_lookup=dict(MOCK_CPIC_LOOKUP),
    )
    load_gwas_tsv(
        mock_gwas_tsv,
        tmp_path / "gwas.sqlite",
        source_url="test://mock-gwas",
    )
    return tmp_path


@pytest.fixture
def cm_stack():
    """Per-test ``contextlib.ExitStack`` for context-managed resources.

    Tests that construct ``Annotator`` instances (or any other context-manager
    resource) should register them with this stack instead of binding the
    instance to a bare local. The stack ``__exit__`` fires at test teardown
    and propagates ``__exit__`` to each registered resource, guaranteeing
    deterministic cleanup. The pattern:

        def test_x(self, cm_stack, clinvar_data_dir):
            clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
            gnomad  = cm_stack.enter_context(GnomadAnnotator(clinvar_data_dir))
            result = run_analysis(...)

    Why this exists (GH allelix-dev #78, blocks #36): the pipeline's own
    ``ExitStack`` closes annotators it receives, but the test still holds the
    closed-instance reference. When ``Annotator.__del__`` is the safety net,
    GC of the test-local reference fires ``__del__`` → ``close()`` again
    (idempotent), no warning. With ``__del__`` removed, GC of the closed
    instance fires sqlite3's own ``Connection`` GC finalizer, which surfaces
    a ``ResourceWarning`` during the next test's startup that pytest's
    ``error::PytestUnraisableExceptionWarning`` filter converts to a hard
    failure. Wrapping construction in this stack puts cleanup back under
    deterministic, ``__del__``-independent control.
    """
    with contextlib.ExitStack() as stack:
        yield stack

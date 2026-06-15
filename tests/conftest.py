# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.databases.gwas_loader import load_gwas_tsv
from allelix.databases.manager import load_clinvar_vcf
from allelix.databases.pharmgkb_loader import (
    FUNCTION_CLASS_DECREASED,
    FUNCTION_CLASS_NO_FUNCTION,
    FUNCTION_CLASS_NORMAL,
    load_pharmgkb_tsv,
)

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
def mock_clinvar_grch37_vcf() -> Path:
    """ADR-0021 fixture: synthetic ClinVar VCF with GRCh37 positions."""
    path = FIXTURES_DIR / "mock_clinvar_grch37.vcf"
    if not path.exists():
        pytest.fail(
            f"Mock GRCh37 ClinVar VCF missing: {path}. Run: "
            "python tests/generate_clinvar_fixture.py"
        )
    return path


@pytest.fixture
def mock_clinvar_grch38_vcf() -> Path:
    """ADR-0021 fixture: synthetic ClinVar VCF with GRCh38 positions and
    build-specific REF/ALT for the NIPA1 rs104894490 strand-inverted case.
    """
    path = FIXTURES_DIR / "mock_clinvar_grch38.vcf"
    if not path.exists():
        pytest.fail(
            f"Mock GRCh38 ClinVar VCF missing: {path}. Run: "
            "python tests/generate_clinvar_fixture.py"
        )
    return path


# Back-compat: many tests reference `mock_clinvar_vcf` as a single path.
# They predate ADR-0021's dual-build split. Most query by rsID (position-
# agnostic), so a single fixture works for them. Point this at the GRCh37
# build to preserve historical semantics. New tests should use the
# build-specific fixtures above.
@pytest.fixture
def mock_clinvar_vcf(mock_clinvar_grch37_vcf: Path) -> Path:
    return mock_clinvar_grch37_vcf


@pytest.fixture
def clinvar_data_dir(
    tmp_path: Path,
    mock_clinvar_grch37_vcf: Path,
    mock_clinvar_grch38_vcf: Path,
) -> Path:
    """Build a fresh data dir with populated per-build ClinVar caches.

    ADR-0021 + ADR-0015: GRCh37 cache loaded from the GRCh37 fixture
    (REF=C ALT=G at NIPA1 etc.), GRCh38 cache from the GRCh38 fixture
    (REF=G ALT=A at NIPA1 etc.). The annotator dispatches by
    `variant.build`; tests that exercise the strand-inverted regression
    case can now observe DIFFERENT results across caches.
    """
    fixture_by_build = {
        "GRCh37": mock_clinvar_grch37_vcf,
        "GRCh38": mock_clinvar_grch38_vcf,
    }
    for build, vcf in fixture_by_build.items():
        load_clinvar_vcf(
            vcf,
            tmp_path / clinvar_db_filename(build),
            source_url=f"test://mock-{build}",
            record_name=clinvar_record_name(build),
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
    mock_clinvar_grch37_vcf: Path,
    mock_clinvar_grch38_vcf: Path,
    mock_pharmgkb_dir: Path,
    mock_gwas_tsv: Path,
) -> Path:
    """Build a fresh data dir with all annotators ready (ClinVar + ClinPGx + GWAS).

    Per-build ClinVar caches from the GRCh37 / GRCh38 fixtures (ADR-0021).
    """
    fixture_by_build = {
        "GRCh37": mock_clinvar_grch37_vcf,
        "GRCh38": mock_clinvar_grch38_vcf,
    }
    for build, vcf in fixture_by_build.items():
        load_clinvar_vcf(
            vcf,
            tmp_path / clinvar_db_filename(build),
            source_url=f"test://mock-{build}",
            record_name=clinvar_record_name(build),
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

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""SQLite schemas for cached reference databases.

Each annotator owns its own SQLite file (e.g. `clinvar.sqlite`, `pharmgkb.sqlite`).
Every per-annotator schema embeds the shared `database_versions` table so that
`get_database_info(db_path, name)` works uniformly across them.
"""

from __future__ import annotations

_DATABASE_VERSIONS_TABLE = """
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

CLINVAR_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS clinvar_variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT NOT NULL,
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    clinical_significance TEXT,
    condition TEXT,
    gene TEXT,
    review_status TEXT,
    allele_id INTEGER
);

CREATE INDEX IF NOT EXISTS idx_clinvar_rsid ON clinvar_variants(rsid);
CREATE INDEX IF NOT EXISTS idx_clinvar_position
    ON clinvar_variants(chromosome, position);
"""
    + _DATABASE_VERSIONS_TABLE
)

PHARMGKB_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS pharmgkb_annotations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT NOT NULL,
    genotype TEXT NOT NULL,
    gene TEXT,
    drugs TEXT,
    phenotype TEXT,
    phenotype_category TEXT,
    annotation_text TEXT,
    level_of_evidence TEXT,
    score REAL,
    pgkb_annotation_id TEXT,
    allele_function TEXT,
    function_class TEXT NOT NULL,
    is_nonfinding INTEGER NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_rsid ON pharmgkb_annotations(rsid);

-- ADR-0018: per-allele function extracted from ClinPGx's canonical CPIC
-- template sentence ("The {allele} allele of {rsid} is assigned {function}
-- function by CPIC."). Populated at load time by a pre-pass over the
-- annotation rows. Drives is_nonfinding classification for SNV rows where
-- the `Allele Function` column is empty (i.e., every in-scope row).
CREATE TABLE IF NOT EXISTS pharmgkb_allele_function (
    rsid TEXT NOT NULL,
    allele TEXT NOT NULL,
    function_class TEXT NOT NULL,
    source TEXT NOT NULL DEFAULT 'cpic_template',
    PRIMARY KEY (rsid, allele)
);
"""
    + _DATABASE_VERSIONS_TABLE
)

GWAS_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS gwas_associations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT NOT NULL,
    risk_allele TEXT,
    trait TEXT NOT NULL,
    p_value REAL,
    or_beta REAL,
    ci_text TEXT,
    gene TEXT,
    study_accession TEXT,
    pubmed_id TEXT,
    risk_allele_frequency REAL,
    context TEXT,
    mapped_trait_uri TEXT,
    trait_category TEXT
);

CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON gwas_associations(rsid);
"""
    + _DATABASE_VERSIONS_TABLE
)

GNOMAD_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS gnomad_frequencies (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    rsid TEXT,
    af REAL,
    af_popmax REAL,
    popmax TEXT,
    af_afr REAL,
    af_amr REAL,
    af_asj REAL,
    af_eas REAL,
    af_fin REAL,
    af_nfe REAL,
    af_sas REAL,
    PRIMARY KEY (chrom, pos, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_gnomad_rsid ON gnomad_frequencies(rsid);
"""
    + _DATABASE_VERSIONS_TABLE
)

ALPHAMISSENSE_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS alphamissense_scores (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    rsid TEXT,
    uniprot_id TEXT,
    transcript_id TEXT,
    protein_variant TEXT,
    am_pathogenicity REAL NOT NULL,
    am_class TEXT NOT NULL,
    PRIMARY KEY (chrom, pos, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_am_rsid ON alphamissense_scores(rsid);
"""
    + _DATABASE_VERSIONS_TABLE
)

CADD_SCHEMA = (
    """
CREATE TABLE IF NOT EXISTS cadd_scores (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    phred REAL NOT NULL,
    PRIMARY KEY (chrom, pos, ref, alt)
);
"""
    + _DATABASE_VERSIONS_TABLE
)

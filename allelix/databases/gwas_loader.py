# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GWAS Catalog association download, parse, and load into SQLite.

EBI/NHGRI publishes the full associations TSV at a stable URL. Public
domain -- no licensing restrictions. One row per study x SNP x trait.
"""

from __future__ import annotations

import contextlib
import csv
import logging
import os
import re
import sqlite3
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix.databases.schema import GWAS_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

GWAS_CATALOG_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip"
)
GWAS_DB_FILENAME = "gwas.sqlite"

INSERT_BATCH_SIZE = 5_000

_CATEGORIZER_VERSION = 3

_RSID_RE = re.compile(r"^rs\d+$")

_REQUIRED_GWAS_COLUMNS = frozenset(
    {
        "rsid",
        "risk_allele",
        "trait",
        "p_value",
        "or_beta",
        "ci_text",
        "gene",
        "study_accession",
        "pubmed_id",
        "risk_allele_frequency",
        "context",
        "mapped_trait_uri",
        "trait_category",
    }
)

_BODY_MEASUREMENT_KW = (
    "body height",
    "body mass",
    "body weight",
    "body fat",
    "body surface",
    "body size",
    "waist circumference",
    "hip circumference",
    "waist-hip",
    "waist hip",
    "trunk fat",
    "arm fat",
    "leg fat",
    "birth weight",
    "birth length",
    "lean body mass",
    "lean mass",
    "bone mineral density",
    "heel bone",
    "grip strength",
    "hand grip",
    "standing height",
    "sitting height",
    "whole body water",
    "body water mass",
    "impedance of whole body",
    "impedance of arm",
    "impedance of leg",
    "impedance of trunk",
    "ukb data field 231",
)

_LIPID_KW = (
    "cholesterol",
    "triglyceride",
    "lipoprotein",
    "apolipoprotein",
    "phospholipids",
    "total lipids",
    "cholesteryl",
    "vldl",
    "fatty acid",
)

_HEMATOLOGICAL_KW = (
    "red blood cell",
    "red cell",
    "white blood cell",
    "hemoglobin",
    "platelet",
    "eosinophil",
    "basophil",
    "lymphocyte",
    "monocyte",
    "neutrophil",
    "hematocrit",
    "mean corpuscular",
    "reticulocyte",
    "blood cell count",
    "blood protein",
    "erythrocyte",
    "granulocyte",
)

_OTHER_MEASUREMENT_KW = (
    "blood pressure",
    "pulse pressure",
    "heart rate",
    "pulse rate",
    "lung function",
    "forced vital capacity",
    "forced expiratory volume",
    "peak expiratory",
    "fev/fvc",
    "telomere length",
    "albumin measurement",
    "creatinine",
    "urate",
    "bilirubin",
    "c-reactive protein",
    "vitamin d",
    "calcium measurement",
    "iron measurement",
    "ferritin",
    "testosterone measurement",
    "cortisol measurement",
    "glomerular filtration",
    "intraocular pressure",
    "corneal",
    "intracranial volume",
    "cortical thickness",
    "homocysteine measurement",
    "pain measurement",
    "electrocardiogra",
    "qt interval",
    "jt interval",
    "aminotransferase",
    "hair color",
    "eye color",
)

_BEHAVIORAL_KW = (
    "educational attainment",
    "intelligence",
    "cognitive function",
    "cognitive performance",
    "smoking",
    "alcohol consumption",
    "coffee consumption",
    "risk taking",
    "well-being",
    "neuroticism",
    "loneliness",
    "subjective well",
    "life satisfaction",
    "number of children",
    "age at first birth",
)

_CANCER_KW = (
    "cancer",
    "carcinoma",
    "neoplasm",
    "tumor",
    "tumour",
    "melanoma",
    "lymphoma",
    "leukemia",
    "leukaemia",
    "myeloma",
    "sarcoma",
    "glioma",
    "glioblastoma",
    "meningioma",
    "neuroblastoma",
)

_DRUG_RESPONSE_KW = (
    "response to",
    "drug response",
    "drug metabolism",
    "warfarin",
    "clopidogrel",
    "metformin",
    "drug-induced",
    "adverse drug",
)

_IMMUNE_KW = (
    "rheumatoid arthritis",
    "lupus",
    "crohn",
    "ulcerative colitis",
    "psoriasis",
    "multiple sclerosis",
    "celiac",
    "coeliac",
    "autoimmune",
    "inflammatory bowel",
    "ankylosing spondylitis",
    "asthma",
    "allergy",
    "allergic",
    "atopic",
)

_CARDIOVASCULAR_KW = (
    "coronary artery disease",
    "coronary heart disease",
    "myocardial infarction",
    "heart failure",
    "atrial fibrillation",
    "stroke",
    "cerebrovascular",
    "aortic",
    "venous thromboembolism",
    "pulmonary embolism",
    "peripheral artery",
    "hypertension",
)

_METABOLIC_KW = (
    "type 2 diabetes",
    "type 1 diabetes",
    "diabetes mellitus",
    "metabolic syndrome",
    "obesity",
    "insulin resistance",
    "glycated hemoglobin",
    "fasting glucose",
    "fasting insulin",
)

_NEUROLOGICAL_KW = (
    "alzheimer",
    "parkinson",
    "epilepsy",
    "seizure",
    "schizophrenia",
    "bipolar disorder",
    "autism",
    "attention deficit",
    "dementia",
    "amyotrophic lateral",
    "huntington",
    "migraine",
    "major depressive disorder",
    "anxiety disorder",
)

_DISEASE_CATCHALL_KW = (
    "disease",
    "disorder",
    "syndrome",
    "deficiency",
    "infection",
    "sepsis",
    "pneumonia",
    "hepatitis",
    "cirrhosis",
    "nephropathy",
    "retinopathy",
    "neuropathy",
    "cardiomyopathy",
    "myopathy",
    "anemia",
    "anaemia",
)


def _is_metabolite_ratio(trait_lc: str) -> bool:
    """True for NMR metabolomics ratio traits (e.g. 'cholesterol-to-phospholipid ratio')."""
    return "-to-" in trait_lc and " ratio" in trait_lc


def _is_uncharacterized_analyte(trait_lc: str) -> bool:
    """True for Nightingale/Metabolon uncharacterized analytes (e.g. 'x-12345 level')."""
    return trait_lc.startswith("x-") and " level" in trait_lc


def classify_gwas_trait(
    mapped_trait: str,
    mapped_trait_uri: str,
    disease_trait: str = "",
) -> str:
    """Classify a GWAS Catalog trait for default report filtering (ADR-0024).

    Keywords are matched against MAPPED_TRAIT + DISEASE/TRAIT combined.
    GWAS Catalog's EFO mapping is inconsistent — UKB data-field rows
    often have empty MAPPED_TRAIT but informative DISEASE/TRAIT (e.g.
    'Impedance of whole body (UKB data field 23106)'). Either field
    populates the substring keyword space.
    """
    trait = f"{mapped_trait} {disease_trait}".lower().strip()
    uri = mapped_trait_uri.lower()

    if "mondo_" in uri:
        return "disease"
    if "oba_" in uri:
        return "other_measurement"

    if _is_metabolite_ratio(trait):
        return "other_measurement"
    if _is_uncharacterized_analyte(trait):
        return "other_measurement"

    if any(kw in trait for kw in _CANCER_KW):
        return "cancer"
    if any(kw in trait for kw in _DRUG_RESPONSE_KW):
        return "drug_response"
    if any(kw in trait for kw in _IMMUNE_KW):
        return "immune"
    if any(kw in trait for kw in _CARDIOVASCULAR_KW):
        return "cardiovascular"
    if any(kw in trait for kw in _METABOLIC_KW):
        return "metabolic"
    if any(kw in trait for kw in _NEUROLOGICAL_KW):
        return "neurological"
    if any(kw in trait for kw in _DISEASE_CATCHALL_KW):
        return "disease"

    if any(kw in trait for kw in _BODY_MEASUREMENT_KW):
        return "body_measurement"
    if any(kw in trait for kw in _LIPID_KW):
        return "lipid_measurement"
    if any(kw in trait for kw in _HEMATOLOGICAL_KW):
        return "hematological_measurement"
    if any(kw in trait for kw in _OTHER_MEASUREMENT_KW):
        return "other_measurement"
    if any(kw in trait for kw in _BEHAVIORAL_KW):
        return "behavioral"

    if "measurement" in trait:
        return "other_measurement"

    return "other"


def schema_is_current(db_path: Path) -> bool:
    """True iff the cache has the expected GWAS schema and categorizer version.

    Checks ``local_version_tag`` for the ``cv:N`` stamp.  Legacy caches
    that still carry the tag inside ``remote_signal`` are migrated by
    ``_stamp_existing_gwas_cache`` before this runs.
    """
    if not db_path.exists():
        return False
    try:
        from allelix.databases.manager import get_database_info

        info = get_database_info(db_path, "gwas")
        if info is None:
            return False
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            try:
                cols = {row[1] for row in conn.execute("PRAGMA table_info(gwas_associations)")}
            except sqlite3.DatabaseError:
                return False
            if not _REQUIRED_GWAS_COLUMNS.issubset(cols):
                return False
        tag = info.get("local_version_tag") or ""
        return tag == f"cv:{_CATEGORIZER_VERSION}"
    except sqlite3.DatabaseError:
        return False


def _parse_risk_allele(field: str) -> str | None:
    """Extract single-base risk allele from 'rs123-A' format, or None."""
    if not field or "-" not in field:
        return None
    allele = field.rsplit("-", 1)[-1].strip().upper()
    if len(allele) == 1 and allele in "ACGT":
        return allele
    return None


def _safe_float(value: str) -> float | None:
    if not value or value.strip().upper() == "NR":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def iter_gwas_records(tsv_path: Path) -> Iterator[dict[str, object]]:
    """Yield one record per qualifying row from the GWAS Catalog TSV.

    Skips multi-SNP haplotypes, rows without traits, and rows where the
    p-value is unparseable. Deduplicates by (rsid, trait), keeping the
    row with the lowest p-value.
    """
    best: dict[tuple[str, str], dict[str, object]] = {}

    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            snp_field = (row.get("SNPS") or "").strip()
            if not _RSID_RE.match(snp_field):
                continue
            trait = (row.get("DISEASE/TRAIT") or "").strip()
            if not trait:
                continue

            p_value = _safe_float(row.get("P-VALUE", ""))
            risk_allele = _parse_risk_allele(row.get("STRONGEST SNP-RISK ALLELE", ""))
            or_beta = _safe_float(row.get("OR or BETA", ""))
            raf = _safe_float(row.get("RISK ALLELE FREQUENCY", ""))
            gene = (row.get("MAPPED_GENE") or "").strip() or None
            study = (row.get("STUDY ACCESSION") or "").strip() or None
            pubmed = (row.get("PUBMEDID") or "").strip() or None
            ci_text = (row.get("95% CI (TEXT)") or "").strip() or None
            context = (row.get("CONTEXT") or "").strip() or None
            mapped_trait = (row.get("MAPPED_TRAIT") or "").strip()
            mapped_trait_uri = (row.get("MAPPED_TRAIT_URI") or "").strip()
            trait_category = classify_gwas_trait(
                mapped_trait,
                mapped_trait_uri,
                disease_trait=trait,
            )

            key = (snp_field, trait)
            record = {
                "rsid": snp_field,
                "risk_allele": risk_allele,
                "trait": trait,
                "p_value": p_value,
                "or_beta": or_beta,
                "ci_text": ci_text,
                "gene": gene,
                "study_accession": study,
                "pubmed_id": pubmed,
                "risk_allele_frequency": raf,
                "context": context,
                "mapped_trait_uri": mapped_trait_uri or None,
                "trait_category": trait_category,
            }

            existing = best.get(key)
            if existing is None:
                best[key] = record
            else:
                ep = existing["p_value"]
                if p_value is not None and (ep is None or p_value < ep):
                    best[key] = record

    yield from best.values()


def load_gwas_tsv(
    tsv_path: Path,
    db_path: Path,
    source_url: str = "",
    remote_signal: str | None = None,
) -> int:
    """Parse a GWAS Catalog TSV into a fresh SQLite cache atomically.

    Writes to a `.tmp` sibling and `os.replace`s onto `db_path` only after a
    successful commit. Returns the number of records loaded.
    """
    tmp_path = db_path.parent / f"{db_path.name}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    version = datetime.now(UTC).strftime("%Y-%m-%d")

    try:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.executescript(GWAS_SCHEMA)
            insert_sql = (
                "INSERT INTO gwas_associations "
                "(rsid, risk_allele, trait, p_value, or_beta, ci_text, "
                "gene, study_accession, pubmed_id, risk_allele_frequency, "
                "context, mapped_trait_uri, trait_category) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            )
            batch: list[tuple] = []
            count = 0
            for record in iter_gwas_records(tsv_path):
                batch.append(
                    (
                        record["rsid"],
                        record["risk_allele"],
                        record["trait"],
                        record["p_value"],
                        record["or_beta"],
                        record["ci_text"],
                        record["gene"],
                        record["study_accession"],
                        record["pubmed_id"],
                        record["risk_allele_frequency"],
                        record["context"],
                        record["mapped_trait_uri"],
                        record["trait_category"],
                    )
                )
                if len(batch) >= INSERT_BATCH_SIZE:
                    conn.executemany(insert_sql, batch)
                    count += len(batch)
                    batch.clear()
            if batch:
                conn.executemany(insert_sql, batch)
                count += len(batch)
            conn.execute(
                "INSERT INTO database_versions "
                "(name, source_url, version, downloaded_at, record_count, "
                "remote_signal, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    "gwas",
                    source_url,
                    version,
                    datetime.now(UTC).isoformat(),
                    count,
                    remote_signal or "",
                    f"cv:{_CATEGORIZER_VERSION}",
                ),
            )
            conn.commit()
        os.replace(tmp_path, db_path)
        return count
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                logger.warning("Could not remove failed temp DB %s", tmp_path)
        raise

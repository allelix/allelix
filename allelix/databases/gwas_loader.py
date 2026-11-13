# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GWAS Catalog association download, parse, and load into SQLite.

EBI/NHGRI publishes the full associations TSV at a stable URL. Public
domain -- no licensing restrictions. One row per study x SNP x trait.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import re
import sqlite3
from collections import Counter
from datetime import UTC, datetime
from typing import TYPE_CHECKING, cast

from allelix.databases.loader_utils import publish_cache
from allelix.databases.schema import GWAS_SCHEMA
from allelix.databases.staging import owned_staging

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


GWAS_CATALOG_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip"
)
GWAS_DB_FILENAME = "gwas.sqlite"

INSERT_BATCH_SIZE = 5_000

_CATEGORIZER_VERSION = 4

# Source rows skipped because their risk allele is not the row's own SNP's:
# the risk field names another rsID, or something that is not an rsID.
# Counted per source release, in ``gwas_source_skips``.
RISK_OTHER_SNP = "other_snp"
RISK_UNATTRIBUTED = "unattributed"
RISK_SKIP_REASONS: tuple[str, ...] = (RISK_OTHER_SNP, RISK_UNATTRIBUTED)

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
        from allelix.databases.cache_state import connect_readonly

        with contextlib.closing(connect_readonly(db_path)) as conn:
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


def _risk_allele_for(snp: str, field: str) -> tuple[str | None, str | None]:
    """``(risk allele, skip reason)`` of a row whose SNP is ``snp``.

    A single-base risk allele belongs to the row only when the text before
    the last ``-`` is the row's own rsID (trimmed, any case). Otherwise the
    row is skipped: ``other_snp`` when that text is another rsID (no merge
    history is available to tell an alias), ``unattributed`` when it is not
    an rsID. A field without a single-base allele keeps no risk allele and
    is not skipped, as before.
    """
    allele = _parse_risk_allele(field)
    if allele is None:
        return None, None
    prefix = field.rsplit("-", 1)[0].strip().lower()
    if prefix == snp.lower():
        return allele, None
    return None, RISK_OTHER_SNP if _RSID_RE.match(prefix) else RISK_UNATTRIBUTED


def _safe_float(value: str) -> float | None:
    if not value or value.strip().upper() == "NR":
        return None
    try:
        return float(value)
    except ValueError:
        return None


class GwasRecords:
    """The records of a GWAS Catalog TSV (``iter_gwas_records``), with its skips.

    ``skipped`` counts, by reason, the source rows left out because their
    risk allele is not their own SNP's; it is complete once iteration ends.
    """

    def __init__(self, tsv_path: Path) -> None:
        """Read ``tsv_path`` lazily, as the records are taken."""
        self.skipped: Counter[str] = Counter()
        self._records = _gwas_records(tsv_path, self.skipped)

    def __iter__(self) -> Iterator[dict[str, object]]:
        """The records, once."""
        return self._records

    def __next__(self) -> dict[str, object]:
        """The next record."""
        return next(self._records)


def iter_gwas_records(tsv_path: Path) -> GwasRecords:
    """One record per qualifying row from the GWAS Catalog TSV (``_gwas_records``)."""
    return GwasRecords(tsv_path)


def _gwas_records(tsv_path: Path, skipped: Counter[str]) -> Iterator[dict[str, object]]:
    """Yield one record per qualifying row from the GWAS Catalog TSV.

    Skips multi-SNP haplotypes, rows without traits, and rows where the
    p-value is unparseable. Skips, before the selection below, rows whose
    risk allele belongs to another SNP or to no SNP (`_risk_allele_for`),
    counting each source row by reason in ``skipped``. Deduplicates by
    (rsid, trait), keeping the row with the lowest p-value.

    Not strictly streaming: the (rsid, trait) dedup-on-min-p-value pass
    has to see every row before any (rsid, trait) is final, so the
    intermediate ``best`` dict holds the post-filter row set in memory
    until the input is exhausted. That is acceptable: the
    GWAS Catalog has ~600k post-filter rows, peak memory stays well
    under a normal-machine budget, and a true streaming pass would
    require either a pre-sort on disk or a two-pass design with
    materially worse wall-clock. Worth revisiting only if the catalog
    grows an order of magnitude.
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
            risk_allele, skip = _risk_allele_for(
                snp_field, row.get("STRONGEST SNP-RISK ALLELE") or ""
            )
            if skip is not None:
                skipped[skip] += 1
                continue
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
                ep = cast("float | None", existing["p_value"])
                if p_value is not None and (ep is None or p_value < ep):
                    best[key] = record

    yield from best.values()


# Truncation sanity floor for production loads. This guards the count
# returned by iter_gwas_records — i.e. rows AFTER haplotype/no-trait
# filtering and (rsid, trait) dedup, not the raw catalog. EBI curates
# ~625K lead associations (GWAS Catalog, 2025); the loaded count is lower
# than that but still far above this floor. 100K only catches gross
# truncation (a mid-stream download committed as "complete") while staying
# permissive against legitimate upstream drift. Set to 0 from tests so
# synthetic fixtures of any size load cleanly.
GWAS_MIN_ROWS = 100_000


def source_skip_counts(conn: sqlite3.Connection) -> dict[str, int] | None:
    """The source rows the cache's release skipped, by reason; None before 2.7.0.

    Read from ``gwas_source_skips``, which ``load_gwas_tsv`` writes with
    every load, so the counts are the cached release's, without its table.
    """
    try:
        rows = conn.execute("SELECT reason, source_rows FROM gwas_source_skips").fetchall()
    except sqlite3.DatabaseError:
        return None
    counts = {reason: int(n) for reason, n in rows}
    return {reason: counts.get(reason, 0) for reason in RISK_SKIP_REASONS}


def gwas_source_hash(tsv_path: Path) -> str:
    """SHA-256 of a GWAS Catalog table's bytes, as stored in ``content_hash``."""
    digest = hashlib.sha256()
    with tsv_path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def retained_source_stamp(db_path: Path) -> tuple[str | None, str | None]:
    """The cache's freshness signal and the hash of the table it was built from (read-only).

    Both come from one read, so they describe the same build. Caches built
    before 2.6.0 have no ``content_hash`` column; it reads as None and the
    cache stays valid.
    """
    from allelix.databases.cache_state import connect_readonly

    try:
        with contextlib.closing(connect_readonly(db_path)) as conn:
            try:
                row = conn.execute(
                    "SELECT remote_signal, content_hash FROM database_versions WHERE name = 'gwas'"
                ).fetchone()
            except sqlite3.OperationalError:
                row = conn.execute(
                    "SELECT remote_signal, NULL FROM database_versions WHERE name = 'gwas'"
                ).fetchone()
    except sqlite3.Error:
        return None, None
    if row is None:
        return None, None
    return row[0] or None, row[1] or None


def load_gwas_tsv(
    tsv_path: Path,
    db_path: Path,
    source_url: str = "",
    remote_signal: str | None = None,
    min_rows: int = 0,
    content_hash: str | None = None,
) -> int:
    """Parse a GWAS Catalog TSV into a fresh SQLite cache atomically.

    Writes to a staging file of its own beside `db_path` (`owned_staging`)
    and `os.replace`s it onto `db_path` only after a successful commit; the
    staging file is removed on any failure or interruption. Returns the number of records loaded.

    ``min_rows`` is a sanity floor checked before the final ``os.replace``.
    Set by production callers (see ``GwasAnnotator.setup``) to
    ``GWAS_MIN_ROWS``; defaults to 0 so test fixtures of any size load.

    ``content_hash`` (``gwas_source_hash`` of ``tsv_path``) is stored so a
    later re-ingest reads only the table this cache was built from.
    """
    version = datetime.now(UTC).strftime("%Y-%m-%d")

    with owned_staging(db_path) as tmp_path:
        with contextlib.closing(sqlite3.connect(tmp_path)) as conn:
            conn.executescript(GWAS_SCHEMA)
            # Optional since 2.6.0; caches without it stay valid.
            conn.execute("ALTER TABLE database_versions ADD COLUMN content_hash TEXT")
            insert_sql = (
                "INSERT INTO gwas_associations "
                "(rsid, risk_allele, trait, p_value, or_beta, ci_text, "
                "gene, study_accession, pubmed_id, risk_allele_frequency, "
                "context, mapped_trait_uri, trait_category) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            )
            batch: list[tuple] = []
            count = 0
            records = iter_gwas_records(tsv_path)
            for record in records:
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
                "remote_signal, local_version_tag, content_hash) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    "gwas",
                    source_url,
                    version,
                    datetime.now(UTC).isoformat(),
                    count,
                    remote_signal or "",
                    f"cv:{_CATEGORIZER_VERSION}",
                    content_hash,
                ),
            )
            conn.execute(
                "CREATE TABLE gwas_source_skips "
                "(reason TEXT PRIMARY KEY, source_rows INTEGER NOT NULL)"
            )
            # A stand-in iterator without counts records none.
            skipped: Counter[str] = (
                records.skipped if isinstance(records, GwasRecords) else Counter()
            )
            for reason in RISK_SKIP_REASONS:
                conn.execute(
                    "INSERT INTO gwas_source_skips (reason, source_rows) VALUES (?, ?)",
                    (reason, skipped[reason]),
                )
            conn.commit()
        if count < min_rows:
            msg = (
                f"GWAS Catalog load aborted: only {count:,} rows ingested "
                f"(floor {min_rows:,}). The download was likely truncated "
                f"in flight (chunked transfer with no Content-Length, or "
                f"connection drop mid-stream). Retry with "
                f"`allelix db update --force`."
            )
            raise OSError(msg)
        publish_cache(tmp_path, db_path)
        return count

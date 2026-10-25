# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ClinVar VCF download, parse, and load into SQLite."""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import logging
import os
import re
import sqlite3
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, TypedDict

from allelix import __version__
from allelix.databases.cache_state import connect_readonly, is_wal
from allelix.databases.schema import CLINVAR_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

# ADR-0021 per-build dispatch is preserved at the SQLite cache
# layer (clinvar.GRCh37.sqlite + clinvar.GRCh38.sqlite), but the source
# is a single TSV pair that carries rows for every supported Assembly.
# variant_summary.txt.gz holds per-(VariationID, Assembly) rows with
# position / ref / alt / rsID / gene; submission_summary.txt.gz holds
# per-SCV (VariationID, ClinicalSignificance, ReportedPhenotypeInfo,
# ReviewStatus, SCV ID) rows that join on VariationID. The loader
# emits one cache record per (variant, SCV) — the per-SCV pairing
# that the prior VCF loader's CLNSIG|CLNDN parse could not.
CLINVAR_VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)
CLINVAR_SUBMISSION_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
)
CLINVAR_MIN_ROWS = {"GRCh37": 5_000_000, "GRCh38": 5_000_000}
INSERT_BATCH_SIZE = 5_000
DOWNLOAD_TIMEOUT_SECONDS = 60
SIGNAL_TIMEOUT_SECONDS = 15
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
USER_AGENT = f"allelix/{__version__} (+https://github.com/allelix/allelix)"


class DatabaseInfo(TypedDict):
    """Cached database version metadata."""

    source_url: str
    version: str
    downloaded_at: str
    record_count: int
    remote_signal: str | None
    local_version_tag: str | None


def fetch_remote_text(
    url: str, timeout: float = SIGNAL_TIMEOUT_SECONDS, *, max_bytes: int | None = None
) -> str | None:
    """Fetch a small text resource (e.g., a `.md5` file) and return its body.

    Returns None on any failure — `db update`'s freshness check treats
    `None` as "can't verify" and falls through to a "skip with notice".
    Never raises.
    """
    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = response.read() if max_bytes is None else response.read(max_bytes + 1)
            if max_bytes is not None and len(body) > max_bytes:
                return None
            return body.decode("utf-8", errors="replace")
    except (OSError, ValueError) as exc:
        if hasattr(exc, "close"):
            exc.close()
        return None


def head_request_headers(
    url: str, timeout: float = SIGNAL_TIMEOUT_SECONDS
) -> dict[str, str] | None:
    """Issue an HTTP HEAD and return the response headers as a plain dict.

    Returns None on any failure. Never raises.
    """
    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD")
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return dict(response.headers.items())
    except (OSError, ValueError) as exc:
        if hasattr(exc, "close"):
            exc.close()
        return None


def download(url: str, dest: Path) -> None:
    """Download `url` to `dest`. Streaming, atomic (.part rename), with timeout.

    - Streams chunks directly to a `.part` sibling file (no in-memory copy).
    - Sets a real User-Agent so CDNs don't reject the default python-urllib UA.
    - `os.replace`s the .part onto `dest` only after a full successful write,
      so a killed mid-download never leaves a half-file at the target name.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part_path = dest.parent / f"{dest.name}.part"
    if part_path.exists():
        part_path.unlink()

    logger.info("Downloading %s -> %s", url, dest)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with (
            urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response,
            part_path.open("wb") as out,
        ):
            expected_size = response.headers.get("Content-Length")
            while chunk := response.read(DOWNLOAD_CHUNK_SIZE):
                out.write(chunk)
            out.flush()
            try:
                os.fsync(out.fileno())
            except OSError:
                logger.debug("fsync unsupported on this filesystem; continuing")

        actual_size = part_path.stat().st_size
        if expected_size is not None:
            expected = int(expected_size)
            if actual_size != expected:
                part_path.unlink(missing_ok=True)
                raise OSError(
                    f"Download truncated: expected {expected:,} bytes, "
                    f"got {actual_size:,} bytes from {url}"
                )
        os.replace(part_path, dest)
    except Exception as exc:
        if hasattr(exc, "close"):
            exc.close()
        if part_path.exists():
            try:
                part_path.unlink()
            except OSError:
                logger.warning("Could not remove stale partial download %s", part_path)
        raise


def verify_file_hash(path: Path, algorithm: str, expected_hex: str) -> None:
    """Verify a file's cryptographic hash and delete it on mismatch.

    Reads the file in streaming chunks to avoid loading multi-GB files
    into memory. On mismatch, deletes the file and raises ``OSError``.
    """
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        while chunk := f.read(DOWNLOAD_CHUNK_SIZE):
            h.update(chunk)
    actual = h.hexdigest()
    if actual != expected_hex:
        path.unlink(missing_ok=True)
        raise OSError(
            f"Integrity check failed for {path.name}: "
            f"expected {algorithm}:{expected_hex}, "
            f"got {algorithm}:{actual}"
        )


# --- per-SCV TSV loader ------------------------------------
#
# Column indices (zero-based) inside the TSV files, taken from the
# verified column-header rows fetched 2026-06-17 from
#   https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
# Both files use `\t` as field separator and `#` as comment-line prefix.

# variant_summary.txt.gz — 43 columns. We only need a handful.
_VS_TYPE = 1  # "Type": "single nucleotide variant", "Deletion", ...
_VS_AGGREGATE = 6  # "ClinicalSignificance": ClinVar's aggregate classification
_VS_NAME = 2  # "Name": HGVS, e.g. "NM_000130.5(F5):c.1601G>A (p.Arg534Gln)"
_VS_RS = 9  # "RS# (dbSNP)"     — number only, no `rs` prefix
_VS_GENE_SYMBOL = 4  # "GeneSymbol"
_VS_REVIEW_STATUS = 24  # "ReviewStatus" — aggregate
_VS_ASSEMBLY = 16  # "Assembly"        — "GRCh37" / "GRCh38" / "NCBI36"
_VS_CHROMOSOME = 18  # "Chromosome"
_VS_VARIATION_ID = 30  # "VariationID"   — join key against submission_summary
_VS_POSITION_VCF = 31  # "PositionVCF"
_VS_REF_ALLELE_VCF = 32  # "ReferenceAlleleVCF"
_VS_ALT_ALLELE_VCF = 33  # "AlternateAlleleVCF"

# submission_summary.txt.gz — 16 columns.
_SS_VARIATION_ID = 0  # "VariationID"   — join key
_SS_CLIN_SIG = 1  # "ClinicalSignificance"  — per-SCV (single token)
_SS_REPORTED_PHENO = 5  # "ReportedPhenotypeInfo" — "MedGenID:Name" tuples
_SS_REVIEW_STATUS = 6  # "ReviewStatus" — per-SCV
_SS_SCV = 10  # "SCV"                   — submission accession
_SS_CONTRIBUTES = 15  # "ContributesToAggregateClassification"

# The header name each column index above must carry. The loader checks
# both files' header rows against these before reading a record, so an
# upstream column change fails the update instead of building a cache
# from shifted columns.
_VS_REQUIRED_COLUMNS: dict[int, str] = {
    _VS_TYPE: "Type",
    _VS_NAME: "Name",
    _VS_AGGREGATE: "ClinicalSignificance",
    _VS_GENE_SYMBOL: "GeneSymbol",
    _VS_RS: "RS# (dbSNP)",
    _VS_ASSEMBLY: "Assembly",
    _VS_CHROMOSOME: "Chromosome",
    _VS_REVIEW_STATUS: "ReviewStatus",
    _VS_VARIATION_ID: "VariationID",
    _VS_POSITION_VCF: "PositionVCF",
    _VS_REF_ALLELE_VCF: "ReferenceAlleleVCF",
    _VS_ALT_ALLELE_VCF: "AlternateAlleleVCF",
}
_SS_REQUIRED_COLUMNS: dict[int, str] = {
    _SS_VARIATION_ID: "VariationID",
    _SS_CLIN_SIG: "ClinicalSignificance",
    _SS_REPORTED_PHENO: "ReportedPhenotypeInfo",
    _SS_REVIEW_STATUS: "ReviewStatus",
    _SS_SCV: "SCV",
    _SS_CONTRIBUTES: "ContributesToAggregateClassification",
}

# Per-SCV review status ClinVar gives a submission it has flagged (for
# example as outdated or in error). Such submissions are not loaded.
_FLAGGED_REVIEW_STATUS = "flagged submission"

# Per-SCV review tiers, strongest first, as ClinVar ranks submissions when
# it computes an aggregate. A status not listed ranks below all of them.
_SUBMISSION_TIERS: dict[str, int] = {
    "practice guideline": 0,
    "reviewed by expert panel": 1,
    "criteria provided, single submitter": 2,
    "no assertion criteria provided": 3,
}
_UNRANKED_TIER = 4

_PATHOGENIC_CLNSIGS: frozenset[str] = frozenset(
    {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "pathogenic, low penetrance",
        "likely pathogenic, low penetrance",
    }
)

# Classifications of disease causation: the terms ClinVar's conflict rule
# compares. Drug response, risk factor, protective and the like are a
# different kind of classification.
_DISEASE_CLNSIGS: frozenset[str] = _PATHOGENIC_CLNSIGS | frozenset(
    {
        "uncertain significance",
        "likely benign",
        "benign",
        "benign/likely benign",
        "vus-high",
        "vus-mid",
        "vus-low",
        "conflicting classifications of pathogenicity",
        "no classifications from unflagged records",
    }
)

# ``outside_aggregate`` values in the cache.
CONTRIBUTES_TO_AGGREGATE = 0
OUTSIDE_AGGREGATE_AGREED = 1
OUTSIDE_AGGREGATE_CONFLICTING = 2


def aggregate_germline_term(aggregate: str) -> str:
    """The germline part of a ClinVar aggregate, lower-cased.

    ClinVar appends other kinds of classification after a semicolon:
    ``Conflicting classifications of pathogenicity; risk factor`` gives
    ``conflicting classifications of pathogenicity``.
    """
    return aggregate.split(";")[0].strip().lower()


_PATHOGENIC_AGGREGATE_PARTS = frozenset(
    {
        "pathogenic",
        "likely pathogenic",
        "pathogenic, low penetrance",
        "likely pathogenic, low penetrance",
    }
)


def is_pathogenic_aggregate(aggregate: str) -> bool:
    """True when ClinVar's aggregate is Pathogenic or Likely pathogenic (ADR-0050).

    Every ``/``-separated part of the germline term is Pathogenic or Likely
    pathogenic (low penetrance included): ``Pathogenic/Likely pathogenic;
    risk factor`` counts, ``Pathogenic/Likely risk allele`` and
    ``Conflicting classifications of pathogenicity`` do not.
    """
    # ClinVar's aggregate uses spaces; the per-submission form uses "_".
    term = aggregate_germline_term(aggregate).replace("_", " ")
    return bool(term) and all(
        part.strip() in _PATHOGENIC_AGGREGATE_PARTS for part in term.split("/")
    )


def is_disease_classification(aggregate: str) -> bool:
    """True when ClinVar's aggregate classifies disease causation.

    False when the top review tier gave another kind of classification
    (``drug response``, ``risk factor``, ``not provided``, ...), so the
    aggregate says nothing about pathogenicity. A combined term such as
    ``Uncertain significance/Uncertain risk allele`` counts as disease.
    """
    term = aggregate_germline_term(aggregate)
    return any(part.strip() in _DISEASE_CLNSIGS for part in term.split("/")) or (
        term in _DISEASE_CLNSIGS
    )


# Exact set of ClinicalSignificance values the TSV loader filters at ingest.
# These are values that ClinVar uses to mark SCVs with no actionable
# classification — placeholders ("-", "", "not specified") AND
# non-classification curatorial categories ("other", "association",
# "association not found") — AND that aren't in
# `allelix.annotators.clinvar._CLNSIG_MAGNITUDE` (so without this
# filter they'd land at the 5.0 default — equal to the analyze
# display floor — and surface to users as bogus annotations).
#
# Values intentionally NOT in this set:
#   - "not provided" / "not_provided" — already maps to 2.0 in
#     _CLNSIG_MAGNITUDE (safe below the floor); dropping at ingest
#     would lose real submitter records.
#   - "no classification for the single variant" — same.
#
# v2.2.1 added the three non-classification terms ("other",
# "association", "association not found"). The per-SCV TSV
# switch surfaced these from the underlying submission_summary
# rows where the old summarized data hid them; they're meaningfully
# distinct from placeholders (a curator chose them) but equally
# unactionable as report content.
#
# The protocol's §7b significance-sentinel ship-gate scans the
# same set to verify the loader's commitment against the live cache.
# Keep both in sync.
_CLINVAR_PLACEHOLDER_CLNSIGS: frozenset[str] = frozenset(
    {
        "",
        "-",
        "not specified",
        "no classification provided",
        # non-classification curatorial terms — kept out of
        # _CLNSIG_MAGNITUDE deliberately (they aren't classifications)
        # and out of _BENIGN_CLNSIGS (they aren't benign; the repute
        # work must not see them as "good").
        "other",
        "association",
        "association not found",
    }
)


# Per-SCV classifications dropped from normal-allele records whose VCF
# alleles differ on the build (see ``describes_reference_allele``).
_NORMAL_ALLELE_DROPPED_CLNSIGS: frozenset[str] = frozenset(
    {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "pathogenic, low penetrance",
        "likely pathogenic, low penetrance",
    }
)

_HGVS_ALLELE = re.compile(r"\b[cgmn]\.([^\s(]+)")


def describes_reference_allele(name: str) -> bool:
    """Return True when a ClinVar record's HGVS name describes the reference sequence.

    ClinVar writes the normal (reference-sequence) allele as ``c.N=``, for
    example ``NM_000130.5(F5):c.1601=``. Only the nucleotide description
    counts: a synonymous protein change such as ``c.123A>G (p.Leu41=)``
    describes a real variant.
    """
    match = _HGVS_ALLELE.search(name)
    return match is not None and match.group(1).endswith("=")


def _open_tsv(path: Path) -> Iterator[list[str]]:
    """Stream a (possibly-gzipped) tab-delimited TSV, skipping comment lines."""
    opener: object = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield line.split("\t")


class ClinVarFormatError(RuntimeError):
    """A ClinVar TSV does not have the columns the loader reads."""


def _read_header(path: Path) -> list[str] | None:
    """Return the column names of a ClinVar TSV, or None when it has no header row.

    The header is the last ``#`` comment line holding a tab before the
    first record. submission_summary precedes it with ``#`` lines that
    explain each column and hold no tab.
    """
    opener: object = gzip.open if path.suffix == ".gz" else open
    header: list[str] | None = None
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if not line.startswith("#"):
                break
            if "\t" in line:
                header = line.lstrip("#").split("\t")
    return header


def check_tsv_header(path: Path, required: dict[int, str]) -> None:
    """Raise ``ClinVarFormatError`` unless each required column has its expected name.

    The loader reads columns by position. A renamed, inserted or removed
    column upstream would otherwise shift every value it reads, skip every
    record, and leave an empty cache that reports ready.
    """
    header = _read_header(path)
    if header is None:
        msg = f"{path.name}: no header row; cannot confirm the ClinVar column layout"
        raise ClinVarFormatError(msg)
    for index, name in sorted(required.items()):
        found = header[index].strip() if index < len(header) else None
        if found != name:
            msg = (
                f"{path.name}: column {index + 1} is {found!r}, expected {name!r}; "
                "the ClinVar file layout has changed"
            )
            raise ClinVarFormatError(msg)


def _decode_reported_phenotype(field: str) -> str:
    """Strip MedGen ID prefixes from a `ReportedPhenotypeInfo` cell.

    submission_summary encodes conditions as ``MedGenID:Name`` tuples
    separated by ``|``. Sometimes ``MedGenID`` is ``na`` when no MedGen
    identifier exists. We want the human-readable name without the
    prefix; multiple conditions stay ``;``-joined to match the VCF
    loader's output convention so the downstream display path is
    unchanged.

    Examples:
        ``"C3150901:Hereditary spastic paraplegia 48"``
            → ``"Hereditary spastic paraplegia 48"``
        ``"na:not provided"``
            → ``"not provided"``
        ``"C1|C2:Cond A|C3:Cond B"``
            → ``"Cond A; Cond B"``  (lone token "C1" dropped — no Name)
    """
    if not field or field == "-":
        return ""
    parts = []
    for token in field.split("|"):
        token = token.strip()
        if not token or token == "-":
            continue
        if ":" not in token:
            # Token without "ID:Name" structure — skip rather than emit
            # a bare MedGen ID that looks like noise.
            continue
        _, _, name = token.partition(":")
        name = name.strip()
        if name and name != "-":
            parts.append(name)
    return "; ".join(parts)


@dataclass(frozen=True)
class _OutrankedDecision:
    """ClinVar's conflict rule applied to a variant's outranked submissions.

    ``tier`` is the strongest review tier among the outranked, non-flagged
    disease classifications; only that tier is judged, as ClinVar judges
    the aggregate. ``counts`` holds each classification's count within it,
    ``total`` their sum, and ``agreed`` is True when every one is P/LP.
    """

    tier: int
    counts: dict[str, int]
    total: int
    agreed: bool


class _Submission(NamedTuple):
    """One submission_summary row the loader keeps (not flagged, not a placeholder)."""

    significance: str
    review_status: str
    contributes: bool
    reported_phenotype: str

    @property
    def sig(self) -> str:
        return self.significance.lower()

    @property
    def tier(self) -> int:
        return _SUBMISSION_TIERS.get(self.review_status.lower(), _UNRANKED_TIER)


def _variation_groups(path: Path) -> Iterator[tuple[int, list[_Submission]]]:
    """Yield each VariationID's kept submissions, one group at a time.

    submission_summary is sorted by VariationID, which lets the split and
    ADR-0045's tier rule be worked out per variant in one pass. A file that
    is not sorted raises ``ClinVarFormatError`` (the old cache is kept)
    rather than giving a variant two partial groups.
    """
    current: int | None = None
    group: list[_Submission] = []
    for cols in _open_tsv(path):
        if len(cols) <= _SS_CONTRIBUTES:
            continue
        try:
            variation_id = int(cols[_SS_VARIATION_ID].strip())
        except ValueError:
            continue
        if variation_id != current:
            if current is not None and variation_id < current:
                msg = (
                    f"{path.name}: VariationID {variation_id} follows {current}; "
                    "the file is no longer sorted by VariationID"
                )
                raise ClinVarFormatError(msg)
            if current is not None and group:
                yield current, group
            current, group = variation_id, []
        review = cols[_SS_REVIEW_STATUS].strip()
        significance = cols[_SS_CLIN_SIG].strip()
        # ClinVar has flagged the submission; it may still count toward the
        # aggregate when every submission on the variant is flagged ("no
        # classifications from unflagged records"). ClinVar's placeholder
        # values ("-", "not specified", ...) record no classification; a
        # case-insensitive match keeps "Not Specified" out too.
        if review.lower() == _FLAGGED_REVIEW_STATUS:
            continue
        if significance.lower() in _CLINVAR_PLACEHOLDER_CLNSIGS:
            continue
        group.append(
            _Submission(
                significance,
                review,
                cols[_SS_CONTRIBUTES].strip().lower() == "yes",
                cols[_SS_REPORTED_PHENO],
            )
        )
    if current is not None and group:
        yield current, group


def _decide_outranked(submissions: list[_Submission]) -> _OutrankedDecision | None:
    """ClinVar's tiered conflict rule over a variant's outranked submissions (ADR-0045).

    Only disease classifications count. Returns None when the strongest
    tier among them holds no P/LP submission.
    """
    disease = [s for s in submissions if not s.contributes and s.sig in _DISEASE_CLNSIGS]
    if not disease:
        return None
    top = min(s.tier for s in disease)
    counts: dict[str, int] = {}
    for s in disease:
        if s.tier == top:
            counts[s.sig] = counts.get(s.sig, 0) + 1
    pathogenic = sum(n for sig, n in counts.items() if sig in _PATHOGENIC_CLNSIGS)
    if pathogenic == 0:
        return None
    total = sum(counts.values())
    return _OutrankedDecision(top, counts, total, pathogenic == total)


def _group_records(
    variation_id: int,
    group: list[_Submission],
    variant: tuple,
    *,
    outranked: _OutrankedDecision | None,
    aggregate_only: bool,
) -> Iterator[dict[str, object]]:
    """The cache records for one variant's kept submissions.

    Each carries ClinVar's aggregate and the submitter split: how many of
    the contributing submissions gave its classification, of how many. An
    outranked submission (ADR-0045) is kept only when ``outranked`` holds a
    decision for the variant, and only from the decided tier; its split
    counts that tier.
    """
    rsid, chromosome, position, ref, alt, gene, review_status_agg, normal, aggregate = variant
    counts: dict[str, int] = {}
    for s in group:
        if s.contributes:
            counts[s.sig] = counts.get(s.sig, 0) + 1
    contributing = sum(counts.values())
    for s in group:
        if normal and s.sig in _NORMAL_ALLELE_DROPPED_CLNSIGS:
            # A normal-allele record names the reference sequence; carrying
            # it is not a pathogenic finding (F5 c.1601=, written T>C on
            # GRCh37). Its other classifications stay as they are.
            continue
        same: int | None
        total: int | None
        if s.contributes:
            same, total, outside = counts[s.sig], contributing, CONTRIBUTES_TO_AGGREGATE
        elif outranked is not None and s.sig in _PATHOGENIC_CLNSIGS and s.tier == outranked.tier:
            same, total = outranked.counts.get(s.sig, 0), outranked.total
            outside = (
                OUTSIDE_AGGREGATE_AGREED if outranked.agreed else OUTSIDE_AGGREGATE_CONFLICTING
            )
        elif aggregate_only:
            continue
        else:
            same, total, outside = None, None, OUTSIDE_AGGREGATE_CONFLICTING
        yield {
            "rsid": rsid,
            "chromosome": chromosome,
            "position": position,
            "ref": ref,
            "alt": alt,
            "clinical_significance": s.significance,
            "condition": _decode_reported_phenotype(s.reported_phenotype),
            "gene": gene,
            # The per-SCV review status is more specific than the aggregate's.
            "review_status": s.review_status or review_status_agg,
            "allele_id": variation_id,
            "aggregate": aggregate,
            "submitters_same": same,
            "submitters_total": total,
            "outside_aggregate": outside,
            "normal_allele": int(normal),
        }


def iter_clinvar_tsv_records(
    variant_summary_path: Path,
    submission_summary_path: Path,
    build: str,
    *,
    aggregate_only: bool = True,
) -> Iterator[dict[str, object]]:
    """Yield one cache record per (variant, SCV submission).

    Joins ``submission_summary.txt.gz`` (per-SCV rows) against
    ``variant_summary.txt.gz`` (per-(VariationID, Assembly) rows) for
    the requested ``build`` and emits one dict per matched submission
    keyed for the downstream cache writer.

    Why two files: ``submission_summary`` carries the true
    (ClinicalSignificance, ReportedPhenotypeInfo) pair for each SCV
    submission. ``variant_summary`` carries the position / ref / alt /
    rsID / gene per (variant, build). The legacy VCF loader's
    CLNSIG|CLNDN Frankenstein pairing (the VCF loader was removed in v2.2.0) is
    structurally absent here because each row is one submission.

    Memory: the variant_summary index uses a temp on-disk SQLite, NOT
    a Python dict. The full dataset is ~50M variants x ~150 bytes per
    dict entry = ~7.5 GB resident, which would OOM small CI runners
    and is wasteful on any host. SQLite's b-tree pages and OS-level
    caching keep peak memory in the low hundreds of MB regardless of
    input size.

    Args:
        variant_summary_path: Path to ``variant_summary.txt.gz`` (or .txt).
        submission_summary_path: Path to ``submission_summary.txt.gz``.
        build: Target genome build — ``"GRCh37"`` or ``"GRCh38"``. Rows
            from other Assembly values are skipped during the
            variant_summary pass.
        aggregate_only: When True (default), skip SCV submissions whose
            ``ContributesToAggregateClassification`` is not ``"yes"``,
            except the outranked P/LP submissions ADR-0045 keeps.
            submission_summary holds current submissions only; ``"no"``
            marks a submission outranked by a higher review tier, which
            ClinVar's aggregate uses alone. Set False to materialize every
            submission (forensic use, not production ingest); outranked
            rows not kept by ADR-0045 then carry
            ``OUTSIDE_AGGREGATE_CONFLICTING`` and no split.

    Yields:
        Dicts shaped for the cache INSERT in ``load_clinvar_tsv``. The
        ``allele_id`` field is repurposed to carry the VariationID
        (per-variant integer, no per-SCV granularity).
    """
    import tempfile

    if build not in ("GRCh37", "GRCh38"):
        msg = f"unsupported build {build!r} — expected 'GRCh37' or 'GRCh38'"
        raise ValueError(msg)
    check_tsv_header(variant_summary_path, _VS_REQUIRED_COLUMNS)
    check_tsv_header(submission_summary_path, _SS_REQUIRED_COLUMNS)

    # Stream variant_summary into a temp SQLite indexed by VariationID
    # for the requested build only. Roughly halves the row count vs.
    # keeping both builds.
    with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        with contextlib.closing(sqlite3.connect(tmp_path)) as tmp_conn:
            tmp_conn.executescript(
                "CREATE TABLE vs ("
                "  variation_id INTEGER PRIMARY KEY, "
                "  rsid TEXT, "
                "  chromosome TEXT, "
                "  position INTEGER, "
                "  ref TEXT, "
                "  alt TEXT, "
                "  gene TEXT, "
                "  review_status TEXT, "
                "  normal_allele INTEGER, "
                "  aggregate TEXT"
                ");"
            )
            # Variants whose aggregate is not a disease classification: their
            # outranked P/LP submissions may still be reported (ADR-0045).
            non_disease: set[int] = set()
            batch: list[tuple[object, ...]] = []
            for cols in _open_tsv(variant_summary_path):
                if len(cols) <= _VS_ALT_ALLELE_VCF:
                    continue
                if cols[_VS_ASSEMBLY] != build:
                    continue
                # A variant without an rsID (RS# -1) is kept with rsid "":
                # VCF input matches it by locus (ADR-0046).
                rs = cols[_VS_RS].strip()
                rsid = "" if not rs or rs == "-1" else f"rs{rs}"
                pos_str = cols[_VS_POSITION_VCF].strip()
                if not pos_str or pos_str == "-1":
                    continue
                try:
                    position = int(pos_str)
                except ValueError:
                    continue
                ref = cols[_VS_REF_ALLELE_VCF].strip()
                alt = cols[_VS_ALT_ALLELE_VCF].strip()
                if not ref or ref in ("-", "na") or not alt or alt in ("-", "na"):
                    # Skip rows without VCF-style ref/alt — these are
                    # complex or copy-number variants the SNV-shaped
                    # cache can't represent. Logged at debug, not warn,
                    # because they're expected (typical clinvar load
                    # drops ~5% of rows here).
                    continue
                try:
                    variation_id = int(cols[_VS_VARIATION_ID].strip())
                except (ValueError, IndexError):
                    continue
                normal_allele = describes_reference_allele(cols[_VS_NAME])
                if ref == alt and cols[_VS_TYPE].strip() != "single nucleotide variant":
                    # A non-SNV with REF == ALT has no alleles on this build.
                    # An SNV with REF == ALT names the allele this build's
                    # reference carries: a variant allele (F5 c.1601G>A on
                    # GRCh37) or, for a normal-allele record, the normal one
                    # (CASP12 rs497116 on GRCh37). Both are kept, so a record
                    # reads the same on both builds; ``normal_allele`` tells
                    # them apart.
                    continue
                gene = cols[_VS_GENE_SYMBOL].strip()
                if " " in gene:
                    # Not a symbol: ClinVar's placeholder for a region,
                    # "covers 10 genes, none of which curated to show
                    # dosage sensitivity". Symbol lists are ";"-joined.
                    gene = ""
                aggregate = cols[_VS_AGGREGATE].strip()
                if not is_disease_classification(aggregate):
                    non_disease.add(variation_id)
                batch.append(
                    (
                        variation_id,
                        rsid,
                        cols[_VS_CHROMOSOME].strip(),
                        position,
                        ref,
                        alt,
                        gene,
                        cols[_VS_REVIEW_STATUS].strip(),
                        int(normal_allele),
                        aggregate,
                    )
                )
                if len(batch) >= INSERT_BATCH_SIZE:
                    tmp_conn.executemany(
                        "INSERT OR IGNORE INTO vs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    batch.clear()
            if batch:
                tmp_conn.executemany(
                    "INSERT OR IGNORE INTO vs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
            tmp_conn.commit()

            # Second pass: stream submission_summary one VariationID at a
            # time, look the variant up once, emit its kept submissions.
            cursor = tmp_conn.cursor()
            for variation_id, group in _variation_groups(submission_summary_path):
                if (
                    aggregate_only
                    and not any(s.contributes for s in group)
                    and (variation_id not in non_disease)
                ):
                    continue
                row = cursor.execute(
                    "SELECT rsid, chromosome, position, ref, alt, gene, review_status, "
                    "normal_allele, aggregate FROM vs WHERE variation_id = ?",
                    (variation_id,),
                ).fetchone()
                if row is None:
                    continue
                yield from _group_records(
                    variation_id,
                    group,
                    row,
                    outranked=_decide_outranked(group) if variation_id in non_disease else None,
                    aggregate_only=aggregate_only,
                )
    finally:
        tmp_path.unlink(missing_ok=True)


def load_clinvar_tsv(
    variant_summary_path: Path,
    submission_summary_path: Path,
    db_path: Path,
    build: str,
    *,
    source_url: str = "",
    remote_signal: str | None = None,
    aggregate_only: bool = True,
    min_rows: int = 0,
    previous_count: int | None = None,
    record_name: str = "clinvar",
) -> None:
    """Build a ClinVar SQLite cache from the per-SCV TSV sources.

    Production ingest path for ClinVar. The annotator returns
    MULTIPLE Annotation objects per variant when multiple SCVs target
    the same (chrom, pos, ref, alt); the annotator-side reconciliation
    landed in v2.2.0.

    The cache is built beside ``db_path`` and replaces it only once it is
    complete, meets ``min_rows`` and at least 90% of ``previous_count``,
    and carries final metadata under ``record_name``. Production setup
    supplies the per-build policy and validates all selected candidates
    before publishing any. Direct synthetic builds can use smaller floors.
    """
    partial = db_path.with_name(db_path.name + ".partial")
    partial.unlink(missing_ok=True)
    try:
        _build_clinvar_cache(
            variant_summary_path,
            submission_summary_path,
            partial,
            build,
            source_url=source_url,
            remote_signal=remote_signal,
            aggregate_only=aggregate_only,
            min_rows=min_rows,
            previous_count=previous_count,
            record_name=record_name,
        )
        os.replace(partial, db_path)
    finally:
        partial.unlink(missing_ok=True)


_CLINVAR_INSERT = (
    "INSERT INTO clinvar_variants "
    "(rsid, chromosome, position, ref, alt, clinical_significance, condition, gene, "
    "review_status, allele_id, aggregate, submitters_same, submitters_total, "
    "outside_aggregate, normal_allele) "
    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
)


def _build_clinvar_cache(
    variant_summary_path: Path,
    submission_summary_path: Path,
    db_path: Path,
    build: str,
    *,
    source_url: str,
    remote_signal: str | None,
    aggregate_only: bool,
    min_rows: int,
    previous_count: int | None,
    record_name: str,
) -> None:
    """Write a new ClinVar cache at ``db_path``; raise when it holds no records."""
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        batch: list[tuple[object, ...]] = []
        for record in iter_clinvar_tsv_records(
            variant_summary_path,
            submission_summary_path,
            build,
            aggregate_only=aggregate_only,
        ):
            batch.append(
                (
                    record["rsid"],
                    record["chromosome"],
                    record["position"],
                    record["ref"],
                    record["alt"],
                    record["clinical_significance"],
                    record["condition"],
                    record["gene"],
                    record["review_status"],
                    record["allele_id"],
                    record["aggregate"],
                    record["submitters_same"],
                    record["submitters_total"],
                    record["outside_aggregate"],
                    record["normal_allele"],
                )
            )
            if len(batch) >= INSERT_BATCH_SIZE:
                conn.executemany(
                    _CLINVAR_INSERT,
                    batch,
                )
                batch.clear()
        if batch:
            conn.executemany(
                _CLINVAR_INSERT,
                batch,
            )
        record_count = conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0]
        relative_floor = (9 * previous_count + 9) // 10 if previous_count is not None else 0
        if record_count == 0 or record_count < min_rows or record_count < relative_floor:
            problem = (
                "the source files gave no records (0 rows)"
                if record_count == 0
                else f"candidate has {record_count:,} rows"
            )
            raise ClinVarFormatError(
                f"ClinVar {build}: {problem}; "
                f"old count {previous_count if previous_count is not None else 'absent'}; "
                f"requires fixed floor {min_rows:,} and 90% floor {relative_floor:,}; "
                "existing caches kept"
            )
        from datetime import UTC, datetime

        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        # Build identity/count metadata is required even when freshness is unknown.
        build_date = datetime.now(UTC).strftime("%Y-%m-%d")
        tag = f"iv:{CLINVAR_INTERPRETER_VERSION}"
        stamp_remote_signal(conn, record_name, remote_signal, source_url)
        conn.execute(
            "UPDATE database_versions SET version = ?, local_version_tag = ?, "
            "record_count = ? WHERE name = ?",
            (build_date, tag, record_count, record_name),
        )
        metadata = conn.execute(
            "SELECT source_url, version, local_version_tag, record_count, remote_signal, "
            "downloaded_at FROM database_versions WHERE name = ?",
            (record_name,),
        ).fetchone()
        if (
            metadata is None
            or metadata[:5] != (source_url, build_date, tag, record_count, remote_signal)
            or not metadata[5]
        ):
            raise ClinVarFormatError(
                f"ClinVar {build}: candidate metadata is invalid; existing caches kept"
            )
        conn.commit()


def get_database_info(db_path: Path, name: str) -> DatabaseInfo | None:
    """Return version metadata for a cached database, or None if not present.

    Tolerates older caches that lack ``remote_signal`` or
    ``local_version_tag`` columns by falling back to progressively
    simpler SELECTs.  Missing columns report as None; the next
    ``db update`` self-heals via the annotator's migration path. The
    file is opened read-only and never changed (ADR-0048). A cache in
    WAL mode is not opened at all, and reads as not present until the
    update step converts it.
    """
    if not db_path.exists() or is_wal(db_path):
        return None
    try:
        with contextlib.closing(connect_readonly(db_path)) as conn:
            remote_signal: str | None = None
            local_version_tag: str | None = None
            try:
                row = conn.execute(
                    "SELECT source_url, version, downloaded_at, record_count, "
                    "remote_signal, local_version_tag "
                    "FROM database_versions WHERE name = ?",
                    (name,),
                ).fetchone()
                if row is None:
                    return None
                (
                    source_url,
                    version,
                    downloaded_at,
                    record_count,
                    remote_signal,
                    local_version_tag,
                ) = row
            except sqlite3.OperationalError:
                try:
                    row = conn.execute(
                        "SELECT source_url, version, downloaded_at, record_count, "
                        "remote_signal FROM database_versions WHERE name = ?",
                        (name,),
                    ).fetchone()
                except sqlite3.OperationalError:
                    try:
                        row = conn.execute(
                            "SELECT source_url, version, downloaded_at, record_count "
                            "FROM database_versions WHERE name = ?",
                            (name,),
                        ).fetchone()
                    except sqlite3.DatabaseError:
                        return None
                    if row is None:
                        return None
                    source_url, version, downloaded_at, record_count = row
                except sqlite3.DatabaseError:
                    return None
                else:
                    if row is None:
                        return None
                    source_url, version, downloaded_at, record_count, remote_signal = row
            except sqlite3.DatabaseError:
                return None
            return DatabaseInfo(
                source_url=source_url,
                version=version,
                downloaded_at=downloaded_at,
                record_count=record_count,
                remote_signal=remote_signal,
                local_version_tag=local_version_tag,
            )
    except sqlite3.DatabaseError:
        return None


def _ensure_local_version_tag_column(conn: sqlite3.Connection) -> None:
    """Add ``local_version_tag`` column if absent (idempotent soft migration)."""
    with contextlib.suppress(sqlite3.OperationalError):
        conn.execute("ALTER TABLE database_versions ADD COLUMN local_version_tag TEXT")


def stamp_remote_signal(
    conn: sqlite3.Connection,
    name: str,
    remote_signal: str | None,
    source_url: str = "",
) -> None:
    """Ensure ``database_versions`` exists and upsert the remote signal.

    Existing rows keep their version / downloaded_at / record_count
    metadata; only ``remote_signal`` is overwritten.  If the table or row
    is missing (pre-built caches shipped without version metadata), both
    are created with placeholder values that ``db status`` degrades
    gracefully on.
    """
    conn.execute(
        "CREATE TABLE IF NOT EXISTS database_versions ("
        "name TEXT PRIMARY KEY, source_url TEXT NOT NULL, version TEXT, "
        "downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
        "remote_signal TEXT, local_version_tag TEXT)"
    )
    _ensure_local_version_tag_column(conn)
    conn.execute(
        "INSERT INTO database_versions"
        " (name, source_url, version, downloaded_at, record_count, remote_signal)"
        " VALUES (?, ?, NULL, datetime('now'), 0, ?)"
        " ON CONFLICT(name) DO UPDATE SET remote_signal = excluded.remote_signal",
        (name, source_url, remote_signal),
    )


def stamp_existing_clinvar_cache(db_path: Path) -> bool:
    """One-shot migration: self-heal a ClinVar cache up to the current iv tag.

    Self-heal is allowed ONLY when there is positive evidence the cache
    was built by the current interpreter version. Absent that evidence
    the function returns False and ``ClinVarAnnotator.is_ready()`` lets
    ``db update`` reingest from upstream. Silent promotion of an
    unknown-version cache is the failure mode this exists to prevent —
    serving pre-format-change data labeled as fresh would be worse than
    paying the redownload cost.

    Decision matrix (for each ``database_versions`` row matching
    ``clinvar%``):

    +-----------------------+---------------------------+-------------------+
    | local_version_tag     | remote_signal `|iv:N`     | Action            |
    +=======================+===========================+===================+
    | == CURRENT            | (any)                     | noop, return True |
    +-----------------------+---------------------------+-------------------+
    | != CURRENT (some N)   | (any)                     | return False      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | baked N == CURRENT        | self-heal: move   |
    |                       |                           | tag to column,    |
    |                       |                           | strip `|iv:`      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | baked N != CURRENT        | return False      |
    +-----------------------+---------------------------+-------------------+
    | NULL                  | no `|iv:` marker          | return False      |
    +-----------------------+---------------------------+-------------------+

    The bottom-right case is the v2.2 safety hole: a
    pre-v2.0.1 cache (iv:1-era, NULL tag, no baked marker) skipping
    straight to v2.2 used to get stamped as iv:CURRENT unconditionally —
    serving old single-row VCF data labeled as fresh per-SCV TSV data.
    Bumping to iv:3 made that promotion newly catastrophic because the
    interpreter changed the data shape, not just the emit rules. NULL
    tag with no marker is now treated as unknown legacy and reingests.

    Returns True only when every clinvar row in the cache now carries
    the current tag.
    """
    if not db_path.exists():
        return False
    import contextlib

    from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

    tag = f"iv:{CLINVAR_INTERPRETER_VERSION}"
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        _ensure_local_version_tag_column(conn)
        try:
            rows = conn.execute(
                "SELECT name, remote_signal, local_version_tag "
                "FROM database_versions WHERE name LIKE 'clinvar%'"
            ).fetchall()
        except sqlite3.OperationalError:
            return False
        if not rows:
            return False
        stamped = False
        for name, sig, existing_tag in rows:
            if existing_tag == tag:
                continue
            if existing_tag is not None:
                # iv:N != iv:CURRENT — interpreter changed, reingest.
                return False
            # NULL tag: self-heal only if remote_signal carries a baked
            # `|iv:N` marker matching CURRENT. No marker or mismatched
            # N → unknown / wrong-version legacy → reingest.
            baked = _parse_baked_iv(sig)
            if baked != CLINVAR_INTERPRETER_VERSION:
                return False
            clean_signal = (sig or "").split("|iv:")[0]
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = ? WHERE name = ?",
                (clean_signal, tag, name),
            )
            stamped = True
        if stamped:
            conn.commit()
        return True


def _parse_baked_iv(remote_signal: str | None) -> int | None:
    """Extract the ``|iv:N`` suffix from a legacy ``remote_signal``, or None.

    Pre-v2.0.1 caches baked the interpreter version into the
    ``remote_signal`` column as ``"<sig>|iv:N"`` before
    ``local_version_tag`` existed. This helper recovers ``N`` so the
    caller can decide between self-heal and reingest. Returns None for
    missing or malformed markers (caller treats as unknown).
    """
    return parse_baked_version(remote_signal, "iv")


def parse_baked_version(remote_signal: str | None, kind: str) -> int | None:
    """Extract ``N`` from a legacy ``|<kind>:N`` suffix in ``remote_signal``, or None.

    Caches written before ``local_version_tag`` existed baked their
    version into ``remote_signal`` (``|iv:N`` for ClinVar and ClinPGx,
    ``|cv:N`` for the GWAS Catalog). None means no marker or a malformed
    one: the cache's version is unknown.
    """
    marker = f"|{kind}:"
    if not remote_signal or marker not in remote_signal:
        return None
    _, _, after = remote_signal.partition(marker)
    head = after.split("|", 1)[0].strip()
    if not head.isdigit():
        return None
    return int(head)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GRCh37 enrichment permits flagged rows with equal reference bases.

All input rows and SQLite caches are synthetic. The persisted flag is varied
independently to preserve existing translations while allowing equal-reference gains.
"""

from __future__ import annotations

import contextlib
import importlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.rsid_identity_loader import (
    RSID_IDENTITY_DB_FILENAME,
    IdentityRow,
    write_identity_db,
)
from allelix.databases.schema import (
    ALPHAMISSENSE_SCHEMA,
    CADD_SCHEMA,
    CLINVAR_SCHEMA,
    GNOMAD_SCHEMA,
)
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import AnalysisResult, run_analysis

if TYPE_CHECKING:
    from pathlib import Path


def _identity(data: Path, rows: list[IdentityRow]) -> None:
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2001-01-01T00:00:00Z",
    )


@pytest.mark.parametrize("flag", [0, 1], ids=["unflagged", "flagged"])
def test_three_allele_same_reference_translates_both_pairs(tmp_path: Path, flag: int) -> None:
    """Both non-palindromic pairs survive an over-flagged three-allele site."""
    _identity(
        tmp_path,
        [("rs900980001", "1", 1000, "C", alt, 2000, "C", alt, flag) for alt in ("A", "T")],
    )
    expected = {("1", 1000, "C", alt): ("1", 2000, "C", alt) for alt in ("A", "T")}
    with RsidIdentityAnnotator(tmp_path) as identity:
        assert identity.is_ready()
        assert identity.translate_to_grch38(set(expected)) == expected


def test_unflagged_three_allele_site_keeps_its_cg_translation(tmp_path: Path) -> None:
    """A C/G pair in an unflagged three-allele site retains its existing translation."""
    _identity(
        tmp_path,
        [("rs900980003", "1", 1000, "C", alt, 2000, "C", alt, 0) for alt in ("A", "G")],
    )
    expected = {("1", 1000, "C", alt): ("1", 2000, "C", alt) for alt in ("A", "G")}
    with RsidIdentityAnnotator(tmp_path) as identity:
        assert identity.is_ready()
        assert identity.translate_to_grch38(set(expected)) == expected


_TRANSLATION_CASES = [
    pytest.param(0, "C", ("A", "T"), "C", ("A", "T"), True, id="unflagged-three-alleles"),
    pytest.param(1, "C", ("A", "T"), "C", ("A", "T"), True, id="flagged-equal-references"),
    pytest.param(1, "C", ("A",), "G", ("T",), False, id="flagged-complementary-references"),
    pytest.param(1, "C", ("A",), "A", ("C",), False, id="flagged-other-different-references"),
    pytest.param(0, "C", ("A",), "G", ("T",), True, id="unflagged-inverted-complement"),
    pytest.param(0, "A", ("T",), "A", ("T",), True, id="unflagged-AT-retained"),
    pytest.param(0, "C", ("G",), "C", ("G",), True, id="unflagged-CG-retained"),
    pytest.param(1, "A", ("T",), "A", ("T",), True, id="flagged-equal-reference-AT-gained"),
    pytest.param(1, "C", ("G",), "C", ("G",), True, id="flagged-equal-reference-CG-gained"),
    pytest.param(0, "C", ("A",), "T", ("A",), True, id="GRCh38-only-AT-retained"),
    pytest.param(1, "C", ("A",), "T", ("A",), False, id="GRCh38-only-AT-flagged-abstains"),
    pytest.param(0, "A", ("A",), "G", ("A",), True, id="reference-carried-retained"),
    pytest.param(0, "A", ("G",), "G", ("G",), False, id="target-reference-still-excluded"),
    pytest.param(1, "C", ("C",), "C", ("C",), False, id="flagged-target-reference-excluded"),
]


@pytest.mark.parametrize(
    ("flag", "ref37", "alts37", "ref38", "alts38", "translate"), _TRANSLATION_CASES
)
def test_translation_changes_only_the_equal_reference_flagged_veto(
    tmp_path: Path,
    flag: int,
    ref37: str,
    alts37: tuple[str, ...],
    ref38: str,
    alts38: tuple[str, ...],
    translate: bool,
) -> None:
    """Keep old translations and abstentions except for the equal-reference flag veto."""
    pairs = list(zip(alts37, alts38, strict=True))
    _identity(
        tmp_path,
        [("rs900980002", "1", 1000, ref37, a, 2000, ref38, b, flag) for a, b in pairs],
    )
    keys = {("1", 1000, ref37, a) for a, _b in pairs}
    expected = (
        {("1", 1000, ref37, a): ("1", 2000, ref38, b) for a, b in pairs} if translate else {}
    )
    with RsidIdentityAnnotator(tmp_path) as identity:
        assert identity.is_ready()
        assert identity.translate_to_grch38(keys) == expected


@pytest.mark.parametrize(
    ("flag", "ref37", "alts37", "ref38", "alts38", "translate"), _TRANSLATION_CASES
)
def test_wrong_allele_checker_uses_the_same_reference_pair_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    pytestconfig: pytest.Config,
    flag: int,
    ref37: str,
    alts37: tuple[str, ...],
    ref38: str,
    alts38: tuple[str, ...],
    translate: bool,
) -> None:
    """The independent report checker follows the same translation-preserving contract."""
    monkeypatch.syspath_prepend(str(pytestconfig.rootpath / "scripts"))
    checker = importlib.import_module("check_wrong_allele")
    pairs = list(zip(alts37, alts38, strict=True))
    _identity(
        tmp_path,
        [("rs900980004", "1", 1000, ref37, a, 2000, ref38, b, flag) for a, b in pairs],
    )
    with contextlib.closing(checker.Caches(tmp_path, "GRCh37", {"rsid_identity"})) as caches:
        for alt37, alt38 in pairs:
            expected = ("1", 2000, ref38, alt38) if translate else None
            assert caches.to38(("1", 1000, ref37, alt37)) == expected


def _stamp(conn: sqlite3.Connection, name: str, count: int, tag: str = "sv:1") -> None:
    conn.execute(
        "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
        "record_count, local_version_tag) VALUES (?, 'test://cache', '1', '2001-01-01', ?, ?)",
        (name, count, tag),
    )


def _run(
    tmp_path: Path,
    rows: list[IdentityRow],
    frequencies: list[float | None],
    *,
    with_identity: bool = True,
) -> AnalysisResult:
    """Analyze homozygous ALT calls, using the first row at each synthetic site."""
    data = tmp_path / "data"
    data.mkdir()
    if with_identity:
        _identity(data, rows)
    sites = list({row[0]: row for row in reversed(rows)}.values())
    sites.sort(key=lambda row: row[2])
    assert len(sites) == len(frequencies)
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total) VALUES "
            "(?, ?, ?, ?, ?, 'Pathogenic', 'Synthetic condition', 'GENE', "
            "'criteria_provided,_single_submitter', ?, 'Pathogenic', 1, 1)",
            [(*row[:5], i) for i, row in enumerate(sites, 1)],
        )
        _stamp(
            conn, clinvar_record_name("GRCh37"), len(sites), f"iv:{CLINVAR_INTERPRETER_VERSION}"
        )
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        # An unrelated row keeps the cache ready even when every requested AF is absent.
        freqs = [("2", 9999, "G", "A", "rs900980099", 0.4)]
        freqs.extend(
            (row[1], row[5], row[6], row[7], row[0], af)
            for row, af in zip(sites, frequencies, strict=True)
            if af is not None
        )
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) VALUES (?,?,?,?,?,?)",
            freqs,
        )
        _stamp(conn, "gnomad", len(freqs))
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "alphamissense.sqlite")) as conn:
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.executemany(
            "INSERT INTO alphamissense_scores VALUES "
            "(?, ?, ?, ?, ?, 'P1', 'T1', 'A1V', 0.91, 'likely_pathogenic')",
            [(row[1], row[5], row[6], row[7], row[0]) for row in sites],
        )
        _stamp(conn, "alphamissense", len(sites))
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / "cadd.sqlite")) as conn:
        conn.executescript(CADD_SCHEMA)
        conn.executemany(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, 23.5)",
            [(row[1], row[5], row[6], row[7]) for row in sites],
        )
        _stamp(conn, "cadd", len(sites))
        conn.commit()
    chip = tmp_path / "synthetic.txt"
    chip.write_text(
        "# MyHappyGenes [TEMPUS]\n# Sample ID\tPLACEHOLDER\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        + "".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[4]}\t{r[4]}\n" for r in sites)
    )
    with contextlib.ExitStack() as stack:
        clinvar = stack.enter_context(ClinVarAnnotator(data, builds=("GRCh37",)))
        gnomad = stack.enter_context(GnomadAnnotator(data))
        am = stack.enter_context(AlphaMissenseAnnotator(data))
        cadd = stack.enter_context(CaddAnnotator(data))
        identity = stack.enter_context(RsidIdentityAnnotator(data)) if with_identity else None
        assert all(a.is_ready() for a in (clinvar, gnomad, am, cadd))
        if identity is not None:
            assert identity.is_ready()
        result = run_analysis(
            chip,
            MyHappyGenesParser(),
            [clinvar],
            build_override="GRCh37",
            gnomad=gnomad,
            alphamissense=am,
            cadd=cadd,
            identity=identity,
        )
    assert len(result.annotations) == len(sites)
    return result


def _three_alleles(flag: int) -> list[IdentityRow]:
    return [("rs900980001", "1", 1000, "C", alt, 2000, "C", alt, flag) for alt in ("A", "T")]


@pytest.mark.parametrize(
    ("field", "expected", "source"),
    [
        ("allele_frequency", 0.2, "gnomad"),
        ("am_pathogenicity", 0.91, "alphamissense"),
        ("cadd_phred", 23.5, "cadd"),
    ],
)
def test_flagged_three_allele_site_reaches_each_enricher(
    tmp_path: Path, field: str, expected: float, source: str
) -> None:
    """The shared translation path reaches each ready GRCh38 enrichment source."""
    result = _run(tmp_path, _three_alleles(1), [0.2])
    assert getattr(result.annotations[0], field) == pytest.approx(expected)
    assert source in dict(result.annotators_used)


@pytest.mark.parametrize("flag", [0, 1], ids=["unflagged", "flagged"])
def test_ready_gnomad_without_the_translated_row_reports_no_frequency(
    tmp_path: Path, flag: int
) -> None:
    """A completed lookup without a row differs from gnomAD not being run."""
    result = _run(tmp_path, _three_alleles(flag), [None])
    annotation = result.annotations[0]
    assert annotation.allele_frequency is None
    assert annotation.plausibility == "not_assessed"
    assert annotation.plausibility_reason == "no_frequency"
    assert "gnomad" in dict(result.annotators_used)


@pytest.mark.parametrize(
    ("frequencies", "expected_count"),
    [
        ([None, None, None], 3),
        ([0.000001, 0.000001, None], 3),
        ([0.000001, 0.000001, 0.000001], 3),
        ([0.2, 0.2, 0.2], 0),
    ],
    ids=[
        "three-looked-up-absent",
        "two-rare-one-looked-up-absent",
        "three-known-rare",
        "three-common",
    ],
)
def test_multiple_homozygous_plp_counts_rare_or_looked_up_absent_alleles(
    tmp_path: Path, frequencies: list[float | None], expected_count: int
) -> None:
    """A translated allele absent from a queried gnomAD cache still counts as rare."""
    rows = [(f"rs90098001{i}", "1", 1000 + i, "G", "A", 2000 + i, "G", "A", 0) for i in range(3)]
    result = _run(tmp_path, rows, frequencies)
    annotations = sorted(result.annotations, key=lambda a: a.rsid)
    assert [a.allele_frequency for a in annotations] == frequencies
    assert all(a.clinvar_aggregate == "Pathogenic" for a in annotations)
    assert all(a.zygosity_code == "homozygous" for a in annotations)
    assert "gnomad" in dict(result.annotators_used)
    for annotation, af in zip(annotations, frequencies, strict=True):
        if af is None:
            assert annotation.plausibility_reason == "no_frequency"
    warnings = [w for w in result.report_warnings if w.code == "multiple_homozygous_plp"]
    assert [w.count for w in warnings] == ([expected_count] if expected_count else [])


@pytest.mark.parametrize(
    ("reason", "mixed"),
    [
        ("complementary-references", False),
        ("different-references", False),
        ("no-identity", False),
        ("complementary-references", True),
        ("different-references", True),
    ],
    ids=[
        "all-complementary",
        "all-different-references",
        "no-identity",
        "two-rare-complementary",
        "two-rare-different-references",
    ],
)
def test_multiple_homozygous_plp_excludes_alleles_not_assessed(
    tmp_path: Path, reason: str, mixed: bool
) -> None:
    """Untranslated calls cannot count, even when other calls were looked up in gnomAD."""
    ref37, alt37, ref38, alt38 = {
        "complementary-references": ("C", "A", "G", "T"),
        "different-references": ("C", "A", "A", "C"),
        "no-identity": ("G", "A", "G", "A"),
    }[reason]
    rows = []
    frequencies = []
    for i in range(3):
        if mixed and i < 2:
            rows.append((f"rs90098001{i}", "1", 1000 + i, "G", "A", 2000 + i, "G", "A", 0))
            frequencies.append(0.000001)
        else:
            rows.append((f"rs90098001{i}", "1", 1000 + i, ref37, alt37, 2000 + i, ref38, alt38, 1))
            # Common in the target cache, but this call cannot safely be translated there.
            frequencies.append(0.2)
    result = _run(tmp_path, rows, frequencies, with_identity=reason != "no-identity")
    annotations = sorted(result.annotations, key=lambda a: a.rsid)
    assert all(a.clinvar_aggregate == "Pathogenic" for a in annotations)
    assert all(a.zygosity_code == "homozygous" for a in annotations)
    assert [a.allele_frequency for a in annotations] == (
        [0.000001, 0.000001, None] if mixed else [None, None, None]
    )
    unassessed = annotations[2:] if mixed else annotations
    assert all(a.plausibility == "not_assessed" for a in unassessed)
    assert all(a.plausibility_reason != "no_frequency" for a in unassessed)
    assert ("gnomad" in dict(result.annotators_used)) == mixed
    assert not any(w.code == "multiple_homozygous_plp" for w in result.report_warnings)

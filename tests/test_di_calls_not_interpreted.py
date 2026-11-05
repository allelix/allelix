# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A chip D/I call is not interpreted by any source (the owner's D/I decision).

An array writes a deletion or insertion as ``D`` or ``I``. No per-probe
specification available to Allelix establishes which allele each letter
stands for, so
ClinVar does not match such a call, and no other source may either:
not SNPedia, whose pages can be written in the same letters (``(D;I)``),
whatever orientation the page states, at an rsID or a vendor i-ID; not the
GWAS Catalog, whose rows without a risk allele fire on the rsID alone; and
not ClinPGx. The call is still counted (``array_quality.di_calls``), and at
a high-value or ClinVar Pathogenic rsID it is named as not checked.

Unchanged: a base call at an rsID whose SNPedia page also has D/I or "-"
rows is matched exactly as before.

The same answer holds across the pairs the mechanism touches: rsID or
i-ID, page orientation plus, minus or unstated, single-call or batch path,
with or without the rsID identity table, chip file or rsID and genotype
list, SNPedia on or off.

All rsIDs, i-IDs, loci and rows are synthetic, apart from the public
high-value rsID BRCA1 rs80357906 (5382insC), which is used only as a
high-value key with synthetic rows.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.cli import main
from allelix.databases._versions import PHARMGKB_INTERPRETER_VERSION
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION, GWAS_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GWAS_SCHEMA, PHARMGKB_SCHEMA
from allelix.databases.snpedia_parser import _PARSER_VERSION
from allelix.models import Variant

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

# rsIDs in the unused 900970xxx range.
DI_PLUS = "rs900970001"  # SNPedia (D;I) rows, page orientation plus
DI_MINUS = "rs900970002"  # SNPedia (D;I) rows, page orientation minus
DI_UNSTATED = "rs900970003"  # SNPedia (D;I) rows, orientation not stated
MIXED = "rs900970004"  # SNPedia D/I rows and base rows on one page
DASH = "rs900970005"  # SNPedia rows written with "-" and bases
GWAS_ONLY = "rs900970006"  # GWAS rows only (one without a risk allele)
IPROBE = "i900970007"  # a vendor i-ID with (D;I) rows
HIGH_VALUE = "rs80357906"  # BRCA1 5382insC: high-value, synthetic rows only
FILLER = "rs900970099"
FILLER2 = "rs900970098"

# (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic)
_IDENTITY = (
    (DI_PLUS, "1", 1000, "G", "A", 101000, "G", "A", 0),
    (DI_MINUS, "1", 2000, "G", "A", 102000, "G", "A", 0),
    (DI_UNSTATED, "1", 3000, "G", "A", 103000, "G", "A", 0),
    (MIXED, "1", 4000, "C", "T", 104000, "C", "T", 0),
    (DASH, "1", 5000, "C", "T", 105000, "C", "T", 0),
    (GWAS_ONLY, "1", 6000, "C", "T", 106000, "C", "T", 0),
    (HIGH_VALUE, "17", 7000, "G", "A", 107000, "G", "A", 0),
    (FILLER, "2", 9000, "A", "G", 109000, "A", "G", 0),
    (FILLER2, "2", 9100, "A", "G", 109100, "A", "G", 0),
)
_POS37 = {row[0]: (row[1], row[2]) for row in _IDENTITY}
_POS37[IPROBE] = ("1", 7700)

# (rsid, allele1, allele2, magnitude, summary, orientation)
_SNPEDIA_ROWS = (
    (DI_PLUS, "D", "D", 3.0, "plus page DD", "plus"),
    (DI_PLUS, "D", "I", 3.0, "plus page DI", "plus"),
    (DI_PLUS, "I", "I", 3.0, "plus page II", "plus"),
    (DI_MINUS, "D", "D", 3.0, "minus page DD", "minus"),
    (DI_MINUS, "D", "I", 3.0, "minus page DI", "minus"),
    (DI_MINUS, "I", "I", 3.0, "minus page II", "minus"),
    (DI_UNSTATED, "D", "D", 3.0, "unstated page DD", None),
    (DI_UNSTATED, "D", "I", 3.0, "unstated page DI", None),
    (DI_UNSTATED, "I", "I", 3.0, "unstated page II", None),
    (MIXED, "D", "I", 3.0, "mixed page DI", "plus"),
    (MIXED, "C", "T", 3.0, "mixed page CT", "plus"),
    (DASH, "-", "C", 3.0, "dash page -C", "plus"),
    (DASH, "C", "C", 3.0, "dash page CC", "plus"),
    (DASH, "C", "T", 3.0, "dash page CT", "plus"),
    (IPROBE, "D", "I", 3.0, "i-ID page DI", "plus"),
    (IPROBE, "D", "D", 3.0, "i-ID page DD", "plus"),
    (HIGH_VALUE, "D", "I", 3.0, "high-value page DI", "plus"),
)

# (rsid, risk_allele, trait, p_value, or_beta, ci_text)
_GWAS_ROWS = (
    (GWAS_ONLY, None, "Trait without a risk allele", 1e-30, 1.5, "[1.2-1.8]"),
    (GWAS_ONLY, "T", "Trait with risk allele T", 1e-30, 1.5, "[1.2-1.8]"),
    (DI_PLUS, None, "Trait without a risk allele at a D/I site", 1e-30, 1.5, "[1.2-1.8]"),
)

_DI_GENOTYPES = ("DD", "DI", "II")
_PAGES = (DI_PLUS, DI_MINUS, DI_UNSTATED)


def _write_snpedia(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.executescript(
            "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
            "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
            "orientation TEXT, scraped_at TEXT);"
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
            "version TEXT, downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
            "remote_signal TEXT, local_version_tag TEXT);"
        )
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, 'Bad', ?, 'GENE', ?, "
            "'2026-10-10T00:00:00')",
            _SNPEDIA_ROWS,
        )
        conn.execute(
            "INSERT INTO database_versions VALUES ('snpedia', 'test://snpedia', 'test', "
            "'2026-10-10T00:00:00Z', ?, '', ?)",
            (len(_SNPEDIA_ROWS), f"pv:{_PARSER_VERSION}"),
        )
        conn.commit()


def _write_gwas(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / GWAS_DB_FILENAME)) as conn:
        conn.executescript(GWAS_SCHEMA)
        conn.executemany(
            "INSERT INTO gwas_associations (rsid, risk_allele, trait, p_value, or_beta, "
            "ci_text, gene, study_accession, pubmed_id, trait_category) "
            "VALUES (?, ?, ?, ?, ?, ?, 'GENE', 'GCST0', '1', 'disease')",
            _GWAS_ROWS,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES ('gwas', 'test://gwas', 'test', '2026-10-16T00:00:00Z', ?, '', ?)",
            (len(_GWAS_ROWS), f"cv:{_CATEGORIZER_VERSION}"),
        )
        conn.commit()


def _write_pharm(data: Path) -> None:
    # Rows a real build never stores (its loader keeps two-base genotypes
    # only), written here so the guard does not lean on the loader alone.
    with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
        conn.executescript(PHARMGKB_SCHEMA)
        for i, genotype in enumerate(("DI", "DD", "II", "CT")):
            conn.execute(
                "INSERT INTO pharmgkb_annotations (rsid, genotype, gene, drugs, phenotype, "
                "phenotype_category, annotation_text, level_of_evidence, score, "
                "pgkb_annotation_id, function_class, is_nonfinding) "
                "VALUES (?, ?, 'GENE', 'Synthetic drug', 'Synthetic response', 'Efficacy', "
                "'Synthetic finding', '1A', 10, ?, 'decreased', 0)",
                (MIXED, genotype, str(900970800 + i)),
            )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES ('pharmgkb', 'test://pharmgkb', 'test', '2026-10-10T00:00:00Z', 4, '', ?)",
            (f"iv:{PHARMGKB_INTERPRETER_VERSION}",),
        )
        conn.commit()


@pytest.fixture
def data(tmp_path: Path) -> Path:
    d = tmp_path / "data"
    d.mkdir()
    write_identity_db(
        _IDENTITY,
        d / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2026-10-16T00:00:00Z",
    )
    _write_snpedia(d)
    _write_gwas(d)
    _write_pharm(d)
    return d


@pytest.fixture
def identity(data: Path) -> Iterator[RsidIdentityAnnotator]:
    with RsidIdentityAnnotator(data) as annotator:
        yield annotator


def _call(rsid: str, genotype: str, build: str = "GRCh37") -> Variant:
    chrom, pos = _POS37[rsid]
    if build == "GRCh38" and rsid.startswith("rs"):
        pos += 100000
    return Variant(rsid, chrom, pos, genotype[0], genotype[1], build=build)


def _snpedia(data: Path, identity: RsidIdentityAnnotator | None) -> SNPediaAnnotator:
    return SNPediaAnnotator(
        data,
        forward_site_provider=identity.forward_sites if identity is not None else None,
    )


def _both_paths(annotator, variant: Variant) -> tuple[list, list]:
    with annotator:
        return annotator.annotate(variant), list(annotator.batch_annotate([variant]))


# ── SNPedia ─────────────────────────────────────────────────────


@pytest.mark.parametrize("with_identity", (True, False))
@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize("genotype", _DI_GENOTYPES)
@pytest.mark.parametrize("rsid", _PAGES)
def test_snpedia_gives_no_row_for_a_di_call_whatever_the_page_orientation(
    data, identity, rsid, genotype, build, with_identity
):
    single, batch = _both_paths(
        _snpedia(data, identity if with_identity else None), _call(rsid, genotype, build)
    )
    assert (single, batch) == ([], [])


@pytest.mark.parametrize("genotype", ("DI", "DD"))
def test_snpedia_gives_no_row_for_a_di_call_at_a_vendor_i_id(data, identity, genotype):
    single, batch = _both_paths(_snpedia(data, identity), _call(IPROBE, genotype))
    assert (single, batch) == ([], [])


def test_snpedia_gives_no_row_for_a_di_call_on_a_page_that_also_has_base_rows(data, identity):
    single, batch = _both_paths(_snpedia(data, identity), _call(MIXED, "DI"))
    assert (single, batch) == ([], [])


@pytest.mark.parametrize("with_identity", (True, False))
def test_a_base_call_on_a_page_with_di_rows_is_matched_as_before(data, identity, with_identity):
    single, batch = _both_paths(
        _snpedia(data, identity if with_identity else None), _call(MIXED, "CT")
    )
    assert [a.description for a in single] == [a.description for a in batch]
    assert len(single) == 1
    assert "mixed page CT" in single[0].description


@pytest.mark.parametrize(("genotype", "summary"), (("CC", "dash page CC"), ("CT", "dash page CT")))
def test_a_base_call_on_a_page_written_with_dashes_is_matched_as_before(
    data, identity, genotype, summary
):
    single, batch = _both_paths(_snpedia(data, identity), _call(DASH, genotype))
    assert [a.description for a in single] == [a.description for a in batch]
    assert len(single) == 1
    assert summary in single[0].description


# ── GWAS Catalog ────────────────────────────────────────────────


@pytest.mark.parametrize("with_identity", (True, False))
@pytest.mark.parametrize("genotype", _DI_GENOTYPES)
@pytest.mark.parametrize("rsid", (GWAS_ONLY, DI_PLUS))
def test_gwas_gives_no_row_for_a_di_call_even_without_a_risk_allele(
    data, identity, rsid, genotype, with_identity
):
    gwas = GWASCatalogAnnotator(
        data, forward_site_provider=identity.forward_sites if with_identity else None
    )
    single, batch = _both_paths(gwas, _call(rsid, genotype))
    assert (single, batch) == ([], [])


def test_gwas_row_without_a_risk_allele_still_fires_for_a_base_call(data, identity):
    gwas = GWASCatalogAnnotator(data, forward_site_provider=identity.forward_sites)
    single, batch = _both_paths(gwas, _call(GWAS_ONLY, "CC"))
    traits = sorted(a.trait for a in single)
    assert traits == sorted(a.trait for a in batch)
    assert traits == ["Trait without a risk allele"]


# ── ClinPGx (class inventory: already unaffected) ───────────────


@pytest.mark.parametrize("genotype", _DI_GENOTYPES)
def test_clinpgx_gives_no_row_for_a_di_call(data, genotype):
    single, batch = _both_paths(PharmGKBAnnotator(data), _call(MIXED, genotype))
    assert (single, batch) == ([], [])


# ── The report: counted, named where high-value, no source row ──


def _array(path: Path, calls: list[tuple[str, str]]) -> Path:
    lines = [
        "# This data file generated by 23andMe at: Mon Jan 01 00:00:00 2001",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    lines += [f"{rsid}\t{_POS37[rsid][0]}\t{_POS37[rsid][1]}\t{gt}" for rsid, gt in calls]
    path.write_text("\n".join(lines) + "\n")
    return path


def _rsid_list(path: Path, calls: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"{rsid}\t{gt}\n" for rsid, gt in calls))
    return path


_FLOORS = ("--min-magnitude", "0", "--gwas-min-magnitude", "0", "--snpedia-min-magnitude", "0")


def _analyze_json(data: Path, input_path: Path, out: Path, *extra: str) -> dict:
    result = CliRunner(env={"COLUMNS": "300"}).invoke(
        main,
        [
            "analyze",
            str(input_path),
            "--no-update",
            "--data-dir",
            str(data),
            "--output",
            str(out),
            *_FLOORS,
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output
    return json.loads(out.read_text())


def _di_rows(report: dict, di_rsids: set[str]) -> list[tuple[str, str]]:
    return sorted((a["source"], a["rsid"]) for a in report["annotations"] if a["rsid"] in di_rsids)


_CHIP_CALLS = [
    (DI_PLUS, "DI"),
    (DI_MINUS, "DD"),
    (DI_UNSTATED, "II"),
    (MIXED, "DI"),
    (GWAS_ONLY, "DI"),
    (IPROBE, "DI"),
    (HIGH_VALUE, "DI"),
    (FILLER, "AG"),
    (FILLER2, "AG"),
]
_DI_RSIDS = {rsid for rsid, gt in _CHIP_CALLS if set(gt) & {"D", "I"}}


def test_chip_report_has_no_source_row_for_any_di_call(tmp_path, data):
    report = _analyze_json(
        data, _array(tmp_path / "a.txt", _CHIP_CALLS), tmp_path / "r.json", "--build", "grch37"
    )
    assert _di_rows(report, _DI_RSIDS) == []


def test_chip_report_still_counts_every_di_call(tmp_path, data):
    report = _analyze_json(
        data, _array(tmp_path / "a.txt", _CHIP_CALLS), tmp_path / "r.json", "--build", "grch37"
    )
    assert report["input"]["array_quality"]["di_calls"] == len(_DI_RSIDS)


def test_a_di_call_at_a_high_value_rsid_is_named_not_checked_and_has_no_snpedia_row(
    tmp_path, data
):
    report = _analyze_json(
        data, _array(tmp_path / "a.txt", _CHIP_CALLS), tmp_path / "r.json", "--build", "grch37"
    )
    named = [s for s in report["high_value_unchecked"] if s["rsid"] == HIGH_VALUE]
    assert len(named) == 1
    assert (named[0]["reason"], named[0]["call"], named[0]["gene"]) == ("di_call", "D/I", "BRCA1")
    assert "cannot check against its sources" in named[0]["message"]
    warnings = [w for w in report["report_warnings"] if w["code"] == "di_calls_not_checked"]
    assert [w["count"] for w in warnings] == [1]
    assert not [a for a in report["annotations"] if a["rsid"] == HIGH_VALUE]


def test_switching_snpedia_off_changes_nothing_for_di_calls(tmp_path, data):
    array = _array(tmp_path / "a.txt", _CHIP_CALLS)
    on = _analyze_json(data, array, tmp_path / "on.json", "--build", "grch37")
    off = _analyze_json(
        data, array, tmp_path / "off.json", "--build", "grch37", "--exclude-snpedia"
    )
    assert _di_rows(on, _DI_RSIDS) == _di_rows(off, _DI_RSIDS) == []
    assert on["input"]["array_quality"] == off["input"]["array_quality"]
    assert on["high_value_unchecked"] == off["high_value_unchecked"]


def test_an_rsid_list_di_row_gets_no_source_row(tmp_path, data):
    calls = [(DI_PLUS, "DI"), (DI_UNSTATED, "DD"), (GWAS_ONLY, "II"), (FILLER, "AG")]
    report = _analyze_json(data, _rsid_list(tmp_path / "list.txt", calls), tmp_path / "r.json")
    assert _di_rows(report, {DI_PLUS, DI_UNSTATED, GWAS_ONLY}) == []
    assert report["input"]["placement"]["unplaced"].get("di_call") == 3


def test_chip_and_rsid_list_give_the_same_answer_for_a_base_call_beside_di_rows(tmp_path, data):
    """Cross-combination: the fix does not move a base call's rows in either input."""
    chip = _analyze_json(
        data,
        _array(tmp_path / "a.txt", [(MIXED, "CT"), (DASH, "CC"), (FILLER, "AG")]),
        tmp_path / "chip.json",
        "--build",
        "grch37",
    )
    listed = _analyze_json(
        data,
        _rsid_list(tmp_path / "l.txt", [(MIXED, "CT"), (DASH, "CC"), (FILLER, "AG")]),
        tmp_path / "list.json",
    )

    def rows(report: dict) -> list[tuple[str, str, str]]:
        return sorted(
            (a["source"], a["rsid"], a["genotype_match"])
            for a in report["annotations"]
            if a["source"] == "snpedia"
        )

    assert rows(chip) == rows(listed)
    assert {(rsid, gt) for _, rsid, gt in rows(chip)} == {(MIXED, "CT"), (DASH, "CC")}

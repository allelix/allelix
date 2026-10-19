# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GWAS and SNPedia alleles read on the forward strand; GWAS effect direction and floor.

Source alleles are read against the rsID identity table's forward alleles
(abstain-first, ADR-0024). Every GWAS row says which way its effect allele
points. The GWAS floor is max(--min-magnitude, --gwas-min-magnitude);
must-include rsIDs keep their exemption.

All rsIDs, loci and rows here are synthetic.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.gwas import GWASCatalogAnnotator, _effect_direction
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.gwas_loader import _CATEGORIZER_VERSION, GWAS_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GWAS_SCHEMA
from allelix.models import Annotation, Variant
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.association_text import (
    effect_note,
    effect_words,
    strand_abstention_lines,
)
from allelix.utils.allele import (
    STRAND_COMPLEMENT,
    STRAND_DIRECT,
    STRAND_NO_FIT,
    STRAND_PALINDROMIC,
    forward_strand_reading,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

# (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic)
_IDENTITY = (
    # G/A site.
    ("rs900159001", "1", 1000, "G", "A", 101000, "G", "A", 0),
    # A/T site.
    ("rs900159002", "1", 2000, "A", "T", 102000, "A", "T", 1),
    # dbSNP lists every base (shaped like F5 rs6025 on GRCh37).
    ("rs900159003", "1", 3000, "T", "A", 103000, "C", "A", 1),
    ("rs900159003", "1", 3000, "T", "C", 103000, "C", "C", 1),
    ("rs900159003", "1", 3000, "T", "G", 103000, "C", "G", 1),
    ("rs900159003", "1", 3000, "T", "T", 103000, "C", "T", 1),
    # C/T site.
    ("rs900159004", "1", 4000, "C", "T", 104000, "C", "T", 0),
    # C/G site.
    ("rs900159005", "1", 5000, "C", "G", 105000, "C", "G", 1),
)

_OR = "[1.2-1.8]"

# (rsid, risk_allele, trait, p_value, or_beta, ci_text)
_GWAS_ROWS = (
    ("rs900159001", "A", "Trait read directly", 1e-30, 1.5, _OR),
    ("rs900159001", "T", "Trait read complemented", 1e-30, 1.5, _OR),
    ("rs900159001", "C", "Trait with C read as G", 1e-30, 1.5, _OR),
    ("rs900159001", None, "Trait without a risk allele", 1e-30, 1.5, _OR),
    ("rs900159002", "A", "Trait at an A/T site", 1e-30, 1.5, _OR),
    ("rs900159003", "T", "Trait at a four-base site", 1e-30, 1.5, _OR),
    ("rs900159099", "A", "Trait without identity", 1e-30, 1.5, _OR),
    ("rs900159004", "C", "Trait with a protective reference allele", 1e-30, 0.30, "[0.2-0.4]"),
)


def _write_gwas(path: Path) -> None:
    """A GWAS cache in the loader's schema, stamped current, so ``analyze`` uses it."""
    with contextlib.closing(sqlite3.connect(path)) as conn:
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


# (rsid, allele1, allele2, magnitude, summary, orientation)
_SNPEDIA_ROWS = (
    # One homozygous row, orientation unknown, at a C/T site: read complemented.
    ("rs900159004", "G", "G", 3.0, "one row read complemented", None),
    # One homozygous row, orientation unknown, at an A/T site: abstain.
    ("rs900159002", "A", "A", 3.0, "one row at an A/T site", None),
    # Stated orientation at a C/G site: abstain.
    ("rs900159005", "C", "C", 3.0, "plus page at a C/G site", "plus"),
    # Minus page at the four-base site: (A;A) is forward TT.
    ("rs900159003", "A", "A", 7.0, "minus page at a four-base site", "minus"),
    # Plus page whose genotype is not among the forward alleles: abstain.
    ("rs900159001", "C", "C", 3.0, "plus page off the forward alleles", "plus"),
    # A 23andMe i-ID: not an rsID, unchanged.
    ("i900159006", "A", "A", 3.0, "i-ID row", "plus"),
)


def _write_snpedia(path: Path) -> None:
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.execute(
            "CREATE TABLE snpedia_genotypes (rsid TEXT, allele1 TEXT, allele2 TEXT, "
            "magnitude REAL, repute TEXT, summary TEXT, gene TEXT, orientation TEXT)"
        )
        conn.executemany(
            "INSERT INTO snpedia_genotypes VALUES (?, ?, ?, ?, 'Bad', ?, 'GENE', ?)",
            _SNPEDIA_ROWS,
        )
        conn.commit()


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    write_identity_db(
        _IDENTITY,
        tmp_path / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2026-10-16T00:00:00Z",
    )
    _write_gwas(tmp_path / GWAS_DB_FILENAME)
    _write_snpedia(tmp_path / "snpedia.sqlite")
    return tmp_path


@pytest.fixture
def identity(data_dir: Path) -> Iterator[RsidIdentityAnnotator]:
    with RsidIdentityAnnotator(data_dir) as annotator:
        yield annotator


def _call(rsid: str, pos: int, a1: str, a2: str, ref: str, build: str = "GRCh37") -> Variant:
    return Variant(rsid, "1", pos, a1, a2, build=build, ref=ref)


_GWAS_CALLS = (
    _call("rs900159001", 1000, "A", "G", "G"),
    _call("rs900159002", 2000, "A", "A", "A"),
    _call("rs900159003", 3000, "T", "T", "T"),
    _call("rs900159099", 9900, "A", "A", "A"),
)


def _gwas(data_dir: Path, identity: RsidIdentityAnnotator | None) -> GWASCatalogAnnotator:
    return GWASCatalogAnnotator(
        data_dir,
        forward_site_provider=identity.forward_sites if identity is not None else None,
    )


# ── The reading ─────────────────────────────────────────────────


@pytest.mark.parametrize(
    ("alleles", "forward", "reading"),
    [
        ({"A"}, {"G", "A"}, STRAND_DIRECT),
        ({"T"}, {"G", "A"}, STRAND_COMPLEMENT),
        # At a two-allele site that is not A/T or C/G every base fits one reading.
        ({"C"}, {"G", "A"}, STRAND_COMPLEMENT),
        ({"A"}, {"A", "T"}, STRAND_PALINDROMIC),
        ({"G"}, {"C", "G"}, STRAND_PALINDROMIC),
        # Rule 1 comes first where the complement is forward too.
        ({"T"}, {"T", "A", "C", "G"}, STRAND_DIRECT),
        ({"G", "A"}, {"C", "T"}, STRAND_COMPLEMENT),
        ({"G", "C"}, {"C", "T"}, STRAND_NO_FIT),
        ({"D"}, {"C", "T"}, STRAND_NO_FIT),
    ],
)
def test_forward_strand_reading(alleles: set[str], forward: set[str], reading: str) -> None:
    assert forward_strand_reading(alleles, frozenset(forward)) == reading


def test_forward_sites_per_build(identity: RsidIdentityAnnotator) -> None:
    sites = identity.forward_sites({"rs900159003", "rs900159001", "rs900159099"}, "GRCh38")
    assert sites == {
        "rs900159001": ("1", 101000, frozenset({"G", "A"})),
        "rs900159003": ("1", 103000, frozenset({"A", "C", "G", "T"})),
    }


# ── GWAS ────────────────────────────────────────────────────────


def test_gwas_reads_risk_alleles_on_the_forward_strand(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    with _gwas(data_dir, identity) as gwas:
        rows = list(gwas.batch_annotate(_GWAS_CALLS))
        counts = gwas.run_counts()
    by_trait = {a.trait: a for a in rows}
    assert set(by_trait) == {
        "Trait read directly",
        "Trait read complemented",
        "Trait with C read as G",
        "Trait without a risk allele",
        "Trait at a four-base site",
    }
    assert by_trait["Trait read complemented"].alt == "A"
    assert by_trait["Trait with C read as G"].alt == "G"
    assert by_trait["Trait at a four-base site"].alt == "T"
    assert counts == {"not_in_identity_table": 1, "palindromic": 1, "no_fit": 0}


def test_gwas_single_and_batch_paths_agree(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    with _gwas(data_dir, identity) as single, _gwas(data_dir, identity) as batch:
        one_by_one = [a for v in _GWAS_CALLS for a in single.annotate(v)]
        batched = list(batch.batch_annotate(_GWAS_CALLS))
        assert one_by_one == batched
        assert single.run_counts() == batch.run_counts()


def test_gwas_abstains_when_the_call_is_elsewhere(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    misplaced = _call("rs900159001", 1001, "A", "G", "G")
    with _gwas(data_dir, identity) as gwas:
        traits = {a.trait for a in gwas.annotate(misplaced)}
        assert traits == {"Trait without a risk allele"}
        assert gwas.run_counts()["not_in_identity_table"] == 3


def test_gwas_without_an_identity_cache_abstains(tmp_path: Path) -> None:
    _write_gwas(tmp_path / GWAS_DB_FILENAME)
    with RsidIdentityAnnotator(tmp_path) as identity, _gwas(tmp_path, identity) as gwas:
        traits = {a.trait for a in gwas.annotate(_GWAS_CALLS[0])}
        assert traits == {"Trait without a risk allele"}
        assert gwas.run_counts()["not_in_identity_table"] == 3


def test_gwas_without_a_provider_compares_as_written(data_dir: Path) -> None:
    with _gwas(data_dir, None) as gwas:
        traits = {a.trait for a in gwas.annotate(_GWAS_CALLS[1])}
        assert traits == {"Trait at an A/T site"}
        assert sum(gwas.run_counts().values()) == 0


# ── Effect direction ────────────────────────────────────────────


@pytest.mark.parametrize(
    ("risk", "value", "ci", "expected"),
    [
        ("A", 1.5, "[1.2-1.8]", ("increases", "odds_ratio")),
        ("A", 0.3, "[0.2-0.4]", ("decreases", "odds_ratio")),
        ("A", 1.0, "[0.9-1.1]", ("unknown", "")),
        ("A", 0.05, "[0.03-0.07] unit increase", ("increases", "beta")),
        ("A", 0.05, "[0.03-0.07] unit decrease", ("decreases", "beta")),
        ("A", 0.05, "unit decrease", ("decreases", "beta")),
        ("A", -0.05, None, ("decreases", "beta")),
        ("A", 1.31, "NR", ("unknown", "")),
        ("A", None, "[1.2-1.8]", ("unknown", "")),
        (None, 1.5, "[1.2-1.8]", ("unknown", "")),
    ],
)
def test_effect_direction(
    risk: str | None, value: float | None, ci: str | None, expected: tuple[str, str]
) -> None:
    assert _effect_direction(risk, value, ci) == expected


def test_protective_reference_allele_is_shown_as_protective(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    """The outside review's case: a hom-ref user, the effect allele is REF, OR 0.30."""
    hom_ref = _call("rs900159004", 4000, "C", "C", "C")
    with _gwas(data_dir, identity) as gwas:
        (row,) = gwas.annotate(hom_ref)
    assert row.alt == "C"
    assert row.effect_direction == "decreases"
    assert effect_words(row) == (
        "C lowers the odds (a protective association). C is the reference allele."
    )
    assert effect_note(row) == " (C lowers odds)"


def test_effect_direction_in_json_and_empty_elsewhere(tmp_path: Path) -> None:
    from allelix.reports.json_report import render_json

    gwas_row = Annotation(
        source="gwas",
        rsid="rs900159004",
        significance="gwas_association",
        category="trait",
        magnitude=8.0,
        description="GWAS Catalog: x",
        attribution="GWAS Catalog",
        genotype_match="CC",
        effect_direction="decreases",
        effect_measure="odds_ratio",
    )
    other = Annotation(
        source="snpedia",
        rsid="rs900159004",
        significance="snpedia_bad",
        category="clinical",
        magnitude=3.0,
        description="SNPedia: x",
        attribution="SNPedia",
        genotype_match="CC",
    )
    result = _result([gwas_row, other])
    result.strand_abstentions = {"gwas": {"palindromic": 2}, "snpedia": {}}
    out = tmp_path / "report.json"
    render_json(result, output_path=out)
    payload = json.loads(out.read_text())
    directions = {a["source"]: a["effect_direction"] for a in payload["annotations"]}
    assert directions == {"gwas": "decreases", "snpedia": ""}
    assert "effect_measure" not in payload["annotations"][0]
    assert payload["strand_abstentions"] == {
        "gwas": {"not_in_identity_table": 0, "palindromic": 2, "no_fit": 0},
        "snpedia": {"not_in_identity_table": 0, "palindromic": 0, "no_fit": 0},
    }


def test_strand_abstention_lines() -> None:
    lines = strand_abstention_lines(
        {"gwas": {"palindromic": 2, "no_fit": 1}, "snpedia": {"not_in_identity_table": 1}}
    )
    assert lines == [
        "GWAS Catalog: 3 association rows not assessed because the strand of the source's "
        "alleles could not be settled (A/T or C/G site: 2, alleles fit neither strand: 1).",
        "SNPedia: 1 site not assessed because the strand of the source's alleles could not "
        "be settled (rsID not placed at the call by the rsID identity table: 1).",
    ]
    assert strand_abstention_lines({"gwas": {}}) == []


# ── SNPedia ─────────────────────────────────────────────────────

_SNPEDIA_CALLS = (
    _call("rs900159004", 4000, "C", "C", "C"),
    _call("rs900159002", 2000, "A", "A", "A"),
    _call("rs900159005", 5000, "C", "C", "C"),
    _call("rs900159003", 3000, "T", "T", "T"),
    _call("rs900159001", 1000, "C", "C", "G"),
    _call("i900159006", 6000, "A", "A", "A"),
)


def _snpedia(data_dir: Path, identity: RsidIdentityAnnotator | None) -> SNPediaAnnotator:
    return SNPediaAnnotator(
        data_dir,
        forward_site_provider=identity.forward_sites if identity is not None else None,
    )


def test_snpedia_reads_genotypes_on_the_forward_strand(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    with _snpedia(data_dir, identity) as snpedia:
        rows = list(snpedia.batch_annotate(_SNPEDIA_CALLS))
        counts = snpedia.run_counts()
    assert {(a.rsid, a.genotype_match) for a in rows} == {
        ("rs900159004", "CC"),
        ("rs900159003", "TT"),
        ("i900159006", "AA"),
    }
    assert counts == {"not_in_identity_table": 0, "palindromic": 2, "no_fit": 1}


def test_snpedia_single_and_batch_paths_agree(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    with _snpedia(data_dir, identity) as single, _snpedia(data_dir, identity) as batch:
        one_by_one = [a for v in _SNPEDIA_CALLS for a in single.annotate(v)]
        batched = list(batch.batch_annotate(_SNPEDIA_CALLS))
        assert one_by_one == batched
        assert single.run_counts() == batch.run_counts()


def test_snpedia_single_row_site_without_orientation_abstains_at_an_at_site(
    data_dir: Path, identity: RsidIdentityAnnotator
) -> None:
    """A one-row homozygous site cannot show it is palindromic; the identity table can."""
    with _snpedia(data_dir, identity) as snpedia:
        assert snpedia.annotate(_SNPEDIA_CALLS[1]) == []
        assert snpedia.run_counts()["palindromic"] == 1


# ── The GWAS floor ──────────────────────────────────────────────


def _result(annotations: list[Annotation]) -> AnalysisResult:
    return AnalysisResult(
        file_path=Path("/dev/null"),
        parser_name="test",
        parser_display_name="Test",
        sample_id="TEST",
        build="GRCh38",
        total_variants=1,
        skipped_count=0,
        annotators_used=[],
        annotations=annotations,
    )


def _gwas_row(rsid: str, magnitude: float, *, must_include: bool = False) -> Annotation:
    return Annotation(
        source="gwas",
        rsid=rsid,
        significance="gwas_association",
        category="trait",
        magnitude=magnitude,
        description="GWAS Catalog: x",
        attribution="GWAS Catalog",
        genotype_match="CC",
        is_must_include=must_include,
    )


def test_gwas_floor_is_the_larger_of_the_two() -> None:
    result = _result([_gwas_row("rs900159001", 6.0)])
    assert result.filter(min_magnitude=7.0, source_min_magnitudes={"gwas": 5.0}) == []
    assert len(result.filter(min_magnitude=5.0, source_min_magnitudes={"gwas": 5.0})) == 1
    assert result.filter(min_magnitude=5.0, source_min_magnitudes={"gwas": 9.0}) == []


def test_must_include_skips_the_gwas_floor_but_not_min_magnitude() -> None:
    result = _result([_gwas_row("rs10737680", 6.0, must_include=True)])
    assert len(result.filter(min_magnitude=5.0, source_min_magnitudes={"gwas": 9.0})) == 1
    assert result.filter(min_magnitude=7.0, source_min_magnitudes={"gwas": 5.0}) == []
    assert result.filter(min_magnitude=7.0, source_min_magnitudes={"gwas": 9.0}) == []


def test_snpedia_floor_still_replaces_min_magnitude() -> None:
    row = Annotation(
        source="snpedia",
        rsid="rs900159001",
        significance="snpedia_bad",
        category="clinical",
        magnitude=3.0,
        description="SNPedia: x",
        attribution="SNPedia",
        genotype_match="CC",
    )
    result = _result([row])
    assert len(result.filter(min_magnitude=5.0, source_min_magnitudes={"snpedia": 2.0})) == 1


# ── Through the CLI ─────────────────────────────────────────────


def test_analyze_wires_the_identity_table(
    tmp_path: Path, mock_mhg_path: Path, all_annotators_data_dir: Path
) -> None:
    """``get_annotators`` gives GWAS the identity table; the JSON carries the counts."""
    from click.testing import CliRunner

    from allelix.cli import main

    out = tmp_path / "report.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(mock_mhg_path),
            "--data-dir",
            str(all_annotators_data_dir),
            "--min-magnitude",
            "0",
            "--gwas-min-magnitude",
            "0",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out.read_text())
    assert set(payload["strand_abstentions"]) == {"gwas"}
    gwas_rows = [a for a in payload["annotations"] if a["source"] == "gwas"]
    assert gwas_rows
    assert {a["effect_direction"] for a in gwas_rows} <= {"increases", "decreases", "unknown"}
    # The mock chip's GWAS rsIDs are all in the identity rows: nothing abstains.
    assert sum(payload["strand_abstentions"]["gwas"].values()) == 0


def _chip(path: Path) -> Path:
    """A 23andMe-style chip of the synthetic calls (GRCh37 loci)."""
    rows = [
        ("rs900159001", 1000, "AG"),
        ("rs900159002", 2000, "AA"),
        ("rs900159004", 4000, "CC"),
    ]
    # A 23andMe banner, so --format 23andme is a choice, not a forced mismatch.
    lines = [
        "# This data file generated by 23andMe at: Mon Jan 01 00:00:00 2001",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    lines += [f"{rsid}\t1\t{pos}\t{gt}" for rsid, pos, gt in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_analyze_reads_gwas_through_the_identity_table(tmp_path: Path, data_dir: Path) -> None:
    """The A/T-site row abstains through ``analyze``; the console says so."""
    from click.testing import CliRunner

    from allelix.cli import main

    out = tmp_path / "report.json"
    result = CliRunner(env={"COLUMNS": "300"}).invoke(
        main,
        [
            "analyze",
            str(_chip(tmp_path / "chip.txt")),
            "--format",
            "23andme",
            "--build",
            "grch37",
            "--no-update",
            "--data-dir",
            str(data_dir),
            "--min-magnitude",
            "0",
            "--gwas-min-magnitude",
            "0",
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(out.read_text())
    traits = {a["condition"] for a in payload["annotations"] if a["source"] == "gwas"}
    assert "Trait at an A/T site" not in traits
    assert "Trait read complemented" in traits
    assert payload["strand_abstentions"]["gwas"]["palindromic"] == 1
    flat = " ".join(result.output.split())
    assert "GWAS Catalog: 1 association row not assessed" in flat


def test_gwas_floor_help_says_larger() -> None:
    from click.testing import CliRunner

    from allelix.cli import main

    result = CliRunner(env={"COLUMNS": "300"}).invoke(main, ["analyze", "--help"])
    flat = " ".join(result.output.split())
    assert "the larger of this and --min-magnitude" in flat


def _protective_row() -> Annotation:
    return Annotation(
        source="gwas",
        rsid="rs900159004",
        significance="gwas_association",
        category="trait",
        magnitude=9.0,
        description="GWAS Catalog: protective trait",
        attribution="GWAS Catalog",
        genotype_match="CC",
        condition="Protective trait",
        alt="C",
        locus_ref="C",
        carried_alleles=("C", "C"),
        effect_direction="decreases",
        effect_measure="odds_ratio",
    )


def test_html_shows_effect_words_and_strand_counts(tmp_path: Path) -> None:
    from allelix.reports.html import render_html

    result = _result([_protective_row()])
    result.strand_abstentions = {"gwas": {"palindromic": 2}}
    out = tmp_path / "report.html"
    render_html(result, output_path=out)
    html = out.read_text()
    assert "C lowers the odds (a protective association). C is the reference allele." in html
    assert "<strong>Strand.</strong>" in html
    assert "GWAS Catalog: 2 association rows not assessed" in html


def test_terminal_notes_the_effect_direction() -> None:
    from rich.console import Console

    from allelix.reports.terminal import render_terminal

    console = Console(record=True, width=300)
    render_terminal(_result([_protective_row()]), console)
    assert "(C lowers odds)" in console.export_text()

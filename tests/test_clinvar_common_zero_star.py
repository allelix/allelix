# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A zero-star ClinVar P/LP row on a common allele leaves the default report (ADR-0052).

The frequency is gnomAD's for ClinVar's own allele, looked up at the
rsID's GRCh38 locus from the identity table, so it works on GRCh37 input
at palindromic sites where the call's own frequency is missing. A row is
scored 4.0, marked and counted when its submission has zero stars and that
frequency is at least 0.05. Every doubt about the allele leaves the row
as it was. All rows, calls and frequencies are synthetic; the HFE and HBB
rsIDs and positions are the public ones, the rows are not ClinVar's.
"""

from __future__ import annotations

import contextlib
import json
import sqlite3
from typing import TYPE_CHECKING, NamedTuple

import pytest
from click.testing import CliRunner
from rich.console import Console

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION, GNOMAD_SCHEMA_VERSION
from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import (
    COMMON_ALLELE_AF,
    COMMON_ALLELE_MAGNITUDE,
    run_analysis,
)
from allelix.reports.clinvar_text import (
    common_allele_count,
    common_allele_line,
    detail_fields,
    terminal_note,
)
from allelix.reports.html import _get_repute, render_html
from allelix.reports.json_report import render_json
from allelix.reports.terminal import render_terminal

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_NONE = "no assertion criteria provided"
_NO_CLASS = "no classification provided"
_ONE = "criteria provided, single submitter"
_PANEL = "reviewed by expert panel"
_CONFLICTING = "Conflicting classifications of pathogenicity"


class Site(NamedTuple):
    rsid: str  # "" for a ClinVar row without an rsID
    chrom: str
    pos37: int
    pos38: int
    sig: str
    review: str
    aggregate: str
    af: float | None  # gnomAD AF of the allele at the GRCh38 key; None: no row
    ref: str = "G"
    alt: str = "A"
    # The identity table's rows as (pos37, pos38, ref, alt), the same REF and ALT on
    # both builds; None: one row at the site's own positions.
    table: tuple[tuple[int, int, str, str], ...] | None = None
    # ClinVar GRCh38's ALT when it differs from GRCh37's.
    alt38: str | None = None
    # More ClinVar records of the rsID, same REF and ALT, at these positions.
    more37: tuple[int, ...] = ()
    more38: tuple[int, ...] = ()

    def table_rows(self) -> tuple[tuple[int, int, str, str], ...]:
        return self.table or ((self.pos37, self.pos38, self.ref, self.alt),)


_SITES: tuple[Site, ...] = (
    # rs1063192-shaped: dbSNP lists G>A, G>C and G>T, so the site is palindromic.
    Site(
        "rs900970001", "2", 1_000, 11_000, "Likely pathogenic", _NONE,
        "Likely pathogenic; protective", 0.05,
        table=((1_000, 11_000, "G", "A"), (1_000, 11_000, "G", "C"), (1_000, 11_000, "G", "T")),
    ),
    Site(
        "rs900970002", "2", 2_000, 12_000, "Likely pathogenic", _NONE, "Likely pathogenic",
        0.0499,
    ),
    Site("rs900970003", "2", 3_000, 13_000, "Pathogenic", _NO_CLASS, "Pathogenic", 0.65),
    Site("rs900970004", "2", 4_000, 14_000, "Pathogenic", _ONE, "Pathogenic", 0.65),
    Site("rs900970005", "2", 5_000, 15_000, "Likely pathogenic", _NONE, "Likely pathogenic", None),
    Site("rs900970006", "2", 6_000, 16_000, "protective", _NONE, "protective", 0.65),
    Site("rs900970007", "2", 7_000, 17_000, "Likely pathogenic", _PANEL, "Likely pathogenic", 0.9),
    Site(
        "rs900970008", "2", 8_000, 18_000, "Pathogenic, low penetrance", _NONE, "Pathogenic",
        0.2,
    ),
    Site("rs900970009", "2", 9_000, 19_000, "Pathogenic", _NONE, _CONFLICTING, 0.3),
    # The table places the rsID at two GRCh38 loci; ClinVar and gnomAD hold both.
    Site(
        "rs900970010", "2", 10_000, 20_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5,
        table=((10_000, 20_000, "G", "A"), (10_000, 20_500, "G", "A")), more38=(20_500,),
    ),
    # The same at two GRCh37 loci.
    Site(
        "rs900970013", "2", 30_000, 40_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5,
        table=((30_000, 40_000, "G", "A"), (30_500, 40_000, "G", "A")), more37=(30_500,),
    ),
    # The call is not at the table's locus; ClinVar also has the rsID there.
    Site(
        "rs900970014", "2", 31_000, 41_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5,
        table=((31_500, 41_500, "G", "A"),), more37=(31_500,), more38=(41_500,),
    ),
    # ClinVar's builds disagree on the ALT; gnomAD has GRCh37's G>A.
    Site(
        "rs900970011", "2", 11_000, 21_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5,
        alt38="T",
    ),
    # The table holds another allele only.
    Site(
        "rs900970012", "2", 12_000, 22_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5,
        table=((12_000, 22_000, "G", "C"),),
    ),
    # No rsID (ADR-0046): no table lookup.
    Site("", "2", 13_000, 23_000, "Likely pathogenic", _NONE, "Likely pathogenic", 0.5),
    # HFE C282Y and HbS: one star or more, untouched whatever the frequency.
    Site(
        "rs1800562", "6", 26_093_141, 26_092_913, "Pathogenic", _ONE,
        "Pathogenic/Pathogenic, low penetrance; risk factor", 0.059,
    ),
    Site(
        "rs334", "11", 5_248_232, 5_227_002, "Pathogenic", _ONE, "Pathogenic", 0.0016,
        ref="T", alt="A",
        table=tuple((5_248_232, 5_227_002, "T", alt) for alt in "ACG"),
    ),
)  # fmt: skip
_HIDDEN = {"rs900970001", "rs900970003", "rs900970008"}
_EXPECTED = {
    # AF exactly 0.05 at a palindromic site: hidden.
    "rs900970001": (4.0, True),
    # AF 0.0499: below BA1, unchanged.
    "rs900970002": (7.0, False),
    # "no classification provided" is zero stars too.
    "rs900970003": (4.0, True),
    # One star: untouched at any frequency.
    "rs900970004": (9.0, False),
    # No gnomAD row: as before.
    "rs900970005": (7.0, False),
    # Not a P/LP classification.
    "rs900970006": (5.0, False),
    # Expert panel: untouched.
    "rs900970007": (7.0, False),
    # Low penetrance is P/LP.
    "rs900970008": (4.0, True),
    # Already capped by a Conflicting aggregate: not marked again.
    "rs900970009": (4.0, False),
    # Doubts: unchanged.
    "rs900970010": (7.0, False),
    "rs900970011": (7.0, False),
    "rs900970012": (7.0, False),
    "rs900970013": (7.0, False),
    "rs900970014": (7.0, False),
    "": (7.0, False),
    "rs1800562": (9.0, False),
    "rs334": (9.0, False),
}


def _write_clinvar(data: Path, build: str) -> None:
    rows = []
    for i, s in enumerate(_SITES, 1):
        pos = s.pos37 if build == "GRCh37" else s.pos38
        alt = s.alt38 if build == "GRCh38" and s.alt38 else s.alt
        more = s.more37 if build == "GRCh37" else s.more38
        rows += [
            (s.rsid, s.chrom, p, s.ref, alt, s.sig, s.review, i, s.aggregate) for p in (pos, *more)
        ]
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total) "
            "VALUES (?, ?, ?, ?, ?, ?, 'Test condition', 'GENE', ?, ?, ?, 1, 1)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://common', '1', '2026-10-10', ?, NULL, ?)",
            (clinvar_record_name(build), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _write_identity(data: Path) -> None:
    rows = []
    for s in _SITES:
        if not s.rsid:
            continue
        table = s.table_rows()
        palindromic = int(len({alt for *_p, alt in table}) > 1)
        rows += [
            (s.rsid, s.chrom, pos37, ref, alt, pos38, ref, alt, palindromic)
            for pos37, pos38, ref, alt in table
        ]
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://common",
        downloaded_at="2026-10-10T00:00:00Z",
    )


def _write_gnomad(data: Path) -> None:
    with contextlib.closing(sqlite3.connect(data / GNOMAD_DB_FILENAME)) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            [
                (s.chrom, pos38, s.ref, s.alt, s.rsid or None, s.af)
                for s in _SITES
                if s.af is not None
                for pos38 in sorted({s.pos38} | {row[1] for row in s.table_rows()})
            ],
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) "
            "VALUES ('gnomad', 'test://common', '4.1', '2026-10-10', ?, ?)",
            (len(_SITES), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()


def _write_vcf(path: Path, build: str) -> Path:
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    for s in _SITES:
        pos = s.pos37 if build == "GRCh37" else s.pos38
        alt = s.alt38 if build == "GRCh38" and s.alt38 else s.alt
        lines.append(f"{s.chrom}\t{pos}\t{s.rsid or '.'}\t{s.ref}\t{alt}\t50\tPASS\t.\tGT\t1/1")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def data(tmp_path: Path) -> Path:
    d = tmp_path / "data"
    d.mkdir()
    for build in ("GRCh37", "GRCh38"):
        _write_clinvar(d, build)
    _write_identity(d)
    _write_gnomad(d)
    return d


def _run(
    tmp_path: Path,
    data: Path,
    stack: contextlib.ExitStack,
    build: str = "GRCh37",
    *,
    gnomad: bool = True,
    identity: bool = True,
):
    return run_analysis(
        _write_vcf(tmp_path / f"s{build}.vcf", build),
        VcfParser(),
        [stack.enter_context(ClinVarAnnotator(data))],
        build_override=build,
        gnomad=stack.enter_context(GnomadAnnotator(data)) if gnomad else None,
        identity=stack.enter_context(RsidIdentityAnnotator(data)) if identity else None,
    )


def _clinvar_rows(result) -> dict[str, Annotation]:
    """ClinVar rows by rsID; the rsID-less row is keyed ""."""
    rows = {}
    for a in result.annotations:
        if a.source == "clinvar":
            rows[a.rsid if a.rsid.startswith("rs") else ""] = a
    assert set(rows) == {s.rsid for s in _SITES}
    return rows


def _scored(result) -> dict[str, tuple[float, bool]]:
    return {k: (a.magnitude, a.common_allele_capped) for k, a in _clinvar_rows(result).items()}


class TestRule:
    def test_the_threshold_is_ba1(self) -> None:
        assert COMMON_ALLELE_AF == 0.05
        assert COMMON_ALLELE_MAGNITUDE < 5.0

    @pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
    def test_rows_scored(self, tmp_path: Path, data: Path, cm_stack, build: str) -> None:
        assert _scored(_run(tmp_path, data, cm_stack, build)) == _EXPECTED

    def test_grch37_palindromic_site_has_no_call_frequency(
        self, tmp_path: Path, data: Path, cm_stack
    ) -> None:
        """The case the rule is for: enrichment abstains, the ClinVar allele lookup does not."""
        row = _clinvar_rows(_run(tmp_path, data, cm_stack))["rs900970001"]
        assert row.allele_frequency is None
        assert row.clinvar_allele_af == 0.05
        assert row.common_allele_capped

    def test_hfe_c282y_and_hbs_are_untouched(self, tmp_path: Path, data: Path, cm_stack) -> None:
        for build in ("GRCh37", "GRCh38"):
            rows = _clinvar_rows(_run(tmp_path, data, cm_stack, build))
            for rsid in ("rs1800562", "rs334"):
                assert rows[rsid].magnitude == 9.0
                assert not rows[rsid].common_allele_capped
                assert rows[rsid].clinvar_allele_af is None

    def test_without_gnomad_nothing_changes(self, tmp_path: Path, data: Path, cm_stack) -> None:
        result = _run(tmp_path, data, cm_stack, gnomad=False)
        assert not any(m for _mag, m in _scored(result).values())
        assert result.clinvar_common_allele == 0

    def test_without_the_identity_table_nothing_changes(
        self, tmp_path: Path, data: Path, cm_stack
    ) -> None:
        result = _run(tmp_path, data, cm_stack, identity=False)
        assert not any(m for _mag, m in _scored(result).values())

    def test_with_one_clinvar_build_nothing_changes(
        self, tmp_path: Path, data: Path, cm_stack
    ) -> None:
        (data / clinvar_db_filename("GRCh38")).unlink()
        result = _run(tmp_path, data, cm_stack)
        assert not any(m for _mag, m in _scored(result).values())
        assert _clinvar_rows(result)["rs900970003"].magnitude == 9.0

    def test_counted_and_out_of_the_default_report(
        self, tmp_path: Path, data: Path, cm_stack
    ) -> None:
        result = _run(tmp_path, data, cm_stack)
        assert result.clinvar_common_allele == 3
        shown = {a.rsid for a in result.filter(min_magnitude=5.0)}
        assert not shown & _HIDDEN
        assert {a.rsid for a in result.filter(min_magnitude=0.0)} >= _HIDDEN


class TestDisplay:
    def test_text(self, tmp_path: Path, data: Path, cm_stack) -> None:
        rows = _clinvar_rows(_run(tmp_path, data, cm_stack))
        hidden = rows["rs900970003"]
        assert terminal_note(hidden) == " (1 of 1; zero-star; common in gnomAD)"
        label, value = detail_fields(hidden)[-1]
        assert label == "Scored"
        assert "gnomAD frequency of 0.650, 5% or more" in value
        assert terminal_note(rows["rs900970002"]) == ""
        assert common_allele_count(list(rows.values())) == 3
        assert common_allele_line(0) == ""
        assert common_allele_line(1).startswith(
            "1 ClinVar Pathogenic or Likely pathogenic row not shown: zero-star "
            "classification on a common allele"
        )
        assert common_allele_line(3).startswith("3 ClinVar Pathogenic or Likely pathogenic rows")

    def test_html_repute_is_neutral(self, tmp_path: Path, data: Path, cm_stack) -> None:
        rows = _clinvar_rows(_run(tmp_path, data, cm_stack))
        assert _get_repute(rows["rs900970001"]) == "neutral"
        assert _get_repute(rows["rs900970002"]) == "bad"

    def test_json(self, tmp_path: Path, data: Path, cm_stack) -> None:
        out = tmp_path / "r.json"
        render_json(_run(tmp_path, data, cm_stack), output_path=out, min_magnitude=0.0)
        payload = json.loads(out.read_text())
        assert payload["clinvar"]["common_allele_rows_scored_down"] == 3
        rows = {a["rsid"]: a for a in payload["annotations"] if a["source"] == "clinvar"}
        assert rows["rs900970001"]["common_allele_capped"] is True
        assert rows["rs900970001"]["clinvar_allele_af"] == 0.05
        assert rows["rs900970001"]["magnitude"] == 4.0
        assert rows["rs900970002"]["clinvar_allele_af"] == 0.0499
        assert rows["rs900970004"]["common_allele_capped"] is False

    def test_html_banner(self, tmp_path: Path, data: Path, cm_stack) -> None:
        out = tmp_path / "r.html"
        render_html(_run(tmp_path, data, cm_stack), output_path=out)
        assert (
            "3 ClinVar Pathogenic or Likely pathogenic rows not shown: zero-star "
            "classification on a common allele"
        ) in out.read_text()

    def test_terminal_note(self, tmp_path: Path, data: Path, cm_stack) -> None:
        console = Console(record=True, width=300)
        render_terminal(_run(tmp_path, data, cm_stack), console, min_magnitude=0.0)
        assert "(1 of 1; zero-star; common in gnomAD)" in console.export_text()

    def test_console_line(self, tmp_path: Path, data: Path) -> None:
        out = tmp_path / "r.json"
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main,
            [
                "analyze",
                str(_write_vcf(tmp_path / "s.vcf", "GRCh37")),
                "--build",
                "grch37",
                "--no-update",
                "--data-dir",
                str(data),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            "3 ClinVar Pathogenic or Likely pathogenic rows not shown: zero-star "
            "classification on a common allele"
        ) in " ".join(result.output.split())
        payload = json.loads(out.read_text())
        assert payload["clinvar"]["common_allele_rows_scored_down"] == 3
        assert not {a["rsid"] for a in payload["annotations"]} & _HIDDEN

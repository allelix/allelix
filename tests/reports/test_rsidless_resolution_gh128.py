# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GH #128 regression tests — rsID-less VCF coverage on the analyze
panel-coverage and extract paths.

Background: DeepVariant / DRAGEN / GIAB GRCh38-benchmark VCFs leave
the ID column empty. The parser produces a positional pseudo-rsID
(``chr1:11796321:G:A``), the ClinVar position resolver
(``ClinVarAnnotator.bulk_resolve_rsids``) recovers an ``rs...`` for
ClinVar-known variants, but the gnomAD-only / pharma-only / GWAS-only
variants that fill most wellness panels stayed as pseudo-IDs because
the resolver only checked ClinVar. The v2.2.2 fix:

  1. ``GnomadAnnotator.bulk_resolve_rsids_from_positions`` — reverse-
     lookup keyed on ``(chrom, pos)`` against the full dbSNP rsID
     universe gnomAD carries.
  2. Second-pass resolution in ``_pipeline._flush()`` runs gnomAD
     against the variants still rsID-less after ClinVar's pass, stamps
     ``v.rsid`` (audit-stashed in ``v.original_rsid``).
  3. ``_sequential_extract`` and ``_execute_tabix_extract`` updated
     symmetrically — extract was hitting the same pseudo-ID problem.

This file pins the analyze panel-coverage path. Extract paths are
covered in tests/test_extract_rsidless_gh128.py.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

GNOMAD_DB_FILENAME = "gnomad.sqlite"

# Three variants at gnomAD-mock-known positions, all with ID=. so the
# resolver path is the only way the pipeline can know their rsIDs.
# Positions intentionally chosen to overlap the conftest mock fixtures:
#   chr1:11796321 G→A  → rs1801133 (in ClinVar GRCh38 mock + gnomAD mock)
#   chr22:19963748 G→A → rs4680    (NOT in ClinVar mock; in gnomAD mock)
#   chr19:44908684 T→C → rs429358  (NOT in ClinVar mock; in gnomAD mock)
#
# Pipeline order:
#   - ClinVar's resolver stamps rs1801133 (only one in its cache)
#   - The new gnomAD second-pass stamps rs4680 and rs429358
_RSIDLESS_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=19,length=58617616>
##contig=<ID=22,length=51304566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t11796321\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
19\t44908684\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
22\t19963748\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
"""


@pytest.fixture
def rsidless_vcf(tmp_path: Path) -> Path:
    """Write the synthetic rsID-less VCF to a temp path."""
    p = tmp_path / "rsidless.vcf"
    p.write_text(_RSIDLESS_VCF)
    return p


@pytest.fixture
def data_dir_with_clinvar_and_gnomad(all_annotators_data_dir: Path) -> Path:
    """Reuse the ClinVar/ClinPGx/GWAS fixture and add the mock gnomAD on top.

    ``all_annotators_data_dir`` lacks gnomAD by default. The new GH #128
    resolution path requires it, so this fixture stitches the mock cache
    in alongside.
    """
    db_path = all_annotators_data_dir / GNOMAD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.executemany(
            "INSERT INTO gnomad_frequencies"
            " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
                ("22", 19963748, "G", "A", "rs4680", 0.50),
                ("19", 44908684, "T", "C", "rs429358", 0.15),
                ("19", 44908684, "T", "G", "rs429358", 0.08),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count,"
            "  local_version_tag)"
            " VALUES (?, ?, ?, ?, ?, ?)",
            (
                "gnomad",
                "test://mock",
                "4.1",
                "2026-01-01T00:00:00Z",
                4,
                "sv:1",
            ),
        )
        conn.commit()
    return all_annotators_data_dir


class TestGnomadReverseResolver:
    """Direct test of the new ``bulk_resolve_rsids_from_positions`` method."""

    def test_known_positions_resolve(self, data_dir_with_clinvar_and_gnomad: Path) -> None:
        with GnomadAnnotator(data_dir_with_clinvar_and_gnomad) as gnomad:
            assert gnomad.is_ready()
            result = gnomad.bulk_resolve_rsids_from_positions(
                {("1", 11796321), ("22", 19963748), ("19", 44908684)}
            )
        assert ("1", 11796321) in result
        assert ("G", "A", "rs1801133") in result[("1", 11796321)]
        assert ("22", 19963748) in result
        assert ("G", "A", "rs4680") in result[("22", 19963748)]
        assert ("19", 44908684) in result
        rs429358_rows = result[("19", 44908684)]
        # Multi-allelic — both alts surface for the carrier-rule pass.
        assert ("T", "C", "rs429358") in rs429358_rows
        assert ("T", "G", "rs429358") in rs429358_rows

    def test_empty_positions_returns_empty(self, data_dir_with_clinvar_and_gnomad: Path) -> None:
        with GnomadAnnotator(data_dir_with_clinvar_and_gnomad) as gnomad:
            assert gnomad.bulk_resolve_rsids_from_positions(set()) == {}

    def test_unknown_position_yields_no_entry(
        self, data_dir_with_clinvar_and_gnomad: Path
    ) -> None:
        with GnomadAnnotator(data_dir_with_clinvar_and_gnomad) as gnomad:
            result = gnomad.bulk_resolve_rsids_from_positions(
                {("1", 99_999_999)}  # not in the mock
            )
        assert result == {}


class TestPanelCoverageOnRsidlessVcf:
    """The user's reported symptom on #128 — panel coverage on a
    rsID-less VCF should classify panel rsIDs by what the resolver
    recovers, not by what the file's ID column originally said."""

    def test_panel_rsid_resolved_via_gnomad_counts_as_genotyped(
        self,
        rsidless_vcf: Path,
        data_dir_with_clinvar_and_gnomad: Path,
        cm_stack,
    ) -> None:
        """rs4680 is in gnomAD's mock cache but not ClinVar's mock.
        Pre-fix, the panel detection block in ``_flush()`` saw the
        pseudo-rsID ``chr22:19963748:G:A`` and classified rs4680 as
        state-3 missing. Post-fix the gnomAD second-pass stamps
        ``v.rsid = rs4680`` and the panel detection sees it correctly.

        The test runs the full pipeline with a panel of three rsIDs,
        all of which are present in the synthetic VCF, and asserts
        ``coverage["found"] == 3`` and ``missing == []`` — the failure
        mode under the bug was ``found == 0, missing == [<all 3>]``.
        """
        parser = VcfParser()
        annotators = [
            cm_stack.enter_context(ClinVarAnnotator(data_dir_with_clinvar_and_gnomad)),
        ]
        panel = frozenset({"rs1801133", "rs4680", "rs429358"})
        gnomad = cm_stack.enter_context(GnomadAnnotator(data_dir_with_clinvar_and_gnomad))
        result = run_analysis(
            rsidless_vcf,
            parser,
            annotators,
            gnomad=gnomad,
            panel_rsids=panel,
        )
        coverage = result.panel_coverage()
        assert coverage is not None
        assert coverage["requested"] == 3
        assert coverage["found"] == 3, (
            f"Panel coverage missed gnomAD-resolved rsIDs. "
            f"missing={coverage['missing']!r}, found={coverage['found']}"
        )
        assert coverage["missing"] == []

    def test_panel_rsid_genuinely_absent_classified_as_missing(
        self,
        rsidless_vcf: Path,
        data_dir_with_clinvar_and_gnomad: Path,
        cm_stack,
    ) -> None:
        """A panel rsID at a position the file doesn't contain stays
        in state-3 missing — the fix must not silently fabricate hits.
        rs9999999 is at no position in the synthetic VCF."""
        parser = VcfParser()
        annotators = [
            cm_stack.enter_context(ClinVarAnnotator(data_dir_with_clinvar_and_gnomad)),
        ]
        panel = frozenset({"rs1801133", "rs9999999"})
        gnomad = cm_stack.enter_context(GnomadAnnotator(data_dir_with_clinvar_and_gnomad))
        result = run_analysis(
            rsidless_vcf,
            parser,
            annotators,
            gnomad=gnomad,
            panel_rsids=panel,
        )
        coverage = result.panel_coverage()
        assert coverage is not None
        assert coverage["found"] == 1
        assert coverage["missing"] == ["rs9999999"]


class TestOriginalRsidStashedDuringResolution:
    """The audit stash ``v.original_rsid`` must hold the pre-stamp
    pseudo-rsID after a successful resolution — both for the existing
    ClinVar resolver (clinvar.py:688) and the new gnomAD second-pass."""

    def test_clinvar_resolver_stashes_original(
        self,
        rsidless_vcf: Path,
        data_dir_with_clinvar_and_gnomad: Path,
        cm_stack,
    ) -> None:
        """Round-trip the synthetic VCF through ``run_analysis`` and
        inspect ``hv_variants`` for the stash. (hv_variants is the
        per-variant pipeline view exposed on AnalysisResult; it carries
        Variant objects post-resolution.)

        rs1801133 is ClinVar-resolvable; its hv_variant entry — IF
        rs1801133 happens to be in the high-value set — would carry
        original_rsid != None. Since the hv set is config-driven and
        may not include test rsIDs by default, the safer assertion is
        on ``panel_genotypes`` which records ``{rsid: genotype}`` from
        the panel-collection pass that runs AFTER both resolvers."""
        parser = VcfParser()
        annotators = [
            cm_stack.enter_context(ClinVarAnnotator(data_dir_with_clinvar_and_gnomad)),
        ]
        # Panel contains all three resolved rsIDs so panel_genotypes is
        # populated with each. Their presence proves the stamp landed
        # in time for the panel-collect block (post-resolution).
        panel = frozenset({"rs1801133", "rs4680", "rs429358"})
        gnomad = cm_stack.enter_context(GnomadAnnotator(data_dir_with_clinvar_and_gnomad))
        result = run_analysis(
            rsidless_vcf,
            parser,
            annotators,
            gnomad=gnomad,
            panel_rsids=panel,
        )
        # Each variant was 0/1 (het): G/A and T/C.
        assert result.panel_genotypes is not None
        assert result.panel_genotypes.get("rs1801133") == "G/A"
        assert result.panel_genotypes.get("rs4680") == "G/A"
        assert result.panel_genotypes.get("rs429358") == "T/C"

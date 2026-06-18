# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the unified analysis pipeline."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.models import Variant
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import (
    _DETECTION_BUFFER_LIMIT,
    AnalysisResult,
    _BuildDetectionState,
    run_analysis,
)
from allelix.utils.build_detect import BUILD_GRCH36, BUILD_GRCH37, KNOWN_SNP_POSITIONS

if TYPE_CHECKING:
    from pathlib import Path


def _ann(**overrides):
    from allelix.models import Annotation

    defaults = {
        "source": "clinvar",
        "rsid": "rs1",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 5.0,
        "description": "x",
        "attribution": "ClinVar",
        "genotype_match": "A",
        "gene": "GENE1",
    }
    defaults.update(overrides)
    return Annotation(**defaults)


class TestAnalysisResultFilter:
    def _result(self, annotations) -> AnalysisResult:
        from pathlib import Path

        return AnalysisResult(
            file_path=Path("dummy.txt"),
            parser_name="x",
            parser_display_name="X",
            sample_id="S",
            build="GRCh37",
            total_variants=0,
            skipped_count=0,
            annotators_used=[],
            annotations=annotations,
        )

    def test_min_magnitude_excludes_low(self):
        r = self._result([_ann(rsid="lo", magnitude=2), _ann(rsid="hi", magnitude=8)])
        kept = r.filter(min_magnitude=5)
        assert [a.rsid for a in kept] == ["hi"]

    def test_category_filter(self):
        r = self._result([_ann(rsid="c", category="clinical"), _ann(rsid="p", category="pharma")])
        assert [a.rsid for a in r.filter(category="pharma")] == ["p"]

    def test_genes_filter_case_insensitive(self):
        r = self._result([_ann(rsid="m", gene="MTHFR"), _ann(rsid="b", gene="BRCA1")])
        kept = r.filter(genes={"mthfr"})
        assert [a.rsid for a in kept] == ["m"]

    def test_sort_is_magnitude_then_rsid(self):
        r = self._result(
            [
                _ann(rsid="rs2", magnitude=5),
                _ann(rsid="rs1", magnitude=5),
                _ann(rsid="rs3", magnitude=8),
            ]
        )
        kept = r.filter()
        assert [a.rsid for a in kept] == ["rs3", "rs1", "rs2"]

    def test_rsids_filter_case_insensitive(self):
        r = self._result([_ann(rsid="rs1801133"), _ann(rsid="rs4680")])
        kept = r.filter(rsids={"RS1801133"})
        assert [a.rsid for a in kept] == ["rs1801133"]

    def test_genes_and_rsids_combine_with_or(self):
        r = self._result(
            [
                _ann(rsid="rs1", gene="MTHFR"),
                _ann(rsid="rs4680", gene="COMT"),
                _ann(rsid="rsXX", gene="OTHER"),
            ]
        )
        kept = r.filter(genes={"MTHFR"}, rsids={"rs4680"})
        assert sorted(a.rsid for a in kept) == ["rs1", "rs4680"]

    def test_empty_genes_set_matches_nothing(self):
        """An empty filter set != None: explicit "match nothing", empty report."""
        r = self._result([_ann(rsid="rs1", gene="MTHFR")])
        kept = r.filter(genes=frozenset())
        assert kept == []

    def test_empty_rsids_set_matches_nothing(self):
        r = self._result([_ann(rsid="rs1", gene="MTHFR")])
        kept = r.filter(rsids=frozenset())
        assert kept == []

    def test_both_empty_sets_match_nothing(self):
        r = self._result([_ann(rsid="rs1", gene="MTHFR")])
        kept = r.filter(genes=frozenset(), rsids=frozenset())
        assert kept == []

    def test_none_means_no_filter(self):
        """None on both means no gene/rsid filter — every annotation passes."""
        r = self._result([_ann(rsid="rs1", gene="MTHFR"), _ann(rsid="rs2", gene="COMT")])
        kept = r.filter(genes=None, rsids=None)
        assert len(kept) == 2

    def test_gene_none_with_rsid_filter_does_not_crash(self):
        """GWAS/intergenic annotations have gene=None; filter must not crash."""
        r = self._result([_ann(rsid="rs1", gene=None), _ann(rsid="rs2", gene="MTHFR")])
        kept = r.filter(rsids={"rs1"})
        assert [a.rsid for a in kept] == ["rs1"]


class TestPanelCoverage:
    """GH #75: panel coverage distinguishes "not on chip" (state 3) from
    "genotyped but no findings" (state 2) from "genotyped with findings"
    (state 1)."""

    def _result(self, **kwargs) -> AnalysisResult:
        from pathlib import Path

        defaults = dict(
            file_path=Path("dummy.txt"),
            parser_name="x",
            parser_display_name="X",
            sample_id="S",
            build="GRCh37",
            total_variants=0,
            skipped_count=0,
            annotators_used=[],
            annotations=[],
        )
        defaults.update(kwargs)
        return AnalysisResult(**defaults)

    def test_no_panel_returns_none(self):
        """No filter-file supplied → panel_coverage() is None and
        nothing should render."""
        r = self._result()
        assert r.panel_coverage() is None

    def test_classifies_three_states(self):
        """A panel of 4 rsIDs where one had findings, one was genotyped
        but produced no annotations, two weren't on the chip."""
        panel = frozenset({"rs1", "rs2", "rs3", "rs4"})
        genotyped = frozenset({"rs1", "rs2"})
        annotated = frozenset({"rs1"})
        r = self._result(
            panel_rsids=panel,
            genotyped_panel_rsids=genotyped,
            panel_annotated_rsids=annotated,
        )
        cov = r.panel_coverage()
        assert cov is not None
        assert cov["requested"] == 4
        assert cov["found"] == 2
        assert cov["missing"] == ["rs3", "rs4"]
        assert cov["no_findings"] == ["rs2"]

    def test_all_genotyped_with_findings_clean_coverage(self):
        """Happy path: every panel rsID was both genotyped and produced
        an annotation. missing + no_findings are both empty."""
        panel = frozenset({"rs1", "rs2"})
        r = self._result(
            panel_rsids=panel,
            genotyped_panel_rsids=frozenset({"rs1", "rs2"}),
            panel_annotated_rsids=frozenset({"rs1", "rs2"}),
        )
        cov = r.panel_coverage()
        assert cov == {"requested": 2, "found": 2, "missing": [], "no_findings": []}

    def test_nothing_genotyped_all_missing(self):
        """The chip didn't carry any of the requested rsIDs."""
        panel = frozenset({"rs1", "rs2"})
        r = self._result(panel_rsids=panel)
        cov = r.panel_coverage()
        assert cov is not None
        assert cov["found"] == 0
        assert cov["missing"] == ["rs1", "rs2"]
        assert cov["no_findings"] == []

    def test_run_analysis_populates_coverage(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """End-to-end: panel_rsids passed into run_analysis are
        recorded on AnalysisResult, and the genotyped/annotated sets
        come from the actual pipeline pass."""
        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        # rs1801133 is in the MHG fixture AND has ClinVar findings;
        # rs999999999 is a synthetic ClinVar rsID never on a real chip.
        panel = frozenset({"rs1801133", "rs999999999"})
        result = run_analysis(mock_mhg_path, parser, annotators, panel_rsids=panel)
        cov = result.panel_coverage()
        assert cov is not None
        assert cov["requested"] == 2
        assert "rs999999999" in cov["missing"]
        assert "rs999999999" not in cov["no_findings"]

    def test_sub_floor_panel_rsid_lands_in_no_findings_not_limbo(self):
        """GH #106 patch: when ``filtered_annotations`` is supplied
        the "annotated" set is derived from THAT list (post-magnitude-
        filter) instead of the raw pre-filter set.

        Repro of the v2.2.0 accounting lie: a panel rsid has an
        annotation in the cache (so it's in
        ``panel_annotated_rsids``), but the annotation's magnitude is
        below the analyze display floor, so the renderer's
        ``filter()`` drops it. Pre-#106, the rsid:
          - counted as "found" (because panel_annotated_rsids has it)
          - dropped from the annotation table (below floor)
          - NOT in no_findings (because the pipeline thought it had
            annotations)
        Pure limbo.

        Post-#106, threading the filtered list in makes "annotated"
        equal "actually rendered" — the sub-floor rsid lands in
        no_findings honestly.

        Pinned invariant: requested == len(missing) + len(no_findings)
        + count of distinct displayed-annotation rsids.
        """
        panel = frozenset({"rs_displayed", "rs_subfloor", "rs_notonchip"})
        # rs_displayed: genotyped, mag 9.0, will survive the filter.
        # rs_subfloor: genotyped, mag 2.0, won't survive the filter.
        # rs_notonchip: never in input.
        genotyped = frozenset({"rs_displayed", "rs_subfloor"})
        # Pre-filter "annotated" set: both genotyped rsids had cache
        # annotations. This is what run_analysis populates today.
        panel_annotated_raw = frozenset({"rs_displayed", "rs_subfloor"})

        r = self._result(
            panel_rsids=panel,
            genotyped_panel_rsids=genotyped,
            panel_annotated_rsids=panel_annotated_raw,
            annotations=[
                _ann(rsid="rs_displayed", magnitude=9.0),
                _ann(rsid="rs_subfloor", magnitude=2.0),
            ],
        )

        # Post-filter list: only rs_displayed survived the 5.0 floor.
        filtered = r.filter(min_magnitude=5.0)
        assert [a.rsid for a in filtered] == ["rs_displayed"]

        # Pre-#106 (back-compat path): the limbo bug — rs_subfloor
        # counted as "found" but rendered nowhere.
        legacy_cov = r.panel_coverage()
        assert legacy_cov is not None
        assert legacy_cov["found"] == 2
        assert "rs_subfloor" not in legacy_cov["no_findings"]

        # Post-#106: thread the filtered list. rs_subfloor falls into
        # no_findings; `found` still means "in input file" (genotyped
        # count) — same v6 contract.
        cov = r.panel_coverage(filtered)
        assert cov is not None
        assert cov["requested"] == 3
        assert cov["found"] == 2  # both genotyped, unchanged
        assert cov["missing"] == ["rs_notonchip"]
        assert cov["no_findings"] == ["rs_subfloor"]  # was [] in legacy

        # Pinned accounting invariant: the math closes against the
        # user-visible report.
        displayed_rsids = {a.rsid for a in filtered} & panel
        assert cov["requested"] == (
            len(cov["missing"]) + len(cov["no_findings"]) + len(displayed_rsids)
        )


class TestRunAnalysis:
    def test_streams_and_collects(self, mock_mhg_path: Path, all_annotators_data_dir: Path):
        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        result = run_analysis(mock_mhg_path, parser, annotators)
        assert result.parser_name == "myhappygenes"
        assert result.sample_id == "MHG000001"
        assert result.total_variants == 2016
        assert any(a.source == "clinvar" for a in result.annotations)
        assert any(a.source == "pharmgkb" for a in result.annotations)
        # ADR-0021: composite version reports both builds when annotator
        # manages both. Single-build instances collapse to a single part.
        clinvar_versions = [v for name, v in result.annotators_used if name == "clinvar"]
        assert clinvar_versions, "ClinVar annotator missing from used set"
        assert "20260101" in clinvar_versions[0]

    def test_annotator_connections_closed_after_run(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ):
        parser = MyHappyGenesParser()
        ann = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        run_analysis(mock_mhg_path, parser, [ann])
        # ExitStack closed every per-build connection.
        assert ann._conns == {}


class TestGRCh36FlushFailSafe:
    """GRCh36 non-confident detection must use GRCh36 as effective build.

    Issue #6: when detection points to GRCh36 but isn't confident
    (matched < inspected), the pipeline was falling back to header_build
    or GRCh37. This bypassed the ClinVar safety guard (no GRCh36 cache)
    and silently annotated GRCh36 data against GRCh37 coordinates.
    """

    def _grch36_variants(self, count=3, discordant=1):
        """Build variants: `count` at GRCh36 positions + `discordant` at junk positions."""
        variants = []
        grch36_rsids = [
            rsid for rsid, builds in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in builds
        ]
        for rsid in grch36_rsids[:count]:
            chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
            variants.append(
                Variant(rsid=rsid, chromosome=chrom, position=pos, allele1="A", allele2="A")
            )
        for rsid in grch36_rsids[count : count + discordant]:
            chrom, _ = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
            variants.append(
                Variant(rsid=rsid, chromosome=chrom, position=99999999, allele1="A", allele2="A")
            )
        return variants

    def test_non_confident_grch36_uses_grch36_effective(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_non_confident_grch36_with_header_grch37_still_uses_grch36(self):
        state = _BuildDetectionState(override=None, header_build=BUILD_GRCH37)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_confident_grch36_uses_grch36(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=0)
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_diagnostics_report_grch36_as_effective(self):
        state = _BuildDetectionState(override=None, header_build=None)
        variants = self._grch36_variants(count=3, discordant=1)
        for v in variants:
            state.feed(v)
        state.flush()
        diag = state.diagnostics()
        assert diag.effective_build == BUILD_GRCH36
        assert diag.detected_build == BUILD_GRCH36

    def test_buffer_limit_with_single_grch36_probe_uses_grch36(self):
        """Buffer-limit path must apply the same GRCh36 safety as flush().

        Real FTDNA GRCh36 files have 687K+ variants but only 1 probe SNP
        in the first 100K lines. The buffer-limit fallback must run
        detect_build and trigger the GRCh36 guard, not hard-fall-back to
        GRCh37.
        """
        state = _BuildDetectionState(override=None, header_build=BUILD_GRCH37)
        grch36_rsids = [
            rsid for rsid, builds in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in builds
        ]
        rsid = grch36_rsids[0]
        chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
        probe = Variant(rsid=rsid, chromosome=chrom, position=pos, allele1="A", allele2="A")
        state.feed(probe)
        filler = [
            Variant(rsid=f"rs9{i:06d}", chromosome="1", position=i, allele1="A", allele2="A")
            for i in range(_DETECTION_BUFFER_LIMIT)
        ]
        for v in filler:
            ready, _batch = state.feed(v)
            if ready:
                break
        assert state.effective_build == BUILD_GRCH36
        diag = state.diagnostics()
        assert diag.detected_build == BUILD_GRCH36


class TestChrPrefixBuildDetection:
    """GH #38: chr-prefix on contigs as a tertiary build-detection signal.

    GRCh38 callers (DeepVariant, DRAGEN, GATK HC) use ``chr1`` /
    ``chrX``; GRCh37 uses bare ``1`` / ``X``. When rsID matching and
    ``##assembly`` both fail, chr-prefix is the strongest remaining
    heuristic.
    """

    def test_chr_prefix_alone_infers_grch38(self) -> None:
        """No override, no header_build, no rsID matches — chr-prefix → GRCh38."""
        state = _BuildDetectionState(
            override=None,
            header_build=None,
            chr_prefix_observed=True,
        )
        state.flush()
        assert state.effective_build == "GRCh38"

    def test_no_chr_prefix_falls_back_to_grch37(self) -> None:
        """No chr-prefix signal — preserves the existing GRCh37 fallback."""
        state = _BuildDetectionState(
            override=None,
            header_build=None,
            chr_prefix_observed=False,
        )
        state.flush()
        assert state.effective_build == "GRCh37"

    def test_override_wins_over_chr_prefix(self) -> None:
        """``--build grch37`` (override) beats chr-prefix-says-GRCh38."""
        state = _BuildDetectionState(
            override="GRCh37",
            header_build=None,
            chr_prefix_observed=True,
        )
        assert state.effective_build == "GRCh37"

    def test_header_build_wins_over_chr_prefix(self) -> None:
        """Explicit header build beats chr-prefix heuristic."""
        state = _BuildDetectionState(
            override=None,
            header_build="GRCh37",
            chr_prefix_observed=True,
        )
        state.flush()
        assert state.effective_build == "GRCh37"

    def test_grch36_detected_still_wins(self) -> None:
        """GH #38 doesn't break the GRCh36 safety guard. Detected GRCh36
        from position data takes priority over the chr-prefix tertiary
        signal — wrong build assignment from chr-prefix bypassing the
        guard would silently mis-annotate."""
        variants = []
        grch36_rsids = [
            rsid for rsid, builds in KNOWN_SNP_POSITIONS.items() if BUILD_GRCH36 in builds
        ]
        for rsid in grch36_rsids[:3]:
            chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36]
            variants.append(
                Variant(rsid=rsid, chromosome=chrom, position=pos, allele1="A", allele2="A")
            )
        state = _BuildDetectionState(
            override=None,
            header_build=None,
            chr_prefix_observed=True,
        )
        for v in variants:
            state.feed(v)
        state.flush()
        assert state.effective_build == BUILD_GRCH36

    def test_diagnostics_flag_set_when_chr_prefix_picked_build(self) -> None:
        """GH #38: chr_prefix_inferred is True exactly when the chr-prefix
        signal is what picked the effective build — no override, no rsID
        detection, no header build, and the signal flipped the fallback
        from GRCh37 to GRCh38. The CLI uses this to print "Inferred
        GRCh38 from chr-prefixed contig names" instead of the
        blind-default warning."""
        state = _BuildDetectionState(
            override=None,
            header_build=None,
            chr_prefix_observed=True,
        )
        state.flush()
        assert state.effective_build == "GRCh38"
        assert state.diagnostics().chr_prefix_inferred is True

    def test_diagnostics_flag_clear_when_other_signal_picked_build(self) -> None:
        """chr_prefix_observed but a higher-priority signal won — flag
        stays False so the CLI doesn't surface a "inferred from chr-prefix"
        message when the build came from rsID detection or header."""
        # Header build wins.
        state = _BuildDetectionState(
            override=None,
            header_build="GRCh38",
            chr_prefix_observed=True,
        )
        state.flush()
        assert state.diagnostics().chr_prefix_inferred is False
        # Override wins.
        state = _BuildDetectionState(
            override="GRCh38",
            header_build=None,
            chr_prefix_observed=True,
        )
        assert state.diagnostics().chr_prefix_inferred is False

    def test_diagnostics_flag_clear_when_no_chr_prefix(self) -> None:
        """No chr-prefix signal → bare GRCh37 fallback → flag False."""
        state = _BuildDetectionState(
            override=None,
            header_build=None,
            chr_prefix_observed=False,
        )
        state.flush()
        assert state.effective_build == "GRCh37"
        assert state.diagnostics().chr_prefix_inferred is False


class TestGnomadEnrichment:
    """gnomAD frequency enrichment stamps allele_frequency on annotations."""

    def test_enrichment_stamps_frequency(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with gnomAD annotator stamps allele_frequency."""
        import sqlite3

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.annotators.gnomad import GnomadAnnotator
        from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
        from allelix.databases.schema import GNOMAD_SCHEMA

        db_path = clinvar_data_dir / GNOMAD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO gnomad_frequencies"
                " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
            )
            from allelix.databases._versions import GNOMAD_SCHEMA_VERSION

            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (
                    "gnomad",
                    "test://mock",
                    "4.1",
                    "2026-01-01T00:00:00Z",
                    1,
                    f"sv:{GNOMAD_SCHEMA_VERSION}",
                ),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        gnomad = cm_stack.enter_context(GnomadAnnotator(clinvar_data_dir))
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            gnomad=gnomad,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.allele_frequency is not None for a in mthfr)
        assert ("gnomad", "4.1") in result.annotators_used

    def test_no_gnomad_no_frequency(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis without gnomAD leaves allele_frequency as None."""
        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.allele_frequency is None for a in result.annotations)
        assert all(name != "gnomad" for name, _ in result.annotators_used)


class TestAlphaMissenseEnrichment:
    """AlphaMissense enrichment stamps am_pathogenicity/am_class on annotations."""

    def test_enrichment_stamps_pathogenicity(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with AlphaMissense stamps am_pathogenicity and am_class."""
        import sqlite3

        from allelix.annotators.alphamissense import AlphaMissenseAnnotator
        from allelix.databases.alphamissense_loader import ALPHAMISSENSE_DB_FILENAME
        from allelix.databases.schema import ALPHAMISSENSE_SCHEMA

        db_path = clinvar_data_dir / ALPHAMISSENSE_DB_FILENAME
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in ALPHAMISSENSE_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO alphamissense_scores"
                " (chrom, pos, ref, alt, rsid, uniprot_id, transcript_id,"
                " protein_variant, am_pathogenicity, am_class)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    "1",
                    11856378,
                    "G",
                    "A",
                    "rs1801133",
                    "P42898",
                    "ENST001",
                    "A222V",
                    0.72,
                    "ambiguous",
                ),
            )
            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count)"
                " VALUES (?, ?, ?, ?, ?)",
                ("alphamissense", "test://mock", "2023.2", "2026-01-01", 1),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        am = cm_stack.enter_context(AlphaMissenseAnnotator(clinvar_data_dir))
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            alphamissense=am,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.am_pathogenicity is not None for a in mthfr)
        assert any(a.am_class == "ambiguous" for a in mthfr)
        assert ("alphamissense", "2023.2") in result.annotators_used

    def test_no_alphamissense_no_pathogenicity(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis without AlphaMissense leaves am_pathogenicity as None."""
        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.am_pathogenicity is None for a in result.annotations)
        assert all(name != "alphamissense" for name, _ in result.annotators_used)


class TestEnrichmentExactVsFallback:
    """Pipeline splits exact (rsid, alt) lookups from MAX-by-rsid fallback."""

    def test_annotation_with_alt_uses_exact_lookup(self) -> None:
        """Annotations with alt set get exact (rsid, alt) enrichment."""
        from allelix.models import Annotation

        a = Annotation(
            source="clinvar",
            rsid="rs429358",
            significance="clinvar_pathogenic",
            category="clinical",
            magnitude=9.0,
            description="test",
            attribution="ClinVar",
            genotype_match="TC",
            alt="C",
        )
        assert a.alt == "C"
        exact_keys = {(x.rsid, x.alt) for x in [a] if x.alt}
        assert ("rs429358", "C") in exact_keys

    def test_annotation_without_alt_uses_max_fallback(self) -> None:
        """Annotations without alt fall back to MAX-by-rsid enrichment."""
        from allelix.models import Annotation

        a = Annotation(
            source="gwas",
            rsid="rs429358",
            significance="gwas_association",
            category="trait",
            magnitude=7.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="TC",
            alt="",
        )
        max_rsids = {x.rsid for x in [a] if not x.alt}
        assert "rs429358" in max_rsids

    def test_gwas_annotations_have_no_alt(self, cm_stack) -> None:
        """GWAS annotations must not set alt (risk allele != VCF ALT)."""
        import sqlite3
        import tempfile
        from pathlib import Path as _Path

        from allelix.annotators.gwas import GWASCatalogAnnotator
        from allelix.databases.gwas_loader import load_gwas_tsv
        from allelix.databases.schema import GWAS_SCHEMA

        # GH #45: these were previously `pytest.skip(...)` guards. The
        # mock fixture is committed (so the path-exists check should
        # never fail); is_ready() is exercised immediately after a
        # successful load (so it should never report not-ready). Silent
        # skips here would hide a real regression — fail loudly instead.
        db_path = _Path(__file__).parent.parent / "fixtures" / "mock_gwas_catalog.tsv"
        assert db_path.exists(), (
            f"committed mock fixture missing at {db_path} — "
            "did a git-clean / merge accidentally drop it?"
        )

        with tempfile.TemporaryDirectory() as td:
            tmp = _Path(td)
            gwas_db = tmp / "gwas.sqlite"
            with contextlib.closing(sqlite3.connect(gwas_db)) as conn:
                for stmt in GWAS_SCHEMA.split(";"):
                    stmt = stmt.strip()
                    if stmt:
                        conn.execute(stmt)
                conn.commit()
            load_gwas_tsv(db_path, gwas_db)
            ann = cm_stack.enter_context(GWASCatalogAnnotator(tmp))
            assert ann.is_ready(), (
                "GWAS annotator reported not-ready immediately after a "
                "successful load_gwas_tsv — loader or is_ready() has regressed"
            )
            v = Variant(
                rsid="rs1801133",
                chromosome="1",
                position=11796321,
                allele1="C",
                allele2="T",
            )
            results = ann.annotate(v)
            for r in results:
                assert r.alt == "", f"GWAS annotation should not set alt, got {r.alt!r}"
            ann.close()


class TestCaddEnrichment:
    """CADD enrichment stamps cadd_phred on annotations via gnomAD coordinate resolution."""

    def test_enrichment_stamps_cadd_phred(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """run_analysis with CADD + gnomAD stamps cadd_phred."""
        import sqlite3

        from allelix.annotators.cadd import CaddAnnotator
        from allelix.annotators.gnomad import GnomadAnnotator
        from allelix.databases.cadd_loader import CADD_DB_FILENAME
        from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
        from allelix.databases.schema import CADD_SCHEMA, GNOMAD_SCHEMA

        gnomad_path = clinvar_data_dir / GNOMAD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(gnomad_path)) as conn:
            for stmt in GNOMAD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT OR REPLACE INTO gnomad_frequencies"
                " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
            )
            from allelix.databases._versions import GNOMAD_SCHEMA_VERSION

            conn.execute(
                "INSERT OR REPLACE INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (
                    "gnomad",
                    "test://mock",
                    "4.1",
                    "2026-01-01T00:00:00Z",
                    1,
                    f"sv:{GNOMAD_SCHEMA_VERSION}",
                ),
            )
            conn.commit()

        cadd_path = clinvar_data_dir / CADD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(cadd_path)) as conn:
            for stmt in CADD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES (?, ?, ?, ?, ?)",
                ("1", 11796321, "G", "A", 24.3),
            )
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                ("cadd", "test://mock", "v1.7", "2026-01-01T00:00:00Z", 1, "sv:1"),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        gnomad = cm_stack.enter_context(GnomadAnnotator(clinvar_data_dir))
        cadd = cm_stack.enter_context(CaddAnnotator(clinvar_data_dir))
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            gnomad=gnomad,
            cadd=cadd,
        )
        mthfr = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert any(a.cadd_phred is not None for a in mthfr)
        assert any(a.cadd_phred == pytest.approx(24.3) for a in mthfr if a.cadd_phred is not None)
        assert ("cadd", "v1.7") in result.annotators_used

    def test_no_cadd_no_phred(self, cm_stack, mock_mhg_path: Path, clinvar_data_dir: Path) -> None:
        """run_analysis without CADD leaves cadd_phred as None."""
        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        result = run_analysis(mock_mhg_path, parser, [clinvar])
        assert all(a.cadd_phred is None for a in result.annotations)
        assert all(name != "cadd" for name, _ in result.annotators_used)

    def test_cadd_without_gnomad_skips_enrichment(
        self,
        cm_stack,
        mock_mhg_path: Path,
        clinvar_data_dir: Path,
    ) -> None:
        """CADD enrichment requires gnomAD for coordinate resolution."""
        import sqlite3

        from allelix.annotators.cadd import CaddAnnotator
        from allelix.databases.cadd_loader import CADD_DB_FILENAME
        from allelix.databases.schema import CADD_SCHEMA

        cadd_path = clinvar_data_dir / CADD_DB_FILENAME
        with contextlib.closing(sqlite3.connect(cadd_path)) as conn:
            for stmt in CADD_SCHEMA.split(";"):
                stmt = stmt.strip()
                if stmt:
                    conn.execute(stmt)
            conn.execute(
                "INSERT INTO database_versions"
                " (name, source_url, version, downloaded_at, record_count,"
                "  local_version_tag)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                ("cadd", "test://mock", "v1.7", "2026-01-01T00:00:00Z", 0, "sv:1"),
            )
            conn.commit()

        parser = MyHappyGenesParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        cadd = cm_stack.enter_context(CaddAnnotator(clinvar_data_dir))
        result = run_analysis(
            mock_mhg_path,
            parser,
            [clinvar],
            cadd=cadd,
        )
        assert all(a.cadd_phred is None for a in result.annotations)


class TestCaddMultiAllelic:
    """CADD enrichment at multi-allelic sites must use the user's allele, not max."""

    def test_multiallelic_uses_user_allele_not_max(self) -> None:
        """At a multi-allelic site, CADD score must match the user's allele."""
        from allelix.reports._pipeline import _lookup_user_allele

        coords = [("1", 100, "A", "C"), ("1", 100, "A", "G")]
        scores = {("1", 100, "A", "C"): 5.0, ("1", 100, "A", "G"): 30.0}

        result = _lookup_user_allele("C", coords, scores)
        assert result == pytest.approx(5.0)

    def test_multiallelic_no_alt_skips_enrichment(self) -> None:
        """GH #23 (suppress-half): annotations without an explicit alt
        (raw GWAS rows; SNPedia drug-response rows in some shapes) used
        to take a ``MAX(phred)`` fallback across all alts at the
        position. At multi-allelic sites that stamps the highest-CADD
        alt next to the annotation as if it described the user's
        variant. Same wrong-allele hazard as the strand path fixed in
        #18. Now: skip enrichment, leave ``cadd_phred=None``. The full
        fix (carrying the user's alt onto every Annotation so GWAS rows
        take the exact-alt path) is tracked for v2.1."""
        from allelix.models import Annotation
        from allelix.reports._pipeline import _enrich_cadd

        coords = [("1", 100, "A", "C"), ("1", 100, "A", "G")]
        scores = {("1", 100, "A", "C"): 5.0, ("1", 100, "A", "G"): 30.0}

        ann = Annotation(
            source="gwas",
            rsid="rs999",
            significance="gwas_association",
            category="trait",
            magnitude=9.0,
            description="test",
            attribution="GWAS Catalog",
            genotype_match="AC",
            alt="",  # GWAS rows always have alt=""
        )

        class MockGnomad:
            def bulk_resolve_coordinates(self, rsids):
                return {"rs999": coords}

        class MockCadd:
            def bulk_lookup(self, keys):
                return scores

        _enrich_cadd([ann], MockGnomad(), MockCadd())
        # Was 30.0 (the MAX at this position). Now None — we don't know
        # which alt the user carries, so we don't claim to.
        assert ann.cadd_phred is None

    def test_no_direct_match_returns_none(self) -> None:
        """GH #18: complement-fallback removed. An allele that doesn't
        directly match any alt at this position returns None — the
        previous behavior accepted ``complement(user_alt) == alt`` and
        coincidentally stamped a wrong-allele CADD score at multi-allelic
        sites."""
        from allelix.reports._pipeline import _lookup_user_allele

        coords = [("1", 200, "C", "A")]
        scores = {("1", 200, "C", "A"): 12.5}

        # "T" does not match alt "A" directly; under the old code,
        # complement("T") = "A" would have returned 12.5. Now: None.
        result = _lookup_user_allele("T", coords, scores)
        assert result is None

    def test_audit_reproduction_skips_enrichment(self) -> None:
        """GH #18 reproduction: a single-alt site (C → T) with the user
        carrying "A". complement("A") = "T" coincidentally matches the
        alt under the old code, stamping the C→T CADD score onto an
        annotation describing the user's "A" carrier. The fix returns
        None so the score is not stamped."""
        from allelix.reports._pipeline import _lookup_user_allele

        coords = [("1", 300, "C", "T")]
        scores = {("1", 300, "C", "T"): 25.0}

        result = _lookup_user_allele("A", coords, scores)
        assert result is None


class TestBatchedPipeline:
    """The two-phase batched pipeline produces the same results as the
    per-variant path and collects hv_variants in a single pass.

    Phase 1 (build detection) + Phase 2 (batched annotation) → results
    must be identical to what the old per-variant `annotate(v)` loop
    produced for the same inputs.
    """

    def test_hv_variants_collected_during_streaming(
        self,
        cm_stack,
        mock_mhg_path: Path,
        all_annotators_data_dir: Path,
    ):
        """High-value rsIDs accumulate on AnalysisResult.hv_variants.

        Eliminates the second parser.parse() pass that _run_analysis_command
        used to do for high-value no-call detection.
        """
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser

        parser = MyHappyGenesParser()
        ann = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        # Pick rsIDs known to exist in the mock_mhg fixture.
        hv_rsids = {"rs1801133", "rs4680"}
        try:
            result = run_analysis(
                mock_mhg_path,
                parser,
                [ann],
                high_value_rsids=hv_rsids,
            )
        finally:
            ann.close()
        collected_rsids = {v.rsid for v in result.hv_variants}
        assert collected_rsids == hv_rsids
        # Each hv rsid appears exactly once (parser yields each variant once).
        assert len(result.hv_variants) == 2

    def test_hv_variants_empty_when_none_requested(
        self,
        cm_stack,
        mock_mhg_path: Path,
        all_annotators_data_dir: Path,
    ):
        """No high_value_rsids → empty hv_variants list (don't collect everything)."""
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser

        parser = MyHappyGenesParser()
        ann = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            result = run_analysis(mock_mhg_path, parser, [ann])
        finally:
            ann.close()
        assert result.hv_variants == []

    def test_hv_variants_empty_for_unmatched_rsids(
        self,
        cm_stack,
        mock_mhg_path: Path,
        all_annotators_data_dir: Path,
    ):
        """Asking for rsIDs not in the file → empty hv_variants, no crash."""
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser

        parser = MyHappyGenesParser()
        ann = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            result = run_analysis(
                mock_mhg_path,
                parser,
                [ann],
                high_value_rsids={"rs_not_present_in_file_99999"},
            )
        finally:
            ann.close()
        assert result.hv_variants == []

    def test_batch_crosses_size_boundary_produces_same_results(
        self,
        cm_stack,
        mock_mhg_path: Path,
        all_annotators_data_dir: Path,
        monkeypatch,
    ):
        """Forcing a tiny _BATCH_SIZE must not change the annotation set.

        Patches the pipeline's batch constant to 7 so a typical fixture
        produces many batch boundaries. The result set must match a
        single-batch run (default size 5000, larger than the fixture).
        """
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser
        from allelix.reports import _pipeline as pipeline_mod

        parser = MyHappyGenesParser()

        # Single-batch baseline
        ann1 = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            baseline = run_analysis(mock_mhg_path, parser, [ann1])
        finally:
            ann1.close()

        # Many-batch run with a tiny batch size
        monkeypatch.setattr(pipeline_mod, "_BATCH_SIZE", 7)
        ann2 = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            tiny_batched = run_analysis(mock_mhg_path, parser, [ann2])
        finally:
            ann2.close()

        # Annotations must match exactly — batching is a performance
        # technique, not a semantic change.
        assert baseline.annotations == tiny_batched.annotations
        assert baseline.total_variants == tiny_batched.total_variants

    def test_batch_size_one_works(
        self,
        cm_stack,
        mock_mhg_path: Path,
        all_annotators_data_dir: Path,
        monkeypatch,
    ):
        """_BATCH_SIZE=1 (each variant flushes immediately) is correct.

        Edge case — pipeline must produce identical results regardless
        of batch size. Pinned because a future change could accidentally
        require size > 1 (e.g., empty-buffer guards).
        """
        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser
        from allelix.reports import _pipeline as pipeline_mod

        parser = MyHappyGenesParser()
        ann1 = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            baseline = run_analysis(mock_mhg_path, parser, [ann1])
        finally:
            ann1.close()

        monkeypatch.setattr(pipeline_mod, "_BATCH_SIZE", 1)
        ann2 = cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir))
        try:
            unit_batched = run_analysis(mock_mhg_path, parser, [ann2])
        finally:
            ann2.close()

        assert baseline.annotations == unit_batched.annotations


class TestRsidlessVcfResolution:
    """End-to-end: rsID-less VCFs from variant callers get annotated. GH #8.

    Real VCFs from GATK HaplotypeCaller, DeepVariant, etc. write `.` to the
    ID column. Before the resolver landed, every annotator returned zero
    hits because the rsID-keyed lookups had nothing to match. This class
    pins the round trip: parse a `.`-ID VCF, resolve through ClinVar by
    position, see annotations from every rsID-keyed source.
    """

    def test_rsidless_vcf_produces_annotations(self, cm_stack, clinvar_data_dir) -> None:
        """rsID-less VCF flows through resolution → ClinVar annotation appears."""
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        fixture_path = Path(__file__).parent.parent / "fixtures" / "mock_vcf_rsidless.vcf"
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                fixture_path,
                parser,
                [clinvar],
                build_override="GRCh37",
            )
        finally:
            clinvar.close()

        # MTHFR rs1801133 at chr1:11856378 G/A is in the fixture and in the
        # test ClinVar cache; resolution should let ClinVar annotate it.
        mthfr_hits = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert mthfr_hits, "Expected MTHFR rs1801133 annotation after rsID resolution"

    def test_unknown_position_leaves_rsid_empty(
        self,
        cm_stack,
        tmp_path,
        clinvar_data_dir,
    ) -> None:
        """Positions not in ClinVar produce no annotations and don't error out."""
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        # VCF containing only a position that doesn't exist in test ClinVar.
        vcf_path = tmp_path / "rsidless_unknown.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621,assembly=GRCh37>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t99999999\t.\tA\tT\t100\tPASS\t.\tGT\t0/1\n"
        )
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                Path(vcf_path),
                parser,
                [clinvar],
                build_override="GRCh37",
            )
        finally:
            clinvar.close()

        assert result.total_variants == 1
        assert result.annotations == []

    def test_positional_synthetic_id_treated_as_rsidless(
        self,
        cm_stack,
        tmp_path,
        clinvar_data_dir,
    ) -> None:
        """1000 Genomes-style positional IDs (`22:10519265:CA:C`) trigger resolution.

        The ID column is non-empty so a naive ``if not v.rsid`` filter would
        skip resolution, but those IDs are not real rsIDs — no rsID-keyed
        annotator can match them. The filter pivots on "doesn't start with
        rs" instead. GH #8.
        """
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        vcf_path = tmp_path / "synthetic_id.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621,assembly=GRCh37>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t11856378\t1:11856378:G:A\tG\tA\t100\tPASS\t.\tGT\t0/1\n"
        )
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                Path(vcf_path),
                parser,
                [clinvar],
                build_override="GRCh37",
            )
        finally:
            clinvar.close()
        assert any(a.rsid == "rs1801133" for a in result.annotations)

    def test_high_value_rsid_matches_after_resolution(self, cm_stack, clinvar_data_dir) -> None:
        """High-value rsID list catches variants resolved via ClinPGx.

        Without this contract: a user passes ``high_value_rsids={"rs1801133"}``,
        the input VCF has the variant as ``ID=.`` at chr1:11856378 (resolves to
        rs1801133), and ``hv_variants`` ends up empty because the legacy
        per-variant match in ``_accept`` ran with ``v.rsid == ""`` before
        resolution had a chance to mutate it. GH #11. The hv_set check runs at
        the tail of ``_flush()`` instead, after ``bulk_resolve_rsids`` has
        settled v.rsid for the whole batch.
        """
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        fixture_path = Path(__file__).parent.parent / "fixtures" / "mock_vcf_rsidless.vcf"
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                fixture_path,
                parser,
                [clinvar],
                build_override="GRCh37",
                high_value_rsids={"rs1801133", "rs80357906"},
            )
        finally:
            clinvar.close()

        collected = {v.rsid for v in result.hv_variants}
        # rs1801133 (chr1:11856378 G→A) and rs80357906 (chr17:41209080 G→A)
        # are both in the rsidless fixture under ID=. and both resolve via
        # the test ClinVar cache.
        assert "rs1801133" in collected
        assert "rs80357906" in collected
        # Variant identity preserved through resolution: the same Variant
        # object that ends up in hv_variants is the one whose rsid was
        # mutated by the resolver — so its rsid field is the post-resolution
        # value, not the empty string it came in with.
        for v in result.hv_variants:
            assert v.rsid.startswith("rs")

    def test_high_value_rsid_no_match_when_unresolved(self, cm_stack, clinvar_data_dir) -> None:
        """Negative case: rsID-less variant that doesn't resolve stays out of hv_variants.

        Pins the symmetry: post-resolution hv_set check must not stamp
        anything when resolution found nothing — v.rsid stays empty,
        ``"" in hv_set`` is False, the variant is not collected.
        """
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        # Position not in test ClinVar — resolution returns nothing, v.rsid
        # stays empty.
        vcf_path = clinvar_data_dir / "rsidless_unknown_pos.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=249250621,assembly=GRCh37>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t99999999\t.\tA\tT\t100\tPASS\t.\tGT\t0/1\n"
        )
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                Path(vcf_path),
                parser,
                [clinvar],
                build_override="GRCh37",
                high_value_rsids={"rs1801133"},
            )
        finally:
            clinvar.close()
        assert result.hv_variants == []

    def test_rsid_bearing_vcf_unchanged(self, cm_stack, clinvar_data_dir) -> None:
        """Regression: VCFs that already have rsIDs are unaffected by resolution.

        Resolution only triggers for variants whose ID column isn't an
        rs-prefixed identifier. Real rsIDs flow straight through.
        """
        from pathlib import Path

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.parsers.vcf import VcfParser

        fixture_path = Path(__file__).parent.parent / "fixtures" / "mock_vcf.vcf"
        parser = VcfParser()
        clinvar = cm_stack.enter_context(ClinVarAnnotator(clinvar_data_dir))
        try:
            result = run_analysis(
                fixture_path,
                parser,
                [clinvar],
                build_override="GRCh37",
            )
        finally:
            clinvar.close()
        # Pre-existing fixture has rsIDs; ClinVar should annotate based on
        # those — no resolution required.
        assert result.annotations
        for a in result.annotations:
            assert a.rsid.startswith("rs"), f"unexpected rsid shape: {a.rsid!r}"

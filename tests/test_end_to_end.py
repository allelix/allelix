# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""End-to-end integration: full analyze pipeline against the mock generators.

ADR-0015: the mock generators are the contract for what real data looks
like. The pipeline (parse → annotate → render) must produce a known,
human-vetted set of annotations against generator output. Any divergence
needs a code+ADR review, not a silent test update.

This is the gate that would have caught the v0.4.2 indel-anchor incident
and the v0.5.0 ClinPGx non-finding / somatic incidents BEFORE shipping,
had it existed earlier. Adding it now so future regressions can't sneak
through.
"""

from __future__ import annotations

import contextlib
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Iterator

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import run_analysis


def _analyze(mhg: Path, data_dir: Path):
    parser = MyHappyGenesParser()
    annotators = [ClinVarAnnotator(data_dir), PharmGKBAnnotator(data_dir)]
    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        return run_analysis(mhg, parser, bound)


class TestEndToEndAgainstMockGenerators:
    """The full analyze pipeline against the canonical mock generator output.

    These assertions are deliberate snapshots. If they fail, do NOT
    blindly bump the numbers — investigate whether a real regression has
    been introduced. Acceptable reasons to update:
      - a mock generator was deliberately extended (document in CHANGELOG)
      - an ADR was added that changes filtering behavior
    Never: "the test was wrong; let me change it to match."
    """

    def test_clinvar_does_not_fire_on_indel_anchor_against_array_data(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """ADR-0011 negative-path gate.

        The mock generator (correctly) reports rs113993960 as no-call
        ("-/-") because real MHG arrays can't call the CFTR ΔF508 indel.
        ClinVar has rs113993960 as REF=CTT ALT=C Pathogenic. The
        analyzer must NOT emit a CFTR finding here. Pre-v0.4.2 it did;
        this test pins that it doesn't anymore.
        """
        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        cftr_hits = [a for a in result.annotations if a.source == "clinvar" and a.gene == "CFTR"]
        assert not cftr_hits, (
            f"ClinVar emitted {len(cftr_hits)} CFTR annotation(s) against "
            f"a no-call MHG genotype. Indel-anchor protection (ADR-0011) "
            f"has regressed. Hits: {cftr_hits!r}"
        )

    def test_pharmgkb_does_not_emit_nonfindings(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """ADR-0013 / ADR-0016 negative-path gate.

        Reads the structured `is_nonfinding` column from the cache (set by
        the loader from ClinPGx's `Allele Function` field per ADR-0016).
        If any annotation surfaces from a row stored with is_nonfinding=1,
        the annotator's SELECT filter has regressed.
        """
        import sqlite3

        db = all_annotators_data_dir / "pharmgkb.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            nonfinding_rows = {
                (rsid, geno)
                for rsid, geno in conn.execute(
                    "SELECT rsid, genotype FROM pharmgkb_annotations WHERE is_nonfinding = 1"
                )
            }

        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        leaked = [
            a
            for a in result.annotations
            if a.source == "pharmgkb" and (a.rsid, a.genotype_match) in nonfinding_rows
        ]
        assert not leaked, (
            f"ClinPGx emitted {len(leaked)} non-finding row(s) at the analyzer "
            f"layer. The structured filter (is_nonfinding column) has regressed. "
            f"Examples: {[(a.rsid, a.genotype_match) for a in leaked[:3]]}"
        )

    def test_known_real_carriers_still_fire(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """Positive-path gate: the carriers the generator deliberately encodes
        must still produce annotations.
        """
        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        by_rsid = {a.rsid for a in result.annotations}

        # MTHFR C677T heterozygous (mock MHG: G/A) — ClinVar Pathogenic + ClinPGx LoE 2A
        assert "rs1801133" in by_rsid, "MTHFR heterozygous carrier vanished"
        # BRCA1 carrier (mock MHG: G/A) — ClinVar Pathogenic
        assert "rs80357906" in by_rsid, "BRCA1 carrier vanished"
        # CYP2C9*2 heterozygous (mock MHG: C/T) — ClinPGx LoE 1A
        assert "rs1799853" in by_rsid, "CYP2C9*2 carrier vanished"

    def test_homozygous_reference_does_not_fire_clinvar(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """TP53 in the mock is homozygous reference (G/G); ClinVar's row is
        REF=G ALT=T. The carrier rule (ADR-0007) requires the user to carry
        the ALT — wild-type homozygotes must not fire.
        """
        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        tp53_hits = [
            a for a in result.annotations if a.source == "clinvar" and a.rsid == "rs121918506"
        ]
        assert not tp53_hits, (
            f"ClinVar fired on TP53 homozygous-reference genotype. "
            f"Carrier rule (ADR-0007) has regressed. Hits: {tp53_hits!r}"
        )

    def test_nipa1_strand_inversion_no_emission_on_grch38_data(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """ADR-0021 regression: NIPA1 rs104894490 G/G on a GRCh38 file
        must NOT emit a pathogenic ClinVar annotation.

        Real-world bug history: user's MHG file labeled "build 37.1" was
        actually GRCh38. The annotator was correctly checking the
        carrier rule, but against the GRCh37 VCF whose REF/ALT are
        strand-inverted at this position (REF=C ALT=G). User's G/G
        matched ALT=G → false pathogenic call.

        The mock MHG default fixture uses GRCh38 positions (matching
        real MHG behavior). Build-detect identifies GRCh38; the
        ClinVar GRCh38 cache has REF=G ALT=A; the user carries zero
        A alleles → no annotation.

        Pinning the dispatch contract: same MHG, same user genotype,
        forcing --build grch37 would reproduce the OLD wrong behavior.
        That branch is exercised by
        `test_nipa1_grch37_dispatch_reproduces_legacy_false_positive`.
        """
        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        nipa1_hits = [
            a for a in result.annotations if a.source == "clinvar" and a.rsid == "rs104894490"
        ]
        assert not nipa1_hits, (
            "NIPA1 rs104894490 G/G must not emit a ClinVar annotation when the "
            "MHG file's positions are GRCh38. The strand-inverted REF/ALT in "
            "the GRCh37 VCF would produce a false positive — ADR-0021's "
            "auto-detection + per-build dispatch is supposed to prevent that. "
            f"Found {len(nipa1_hits)} unwanted annotation(s): {nipa1_hits!r}"
        )
        # Confirm detection actually picked GRCh38 (otherwise we'd be
        # testing the wrong branch).
        assert result.build_diagnostics is not None
        assert result.build_diagnostics.effective_build == "GRCh38"

    def test_nipa1_grch37_dispatch_reproduces_legacy_false_positive(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ):
        """Pinning the OTHER direction of the dispatch contract.

        With `--build grch37` forced on the (actually GRCh38) MHG file,
        the annotator dispatches to the GRCh37 ClinVar cache which has
        rs104894490 REF=C ALT=G. The user's G/G matches ALT=G → the
        pathogenic annotation DOES emit. This is the OLD wrong behavior;
        it's pinned here so a future "let's always default to GRCh37"
        regression would visibly flip this assertion.
        """
        from contextlib import ExitStack

        from allelix.annotators.clinvar import ClinVarAnnotator
        from allelix.annotators.pharmgkb import PharmGKBAnnotator
        from allelix.parsers.myhappygenes import MyHappyGenesParser
        from allelix.reports._pipeline import run_analysis

        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        with ExitStack() as stack:
            bound = [stack.enter_context(a) for a in annotators]
            result = run_analysis(mock_mhg_path, parser, bound, build_override="GRCh37")

        nipa1_hits = [
            a for a in result.annotations if a.source == "clinvar" and a.rsid == "rs104894490"
        ]
        assert nipa1_hits, (
            "NIPA1 rs104894490 G/G under --build grch37 should reproduce the "
            "legacy cross-build false positive (REF=C ALT=G in GRCh37; user's "
            "G matches ALT=G). This pin guards against accidentally changing "
            "the per-build dispatch contract."
        )
        assert result.build_diagnostics is not None
        assert result.build_diagnostics.effective_build == "GRCh37"
        assert result.build_diagnostics.override is True

    def test_annotation_count_snapshot(self, mock_mhg_path: Path, all_annotators_data_dir: Path):
        """Snapshot the exact count of emitted annotations.

        Any change here demands code review — either a new fixture row
        was added (update with explanation) or a regression introduced
        an extra hit (fix the code).
        """
        result = _analyze(mock_mhg_path, all_annotators_data_dir)
        clinvar_count = sum(1 for a in result.annotations if a.source == "clinvar")
        pharmgkb_count = sum(1 for a in result.annotations if a.source == "pharmgkb")
        # If these numbers change, investigate WHY before updating.
        # Verified 2026-05-19 after benign-filter default (ADR-0008 amendment):
        #   clinvar (6): rs1065852, rs1799853, rs1801133, rs4149056, rs4680,
        #     rs80357906. rs1801394 (Likely_benign) now suppressed by default.
        #   pharmgkb (5): rs1799853, rs1801133, rs4149056, rs4680,
        #     rs900000020 (PA-010: non-CPIC gene, no ClinVar/CPIC data → emits
        #     per ADR-0022; added for ClinVar REF regression coverage)
        # Notably absent:
        #   rs113993960 CFTR (no-call in MHG fixture — ADR-0011 pin)
        #   rs121918506 TP53 (hom-ref in MHG fixture — ADR-0007 pin)
        #   rs1801265 DPYD (not in MHG fixture; regression-only in unit tests)
        #   rs1801394 MTRR (Likely_benign — filtered by benign suppression)
        assert (clinvar_count, pharmgkb_count) == (6, 5), (
            f"End-to-end annotation count drifted: clinvar={clinvar_count}, "
            f"pharmgkb={pharmgkb_count}. Expected (6, 5). Investigate whether "
            f"a fixture row was added (update the snapshot with explanation) "
            f"or a regression introduced an extra hit (fix the code, don't "
            f"update the snapshot)."
        )


def _analyze_all(mhg: Path, data_dir: Path, *, gwas_filter_traits: bool = False):
    """Full pipeline with all three annotators (ClinVar + ClinPGx + GWAS)."""
    parser = MyHappyGenesParser()
    annotators = [
        ClinVarAnnotator(data_dir),
        PharmGKBAnnotator(data_dir),
        GWASCatalogAnnotator(data_dir, filter_traits=gwas_filter_traits),
    ]
    with contextlib.ExitStack() as stack:
        bound = [stack.enter_context(a) for a in annotators]
        return run_analysis(mhg, parser, bound)


class TestGWASEndToEnd:
    """End-to-end GWAS Catalog annotation against the canonical mock fixture.

    Pinned snapshots like TestEndToEndAgainstMockGenerators — investigate
    before updating numbers.
    """

    def test_gwas_carrier_fires(self, mock_mhg_path: Path, all_annotators_data_dir: Path) -> None:
        """Known carriers in MHG fixture produce GWAS annotations."""
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir)
        gwas_rsids = {a.rsid for a in result.annotations if a.source == "gwas"}
        assert "rs1801133" in gwas_rsids, "MTHFR GWAS carrier vanished"
        assert "rs4680" in gwas_rsids, "COMT GWAS carrier vanished"
        assert "rs1799853" in gwas_rsids, "CYP2C9 GWAS carrier vanished"
        assert "rs80357906" in gwas_rsids, "BRCA1 GWAS carrier vanished"
        assert "rs4149056" in gwas_rsids, "SLCO1B1 GWAS carrier vanished"

    def test_gwas_hom_ref_does_not_fire(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """TP53 rs121918506 G/G — risk allele T, user doesn't carry it."""
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir)
        tp53_gwas = [
            a for a in result.annotations if a.source == "gwas" and a.rsid == "rs121918506"
        ]
        assert not tp53_gwas, (
            f"GWAS fired on TP53 homozygous-reference genotype. "
            f"Carrier rule (ADR-0007) has regressed. Hits: {tp53_gwas!r}"
        )

    def test_gwas_annotation_count_snapshot(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """Pin GWAS annotation count. See test_annotation_count_snapshot
        for the investigation protocol.
        """
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir)
        gwas_count = sum(1 for a in result.annotations if a.source == "gwas")
        # Verified 2026-05-19 against mock_gwas_catalog.tsv (8 records
        # after dedup) x mock_myhappygenes.txt genotypes:
        #   rs1801133 G/A: Homocysteine levels (risk A, fires) +
        #                  Height (risk ?, fires capped at 3.0) = 2
        #   rs4680 A/A:    Pain sensitivity (risk A, fires) +
        #                  Blood pressure (risk A, fires) = 2
        #   rs1799853 C/T: Warfarin dose (risk T, fires) = 1
        #   rs80357906 G/A: Breast cancer (risk A, fires) = 1
        #   rs4149056 T/C: Statin myopathy (risk C, fires) = 1
        #   rs121918506 G/G: Synthetic cancer (risk T, does NOT fire) = 0
        # Total: 7
        assert gwas_count == 7, (
            f"GWAS annotation count drifted: got {gwas_count}, expected 7. "
            f"Investigate whether a fixture row was added or a regression "
            f"introduced an extra hit."
        )

    def test_gwas_filtered_count_snapshot(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """Pin GWAS count with default trait filtering (ADR-0024).

        Trait filtering excludes measurement/behavioral categories:
          - Homocysteine levels (other_measurement) -> excluded
          - Height (body_measurement) -> excluded
          - Pain sensitivity (other_measurement) -> excluded
          - Blood pressure (other_measurement) -> excluded
        Remaining: Warfarin dose (drug_response), Breast cancer (cancer),
          Statin-induced myopathy (disease via MONDO_) = 3
        """
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir, gwas_filter_traits=True)
        gwas_count = sum(1 for a in result.annotations if a.source == "gwas")
        assert gwas_count == 3, f"Filtered GWAS count drifted: got {gwas_count}, expected 3."

    def test_clinvar_pharmgkb_counts_unchanged_with_gwas(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """Adding GWAS must not perturb ClinVar or ClinPGx counts."""
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir)
        clinvar_count = sum(1 for a in result.annotations if a.source == "clinvar")
        pharmgkb_count = sum(1 for a in result.annotations if a.source == "pharmgkb")
        assert (clinvar_count, pharmgkb_count) == (6, 5), (
            f"ClinVar/ClinPGx counts changed when GWAS annotator was added: "
            f"clinvar={clinvar_count}, pharmgkb={pharmgkb_count}. Expected (6, 5)."
        )


class TestBenignSuppressionEndToEnd:
    """ADR-0008 amendment: benign annotations suppressed by default."""

    def test_include_benign_restores_likely_benign(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """include_benign=True restores the rs1801394 Likely_benign row."""
        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir, include_benign=True),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        with contextlib.ExitStack() as stack:
            bound = [stack.enter_context(a) for a in annotators]
            result = run_analysis(mock_mhg_path, parser, bound)
        clinvar_count = sum(1 for a in result.annotations if a.source == "clinvar")
        assert clinvar_count == 7, (
            f"include_benign=True should restore rs1801394 Likely_benign. "
            f"Got clinvar_count={clinvar_count}, expected 7."
        )


class TestDefaultReportSanity:
    """Pin that default-invocation annotation count stays below a threshold.

    The mock fixture is small (8 GWAS records, 13 ClinVar records per build,
    16 ClinPGx records), but this test catches regressions where defaults
    let too many annotations through — the same class of defect that let
    95,509 rows through on real data.
    """

    def test_default_rendered_count_below_threshold(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """Default analyze (no flags) on canonical mocks stays under 20."""
        result = _analyze_all(mock_mhg_path, all_annotators_data_dir)
        rendered = result.filter(
            min_magnitude=5.0,
            source_min_magnitudes={"gwas": 9.0},
        )
        # With defaults (mag >= 5.0, gwas mag >= 9.0, no benign):
        #   ClinVar: 6 (see annotation_count_snapshot)
        #   minus those < 5.0: keeps only Pathogenic/Likely_path/Drug_response/Risk_factor
        #   ClinPGx: 5 minus those < 5.0: keeps LoE 1A/2A/2B
        #   GWAS: 7 minus those < 9.0: only hyper-significant + large OR pass
        # Total rendered should be well under 20.
        assert len(rendered) <= 20, (
            f"Default-invocation rendered {len(rendered)} annotations, "
            f"exceeding the 20-row sanity threshold. If a fixture was "
            f"added, update this threshold with justification."
        )


class TestMethylationSanity:
    """Methylation command produces a focused, bounded output."""

    def test_methylation_without_gwas_below_threshold(
        self, mock_mhg_path: Path, all_annotators_data_dir: Path
    ) -> None:
        """Methylation report (ClinVar + ClinPGx only) stays under 20 rows.

        GWAS is excluded from methylation by default. The remaining
        sources (ClinVar clinical + ClinPGx drug response) should produce
        a focused, interpretable set of methylation-pathway annotations.
        """
        parser = MyHappyGenesParser()
        annotators = [
            ClinVarAnnotator(all_annotators_data_dir),
            PharmGKBAnnotator(all_annotators_data_dir),
        ]
        from allelix.reports.methylation import METHYLATION_PANEL_GENES

        with contextlib.ExitStack() as stack:
            bound = [stack.enter_context(a) for a in annotators]
            result = run_analysis(mock_mhg_path, parser, bound)

        rendered = result.filter(min_magnitude=5.0, genes=METHYLATION_PANEL_GENES)
        assert len(rendered) < 20, (
            f"Methylation report (no GWAS) rendered {len(rendered)} rows, "
            f"exceeding the 20-row threshold. Investigate whether new fixture "
            f"rows were added or the gene panel / magnitude filter drifted."
        )


_REAL_GWAS_ZIP = Path(__file__).resolve().parent.parent / "test_data" / "gwas_catalog.zip"


@pytest.mark.slow
class TestRealDataGwasSanity:
    """Sanity checks against the real GWAS Catalog (test_data/, gitignored).

    Skipped when the real data hasn't been downloaded. To populate:
        curl -L -o test_data/gwas_catalog.zip \\
          "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/\\
          gwas-catalog-associations_ontology-annotated-full.zip"
    """

    @pytest.fixture(scope="class")
    def real_gwas_data_dir(self, tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
        """Build the real GWAS Catalog once, in temp, and delete it when done."""
        if not _REAL_GWAS_ZIP.exists():
            pytest.skip("Real GWAS Catalog not downloaded (test_data/gwas_catalog.zip)")

        import zipfile

        from allelix.databases.gwas_loader import load_gwas_tsv

        tmp_path = tmp_path_factory.mktemp("real_gwas")
        with zipfile.ZipFile(_REAL_GWAS_ZIP) as zf:
            tsv_names = [n for n in zf.namelist() if n.endswith(".tsv")]
            assert tsv_names, "No TSV found in GWAS ZIP"
            zf.extract(tsv_names[0], tmp_path)
        tsv_path = tmp_path / tsv_names[0]
        db_path = tmp_path / "gwas.sqlite"
        try:
            load_gwas_tsv(tsv_path, db_path, source_url="test://real-gwas")
        finally:
            tsv_path.unlink(missing_ok=True)
        yield tmp_path
        shutil.rmtree(tmp_path, ignore_errors=True)

    def test_default_gwas_floor_keeps_output_bounded(
        self, mock_mhg_path: Path, real_gwas_data_dir: Path
    ) -> None:
        """Default filters (mag>=5, gwas>=9) on real GWAS data stay under 50 rows.

        Verified 2026-05-19: mock MHG x full GWAS Catalog (795k records)
        produces 331 raw annotations. With defaults, only 7 pass (all mag 9).
        Threshold set at 50 to allow for GWAS Catalog growth.
        """
        parser = MyHappyGenesParser()
        ann = GWASCatalogAnnotator(real_gwas_data_dir)
        with contextlib.ExitStack() as stack:
            bound = [stack.enter_context(ann)]
            result = run_analysis(mock_mhg_path, parser, bound)

        rendered = result.filter(
            min_magnitude=5.0,
            source_min_magnitudes={"gwas": 9.0},
        )
        assert len(rendered) <= 50, (
            f"Default GWAS floor (9.0) let {len(rendered)} annotations through "
            f"on real data. Expected ≤50. The floor may need raising, or the "
            f"magnitude formula produces too many high scores."
        )

    def test_old_floor_would_have_been_unmanageable(
        self, mock_mhg_path: Path, real_gwas_data_dir: Path
    ) -> None:
        """The old floor (7.0) produces significantly more output — regression guard."""
        parser = MyHappyGenesParser()
        ann = GWASCatalogAnnotator(real_gwas_data_dir)
        with contextlib.ExitStack() as stack:
            bound = [stack.enter_context(ann)]
            result = run_analysis(mock_mhg_path, parser, bound)

        old_floor = result.filter(
            min_magnitude=5.0,
            source_min_magnitudes={"gwas": 7.0},
        )
        new_floor = result.filter(
            min_magnitude=5.0,
            source_min_magnitudes={"gwas": 9.0},
        )
        assert len(old_floor) > len(new_floor), (
            "Expected old floor (7.0) to produce more annotations than new (9.0)"
        )

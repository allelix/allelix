# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the per-SCV ClinVar TSV loader (#42 stage A).

The loader joins ``submission_summary.txt.gz`` (per-SCV rows) against
``variant_summary.txt.gz`` (per-(VariationID, Assembly) rows) and emits
one cache record per submission. Stage A added the loader, stage B
wired it into ``db update`` with an interpreter-version bump, stage C
(#42) removed the legacy VCF path entirely — this is now the only
production ingest route.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.databases.manager import (
    _decode_reported_phenotype,
    iter_clinvar_tsv_records,
    load_clinvar_tsv,
)

if TYPE_CHECKING:
    from pathlib import Path


# Synthetic variant_summary.txt header + rows. Column count matches the
# real file (43); only the columns the loader reads need to carry
# meaningful values. Tabs are mandatory.
_VS_HEADER = (
    "#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance"
    "\tClinSigSimple\tLastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)"
    "\tRCVaccession\tPhenotypeIDS\tPhenotypeList\tOrigin\tOriginSimple"
    "\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele"
    "\tAlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters"
    "\tGuidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\tVariationID"
    "\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF"
    "\tSomaticClinicalImpact\tSomaticClinicalImpactLastEvaluated"
    "\tReviewStatusClinicalImpact\tOncogenicity\tOncogenicityLastEvaluated"
    "\tReviewStatusOncogenicity\tSCVsForAggregateGermlineClassification"
    "\tSCVsForAggregateSomaticClinicalImpact"
    "\tSCVsForAggregateOncogenicityClassification"
)


def _vs_row(
    *,
    allele_id: str = "1000",
    gene: str = "BRCA1",
    rs: str = "12345",
    assembly: str = "GRCh38",
    chrom: str = "17",
    position_vcf: str = "41197737",
    ref: str = "G",
    alt: str = "A",
    variation_id: str = "100",
    review_status: str = "criteria provided, single submitter",
) -> str:
    """Build one tab-separated variant_summary row at the right column positions."""
    cols = [""] * 43
    cols[0] = allele_id
    cols[1] = "single nucleotide variant"
    cols[2] = f"NM_007294.3({gene}):c.5266dupC"
    cols[3] = "672"
    cols[4] = gene
    cols[5] = "HGNC:1100"
    cols[6] = "Pathogenic"
    cols[7] = "1"
    cols[8] = "Jan 01, 2024"
    cols[9] = rs
    cols[10] = "-"
    cols[11] = "RCV00001"
    cols[12] = "MedGen:C0006142"
    cols[13] = "Hereditary breast and ovarian cancer syndrome"
    cols[14] = "germline"
    cols[15] = "germline"
    cols[16] = assembly
    cols[17] = "NC_000017.10"
    cols[18] = chrom
    cols[19] = position_vcf
    cols[20] = position_vcf
    cols[21] = ref
    cols[22] = alt
    cols[23] = "17q21.31"
    cols[24] = review_status
    cols[25] = "3"
    cols[26] = "-"
    cols[27] = "N"
    cols[28] = "ClinGen:CA000001"
    cols[29] = "3"
    cols[30] = variation_id
    cols[31] = position_vcf
    cols[32] = ref
    cols[33] = alt
    return "\t".join(cols)


_SS_HEADER = (
    "#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription"
    "\tSubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus"
    "\tCollectionMethod\tOriginCounts\tSubmitter\tSCV"
    "\tSubmittedGeneSymbol\tExplanationOfInterpretation"
    "\tSomaticClinicalImpact\tOncogenicity"
    "\tContributesToAggregateClassification"
)


def _ss_row(
    *,
    variation_id: str = "100",
    significance: str = "Pathogenic",
    reported_phenotype: str = "C0006142:Hereditary breast and ovarian cancer syndrome",
    review_status: str = "criteria provided, single submitter",
    scv: str = "SCV000001.1",
    contributes: str = "yes",
) -> str:
    cols = [""] * 16
    cols[0] = variation_id
    cols[1] = significance
    cols[2] = "Dec 17, 2024"
    cols[3] = "-"
    cols[4] = "BRCA1 c.5266dupC"
    cols[5] = reported_phenotype
    cols[6] = review_status
    cols[7] = "clinical testing"
    cols[8] = "germline:1"
    cols[9] = "Test Submitter"
    cols[10] = scv
    cols[11] = "BRCA1"
    cols[12] = "-"
    cols[13] = "-"
    cols[14] = "-"
    cols[15] = contributes
    return "\t".join(cols)


def _write_tsv(path: Path, header: str, rows: list[str]) -> Path:
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


# ---------------- _decode_reported_phenotype unit tests ----------------


class TestDecodeReportedPhenotype:
    def test_strips_single_medgen_prefix(self):
        assert (
            _decode_reported_phenotype("C3150901:Hereditary spastic paraplegia 48")
            == "Hereditary spastic paraplegia 48"
        )

    def test_strips_na_prefix(self):
        assert _decode_reported_phenotype("na:not provided") == "not provided"

    def test_joins_multiple_conditions_with_semicolon(self):
        assert (
            _decode_reported_phenotype("C1:Condition A|C2:Condition B|C3:Condition C")
            == "Condition A; Condition B; Condition C"
        )

    def test_skips_bare_medgen_id_without_colon(self):
        # "C1" alone (no Name half) is dropped — emitting it would surface
        # a bare ID that looks like noise in the report.
        assert _decode_reported_phenotype("C1|C2:Cond A|C3:Cond B") == "Cond A; Cond B"

    def test_empty_input_returns_empty(self):
        assert _decode_reported_phenotype("") == ""

    def test_dash_input_returns_empty(self):
        assert _decode_reported_phenotype("-") == ""

    def test_skips_dash_name_half(self):
        assert _decode_reported_phenotype("C1:-") == ""

    def test_strips_whitespace(self):
        assert (
            _decode_reported_phenotype("  C1:  Padded condition  | C2: Other ")
            == "Padded condition; Other"
        )


# ---------------- iter_clinvar_tsv_records tests ----------------


class TestIterTsvRecords:
    def test_single_variant_single_scv_emits_one_record(self, tmp_path):
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        r = records[0]
        assert r["rsid"] == "rs12345"
        assert r["chromosome"] == "17"
        assert r["position"] == 41197737
        assert r["ref"] == "G"
        assert r["alt"] == "A"
        assert r["clinical_significance"] == "Pathogenic"
        assert r["condition"] == "Hereditary breast and ovarian cancer syndrome"
        assert r["gene"] == "BRCA1"
        assert r["allele_id"] == 100  # variation_id repurposed

    def test_multiple_scvs_per_variant_emit_multiple_records(self, tmp_path):
        """The whole point of #42: per-SCV pairing instead of one row per variant."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(
                    significance="Pathogenic",
                    reported_phenotype="C1:Hereditary BRCA1 syndrome",
                    scv="SCV000001.1",
                ),
                _ss_row(
                    significance="Likely pathogenic",
                    reported_phenotype="C2:Macular dystrophy",
                    scv="SCV000002.1",
                ),
                _ss_row(
                    significance="risk_factor",
                    reported_phenotype="C3:Increased cancer risk",
                    scv="SCV000003.1",
                ),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 3
        # Per-SCV pairing — each (significance, condition) pair came from
        # one submission, not Frankensteined across submissions.
        pairs = {(r["clinical_significance"], r["condition"]) for r in records}
        assert ("Pathogenic", "Hereditary BRCA1 syndrome") in pairs
        assert ("Likely pathogenic", "Macular dystrophy") in pairs
        assert ("risk_factor", "Increased cancer risk") in pairs

    def test_build_filter_grch37_drops_grch38_rows(self, tmp_path):
        """variant_summary has per-(VariationID, Assembly) rows. Build filter
        selects the requested Assembly only."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(assembly="GRCh37", chrom="17", position_vcf="41245466"),
                _vs_row(assembly="GRCh38", chrom="17", position_vcf="41197737"),
            ],
        )
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])

        records_37 = list(iter_clinvar_tsv_records(vs, ss, "GRCh37"))
        records_38 = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records_37) == 1
        assert records_37[0]["position"] == 41245466
        assert len(records_38) == 1
        assert records_38[0]["position"] == 41197737

    def test_invalid_build_raises(self, tmp_path):
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        with pytest.raises(ValueError, match="unsupported build"):
            list(iter_clinvar_tsv_records(vs, ss, "GRCh99"))

    def test_skips_rows_without_rsid(self, tmp_path):
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(rs="-1", variation_id="100"),  # no dbSNP id
                _vs_row(rs="500", variation_id="200"),  # valid
            ],
        )
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [_ss_row(variation_id="100"), _ss_row(variation_id="200")],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["rsid"] == "rs500"

    def test_skips_rows_with_na_ref_alt(self, tmp_path):
        """Complex / copy-number variants encode ref/alt as 'na' — the
        SNV-shaped cache can't represent them. Skip without warning."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(ref="na", alt="na", variation_id="100"),
                _vs_row(ref="G", alt="A", variation_id="200"),
            ],
        )
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [_ss_row(variation_id="100"), _ss_row(variation_id="200")],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["allele_id"] == 200

    def test_skips_scv_without_variant(self, tmp_path):
        """SCV references a VariationID with no corresponding row in
        variant_summary (different assembly, or filtered upstream).
        The submission_summary pass silently skips."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row(variation_id="100")])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [_ss_row(variation_id="100"), _ss_row(variation_id="9999")],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["allele_id"] == 100

    def test_aggregate_only_filters_contributes_no(self, tmp_path):
        """Default `aggregate_only=True` drops SCVs flagged
        ContributesToAggregateClassification=no."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(scv="SCV001.1", contributes="yes"),
                _ss_row(scv="SCV002.1", contributes="no"),
                _ss_row(scv="SCV003.1", contributes="yes"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 2

    def test_aggregate_only_false_keeps_everything(self, tmp_path):
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(contributes="yes"),
                _ss_row(contributes="no"),
                _ss_row(contributes="yes"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38", aggregate_only=False))
        assert len(records) == 3

    def test_skips_scv_without_significance(self, tmp_path):
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(significance=""),
                _ss_row(significance="Pathogenic"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["clinical_significance"] == "Pathogenic"

    def test_skips_all_placeholder_clnsigs(self, tmp_path):
        """GH #42 follow-up (evaluator defect 5 + cross-PR review):
        ClinVar's submission_summary carries multiple placeholder
        values meaning "no classification recorded." The loader must
        drop every placeholder NOT in ``_CLNSIG_MAGNITUDE`` (those
        would otherwise default to 5.0, the analyze display floor,
        and surface as bogus annotations).

        Pin every value in ``_CLINVAR_PLACEHOLDER_CLNSIGS`` individually
        so if any one is removed from the skip set this test fires;
        and pin that ``not provided`` (which IS in _CLNSIG_MAGNITUDE
        at 2.0) is *kept*, so we don't lose real submitter records
        that just happened to look placeholder-shaped."""
        from allelix.databases.manager import _CLINVAR_PLACEHOLDER_CLNSIGS

        # Every placeholder gets dropped; only the real Pathogenic row
        # and the safe "not provided" row survive.
        placeholder_rows = [_ss_row(significance=p) for p in _CLINVAR_PLACEHOLDER_CLNSIGS]
        keep_rows = [
            _ss_row(significance="Pathogenic"),
            _ss_row(significance="not provided"),  # maps to 2.0, kept
        ]

        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            placeholder_rows + keep_rows,
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        sigs = sorted(r["clinical_significance"] for r in records)
        assert sigs == ["Pathogenic", "not provided"]
        # And the placeholder set must be the exact set used by the
        # protocol §7b ship-gate. If you broaden one, broaden both
        # together — that's the cross-PR-review lesson.
        # GH #116: the set grew with "other", "association",
        # "association not found" — non-classification curatorial
        # terms surfaced by the #42 per-SCV switch that fell to the
        # 5.0 default before this filter.
        expected = frozenset(
            {
                "",
                "-",
                "not specified",
                "no classification provided",
                "other",
                "association",
                "association not found",
            }
        )
        assert expected == _CLINVAR_PLACEHOLDER_CLNSIGS

    def test_skips_junk_clnsig_curatorial_terms(self, tmp_path):
        """GH #116: ClinVar's per-SCV TSV surfaces non-classification
        curatorial terms — "other", "association", "association not
        found" — that the old summarized data hid. Each is a real
        ClinVar value (a curator chose it) but is not a classification,
        is not in _CLNSIG_MAGNITUDE, and would surface at the 5.0
        default as a fake finding.

        Pinned: ingest drops these three terms while real
        classifications survive."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(significance="other"),
                _ss_row(significance="association"),
                _ss_row(significance="association not found"),
                _ss_row(significance="Pathogenic"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        sigs = sorted(r["clinical_significance"] for r in records)
        assert sigs == ["Pathogenic"]

    def test_placeholder_skip_is_case_insensitive(self, tmp_path):
        """The skip filter must match regardless of casing because
        `_magnitude()` normalizes via `_normalize_clnsig()` (which
        lowercases). A case-sensitive filter would let "Not Specified"
        slip past the skip and still land at the 5.0 default — silent
        regression of the same defect 5 surface bug.

        `-` (the confirmed real defect) is moot here, but the prose
        placeholders ("not specified", "no classification provided")
        were added by reasoning, not confirmed casing in real data —
        so this pins the case-insensitive contract before any
        real-data variation surfaces."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(significance="Not Specified"),
                _ss_row(significance="NOT SPECIFIED"),
                _ss_row(significance="No Classification Provided"),
                _ss_row(significance="NO CLASSIFICATION PROVIDED"),
                _ss_row(significance="Pathogenic"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        sigs = [r["clinical_significance"] for r in records]
        # All four casing variants of the placeholders dropped;
        # only the real classification survives.
        assert sigs == ["Pathogenic"]

    def test_per_scv_review_status_preferred_over_aggregate(self, tmp_path):
        """submission_summary's ReviewStatus is per-SCV — more specific
        than variant_summary's aggregate. The loader prefers it."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [_vs_row(review_status="criteria provided, multiple submitters")],
        )
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [_ss_row(review_status="reviewed by expert panel")],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert records[0]["review_status"] == "reviewed by expert panel"

    def test_short_variant_row_skipped_gracefully(self, tmp_path):
        """A truncated variant_summary row (fewer columns than expected)
        is dropped without crashing."""
        vs = tmp_path / "vs.tsv"
        vs.write_text(_VS_HEADER + "\nfoo\tbar\n" + _vs_row() + "\n", encoding="utf-8")
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1

    def test_comment_lines_ignored(self, tmp_path):
        """Lines starting with `#` (the real-file header explanations) are
        skipped by `_open_tsv`."""
        vs = tmp_path / "vs.tsv"
        body = "\n".join(
            [
                "##explanatory comment 1",
                "##explanatory comment 2",
                _VS_HEADER,
                _vs_row(),
            ]
        )
        vs.write_text(body + "\n", encoding="utf-8")
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1

    def test_skips_variant_with_invalid_position(self, tmp_path):
        """variant_summary PositionVCF not parseable as int → row dropped."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(position_vcf="not-a-number", variation_id="100"),
                _vs_row(position_vcf="200", variation_id="200"),
            ],
        )
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [_ss_row(variation_id="100"), _ss_row(variation_id="200")],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["allele_id"] == 200

    def test_skips_variant_with_empty_position(self, tmp_path):
        """Empty PositionVCF (e.g. structural variants) → row dropped before
        the int() parse path."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [_vs_row(position_vcf="", variation_id="100")],
        )
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row(variation_id="100")])
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert records == []

    def test_skips_variant_with_invalid_variation_id(self, tmp_path):
        """variant_summary VariationID not parseable → row dropped after
        ref/alt validation but before joining."""
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(variation_id="not-an-int"),
                _vs_row(variation_id="200"),
            ],
        )
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row(variation_id="200")])
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1
        assert records[0]["allele_id"] == 200

    def test_skips_short_submission_row(self, tmp_path):
        """submission_summary row truncated below the column the loader reads
        (e.g. _SS_CONTRIBUTES at col 15) — dropped without crashing."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = tmp_path / "ss.tsv"
        # First row has only 5 cols — way short of the 16-col schema.
        ss.write_text(
            _SS_HEADER + "\nfoo\tbar\tbaz\tqux\tquux\n" + _ss_row() + "\n", encoding="utf-8"
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1  # only the well-formed row survived

    def test_skips_submission_with_invalid_variation_id(self, tmp_path):
        """submission_summary VariationID not parseable → SCV row dropped
        before the temp-cache lookup."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row(variation_id="100")])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(variation_id="not-an-int"),
                _ss_row(variation_id="100"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert len(records) == 1

    def test_batch_flush_triggers_at_threshold(self, tmp_path, monkeypatch):
        """variant_summary loop flushes a batch every INSERT_BATCH_SIZE rows.
        Mock a low threshold so a small fixture exercises the mid-loop flush
        path (not just the end-of-stream flush)."""
        import allelix.databases.manager as manager_mod

        monkeypatch.setattr(manager_mod, "INSERT_BATCH_SIZE", 2)

        # 3 valid VS rows; with batch size 2 the first 2 flush mid-loop and
        # the third flushes at end-of-stream.
        vs = _write_tsv(
            tmp_path / "vs.tsv",
            _VS_HEADER,
            [
                _vs_row(variation_id="100", rs="100"),
                _vs_row(variation_id="200", rs="200"),
                _vs_row(variation_id="300", rs="300"),
            ],
        )
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(variation_id="100"),
                _ss_row(variation_id="200"),
                _ss_row(variation_id="300"),
            ],
        )
        records = list(iter_clinvar_tsv_records(vs, ss, "GRCh38"))
        assert {r["allele_id"] for r in records} == {100, 200, 300}


# ---------------- load_clinvar_tsv end-to-end tests ----------------


class TestLoadClinvarTsv:
    def test_writes_cache_rows(self, tmp_path):
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(
            tmp_path / "ss.tsv",
            _SS_HEADER,
            [
                _ss_row(scv="SCV000001.1"),
                _ss_row(
                    significance="Likely pathogenic",
                    reported_phenotype="C2:Other condition",
                    scv="SCV000002.1",
                ),
            ],
        )
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_tsv(vs, ss, db, "GRCh38")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            rows = conn.execute(
                "SELECT rsid, clinical_significance, condition FROM clinvar_variants ORDER BY id"
            ).fetchall()
        assert len(rows) == 2
        assert rows[0] == (
            "rs12345",
            "Pathogenic",
            "Hereditary breast and ovarian cancer syndrome",
        )
        assert rows[1] == ("rs12345", "Likely pathogenic", "Other condition")

    def test_stamps_remote_signal_and_interpreter_version(self, tmp_path):
        from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION

        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        db = tmp_path / "clinvar.sqlite"
        load_clinvar_tsv(
            vs,
            ss,
            db,
            "GRCh38",
            source_url="https://test/url",
            remote_signal="md5:abcdef",
        )
        with contextlib.closing(sqlite3.connect(db)) as conn:
            row = conn.execute(
                "SELECT source_url, remote_signal, local_version_tag, "
                "record_count, version FROM database_versions "
                "WHERE name = 'clinvar'"
            ).fetchone()
        assert row[0] == "https://test/url"
        assert row[1] == "md5:abcdef"
        assert row[2] == f"iv:{CLINVAR_INTERPRETER_VERSION}"
        assert row[3] == 1
        # GH #42 follow-up (evaluator defect 3): the TSV loader stamps
        # the cache build date as version so `db status` doesn't show
        # "ClinVar version: None". The build date is YYYY-MM-DD UTC.
        assert row[4] is not None
        assert len(row[4]) == 10  # YYYY-MM-DD
        assert row[4][:4].isdigit()
        assert row[4][4] == "-"

    def test_overwrites_existing_dest(self, tmp_path):
        """Re-running the loader against an existing cache file replaces it.
        The unlink-then-rebuild path keeps the loader idempotent."""
        vs = _write_tsv(tmp_path / "vs.tsv", _VS_HEADER, [_vs_row()])
        ss = _write_tsv(tmp_path / "ss.tsv", _SS_HEADER, [_ss_row()])
        db = tmp_path / "clinvar.sqlite"
        db.write_bytes(b"stale-cache-bytes")
        load_clinvar_tsv(vs, ss, db, "GRCh38")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0]
        assert count == 1


# ---------------- real-data integration test (@slow) ----------------


def _real_data_cache_dir() -> Path:
    """Persistent cache dir for the real TSV files so re-runs don't
    re-download 800 MB. Lives under the user's home, gitignored by
    construction."""
    from pathlib import Path

    p = Path.home() / ".cache" / "allelix-test-data" / "clinvar-tsv"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _ensure_real_tsv(name: str, url: str) -> Path:
    """Download the named TSV to the persistent cache if missing.

    Returns the path. Caller should pytest.skip if the download fails;
    we don't want CI matrix red on transient network issues."""
    from allelix.databases.manager import download

    cache = _real_data_cache_dir()
    dest = cache / name
    if dest.exists() and dest.stat().st_size > 1_000_000:  # sanity bound
        return dest
    try:
        download(url, dest)
    except OSError as exc:
        pytest.skip(f"Could not fetch real {name} ({exc}); skipping integration test")
    return dest


@pytest.mark.slow
@pytest.mark.integration
class TestRealClinvarTsvIngest:
    """End-to-end integration against the actual NCBI-hosted ClinVar TSVs.

    These tests auto-fetch ~800 MB compressed (~5 GB uncompressed) and
    run the full loader. They take 5-15 minutes depending on hardware
    and network. CI skips them via the @slow + @integration markers
    (per CONTRIBUTING.md "Hooks and CI"); they run at ship time.

    The test's job is to verify that the loader, as written, handles
    real-shaped data — not synthetic fixtures — without crashing,
    produces a sensible row count, surfaces per-SCV multi-row variants
    (the whole point of #42), and writes a queryable cache.
    """

    def test_full_ingest_to_cache(self, tmp_path):
        """Build a real GRCh38 cache from the real TSVs and check invariants."""
        from allelix.databases.manager import (
            _CLINVAR_PLACEHOLDER_CLNSIGS,
            CLINVAR_SUBMISSION_SUMMARY_URL,
            CLINVAR_VARIANT_SUMMARY_URL,
        )

        vs_path = _ensure_real_tsv("variant_summary.txt.gz", CLINVAR_VARIANT_SUMMARY_URL)
        ss_path = _ensure_real_tsv("submission_summary.txt.gz", CLINVAR_SUBMISSION_SUMMARY_URL)

        db = tmp_path / "clinvar.GRCh38.sqlite"
        load_clinvar_tsv(vs_path, ss_path, db, "GRCh38", remote_signal="real-data-test")

        with contextlib.closing(sqlite3.connect(db)) as conn:
            row_count = conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0]
            distinct_rsids = conn.execute(
                "SELECT COUNT(DISTINCT rsid) FROM clinvar_variants"
            ).fetchone()[0]
            multi_scv_count = conn.execute(
                "SELECT COUNT(*) FROM ("
                "  SELECT chromosome, position, ref, alt FROM clinvar_variants "
                "  GROUP BY chromosome, position, ref, alt HAVING COUNT(*) > 1"
                ")"
            ).fetchone()[0]
            distinct_sigs = conn.execute(
                "SELECT COUNT(DISTINCT clinical_significance) FROM clinvar_variants"
            ).fetchone()[0]
            # Evaluator defect 5 follow-up: the old spot-check below
            # used `assert sig` (non-empty), which passed on the
            # literal "-" placeholder because "-" is non-empty. Pin
            # the sentinel guard against the whole cache, not a tiny
            # sample, so this test catches any placeholder regression
            # at fast-tier speed against real data — without waiting
            # for the live-cache drift guard to fire.
            #
            # Round-2 evaluator catch: the SQL set MUST match
            # _CLINVAR_PLACEHOLDER_CLNSIGS (the loader's actual skip
            # set). Hand-rolling this set diverged once already — the
            # initial version included "not provided" (which is in
            # _CLNSIG_MAGNITUDE at 2.0, safe, and kept by the loader
            # by design) and omitted "no classification provided"
            # (which the loader DOES filter). The protocol §7b SQL
            # ended up correct; this one didn't — pure copy-paste
            # divergence the cross-PR review lesson should have
            # already retired. Build the SQL straight from the
            # constant so they literally cannot diverge again.
            # LOWER() mirrors the loader's case-insensitive
            # membership test.
            placeholders = sorted(_CLINVAR_PLACEHOLDER_CLNSIGS)
            placeholder_sql = ", ".join("?" * len(placeholders))
            sentinel_count = conn.execute(
                f"SELECT COUNT(*) FROM clinvar_variants "
                f"WHERE LOWER(clinical_significance) IN ({placeholder_sql})",
                placeholders,
            ).fetchone()[0]
            sample = conn.execute(
                "SELECT rsid, clinical_significance, condition FROM clinvar_variants LIMIT 5"
            ).fetchall()

        # Sanity bounds — these aren't trying to assert ClinVar's exact
        # current row count (which drifts weekly with new submissions);
        # they assert "the loader produced realistic data."
        assert row_count > 1_000_000, (
            f"Expected millions of SCV-level rows, got {row_count}. "
            "Either the loader's filtering is too aggressive or the "
            "real ClinVar shrank by an order of magnitude."
        )
        assert distinct_rsids > 100_000, (
            f"Expected 100K+ distinct rsIDs in the cache; got {distinct_rsids}."
        )
        assert multi_scv_count > 1_000, (
            "The whole point of #42 is per-SCV multi-row variants. Got "
            f"{multi_scv_count} variants with >1 SCV row — should be thousands."
        )
        assert distinct_sigs >= 10, (
            f"Expected diverse classification vocabulary; got only "
            f"{distinct_sigs} distinct values. Filtering bug?"
        )
        # Evaluator defect 5: the `-` placeholder (and prose variants)
        # must NEVER land as a literal clinical_significance value.
        # If this ever fires the loader's sentinel filter regressed —
        # CI is now the primary gate, the slow drift-guard against the
        # live cache is backup, the protocol §7 sanity is the final
        # human-eye gate.
        assert sentinel_count == 0, (
            f"Found {sentinel_count} rows where clinical_significance is "
            f"in _CLINVAR_PLACEHOLDER_CLNSIGS ({sorted(_CLINVAR_PLACEHOLDER_CLNSIGS)!r}). "
            f"These would surface in analyze output as bogus annotations "
            f"passing the default magnitude floor. Loader's sentinel filter "
            f"regressed (see iter_clinvar_tsv_records)."
        )
        # Spot-check: every sampled row has the shape we expect.
        # `assert sig` here is a fast-path sanity (every row must
        # carry a value at all) — the placeholder catch is the
        # cache-wide sentinel_count assertion above, which closes
        # the "non-empty != non-placeholder" blind spot the evaluator
        # caught.
        for rsid, sig, condition in sample:
            assert rsid.startswith("rs") and rsid[2:].isdigit(), rsid
            assert sig
            # condition can be empty (some SCVs have no MedGen-mapped phenotype)
            assert isinstance(condition, str)

    # NOTE: test_aggregate_only_drops_supersedeed_submissions used to live
    # here. It ran two more full 50M-variant ingests purely to assert
    # n_default < n_all on real data. The same invariant is already pinned
    # against mock data by TestIterTsvRecords.test_aggregate_only_filters_*
    # (line ~316 and ~332), and the kept test_full_ingest_to_cache still
    # confirms the default mode against real data — so the only thing the
    # deleted test added was real-data confirmation that the filter "fires"
    # (i.e. real ClinVar still contains rows aggregate_only drops). The
    # mock tests prove the logic; R-4 covers vocabulary drift separately.
    # Net trim: slow class 3 passes → 1, ~15 min → ~5 min on slow.yml's
    # weekly cron and on local run-tests.sh. See PR #88 review thread.

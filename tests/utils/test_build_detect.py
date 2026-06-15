# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for genome build detection (ADR-0021)."""

from __future__ import annotations

from allelix.models import Variant
from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    KNOWN_SNP_POSITIONS,
    detect_build,
    normalize_build_label,
)


def _v(rsid: str, chrom: str, position: int) -> Variant:
    return Variant(rsid=rsid, chromosome=chrom, position=position, allele1="A", allele2="A")


class TestDetectBuild:
    def test_grch37_positions_detect_grch37(self):
        variants = [
            _v(rsid, *KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH37])
            for rsid in ("rs1801133", "rs4680", "rs4149056")
        ]
        result = detect_build(variants)
        assert result.build == BUILD_GRCH37
        assert result.matched == 3
        assert result.inspected == 3
        assert result.is_confident

    def test_grch38_positions_detect_grch38(self):
        variants = [
            _v(rsid, *KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH38])
            for rsid in ("rs1801133", "rs4680", "rs4149056")
        ]
        result = detect_build(variants)
        assert result.build == BUILD_GRCH38
        assert result.matched == 3
        assert result.inspected == 3

    def test_single_match_detects_but_not_confident(self):
        """One rsID returns a build but is_confident requires >= 3 matches."""
        variants = [_v("rs1801133", "1", 11796321)]  # GRCh38
        result = detect_build(variants)
        assert result.build == BUILD_GRCH38
        assert result.matched == 1
        assert not result.is_confident

    def test_no_known_snps_returns_none(self):
        variants = [_v("rs999000111", "5", 12345)]
        result = detect_build(variants)
        assert result.build is None
        assert result.matched == 0
        assert result.inspected == 0
        assert not result.is_confident

    def test_known_rsid_at_wrong_position_counted_as_inspected(self):
        """A known rsID at a non-canonical position is inspected but not matched."""
        variants = [_v("rs1801133", "1", 99999999)]
        result = detect_build(variants)
        assert result.build is None
        assert result.inspected == 1
        assert result.matched == 0

    def test_mixed_grch37_grch38_inconsistent_returns_none(self):
        """If half the rsIDs say GRCh37 and half say GRCh38, the file is
        internally inconsistent and detection refuses to pick.
        """
        v_g37 = _v("rs1801133", *KNOWN_SNP_POSITIONS["rs1801133"][BUILD_GRCH37])
        v_g38 = _v("rs4680", *KNOWN_SNP_POSITIONS["rs4680"][BUILD_GRCH38])
        result = detect_build([v_g37, v_g38])
        assert result.build is None
        assert result.inspected == 2

    def test_streaming_completes_early_when_table_exhausted(self):
        """detect_build doesn't drain the iterator — it bails once every
        entry in the table has been seen.
        """
        seen = []

        def gen():
            for rsid in KNOWN_SNP_POSITIONS:
                chrom, pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH38]
                seen.append(rsid)
                yield _v(rsid, chrom, pos)
            # Sentinel that should NOT be yielded if early-exit works.
            seen.append("SHOULD_NOT_BE_YIELDED")
            yield _v("rs999", "1", 1)

        result = detect_build(gen())
        assert result.build == BUILD_GRCH38
        assert "SHOULD_NOT_BE_YIELDED" not in seen

    def test_unknown_rsids_ignored(self):
        """Non-table rsIDs don't interfere with detection."""
        variants = [
            _v("rs999000001", "5", 12345),
            _v("rs999000002", "5", 67890),
            _v("rs1801133", *KNOWN_SNP_POSITIONS["rs1801133"][BUILD_GRCH37]),
            _v("rs4680", *KNOWN_SNP_POSITIONS["rs4680"][BUILD_GRCH37]),
            _v("rs4149056", *KNOWN_SNP_POSITIONS["rs4149056"][BUILD_GRCH37]),
        ]
        result = detect_build(variants)
        assert result.build == BUILD_GRCH37
        assert result.is_confident

    def test_grch36_positions_detect_grch36(self):
        variants = [
            _v(rsid, *KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH36])
            for rsid in ("rs1801133", "rs4680", "rs4149056")
        ]
        result = detect_build(variants)
        assert result.build == BUILD_GRCH36
        assert result.matched == 3
        assert result.inspected == 3
        assert result.is_confident

    def test_grch36_single_match_not_confident(self):
        variants = [_v("rs4680", "22", 18331271)]  # GRCh36
        result = detect_build(variants)
        assert result.build == BUILD_GRCH36
        assert result.matched == 1
        assert not result.is_confident

    def test_grch36_majority_wins_over_minority(self):
        """Three GRCh36 + one GRCh37 → GRCh36 wins by majority."""
        variants = [
            _v("rs1801133", *KNOWN_SNP_POSITIONS["rs1801133"][BUILD_GRCH36]),
            _v("rs4680", *KNOWN_SNP_POSITIONS["rs4680"][BUILD_GRCH36]),
            _v("rs4149056", *KNOWN_SNP_POSITIONS["rs4149056"][BUILD_GRCH36]),
            _v("rs1800497", *KNOWN_SNP_POSITIONS["rs1800497"][BUILD_GRCH37]),
        ]
        result = detect_build(variants)
        assert result.build == BUILD_GRCH36
        assert result.matched == 3
        assert result.inspected == 4


class TestNormalizeBuildLabel:
    def test_grch37_canonical(self):
        assert normalize_build_label("GRCh37") == BUILD_GRCH37
        assert normalize_build_label("grch37") == BUILD_GRCH37

    def test_grch38_canonical(self):
        assert normalize_build_label("GRCh38") == BUILD_GRCH38
        assert normalize_build_label("grch38") == BUILD_GRCH38

    def test_grch36_canonical(self):
        assert normalize_build_label("GRCh36") == BUILD_GRCH36
        assert normalize_build_label("grch36") == BUILD_GRCH36

    def test_hg_aliases(self):
        assert normalize_build_label("hg18") == BUILD_GRCH36
        assert normalize_build_label("HG18") == BUILD_GRCH36
        assert normalize_build_label("hg19") == BUILD_GRCH37
        assert normalize_build_label("HG19") == BUILD_GRCH37
        assert normalize_build_label("hg38") == BUILD_GRCH38

    def test_numeric_build(self):
        assert normalize_build_label("36") == BUILD_GRCH36
        assert normalize_build_label("build 36") == BUILD_GRCH36
        assert normalize_build_label("37") == BUILD_GRCH37
        assert normalize_build_label("38") == BUILD_GRCH38
        assert normalize_build_label("Build 37.1") == BUILD_GRCH37
        assert normalize_build_label("build 38") == BUILD_GRCH38

    def test_unknown_returns_none(self):
        assert normalize_build_label("CHM13") is None
        assert normalize_build_label("") is None
        assert normalize_build_label(None) is None
        assert normalize_build_label("   ") is None


class TestNormalizeBuildLabelTokenBoundaries:
    """GH #16: anchored tokenization, ambiguity rejection.

    The previous implementation was first-substring-wins over a bare
    `"3X" in s` chain. A date like `2038-01-01` would match `38` and
    silently retag a file as GRCh38; `hg19/hg38` would match `hg19`
    first and pick GRCh37 even though hg38 was also present.
    """

    def test_date_does_not_match(self) -> None:
        # `2038` has no word boundary before `38` (digit-to-digit), so
        # the standalone-number token cannot fire from inside a year.
        assert normalize_build_label("2038-01-01") is None
        assert normalize_build_label("# downloaded 2037-12-31") is None
        assert normalize_build_label("v37 of platform") is None

    def test_ambiguous_returns_none(self) -> None:
        # Mentioning multiple distinct builds → no safe answer; the
        # previous first-substring-wins picked one arbitrarily.
        assert normalize_build_label("hg19/hg38") is None
        assert normalize_build_label("38 (37 liftover)") is None
        assert normalize_build_label("GRCh37 lifted to GRCh38") is None

    def test_unambiguous_repeated_token_still_matches(self) -> None:
        # Repeating the SAME build is unambiguous.
        assert normalize_build_label("GRCh37 (build 37, hg19)") == BUILD_GRCH37

    def test_grch_with_patch_suffix(self) -> None:
        # `GRCh37.p13` and `Build 37.1` — boundaries on `.` are fine.
        assert normalize_build_label("GRCh37.p13") == BUILD_GRCH37
        assert normalize_build_label("Build 37.1") == BUILD_GRCH37

    def test_ncbi_variants(self) -> None:
        # NCBI 36 / NCBI36 / NCBI 037 all canonicalize.
        assert normalize_build_label("NCBI 36") == BUILD_GRCH36
        assert normalize_build_label("NCBI36") == BUILD_GRCH36
        assert normalize_build_label("NCBI 037") == BUILD_GRCH37


class TestKnownSnpTable:
    """ADR-0021: the canonical SNP position table is authoritative.

    These tests pin invariants the table MUST satisfy. If any fails,
    the table is wrong — fix it against NCBI dbSNP before relaxing the
    assertion.
    """

    def test_every_entry_has_all_three_builds(self):
        for rsid, entry in KNOWN_SNP_POSITIONS.items():
            assert BUILD_GRCH36 in entry, f"{rsid} missing GRCh36 position"
            assert BUILD_GRCH37 in entry, f"{rsid} missing GRCh37 position"
            assert BUILD_GRCH38 in entry, f"{rsid} missing GRCh38 position"

    def test_chromosomes_agree_across_builds(self):
        """A SNP's chromosome doesn't change between builds."""
        for rsid, entry in KNOWN_SNP_POSITIONS.items():
            chroms = {entry[b][0] for b in (BUILD_GRCH36, BUILD_GRCH37, BUILD_GRCH38)}
            assert len(chroms) == 1, f"{rsid}: chr differs across builds: {chroms}"

    def test_positions_differ_between_builds(self):
        """Positions must be unique across all three builds for each SNP."""
        for rsid, entry in KNOWN_SNP_POSITIONS.items():
            positions = [entry[b][1] for b in (BUILD_GRCH36, BUILD_GRCH37, BUILD_GRCH38)]
            assert len(set(positions)) == 3, (
                f"{rsid}: duplicate position across builds — useless for detection"
            )

    def test_positions_are_positive_integers(self):
        for rsid, entry in KNOWN_SNP_POSITIONS.items():
            for build, (_, pos) in entry.items():
                assert isinstance(pos, int) and pos > 0, f"{rsid}/{build}: bad position {pos!r}"

    def test_table_covers_multiple_chromosomes(self):
        """A single-chromosome table would fail to detect builds on files
        that happen to lack that chromosome's probes.
        """
        chroms = {entry[BUILD_GRCH37][0] for entry in KNOWN_SNP_POSITIONS.values()}
        assert len(chroms) >= 5, f"Table covers only {chroms} — need broader spread"

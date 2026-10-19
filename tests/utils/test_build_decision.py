# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The build decision rule (ADR-0044): evidence before labels.

Unit tests of ``contig_length_verdict``, ``decide_build``,
``check_late_seeds`` and ``needs_late_seeds``. Every input is synthetic.
"""

from __future__ import annotations

import pytest

from allelix.utils.build_detect import (
    _CONTIG_LENGTHS,
    BUILD_GRCH36,
    BUILD_GRCH37,
    BUILD_GRCH38,
    BUILD_SOURCE_LABELS,
    BUILD_SOURCES,
    HG19_CHRM_LENGTH,
    KIND_ARRAY,
    KIND_VCF,
    NO_SEEDS,
    SOURCE_ASSEMBLY_HEADER,
    SOURCE_CONTIG_LENGTHS,
    SOURCE_DEFAULT,
    SOURCE_FILE_LABEL,
    SOURCE_OVERRIDE,
    SOURCE_SEEDS,
    BuildDecisionError,
    BuildDetectionResult,
    LengthVerdict,
    check_late_seeds,
    contig_length_verdict,
    decide_build,
    needs_late_seeds,
)

# T2T-CHM13v2.0 chr1, an assembly Allelix does not support.
_T2T_CHR1 = 248_387_328


def _lengths(build: str, *names: str) -> dict[str, int]:
    return {name: _CONTIG_LENGTHS[build][name.removeprefix("chr")] for name in names}


def _seeds(build: str | None, matched: int, inspected: int) -> BuildDetectionResult:
    return BuildDetectionResult(build=build, matched=matched, inspected=inspected)


def _decide(
    *,
    kind: str = KIND_ARRAY,
    override: str | None = None,
    file_label: str | None = None,
    lengths: LengthVerdict | None = None,
    seeds: BuildDetectionResult = NO_SEEDS,
):
    return decide_build(
        file_name="sample.txt",
        kind=kind,
        override=override,
        file_label=file_label,
        default=BUILD_GRCH37 if kind == KIND_ARRAY else None,
        lengths=lengths,
        seeds=seeds,
    )


CONFIDENT_37 = _seeds(BUILD_GRCH37, 10, 10)
CONFIDENT_38 = _seeds(BUILD_GRCH38, 9, 10)
TENTATIVE_37 = _seeds(BUILD_GRCH37, 3, 4)
TENTATIVE_38 = _seeds(BUILD_GRCH38, 1, 1)
TENTATIVE_36 = _seeds(BUILD_GRCH36, 3, 4)
ZERO_VOTES = _seeds(None, 0, 2)  # two seeds found, neither at any build's coordinate
TIE_37_38 = BuildDetectionResult(None, 0, 4, (BUILD_GRCH37, BUILD_GRCH38), 2)
TIE_36_38 = BuildDetectionResult(None, 0, 4, (BUILD_GRCH36, BUILD_GRCH38), 2)
TIE_36_37 = BuildDetectionResult(None, 0, 4, (BUILD_GRCH36, BUILD_GRCH37), 2)


class TestContigLengthVerdict:
    @pytest.mark.parametrize("build", [BUILD_GRCH37, BUILD_GRCH38])
    def test_any_nuclear_contig_decides(self, build: str) -> None:
        for name in ("1", "chr7", "X", "chrY"):
            verdict = contig_length_verdict(_lengths(build, name))
            assert (verdict.status, verdict.build) == ("match", build), name

    def test_every_nuclear_contig_differs_between_builds(self) -> None:
        for name, length in _CONTIG_LENGTHS[BUILD_GRCH37].items():
            assert _CONTIG_LENGTHS[BUILD_GRCH38][name] != length, name

    def test_naming_is_not_evidence(self) -> None:
        """hg19 chr-prefixed lengths are GRCh37; bare GRCh38 lengths are GRCh38."""
        assert contig_length_verdict(_lengths(BUILD_GRCH37, "chr1", "chr22")).build == (
            BUILD_GRCH37
        )
        assert contig_length_verdict(_lengths(BUILD_GRCH38, "1", "22")).build == BUILD_GRCH38

    def test_mt_alone_is_not_evidence(self) -> None:
        assert contig_length_verdict({"MT": 16_569}).status == "none"
        assert contig_length_verdict({"chrM": 16_569}).status == "none"

    def test_alt_decoy_and_random_contigs_are_ignored(self) -> None:
        verdict = contig_length_verdict(
            {"chr1_KI270706v1_random": 175_055, "hs37d5": 35_477_943, "GL000192.1": 547_496}
        )
        assert verdict.status == "none"

    def test_mixed_lengths(self) -> None:
        lengths = {**_lengths(BUILD_GRCH37, "1"), **_lengths(BUILD_GRCH38, "22")}
        verdict = contig_length_verdict(lengths)
        assert verdict.status == "mixed"
        assert verdict.build is None

    def test_unsupported_assembly(self) -> None:
        verdict = contig_length_verdict({"chr1": _T2T_CHR1})
        assert verdict.status == "unsupported"
        assert "248,387,328" in verdict.detail

    def test_hg19_chrm_length_marks_mt_for_skipping(self) -> None:
        lengths = {**_lengths(BUILD_GRCH37, "chr1"), "chrM": HG19_CHRM_LENGTH}
        verdict = contig_length_verdict(lengths)
        assert (verdict.build, verdict.skip_mt) == (BUILD_GRCH37, True)
        assert not contig_length_verdict({**lengths, "chrM": 16_569}).skip_mt


class TestDecideBuild:
    def test_vocabulary(self) -> None:
        assert BUILD_SOURCES == (
            SOURCE_OVERRIDE,
            SOURCE_CONTIG_LENGTHS,
            SOURCE_SEEDS,
            SOURCE_ASSEMBLY_HEADER,
            SOURCE_FILE_LABEL,
            SOURCE_DEFAULT,
        )
        assert set(BUILD_SOURCE_LABELS) == set(BUILD_SOURCES)

    def test_override_beats_everything_and_records_the_label(self) -> None:
        mixed = contig_length_verdict(
            {**_lengths(BUILD_GRCH37, "1"), **_lengths(BUILD_GRCH38, "22")}
        )
        decision = _decide(
            kind=KIND_VCF, override=BUILD_GRCH38, file_label=BUILD_GRCH37, lengths=mixed
        )
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_OVERRIDE)
        assert decision.file_label == BUILD_GRCH37

    # ── VCF contig lengths ──

    def test_lengths_decide_without_a_label(self) -> None:
        lengths = contig_length_verdict(_lengths(BUILD_GRCH37, "chr1"))
        decision = _decide(kind=KIND_VCF, lengths=lengths)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_CONTIG_LENGTHS)

    def test_lengths_agreeing_with_the_label(self) -> None:
        lengths = contig_length_verdict(_lengths(BUILD_GRCH38, "1"))
        decision = _decide(kind=KIND_VCF, file_label=BUILD_GRCH38, lengths=lengths)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_CONTIG_LENGTHS)

    def test_lengths_contradicting_the_assembly_stop(self) -> None:
        lengths = contig_length_verdict(_lengths(BUILD_GRCH37, "1"))
        with pytest.raises(BuildDecisionError, match="header says GRCh38 but its contig"):
            _decide(kind=KIND_VCF, file_label=BUILD_GRCH38, lengths=lengths)

    def test_mixed_lengths_stop(self) -> None:
        lengths = contig_length_verdict(
            {**_lengths(BUILD_GRCH37, "1"), **_lengths(BUILD_GRCH38, "22")}
        )
        with pytest.raises(BuildDecisionError, match="contradict each other"):
            _decide(kind=KIND_VCF, lengths=lengths)

    def test_unsupported_lengths_stop_even_with_a_label(self) -> None:
        lengths = contig_length_verdict({"chr1": _T2T_CHR1})
        with pytest.raises(BuildDecisionError, match="unsupported assembly"):
            _decide(kind=KIND_VCF, file_label=BUILD_GRCH38, lengths=lengths, seeds=CONFIDENT_38)

    def test_hg19_chrm_passes_skip_mt_on(self) -> None:
        lengths = contig_length_verdict(
            {**_lengths(BUILD_GRCH37, "chr1"), "chrM": HG19_CHRM_LENGTH}
        )
        assert _decide(kind=KIND_VCF, lengths=lengths).skip_mt

    # chrM alone declares a length: the lengths do not decide the build,
    # but the hg19 chrM still holds whatever does.
    _HG19_CHRM_ONLY = contig_length_verdict({"chrM": HG19_CHRM_LENGTH})

    def test_hg19_chrm_skipped_under_build_override(self) -> None:
        decision = _decide(kind=KIND_VCF, override=BUILD_GRCH37, lengths=self._HG19_CHRM_ONLY)
        assert decision.source == SOURCE_OVERRIDE
        assert decision.skip_mt

    def test_hg19_chrm_skipped_under_assembly_label(self) -> None:
        decision = _decide(
            kind=KIND_VCF,
            file_label=BUILD_GRCH37,
            lengths=self._HG19_CHRM_ONLY,
            seeds=TENTATIVE_37,
        )
        assert decision.source == SOURCE_ASSEMBLY_HEADER
        assert decision.skip_mt

    def test_hg19_chrm_skipped_when_seeds_decide(self) -> None:
        decision = _decide(kind=KIND_VCF, lengths=self._HG19_CHRM_ONLY, seeds=CONFIDENT_37)
        assert decision.source == SOURCE_SEEDS
        assert decision.skip_mt

    @pytest.mark.parametrize(
        "how",
        [
            {"override": BUILD_GRCH37},
            {"file_label": BUILD_GRCH37, "seeds": TENTATIVE_37},
            {"seeds": CONFIDENT_37},
        ],
    )
    def test_rcrs_chrm_is_kept_on_every_path(self, how: dict[str, object]) -> None:
        lengths = contig_length_verdict({"chrM": 16_569})
        assert not _decide(kind=KIND_VCF, lengths=lengths, **how).skip_mt  # type: ignore[arg-type]

    def test_arrays_never_skip_mt(self) -> None:
        assert not _decide(kind=KIND_ARRAY, override=BUILD_GRCH37).skip_mt

    # ── known-SNP positions ──

    def test_confident_seeds_override_a_wrong_label(self) -> None:
        decision = _decide(file_label=BUILD_GRCH37, seeds=CONFIDENT_38)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_SEEDS)
        assert decision.label_overridden

    def test_confident_seeds_agreeing_with_the_label(self) -> None:
        decision = _decide(file_label=BUILD_GRCH37, seeds=CONFIDENT_37)
        assert decision.source == SOURCE_SEEDS
        assert not decision.label_overridden

    def test_tentative_contradiction_of_a_label_stops(self) -> None:
        with pytest.raises(BuildDecisionError) as exc:
            _decide(file_label=BUILD_GRCH37, seeds=_seeds(BUILD_GRCH38, 3, 4))
        message = str(exc.value)
        assert "sample.txt" in message
        assert "point to GRCh38 but the file says GRCh37" in message
        assert "3 of 4 known-SNP positions match GRCh38 (tentative)" in message
        assert "--build grch37 or --build grch38" in message

    def test_tentative_agreement_keeps_the_label(self) -> None:
        decision = _decide(file_label=BUILD_GRCH37, seeds=TENTATIVE_37)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_FILE_LABEL)
        decision = _decide(kind=KIND_VCF, file_label=BUILD_GRCH37, seeds=TENTATIVE_37)
        assert decision.source == SOURCE_ASSEMBLY_HEADER

    # A tie (non-zero votes, no single winner) lets the run proceed
    # on the build that would otherwise decide only when that build is one of
    # the tied builds; otherwise it stops.

    def test_tie_including_the_label_keeps_it(self) -> None:
        for kind, source in ((KIND_ARRAY, SOURCE_FILE_LABEL), (KIND_VCF, SOURCE_ASSEMBLY_HEADER)):
            decision = _decide(kind=kind, file_label=BUILD_GRCH37, seeds=TIE_37_38)
            assert (decision.build, decision.source) == (BUILD_GRCH37, source), kind
            assert decision.seeds.tied == (BUILD_GRCH37, BUILD_GRCH38)

    def test_tie_excluding_the_label_stops(self) -> None:
        """Without GRCh36 among the tied builds, a tie that excludes the
        label stops."""
        for kind in (KIND_ARRAY, KIND_VCF):
            with pytest.raises(BuildDecisionError) as exc:
                _decide(kind=kind, file_label=BUILD_GRCH36, seeds=TIE_37_38)
            message = str(exc.value)
            assert "tied between GRCh37 and GRCh38 but the file says GRCh36" in message
            assert (
                "4 known-SNP positions found, tied between GRCh37 and GRCh38 (2 each)" in message
            )
            assert "--build grch37 or --build grch38" in message

    @pytest.mark.parametrize("tie", ["36/37", "36/38"])
    @pytest.mark.parametrize("label", [None, BUILD_GRCH36, BUILD_GRCH37, BUILD_GRCH38])
    @pytest.mark.parametrize("kind", [KIND_ARRAY, KIND_VCF])
    def test_tie_including_grch36_takes_the_fail_safe(
        self, kind: str, label: str | None, tie: str
    ) -> None:
        """Without a contig-length decision, a tie that includes GRCh36
        means GRCh36, whether or not the label is among the tied builds."""
        seeds = TIE_36_37 if tie == "36/37" else TIE_36_38
        lengths = contig_length_verdict({}) if kind == KIND_VCF else None
        decision = _decide(kind=kind, file_label=label, lengths=lengths, seeds=seeds)
        assert (decision.build, decision.source) == (BUILD_GRCH36, SOURCE_SEEDS)
        assert decision.tentative_grch36
        assert decision.label_overridden is (label not in (None, BUILD_GRCH36))

    def test_override_beats_a_grch36_tie(self) -> None:
        decision = _decide(override=BUILD_GRCH38, file_label=BUILD_GRCH37, seeds=TIE_36_38)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_OVERRIDE)

    def test_unlabelled_chip_tie_including_the_default_uses_it(self) -> None:
        decision = _decide(seeds=TIE_37_38)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_DEFAULT)

    def test_unlabelled_chip_tie_excluding_the_default_takes_the_fail_safe(self) -> None:
        """With GRCh37 the default, a tie that excludes it includes GRCh36: fail-safe."""
        decision = _decide(seeds=TIE_36_38)
        assert (decision.build, decision.tentative_grch36) == (BUILD_GRCH36, True)

    def test_unlabelled_vcf_tie_stops(self) -> None:
        """No lengths and no label: nothing would otherwise decide."""
        with pytest.raises(BuildDecisionError, match="no contig lengths, no assembly"):
            _decide(kind=KIND_VCF, lengths=contig_length_verdict({}), seeds=TIE_37_38)

    def test_prescan_tie_excluding_the_contig_length_build_stops(self) -> None:
        lengths = contig_length_verdict(_lengths(BUILD_GRCH37, "1"))
        with pytest.raises(BuildDecisionError, match="tied between GRCh36 and GRCh38"):
            _decide(kind=KIND_VCF, lengths=lengths, seeds=TIE_36_38)
        decision = _decide(kind=KIND_VCF, lengths=lengths, seeds=TIE_37_38)
        assert decision.source == SOURCE_CONTIG_LENGTHS

    def test_contig_lengths_decide_a_grch36_tie(self) -> None:
        """The GRCh36 tie fail-safe does not apply when contig lengths decided."""
        lengths = contig_length_verdict(_lengths(BUILD_GRCH37, "1"))
        decision = _decide(kind=KIND_VCF, lengths=lengths, seeds=TIE_36_37)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_CONTIG_LENGTHS)
        assert not decision.tentative_grch36

    def test_confident_grch36_decides(self) -> None:
        decision = _decide(file_label=BUILD_GRCH37, seeds=_seeds(BUILD_GRCH36, 10, 10))
        assert (decision.build, decision.source) == (BUILD_GRCH36, SOURCE_SEEDS)
        assert decision.label_overridden
        assert not decision.tentative_grch36

    @pytest.mark.parametrize("kind", [KIND_ARRAY, KIND_VCF])
    def test_label_without_seeds_stops(self, kind: str) -> None:
        """a label with nothing to check it against does not decide alone."""
        lengths = contig_length_verdict({}) if kind == KIND_VCF else None
        with pytest.raises(BuildDecisionError) as exc:
            _decide(kind=kind, file_label=BUILD_GRCH38, lengths=lengths)
        message = str(exc.value)
        assert "the file says GRCh38 but carries no known-SNP position to check it" in message
        assert "--build grch37 or --build grch38" in message

    @pytest.mark.parametrize("kind", [KIND_ARRAY, KIND_VCF])
    def test_label_with_one_agreeing_seed_decides(self, kind: str) -> None:
        lengths = contig_length_verdict({}) if kind == KIND_VCF else None
        decision = _decide(kind=kind, file_label=BUILD_GRCH38, lengths=lengths, seeds=TENTATIVE_38)
        source = SOURCE_ASSEMBLY_HEADER if kind == KIND_VCF else SOURCE_FILE_LABEL
        assert (decision.build, decision.source) == (BUILD_GRCH38, source)

    def test_grch36_label_without_seeds_runs_as_grch36(self) -> None:
        """A vendor GRCh36 label is exempt from the stop, like --build grch36:
        GRCh36 turns position lookups off (ADR-0025)."""
        decision = _decide(file_label=BUILD_GRCH36)
        assert (decision.build, decision.source) == (BUILD_GRCH36, SOURCE_FILE_LABEL)

    def test_vcf_grch36_assembly_without_seeds_still_stops(self) -> None:
        with pytest.raises(BuildDecisionError, match="carries no known-SNP position"):
            _decide(kind=KIND_VCF, file_label=BUILD_GRCH36, lengths=contig_length_verdict({}))

    def test_label_without_seeds_passes_with_build(self) -> None:
        decision = _decide(file_label=BUILD_GRCH38, override=BUILD_GRCH38)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_OVERRIDE)

    @pytest.mark.parametrize("label", [None, BUILD_GRCH37])
    def test_tentative_grch36_fails_safe_to_grch36(self, label: str | None) -> None:
        """ADR-0025: GRCh36 turns position lookups off, so it is the safe answer."""
        decision = _decide(file_label=label, seeds=TENTATIVE_36)
        assert decision.build == BUILD_GRCH36
        # A single GRCh36 probe is tentative too, and still means GRCh36.
        assert _decide(file_label=label, seeds=_seeds(BUILD_GRCH36, 1, 1)).build == BUILD_GRCH36
        assert decision.build == BUILD_GRCH36
        assert decision.tentative_grch36
        assert decision.label_overridden is (label is not None)

    # ── no label ──

    def test_unlabelled_chip_without_evidence_stops(self) -> None:
        """No label and no seed inspected at all stops."""
        with pytest.raises(BuildDecisionError) as exc:
            _decide()
        message = str(exc.value)
        assert "sample.txt" in message
        assert "states no genome build and carries no known-SNP positions" in message
        assert "Evidence: no known-SNP positions in the file; the file states no build" in message
        assert "--build grch37 or --build grch38" in message

    def test_unlabelled_chip_without_evidence_with_override(self) -> None:
        decision = _decide(override=BUILD_GRCH38)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_OVERRIDE)

    # Seeds inspected but placed on no build (zero votes) stop the
    # run unless contig lengths decided.

    @pytest.mark.parametrize("label", [None, BUILD_GRCH37, BUILD_GRCH38])
    @pytest.mark.parametrize("kind", [KIND_ARRAY, KIND_VCF])
    def test_zero_vote_seeds_stop(self, kind: str, label: str | None) -> None:
        lengths = contig_length_verdict({}) if kind == KIND_VCF else None
        with pytest.raises(BuildDecisionError) as exc:
            _decide(kind=kind, file_label=label, lengths=lengths, seeds=ZERO_VOTES)
        message = str(exc.value)
        assert "its known-SNP positions match none of the supported builds" in message
        assert "2 known-SNP position(s) found, none at any supported build's coordinate" in message
        assert "--build grch37 or --build grch38" in message

    def test_zero_vote_seeds_with_override(self) -> None:
        decision = _decide(file_label=BUILD_GRCH37, override=BUILD_GRCH37, seeds=ZERO_VOTES)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_OVERRIDE)

    def test_zero_vote_seeds_do_not_stop_a_contig_length_decision(self) -> None:
        lengths = contig_length_verdict(_lengths(BUILD_GRCH37, "1"))
        decision = _decide(kind=KIND_VCF, lengths=lengths, seeds=ZERO_VOTES)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_CONTIG_LENGTHS)
        check_late_seeds(decision, ZERO_VOTES, "x.vcf")

    def test_unlabelled_chip_with_tentative_default_evidence_uses_the_default(self) -> None:
        assert _decide(seeds=TENTATIVE_37).source == SOURCE_DEFAULT

    @pytest.mark.parametrize(
        ("votes", "matched", "inspected"),
        [
            ("2 GRCh37, 1 GRCh38", 2, 3),
            ("3 GRCh37, 1 GRCh38", 3, 4),
            ("1 GRCh37, 2 at no build's coordinate", 1, 3),
        ],
    )
    def test_default_follows_the_tentative_winner(
        self, votes: str, matched: int, inspected: int
    ) -> None:
        """The default applies when the seeds' tentative winner is the default; a
        minority for another build, or at no build's coordinate, does not stop it."""
        decision = _decide(seeds=_seeds(BUILD_GRCH37, matched, inspected))
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_DEFAULT), votes

    def test_unlabelled_chip_tentatively_another_build_stops(self) -> None:
        with pytest.raises(BuildDecisionError, match="not the GRCh37 default"):
            _decide(seeds=TENTATIVE_38)

    def test_vcf_without_lengths_assembly_or_confident_seeds_stops(self) -> None:
        none = contig_length_verdict({})
        for seeds in (NO_SEEDS, TENTATIVE_38, TENTATIVE_37):
            with pytest.raises(BuildDecisionError, match="no contig lengths, no assembly"):
                _decide(kind=KIND_VCF, lengths=none, seeds=seeds)

    def test_vcf_with_confident_seeds_and_nothing_else(self) -> None:
        decision = _decide(kind=KIND_VCF, lengths=contig_length_verdict({}), seeds=CONFIDENT_38)
        assert (decision.build, decision.source) == (BUILD_GRCH38, SOURCE_SEEDS)


class TestLateSeedCheck:
    def _by_lengths(self, build: str):
        return _decide(kind=KIND_VCF, lengths=contig_length_verdict(_lengths(build, "1")))

    def test_agreeing_seeds_pass(self) -> None:
        check_late_seeds(self._by_lengths(BUILD_GRCH37), CONFIDENT_37, "x.vcf")

    def test_late_tie_including_the_lengths_build_passes(self) -> None:
        check_late_seeds(self._by_lengths(BUILD_GRCH37), TIE_37_38, "x.vcf")

    def test_late_grch36_tie_including_the_lengths_build_passes(self) -> None:
        """After a contig-length decision a GRCh36 tie that includes the lengths build passes."""
        check_late_seeds(self._by_lengths(BUILD_GRCH37), TIE_36_37, "x.vcf")

    def test_late_tie_excluding_the_lengths_build_stops(self) -> None:
        with pytest.raises(BuildDecisionError) as exc:
            check_late_seeds(self._by_lengths(BUILD_GRCH37), TIE_36_38, "x.vcf")
        message = str(exc.value)
        assert "contig lengths say GRCh37 but its known-SNP positions are tied" in message
        assert "annotation pass" in message
        assert "nothing was written" in message

    def test_no_seeds_pass(self) -> None:
        check_late_seeds(self._by_lengths(BUILD_GRCH37), NO_SEEDS, "x.vcf")

    @pytest.mark.parametrize("seeds", [CONFIDENT_38, TENTATIVE_38])
    def test_any_majority_against_the_lengths_stops(self, seeds) -> None:
        with pytest.raises(BuildDecisionError) as exc:
            check_late_seeds(self._by_lengths(BUILD_GRCH37), seeds, "x.vcf")
        message = str(exc.value)
        assert "contig lengths say GRCh37" in message
        assert "annotation pass" in message
        assert "nothing was written" in message

    def test_only_a_lengths_or_override_decision_is_checked(self) -> None:
        check_late_seeds(
            _decide(file_label=BUILD_GRCH37, seeds=TENTATIVE_37), CONFIDENT_38, "x.txt"
        )
        check_late_seeds(_decide(seeds=CONFIDENT_37), CONFIDENT_38, "x.txt")


class TestBuildOverrideChecked:
    """--build against the file's own evidence."""

    _LENGTHS_37 = contig_length_verdict(_lengths(BUILD_GRCH37, "1"))

    @pytest.mark.parametrize("kind", [KIND_ARRAY, KIND_VCF])
    def test_confident_contrary_seeds_stop_naming_both(self, kind: str) -> None:
        with pytest.raises(BuildDecisionError) as exc:
            _decide(kind=kind, override=BUILD_GRCH37, seeds=CONFIDENT_38)
        message = str(exc.value)
        assert "--build grch37 was given but the file's known-SNP positions say GRCh38" in message
        assert "9 of 10 known-SNP positions match GRCh38 (confident)" in message

    def test_contrary_contig_lengths_stop_naming_both(self) -> None:
        with pytest.raises(BuildDecisionError) as exc:
            _decide(kind=KIND_VCF, override=BUILD_GRCH38, lengths=self._LENGTHS_37)
        assert "--build grch38 was given but the file's contig lengths are GRCh37" in str(
            exc.value
        )

    @pytest.mark.parametrize(
        "seeds",
        # 8 of 11 is many contrary positions, but under 80%: not confident.
        [TENTATIVE_38, _seeds(BUILD_GRCH38, 8, 11), TIE_37_38, TIE_36_38, ZERO_VOTES, NO_SEEDS],
    )
    def test_tentative_or_no_evidence_runs_as_forced(self, seeds: BuildDetectionResult) -> None:
        decision = _decide(override=BUILD_GRCH37, seeds=seeds)
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_OVERRIDE)
        assert decision.seeds == seeds  # kept for the warning notice

    def test_agreeing_evidence_runs(self) -> None:
        decision = _decide(
            kind=KIND_VCF, override=BUILD_GRCH37, lengths=self._LENGTHS_37, seeds=CONFIDENT_37
        )
        assert (decision.build, decision.source) == (BUILD_GRCH37, SOURCE_OVERRIDE)

    def test_grch36_is_exempt(self) -> None:
        """ADR-0025: GRCh36 turns position lookups off, so it cannot attach wrong records."""
        decision = _decide(override=BUILD_GRCH36, seeds=CONFIDENT_38)
        assert decision.build == BUILD_GRCH36
        decision = _decide(kind=KIND_VCF, override=BUILD_GRCH36, lengths=self._LENGTHS_37)
        assert decision.build == BUILD_GRCH36
        assert not needs_late_seeds(decision)
        check_late_seeds(decision, CONFIDENT_38, "x.vcf")

    def test_needs_late_seeds(self) -> None:
        by_lengths = _decide(kind=KIND_VCF, lengths=self._LENGTHS_37)
        forced_with_lengths = _decide(
            kind=KIND_VCF, override=BUILD_GRCH37, lengths=self._LENGTHS_37
        )
        forced_without = _decide(
            kind=KIND_VCF, override=BUILD_GRCH37, lengths=contig_length_verdict({})
        )
        forced_array = _decide(override=BUILD_GRCH37)
        assert needs_late_seeds(by_lengths)
        assert needs_late_seeds(forced_with_lengths)
        assert not needs_late_seeds(forced_without)
        assert not needs_late_seeds(forced_array)

    def test_late_confident_contrary_seeds_stop(self) -> None:
        decision = _decide(kind=KIND_VCF, override=BUILD_GRCH37, lengths=self._LENGTHS_37)
        with pytest.raises(BuildDecisionError) as exc:
            check_late_seeds(decision, CONFIDENT_38, "x.vcf")
        message = str(exc.value)
        assert "--build grch37 was given but the file's known-SNP positions say GRCh38" in message
        assert "nothing was written" in message

    @pytest.mark.parametrize("seeds", [TENTATIVE_38, TIE_36_38, CONFIDENT_37])
    def test_late_tentative_or_agreeing_seeds_pass(self, seeds: BuildDetectionResult) -> None:
        decision = _decide(kind=KIND_VCF, override=BUILD_GRCH37, lengths=self._LENGTHS_37)
        check_late_seeds(decision, seeds, "x.vcf")

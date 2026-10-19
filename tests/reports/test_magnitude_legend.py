# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The HTML magnitude legend states the scores the code assigns.

Renders a report, takes its "Understanding Magnitude Scores" legend, and
checks every ClinVar, ClinPGx and GWAS value in it against the scoring
tables and functions: a score that changes in the code without the legend
following fails here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from allelix.annotators.clinvar import (
    _AGGREGATE_CAP,
    _BENIGN_CLNSIGS,
    _CLNSIG_MAGNITUDE,
    _OUTSIDE_CONFLICT_CAP,
)
from allelix.annotators.gwas import _UNKNOWN_RISK_ALLELE_MAG_CAP
from allelix.annotators.gwas import _magnitude as gwas_magnitude
from allelix.annotators.pharmgkb import _LOE_MAGNITUDE
from allelix.annotators.pharmgkb import _magnitude as clinpgx_magnitude
from allelix.reports._pipeline import COMMON_ALLELE_AF, COMMON_ALLELE_MAGNITUDE, AnalysisResult
from allelix.reports.html import render_html

_ROW = re.compile(r"<tr><td><strong>([^<]+)</strong></td><td>(.*?)</td></tr>")
_OR_CI = "[1.10-1.30]"  # the bare ratio interval that marks an odds ratio


@pytest.fixture(scope="module")
def legend(tmp_path_factory) -> str:
    out = tmp_path_factory.mktemp("legend") / "r.html"
    render_html(
        AnalysisResult(
            file_path=Path("genotype.txt"),
            parser_name="vcf",
            parser_display_name="VCF",
            sample_id="S1",
            build="GRCh38",
            total_variants=0,
            skipped_count=0,
            annotators_used=[],
            annotations=[],
        ),
        output_path=out,
    )
    page = out.read_text()
    start = page.index("<strong>Understanding Magnitude Scores</strong>")
    return page[start : page.index("</details>", start)]


def _section(legend: str, heading: str) -> str:
    start = legend.index(f"<h3>{heading}")
    end = legend.find("<h3>", start + 1)
    return legend[start : end if end != -1 else None]


def _rows(section: str) -> list[tuple[float, str]]:
    return [(float(score), text) for score, text in _ROW.findall(section)]


def _label(term: str) -> str:
    text = term.replace("_", " ").replace("vus-", "VUS-")
    return text[:1].upper() + text[1:]


def _num(value: float) -> str:
    return f"{value:g}"


class TestClinVar:
    def test_every_classification_with_its_score(self, legend: str) -> None:
        listed: dict[str, float] = {}
        for score, text in _rows(_section(legend, "ClinVar")):
            for label in text.split(" / "):
                listed[label.removesuffix(" (shown with --include-benign)")] = score
        assert listed == {_label(term): score for term, score in _CLNSIG_MAGNITUDE.items()}

    def test_scores_are_listed_highest_first_once_each(self, legend: str) -> None:
        scores = [score for score, _text in _rows(_section(legend, "ClinVar"))]
        assert scores == sorted(set(_CLNSIG_MAGNITUDE.values()), reverse=True)

    def test_benign_rows_say_how_to_show_them(self, legend: str) -> None:
        section = _section(legend, "ClinVar")
        for term in _BENIGN_CLNSIGS:
            assert f"{_label(term)} (shown with --include-benign)" in section

    def test_review_status_row_is_gone(self, legend: str) -> None:
        assert "single submitter" not in legend
        assert "no assertion criteria)" not in _section(legend, "ClinVar").split("<p>")[1]

    def test_aggregate_caps(self, legend: str) -> None:
        section = _section(legend, "ClinVar")
        for term, cap in _AGGREGATE_CAP.items():
            assert f"{_label(term)} {_num(cap)}" in section
        assert f"are scored {_num(_OUTSIDE_CONFLICT_CAP)} when their review tier" in section

    def test_common_zero_star_rule(self, legend: str) -> None:
        section = _section(legend, "ClinVar")
        assert (
            f"is scored {_num(COMMON_ALLELE_MAGNITUDE)} when gnomAD gives ClinVar&rsquo;s "
            f"allele a frequency of {_num(COMMON_ALLELE_AF * 100)}% or more"
        ) in section

    def test_unscored_labels_are_not_shown(self, legend: str) -> None:
        assert "A classification without a score here is not shown" in legend


class TestClinPGx:
    def test_every_level_with_its_score(self, legend: str) -> None:
        listed = {}
        for score, text in _rows(_section(legend, "ClinPGx")):
            level = re.match(r"Level (\w+) &mdash;", text)
            assert level is not None, text
            listed[level.group(1).lower()] = score
        assert listed == _LOE_MAGNITUDE
        assert all(clinpgx_magnitude(level) == score for level, score in listed.items())


class TestGwas:
    def _tiers(self, legend: str) -> tuple[list[tuple[float, float]], float]:
        tiers = []
        rest = None
        for score, text in _rows(_section(legend, "GWAS")):
            p = re.match(r"p &le; (\d)&times;10<sup>(-?\d+)</sup>", text)
            if p is None:
                assert text == "weaker or no p-value"
                rest = score
            else:
                tiers.append((float(f"{p.group(1)}e{p.group(2)}"), score))
        assert rest is not None
        return tiers, rest

    def test_each_tier_matches_the_code(self, legend: str) -> None:
        tiers, rest = self._tiers(legend)
        assert [p for p, _s in tiers] == sorted(p for p, _s in tiers)
        for p, score in tiers:
            assert gwas_magnitude(p, None) == score  # the bound is inclusive
        assert gwas_magnitude(None, None) == rest

    def test_every_p_value_scores_as_the_legend_says(self, legend: str) -> None:
        """p from 1e-325 to 1 in quarter decades, and each tier bound and just above it."""
        tiers, rest = self._tiers(legend)
        probes = [10 ** (e / 4) for e in range(-1300, 1)]
        probes += [p * f for p, _s in tiers for f in (1.0, 1.0001)]
        for p in probes:
            expected = next((score for bound, score in tiers if p <= bound), rest)
            assert gwas_magnitude(p, None) == expected, p

    def test_odds_ratio_steps(self, legend: str) -> None:
        section = _section(legend, "GWAS")
        steps = re.findall(r"\+([\d.]+) at ([\d.]+) or more \(or ([\d.]+) or less\)", section)
        assert len(steps) == 2
        base = gwas_magnitude(5e-8, None)
        previous_high = None
        for added, high, low in steps:
            added_f, high_f, low_f = float(added), float(high), float(low)
            assert gwas_magnitude(5e-8, high_f, _OR_CI) == base + added_f
            assert gwas_magnitude(5e-8, low_f, _OR_CI) == base + added_f
            if previous_high is not None:
                just_below = (high_f + previous_high) / 2
                assert gwas_magnitude(5e-8, just_below, _OR_CI) == base + added_f
            previous_high = high_f
        # Below the last step: no change; without a confirmed odds ratio: no change.
        assert gwas_magnitude(5e-8, 1.5, _OR_CI) == base
        assert gwas_magnitude(5e-8, 5.0, None) == base

    def test_ceiling_and_unknown_allele_cap(self, legend: str) -> None:
        section = _section(legend, "GWAS")
        ceiling = float(re.search(r"at most ([\d.]+)\.", section).group(1))
        assert gwas_magnitude(5e-100, 5.0, _OR_CI) == min(
            gwas_magnitude(5e-100, None) + 1, ceiling
        )
        assert gwas_magnitude(1e-300, 5.0, _OR_CI) <= ceiling
        assert f"the score is at most {_num(_UNKNOWN_RISK_ALLELE_MAG_CAP)}." in section

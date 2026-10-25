# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Self-contained HTML report renderer.

The output is a single ``.html`` file with inline CSS and JS. No external
assets, works from ``file://``. Per ADR-0003, every annotation carries its
source attribution and the page header restates the informational posture.

v1.8.0 redesign: 5-column table (Magnitude, Gene, Genotype, Repute,
Summary) with annotations grouped by variant and ``genotype_match``. Clicking a
row opens a sliding detail sidebar showing all source annotations vertically.
"""

from __future__ import annotations

import html
import json
from collections import defaultdict
from datetime import UTC, datetime
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.annotators.clinvar import (
    _AGGREGATE_CAP,
    _BENIGN_CLNSIGS,
    _CLNSIG_MAGNITUDE,
    _OUTSIDE_CONFLICT_CAP,
)
from allelix.annotators.gwas import _UNKNOWN_RISK_ALLELE_MAG_CAP
from allelix.annotators.gwas import _magnitude as _gwas_magnitude
from allelix.annotators.pharmgkb import _LOE_MAGNITUDE
from allelix.models import annotation_label
from allelix.reports import (
    REGULATORY_NOTICE,
    array_quality_lines,
    atomic_write_text,
    input_left_unchecked,
    input_method_label,
    is_imputed,
    vcf_quality_lines,
)
from allelix.reports._pipeline import (
    COMMON_ALLELE_AF,
    COMMON_ALLELE_MAGNITUDE,
    PANEL_NOT_ON_CHIP,
    PANEL_STATES,
    _enrichment_allele,
    is_vcf_parser,
    panel_reason_label,
    panel_state_label,
    rollup_gwas_duplicates,
)
from allelix.reports._sources import SKIPPED
from allelix.reports.association_text import effect_words, strand_abstention_lines
from allelix.reports.build_text import INCONCLUSIVE, MISMATCH, build_notices, short_source
from allelix.reports.clinvar_text import common_allele_line, detail_fields, scored_down_line
from allelix.utils.plausibility import (
    CHROMOSOME,
    GNOMAD_NOT_RUN,
    IMPLAUSIBLE,
    NO_ALLELE,
    NO_FREQUENCY,
    NOT_ASSESSED,
    SUSPICIOUS,
    most_severe,
)
from allelix.utils.plausibility import HALF_CALL as REASON_HALF_CALL

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult, PanelCoverage
    from allelix.reports._sources import ReportWarning, SourceRecord
    from allelix.reports.diff import DiffResult


# ---------------------------------------------------------------------------
# Repute classification
# ---------------------------------------------------------------------------

_BAD_SIGNIFICANCE = frozenset(
    {
        "clinvar_pathogenic",
        "clinvar_pathogenic/likely_pathogenic",
        "clinvar_likely_pathogenic",
        "clinvar_risk_factor",
        "clinvar_pathogenic,_low_penetrance",
        "clinvar_likely_pathogenic,_low_penetrance",
        "clinvar_established_risk_allele",
        "clinvar_confers_sensitivity",
        "clinvar_likely_risk_allele",
        "snpedia_bad",
    }
)

_GOOD_SIGNIFICANCE = frozenset(
    {
        "clinvar_benign",
        "clinvar_benign/likely_benign",
        "clinvar_likely_benign",
        "clinvar_protective",
        "snpedia_good",
    }
)


def _classify_repute(significance: str) -> str:
    """Derive CSS class from the significance field."""
    sig = significance.lower()
    if sig in _BAD_SIGNIFICANCE:
        return "repute-bad"
    if sig in _GOOD_SIGNIFICANCE:
        return "repute-good"
    return "repute-neutral"


def _get_repute(ann: Annotation) -> str:
    """Return ``'bad'``, ``'good'``, or ``'neutral'`` for an annotation.

    A ClinVar row scored down to ClinVar's own Conflicting or Benign view
    (ADR-0045), or as a zero-star submission on a common allele (ADR-0052),
    is neutral, whatever its submitter said.
    """
    if ann.aggregate_capped or ann.common_allele_capped:
        return "neutral"
    return _classify_repute(ann.significance).removeprefix("repute-")


# ---------------------------------------------------------------------------
# HTML escaping
# ---------------------------------------------------------------------------


def _escape(value: str) -> str:
    return html.escape(value or "", quote=True)


# ---------------------------------------------------------------------------
# Static assets
# ---------------------------------------------------------------------------

_CSS = """\
:root {
  --bg: #fafafa;
  --bg-surface: #fff;
  --text: #212121;
  --text-muted: #757575;
  --border: #e0e0e0;
  --border-light: #f0f0f0;
  --hover: rgba(0, 0, 0, 0.04);
  --selected: rgba(25, 118, 210, 0.08);
  --backdrop: rgba(0, 0, 0, 0.3);
  --panel-bg: #fff;
  --panel-shadow: rgba(0, 0, 0, 0.15);
  --notice-bg: #fff8e1;
  --notice-border: #f9a825;
  --notice-warn-bg: #fff3e0;
  --notice-warn-border: #e65100;
  --link: #1976d2;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #121212;
    --bg-surface: #1e1e1e;
    --text: #e0e0e0;
    --text-muted: #9e9e9e;
    --border: #333;
    --border-light: #2a2a2a;
    --hover: rgba(255, 255, 255, 0.06);
    --selected: rgba(100, 181, 246, 0.12);
    --backdrop: rgba(0, 0, 0, 0.5);
    --panel-bg: #1e1e1e;
    --panel-shadow: rgba(0, 0, 0, 0.4);
    --notice-bg: #332b00;
    --notice-border: #f9a825;
    --notice-warn-bg: #331a00;
    --notice-warn-border: #e65100;
    --link: #93c5fd;
  }
}
[data-theme="dark"] {
  --bg: #121212;
  --bg-surface: #1e1e1e;
  --text: #e0e0e0;
  --text-muted: #9e9e9e;
  --border: #333;
  --border-light: #2a2a2a;
  --hover: rgba(255, 255, 255, 0.06);
  --selected: rgba(100, 181, 246, 0.12);
  --backdrop: rgba(0, 0, 0, 0.5);
  --panel-bg: #1e1e1e;
  --panel-shadow: rgba(0, 0, 0, 0.4);
  --notice-bg: #332b00;
  --notice-border: #f9a825;
  --notice-warn-bg: #331a00;
  --notice-warn-border: #e65100;
  --link: #93c5fd;
}

*, *::before, *::after { box-sizing: border-box; }

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  font-size: 14px;
  line-height: 1.5;
  color: var(--text);
  background: var(--bg);
  padding: 24px;
  margin: 0;
}
h1 { margin-bottom: .25rem; }
.subtitle { color: var(--text-muted); margin-top: 0; }
a { color: var(--link); }

.notice {
  background: var(--notice-bg, #fff8e1); border-left: 4px solid var(--notice-border, #f9a825);
  padding: 1rem; margin: 1.5rem 0; border-radius: 4px; font-size: .95rem;
}
.sources-off { color: var(--text-muted); font-size: 0.9em; margin: 0.5em 0; }
.details-summary { margin: .5rem 0; }
.tech-details { margin: 0 0 1.5rem; }
.tech-details > summary { cursor: pointer; font-weight: 600; }
.notice-warn {
  background: var(--notice-warn-bg, #fff3e0);
  border-left: 4px solid var(--notice-warn-border, #e65100);
  padding: 1rem; margin: 1.5rem 0; border-radius: 4px; font-size: .95rem;
}
.notice table.panel-coverage {
  margin-top: .75rem; border-collapse: collapse; font-size: .85rem;
}
.notice table.panel-coverage th,
.notice table.panel-coverage td {
  text-align: left; padding: .25rem .75rem .25rem 0;
  border-bottom: 1px solid var(--border);
}
.notice table.panel-coverage th { font-weight: 600; }
.education {
  background: var(--bg-surface); border-left: 4px solid var(--border);
  padding: 1rem 1.25rem; margin: 1.5rem 0; border-radius: 4px; font-size: .9rem;
}
.education h2 { margin-top: 0; font-size: 1.05rem; }
.education h3 { font-size: .95rem; margin: .75rem 0 .25rem; }
.education p { margin: .35rem 0; }
details.education summary { cursor: pointer; font-size: 1.05rem; }
details.education[open] summary { margin-bottom: .5rem; }

/* Summary cards */
.summary {
  display: flex; flex-wrap: wrap; gap: .75rem; margin-bottom: 1.5rem;
}
.card {
  background: var(--bg-surface); padding: .75rem 1rem; border-radius: 8px;
  flex: 1 1 140px; min-width: 110px;
  border: 1px solid var(--border);
}
.card .label {
  font-size: .75rem; color: var(--text-muted);
  text-transform: uppercase; letter-spacing: .05em;
}
.card .value { font-size: 1.2rem; font-weight: 600; }
.card .build-source { font-size: .75rem; font-weight: 400; color: var(--text-muted); }
.card-bad .value { color: #c62828; }
.card-good .value { color: #2e7d32; }

/* Controls */
.controls {
  display: flex; flex-wrap: wrap; gap: .75rem;
  align-items: center; margin-bottom: 1rem;
}
#search {
  flex: 1 1 250px; padding: 8px 12px; border: 1px solid var(--border);
  border-radius: 6px; font-size: 14px; outline: none;
  background: var(--bg-surface); color: var(--text);
}
#search:focus { border-color: #1976d2; box-shadow: 0 0 0 2px rgba(25,118,210,.15); }
.filters { display: flex; gap: 4px; }
.filter-btn {
  padding: 6px 14px; border: 1px solid var(--border); border-radius: 6px;
  background: var(--bg-surface); cursor: pointer; font-size: 13px; color: var(--text);
}
.filter-btn:hover { background: var(--hover); }
.filter-btn.active {
  background: #1976d2; color: #fff; border-color: #1976d2;
}

/* Table */
table {
  width: 100%; table-layout: fixed; border-collapse: collapse;
  background: var(--bg-surface); border: 1px solid var(--border);
  border-radius: 8px; overflow: hidden;
}
thead th {
  text-align: left; padding: 10px 12px;
  font-size: 11px; font-weight: 600;
  text-transform: uppercase; letter-spacing: 0.5px;
  color: var(--text-muted); background: var(--bg);
  border-bottom: 2px solid var(--border); user-select: none;
}
th.sortable { cursor: pointer; }
th.sortable:hover { color: var(--text); }
th .sort-arrow { font-size: .7rem; margin-left: .25rem; color: var(--text-muted); }

th:nth-child(1) { width: 50px; }
th:nth-child(2) { width: 100px; }
th:nth-child(3) { width: 80px; }
th:nth-child(4) { width: 80px; }

tbody td {
  padding: 8px 12px; font-size: 13px;
  border-bottom: 1px solid var(--border-light); vertical-align: middle;
}
tbody tr { cursor: pointer; transition: background-color 0.15s; }
tbody tr:hover { background: var(--hover); }
tbody tr.selected {
  background: var(--selected);
  box-shadow: inset 3px 0 0 #1976d2;
}

.gene-cell, .gt-cell {
  font-family: 'SF Mono', 'Cascadia Code', Consolas, monospace;
  font-size: 12px;
}
.plaus { margin-left: 6px; cursor: help; font-family: sans-serif; }
.plaus-implausible { color: var(--notice-warn-border); }
.plaus-suspicious  { color: var(--notice-border); }
.plaus-not_assessed { color: var(--text-muted); }
.sum-cell {
  overflow: hidden; text-overflow: ellipsis;
  white-space: nowrap; max-width: 0;
}

/* Badges and pills */
.badge {
  display: inline-flex; align-items: center; justify-content: center;
  width: 28px; height: 28px; border-radius: 50%;
  font-weight: 700; font-size: 13px; line-height: 1;
}
.badge-bad     { background: #ef5350; color: #fff; }
.badge-good    { background: #66bb6a; color: #fff; }
.badge-neutral { background: #bdbdbd; color: #424242; }

.pill {
  display: inline-block; padding: 2px 10px; border-radius: 12px;
  font-size: 12px; font-weight: 600; white-space: nowrap;
}
.pill-bad     { background: #ffebee; color: #c62828; border: 1px solid #ef9a9a; }
.pill-good    { background: #e8f5e9; color: #2e7d32; border: 1px solid #a5d6a7; }
.pill-neutral { background: var(--bg); color: var(--text-muted); border: 1px solid var(--border); }

/* Sidebar */
.sidebar-backdrop {
  display: none; position: fixed;
  top: 0; left: 0; right: 0; bottom: 0;
  background: var(--backdrop); z-index: 999;
}
.sidebar-backdrop.open { display: block; }

.detail-panel {
  position: fixed; top: 0; right: 0;
  width: 480px; height: 100vh;
  background: var(--panel-bg);
  box-shadow: -2px 0 12px var(--panel-shadow);
  transform: translateX(100%);
  transition: transform 0.25s ease;
  z-index: 1000; overflow-y: auto;
  display: flex; flex-direction: column;
}
.detail-panel.open { transform: translateX(0); }

.panel-header {
  position: sticky; top: 0; background: var(--panel-bg); z-index: 1;
  display: flex; justify-content: space-between; align-items: center;
  padding: 12px 16px; border-bottom: 1px solid var(--border);
}
.panel-close {
  font-size: 24px; background: none; border: none;
  cursor: pointer; color: var(--text-muted); padding: 0 4px; line-height: 1;
}
.panel-close:hover { color: var(--text); }
.panel-nav {
  display: flex; align-items: center; gap: 8px;
  font-size: 13px; color: var(--text-muted);
}
.panel-nav button {
  background: none; border: 1px solid var(--border); border-radius: 4px;
  padding: 2px 8px; cursor: pointer; font-size: 13px; color: var(--text);
}
.panel-nav button:hover { background: var(--hover); }

.panel-body { padding: 16px; flex: 1; }
.panel-gene { font-size: 20px; font-weight: 700; margin: 0; }
.panel-rsid { font-size: 14px; color: var(--text-muted); margin: 2px 0 8px; }
.panel-meta { display: flex; gap: 8px; align-items: center; margin-bottom: 16px; }

.field-row {
  display: grid; grid-template-columns: 120px 1fr;
  gap: 4px 12px; padding: 4px 0; font-size: 13px; line-height: 1.5;
}
.field-label { color: var(--text-muted); font-size: 12px; }
.field-value { color: var(--text); word-break: break-word; }

.source-header {
  font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px;
  color: var(--text-muted); border-bottom: 1px solid var(--border);
  padding: 16px 0 4px; margin-top: 8px;
}
.ann-entry { padding-bottom: 8px; margin-bottom: 8px; }
.ann-entry:not(:last-child) {
  border-bottom: 1px solid var(--border-light);
}

.empty { padding: 2rem; text-align: center; color: var(--text-muted); font-style: italic; }

footer {
  margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border);
  font-size: .8rem; color: var(--text-muted);
}

/* Theme toggle */
.theme-toggle {
  background: none; border: 1px solid var(--border); border-radius: 6px;
  padding: 6px 10px; cursor: pointer; font-size: 16px; line-height: 1;
  color: var(--text);
}
.theme-toggle:hover { background: var(--hover); }

@media (max-width: 1023px) {
  .detail-panel { width: 100%; }
  th:nth-child(2), td:nth-child(2) { width: 80px; }
  th:nth-child(3), td:nth-child(3) { width: 60px; }
  th:nth-child(4), td:nth-child(4) { width: 70px; }
}
"""

# The array caveat in "Reading This Report", and its replacement for input
# whose genotypes are imputed (MyHeritage low-pass WGS).
_PSEUDOGENE_NOTE = """\
<h3>Pseudogene cross-hybridization</h3>
<p>Some genes have known pseudogenes with high sequence similarity
(e.g., PKD1 has six pseudogenes sharing &gt;97% identity). Array-based
genotyping probes can cross-hybridize, producing false genotype calls.
If a result seems inconsistent with your health history, confirmatory
testing by a different method is recommended.</p>
"""

# How each imputed input was produced, for the "Imputed genotypes" note.
_IMPUTATION_ORIGIN: dict[str, str] = {
    "low-pass-wgs-imputed": (
        "This file comes from low-pass whole genome sequencing: each position was\n"
        "read only about twice, and most genotypes were filled in (imputed) from a\n"
        "reference panel of other people."
    ),
    "array-imputed": (
        "This file extends a genotyping array by imputation: most of its genotypes\n"
        "were not measured but filled in (imputed) from a reference panel of other\n"
        "people."
    ),
}
_IMPUTATION_ORIGIN_DEFAULT = (
    "Most genotypes in this file were filled in (imputed) from a reference panel\n"
    "of other people rather than measured."
)
_IMPUTATION_NOTE = """\
<h3>Imputed genotypes</h3>
<p>{origin} Imputation is accurate for common variants
and least reliable for rare ones, which are the variants most findings in this
report are about. Confirm any finding by direct genotyping before acting on
it.</p>
"""


def _education_section(input_method: str) -> str:
    """The "Reading This Report" section, with the caveat that fits the input."""
    if is_imputed(input_method):
        origin = _IMPUTATION_ORIGIN.get(input_method, _IMPUTATION_ORIGIN_DEFAULT)
        return _EDUCATION_SECTION.replace(_PSEUDOGENE_NOTE, _IMPUTATION_NOTE.format(origin=origin))
    return _EDUCATION_SECTION


_EDUCATION_SECTION = """\
<details class="education">
<summary><strong>Reading This Report</strong></summary>

<h3>Pseudogene cross-hybridization</h3>
<p>Some genes have known pseudogenes with high sequence similarity
(e.g., PKD1 has six pseudogenes sharing &gt;97% identity). Array-based
genotyping probes can cross-hybridize, producing false genotype calls.
If a result seems inconsistent with your health history, confirmatory
testing by a different method is recommended.</p>

<h3>ClinVar aggregation</h3>
<p>ClinVar aggregates submissions from multiple sources. Different submitters
may classify the same variant differently, and significance labels may be
paired with conditions from a different submitter&rsquo;s entry. Single-submitter
entries carry less weight than expert-panel-reviewed classifications.</p>

<h3>Carrier vs. affected</h3>
<p>A variant classified as pathogenic in a recessive condition requires two
copies to cause disease. If you are heterozygous (one copy), you are a carrier.
Carrier status does not typically cause symptoms but may be relevant for family
planning.</p>
<p>A call marked <strong>Hemizygous</strong> is one the file declares as a single
copy, such as an X chromosome call in a male. For an X-linked recessive
condition, one copy of a pathogenic allele on a male's single X can cause the
condition; it does not make him a carrier in the sense above.</p>

<h3>Confirmatory testing</h3>
<p>No genotyping platform is 100% accurate. Clinically significant findings
should be confirmed with an independent method before making medical
decisions.</p>
</details>
"""

_MAGNITUDE_LEGEND_TEMPLATE = """\
<details class="education">
<summary><strong>Understanding Magnitude Scores</strong></summary>

<p>Each annotation source uses its own criteria to assign a magnitude score
(0&ndash;10) that reflects clinical importance. Higher scores warrant more
attention. The score shown in the table is the maximum across all source
annotations for that variant.</p>

{scored_sources}
<h3>SNPedia (community-curated)</h3>
<p>Magnitude is assigned directly by community editors (0&ndash;10 scale).
Higher scores indicate greater clinical or personal relevance as judged by
contributors. See <a href="https://www.snpedia.com/index.php/Magnitude">\
SNPedia&rsquo;s magnitude documentation</a> for details.</p>

<h3>CADD (variant deleteriousness)</h3>
<p>CADD PHRED scores rank how deleterious a variant is relative to all
possible human SNVs. Higher scores = more likely to be deleterious.</p>
<table style="width:auto; font-size:.85rem; margin:.5rem 0;">
<tr><td><strong>&ge; 30</strong></td><td>Top 0.1% most deleterious</td></tr>
<tr><td><strong>&ge; 20</strong></td><td>Top 1% most deleterious</td></tr>
<tr><td><strong>&ge; 10</strong></td><td>Top 10% most deleterious</td></tr>
<tr><td><strong>&lt; 10</strong></td><td>Below top 10%</td></tr>
</table>

<h3>AlphaMissense (missense pathogenicity)</h3>
<p>DeepMind&rsquo;s protein-structure-based pathogenicity prediction for
missense variants. Score 0&ndash;1; higher = more likely pathogenic.</p>
<table style="width:auto; font-size:.85rem; margin:.5rem 0;">
<tr><td><strong>&ge; 0.564</strong></td><td>Likely pathogenic</td></tr>
<tr><td><strong>0.340 &ndash; 0.564</strong></td><td>Ambiguous</td></tr>
<tr><td><strong>&lt; 0.340</strong></td><td>Likely benign</td></tr>
</table>

</details>
"""

# The ClinVar, ClinPGx and GWAS sections of the legend are generated from the
# scoring code, so the legend cannot drift from the scores again.
_LEGEND_TABLE = '<table style="width:auto; font-size:.85rem; margin:.5rem 0;">'

_LOE_TEXT: dict[str, str] = {
    "1a": "CPIC guideline or FDA label",
    "1b": "strong clinical evidence",
    "2a": "moderate evidence",
    "2b": "moderate evidence, weaker replication",
    "3": "a single study or unreplicated evidence",
    "4": "case reports, non-significant studies or in vitro evidence",
}

# GWAS p-value tiers as ``_magnitude`` tests them (``p <= threshold``),
# strongest first, with a note; the score shown is the code's own.
_GWAS_TIERS: tuple[tuple[float, str], ...] = (
    (5e-100, ""),
    (5e-20, ""),
    (5e-8, " (genome-wide significance)"),
    (5e-6, " (suggestive)"),
    (5e-4, " (nominal)"),
)
# ``_magnitude``'s odds-ratio steps: (at or above, at or below, added).
_GWAS_ODDS_RATIO_STEPS: tuple[tuple[float, float, float], ...] = (
    (3.0, 0.33, 1.0),
    (2.0, 0.5, 0.5),
)
_GWAS_MAGNITUDE_CEILING = 9.0


def _num(value: float) -> str:
    """``9.0`` as ``9``, ``8.5`` as ``8.5``."""
    return f"{value:g}"


def _p_html(p: float) -> str:
    """``5e-08`` as ``5&times;10<sup>-8</sup>``."""
    mantissa, exponent = f"{p:.0e}".split("e")
    return f"{mantissa}&times;10<sup>{int(exponent)}</sup>"


def _clinvar_label(term: str) -> str:
    """``likely_pathogenic,_low_penetrance`` as ``Likely pathogenic, low penetrance``."""
    text = term.replace("_", " ").replace("vus-", "VUS-")
    return text[:1].upper() + text[1:]


def _legend_rows(rows: list[tuple[str, str]]) -> str:
    cells = "".join(
        f"<tr><td><strong>{score}</strong></td><td>{text}</td></tr>" for score, text in rows
    )
    return f"{_LEGEND_TABLE}{cells}</table>"


def _clinvar_legend() -> str:
    by_score: dict[float, list[str]] = defaultdict(list)
    for term, score in _CLNSIG_MAGNITUDE.items():
        label = _clinvar_label(term)
        if term in _BENIGN_CLNSIGS:
            label += " (shown with --include-benign)"
        by_score[score].append(label)
    rows = [(_num(score), " / ".join(by_score[score])) for score in sorted(by_score, reverse=True)]
    caps = ", ".join(
        f"{_escape(_clinvar_label(term))} {_num(cap)}"
        for term, cap in sorted(_AGGREGATE_CAP.items(), key=lambda item: -item[1])
    )
    return (
        "<h3>ClinVar (clinical significance)</h3>\n"
        "<p>Each submission is scored by its own classification:</p>\n"
        f"{_legend_rows(rows)}\n"
        "<p>A disease classification is never scored above ClinVar&rsquo;s own "
        f"classification of the variant: {caps}. Outranked Pathogenic or Likely pathogenic "
        "submissions kept beside a classification of another kind are scored "
        f"{_num(_OUTSIDE_CONFLICT_CAP)} when their review tier also holds other disease "
        "classifications.</p>\n"
        "<p>A Pathogenic or Likely pathogenic submission with zero stars (no assertion "
        "criteria, or no classification provided) is scored "
        f"{_num(COMMON_ALLELE_MAGNITUDE)} when gnomAD gives ClinVar&rsquo;s allele a "
        f"frequency of {_num(COMMON_ALLELE_AF * 100)}% or more.</p>\n"
        "<p>A classification without a score here is not shown; the report counts it.</p>\n"
    )


def _clinpgx_legend() -> str:
    rows = [
        (_num(score), f"Level {level.upper()} &mdash; {_LOE_TEXT[level]}")
        for level, score in sorted(_LOE_MAGNITUDE.items(), key=lambda item: -item[1])
    ]
    return f"<h3>ClinPGx (pharmacogenomic evidence)</h3>\n{_legend_rows(rows)}\n"


def _gwas_legend() -> str:
    rows = [
        (_num(_gwas_magnitude(threshold, None)), f"p &le; {_p_html(threshold)}{note}")
        for threshold, note in _GWAS_TIERS
    ]
    rows.append((_num(_gwas_magnitude(None, None)), "weaker or no p-value"))
    steps = "; ".join(
        f"+{_num(added)} at {_num(high)} or more (or {_num(low)} or less)"
        for high, low, added in _GWAS_ODDS_RATIO_STEPS
    )
    return (
        "<h3>GWAS Catalog (trait associations)</h3>\n"
        f"{_legend_rows(rows)}\n"
        "<p>When the study reports a confirmed odds ratio, the score rises: "
        f"{steps}; at most {_num(_GWAS_MAGNITUDE_CEILING)}. When the study does not name "
        f"its risk allele, the score is at most {_num(_UNKNOWN_RISK_ALLELE_MAG_CAP)}.</p>\n"
    )


_MAGNITUDE_LEGEND = _MAGNITUDE_LEGEND_TEMPLATE.replace(
    "{scored_sources}", _clinvar_legend() + "\n" + _clinpgx_legend() + "\n" + _gwas_legend()
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _input_method_banner(method: str) -> str:
    """The input's declared genotyping method, e.g. MyHeritage low-pass WGS."""
    label = input_method_label(method)
    if not label:
        return ""
    caution = (
        " Rare variants are the least reliable imputed calls; confirm any finding by"
        " direct genotyping."
        if is_imputed(method)
        else ""
    )
    return (
        '<div class="notice"><strong>Input method.</strong> '
        f"{_escape(label[0].upper() + label[1:])}.{_escape(caution)}</div>"
    )


def _vcf_quality_banner(lines: list[str]) -> str:
    """What the input's quality checks found: VCF gating, array calls ClinVar cannot match."""
    if not lines:
        return ""
    items = "".join(f"<li>{_escape(line)}</li>" for line in lines)
    return f'<div class="notice"><strong>Input quality.</strong><ul>{items}</ul></div>'


def _report_warnings_banner(warnings: list[ReportWarning]) -> str:
    """Report-level warnings (``report_warnings``), one ``notice-warn`` each."""
    return "".join(
        f'<div class="notice-warn"><strong>Warning.</strong> {_escape(w.message)}</div>'
        for w in warnings
    )


# A report warning that already names a source's failure, by warning code:
# that source is left out of "Sources not used", so the top says it once.
_WARNING_NAMES_SOURCE = {"clinvar_skipped": "clinvar"}


def _sources_failed_banner(sources: list[SourceRecord], warnings: list[ReportWarning]) -> str:
    """Sources whose cache could not be used: the user can act on these.

    A source a report warning already names (``_WARNING_NAMES_SOURCE``) is
    left out; with nothing else, there is no box.
    """
    warned = {_WARNING_NAMES_SOURCE[w.code] for w in warnings if w.code in _WARNING_NAMES_SOURCE}
    failed = [s for s in sources if s.failed and s.name not in warned]
    if not failed:
        return ""
    items = "".join(f"<li>{_escape(s.message)}</li>" for s in failed)
    return f'<div class="notice"><strong>Sources not used.</strong><ul>{items}</ul></div>'


def _sources_off_line(sources: list[SourceRecord]) -> str:
    """Sources the user turned off."""
    chosen = [s for s in sources if s.state == SKIPPED and not s.failed]
    if not chosen:
        return ""
    names = ", ".join(f"{_escape(s.display_name)} ({_escape(s.choice_text)})" for s in chosen)
    return f'<p class="sources-off">Not run by choice: {names}.</p>'


def _clinvar_banner(scored_down: int, unrecognised: int, common_allele: int) -> str:
    """ADR-0045, ADR-0052: ClinVar rows that left the default report, and why."""
    lines = [
        line for line in (scored_down_line(scored_down), common_allele_line(common_allele)) if line
    ]
    if unrecognised:
        lines.append(
            f"{unrecognised} ClinVar row(s) not shown: their classification is not one "
            "Allelix scores yet."
        )
    if not lines:
        return ""
    items = "".join(f"<li>{_escape(line)}</li>" for line in lines)
    return f'<div class="notice"><strong>ClinVar.</strong><ul>{items}</ul></div>'


def _strand_banner(abstentions: dict[str, dict[str, int]]) -> str:
    """GWAS and SNPedia rows not assessed because their strand was not settled."""
    lines = strand_abstention_lines(abstentions)
    if not lines:
        return ""
    items = "".join(f"<li>{_escape(line)}</li>" for line in lines)
    return f'<div class="notice"><strong>Strand.</strong><ul>{items}</ul></div>'


# The top of the HTML report: a short disclaimer and the warnings that change
# what a reader should do stay visible; the other notices, unchanged, sit in a
# collapsed "Technical details" section. JSON and terminal output keep them all.
_SHORT_DISCLAIMER = (
    "This report shows what public databases say about the variants in your file. "
    "It is for information only, not medical advice or a diagnosis."
)
_DETAILS_UNCHECKED = "Some calls could not be checked."
_DETAILS_BELOW = "Details below."
_DETAILS_SUMMARY_PLAIN = "More about how this report was made is in the details below."


def _rows(count: int, verb: str) -> str:
    """``1 ClinVar row was <verb>`` or ``N ClinVar rows were <verb>``."""
    if count == 1:
        return f"1 ClinVar row was {verb}."
    return f"{count:,} ClinVar rows were {verb}."


def _details_summary(*, unchecked: bool, scored_down: int, not_shown: int) -> str:
    """The line above "Technical details", naming only what the section holds.

    ``unchecked``: calls not checked (input quality or strand);
    ``scored_down``: ClinVar rows scored down (aggregate cap and the
    zero-star common-allele rule); ``not_shown``: ClinVar rows with a
    classification Allelix does not score. With none of them, the
    section holds only general notes and the line says so.
    """
    parts = [_DETAILS_UNCHECKED] if unchecked else []
    if scored_down:
        parts.append(_rows(scored_down, "scored down"))
    if not_shown:
        parts.append(_rows(not_shown, "not shown"))
    if not parts:
        return _DETAILS_SUMMARY_PLAIN
    return " ".join([*parts, _DETAILS_BELOW])


def _technical_details(notices: list[str], *, summary: str) -> str:
    """The summary line and the collapsed section holding ``notices`` (HTML)."""
    line = _escape(summary)
    return (
        f'<p class="details-summary">{line}</p>'
        '<details class="tech-details"><summary>Technical details</summary>'
        f"{''.join(notices)}</details>"
    )


def _hv_nocall_banner(warnings: list[str] | None) -> str:
    if not warnings:
        return ""
    items = "".join(f"<li>{_escape(w)}</li>" for w in warnings)
    return (
        '<div class="notice-warn"><strong>High-value no-calls.</strong> '
        "The following clinically important SNPs returned no genotype call:"
        f"<ul>{items}</ul></div>"
    )


def _panel_coverage_section(coverage: PanelCoverage | None, parser_name: str = "") -> str:
    """Render the Panel coverage section (ADR-0051); quiet when no panel was supplied.

    One row per requested rsID with its gene, genotype, state and number
    of shown annotations. The four state counts sum to the number
    requested. VCF input has no chip: an rsID absent from it is "not in
    file", which may mean reference at that position or no variant call.
    """
    if coverage is None:
        return ""
    vcf = is_vcf_parser(parser_name)
    counts = ", ".join(
        f"{coverage['counts'][state]} {panel_state_label(state, vcf=vcf)}"
        for state in PANEL_STATES
    )
    note = (
        " An rsID not in the file may be reference at that position, or have no variant call."
        if vcf and coverage["counts"][PANEL_NOT_ON_CHIP]
        else ""
    )
    rows = "".join(
        f"<tr><td><code>{_escape(v['rsid'])}</code></td>"
        f"<td>{_escape(v['gene'] or '')}</td>"
        f"<td>{'<code>' + _escape(v['genotype']) + '</code>' if v['genotype'] else ''}</td>"
        f"<td>{panel_state_label(v['state'], vcf=vcf)}</td>"
        + (f"<td>{_escape(panel_reason_label(v['reason']))}</td>" if vcf else "")
        + f"<td>{v['annotation_count']}</td></tr>"
        for v in coverage["variants"]
    )
    table = (
        "<table class='panel-coverage'><thead><tr><th>rsID</th><th>Gene</th>"
        "<th>Genotype</th><th>State</th>"
        + ("<th>Reason</th>" if vcf else "")
        + "<th>Annotations</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
        if rows
        else ""
    )
    return (
        '<div class="notice"><strong>Panel coverage.</strong> '
        f"{coverage['requested']} requested: {counts}.{note}{table}</div>"
    )


def _license_attributions(annotators_used: list[tuple[str, str | None]]) -> str:
    """Build license attribution HTML from annotator LicenseDescriptors."""
    import logging

    from allelix.annotators import get_annotator_class

    logger = logging.getLogger(__name__)
    parts: list[str] = []
    for name, _version in annotators_used:
        cls = get_annotator_class(name)
        if cls is None:
            logger.warning("No annotator class found for '%s' — attribution omitted", name)
            continue
        desc = cls.license
        source_link = desc.source_url or desc.license_url
        parts.append(
            f" <a href='{html.escape(source_link)}'>"
            f"{html.escape(cls.display_name)}</a>: "
            f"{html.escape(desc.attribution_text)}"
            f" (<a href='{html.escape(desc.license_url)}'>license</a>)"
        )
    return "".join(parts)


def _review_stars(status: int | str | None) -> str:
    """Render ClinVar review status as filled/empty stars."""
    if status is None or status == "" or status == "—":
        return ""
    if isinstance(status, str):
        star_count = status.count("_") if "_" in status else 0
        if "expert" in status.lower() or "practice" in status.lower():
            star_count = 4
        elif "multiple" in status.lower():
            star_count = 3
        elif "single" in status.lower():
            star_count = 1
        elif "no_assertion" in status.lower() or "no assertion" in status.lower():
            star_count = 0
    else:
        star_count = int(status)
    star_count = max(0, min(star_count, 4))
    return "★" * star_count + "☆" * (4 - star_count)


# ---------------------------------------------------------------------------
# Variant grouping
# ---------------------------------------------------------------------------


def _group_annotations(
    annotations: list[Annotation],
) -> list[list[Annotation]]:
    """Group annotations by variant and ``genotype_match``; sort by max magnitude.

    The variant is its rsID, or its locus when it has none (ADR-0046).
    """
    groups: dict[tuple[str, str], list[Annotation]] = defaultdict(list)
    for ann in annotations:
        key = (annotation_label(ann), ann.genotype_match)
        groups[key].append(ann)
    return sorted(
        groups.values(),
        key=lambda g: max(a.magnitude for a in g),
        reverse=True,
    )


# ---------------------------------------------------------------------------
# Variant data JSON
# ---------------------------------------------------------------------------


def _add_metrics(target: dict, ann: Annotation) -> None:
    """Copy the annotation's enrichment values into ``target`` where it has none yet."""
    if ann.allele_frequency is not None and "allele_frequency" not in target:
        target["allele_frequency"] = round(ann.allele_frequency, 6)
    if ann.am_pathogenicity is not None and "am_pathogenicity" not in target:
        target["am_pathogenicity"] = round(ann.am_pathogenicity, 4)
    if ann.am_class and "am_class" not in target:
        target["am_class"] = ann.am_class
    if ann.cadd_phred is not None and "cadd_phred" not in target:
        target["cadd_phred"] = round(ann.cadd_phred, 1)


def _build_variant_data(sorted_groups: list[list[Annotation]]) -> str:
    """Serialize grouped annotations to JSON for the detail sidebar."""
    variant_data = []
    for group in sorted_groups:
        display = max(group, key=lambda a: a.magnitude)
        entry: dict = {
            "rsid": annotation_label(display),
            "gene": display.gene or "",
            "genotype": display.genotype_match,
            "zygosity": display.zygosity,
            "plausibility": most_severe([a.plausibility for a in group]),
            "plausibility_reason": _not_assessed_reason(group),
            "annotations": [],
        }
        ordered = sorted(group, key=lambda a: -a.magnitude)
        # gnomAD, AlphaMissense and CADD describe one allele each. A variant
        # with more than one allele lists them per allele, never on one line.
        alleles = [_enrichment_allele(ann) or ann.alt for ann in ordered]
        several = len(set(alleles) - {""}) > 1
        per_allele: dict[str, dict] = {}
        for ann, allele in zip(ordered, alleles, strict=True):
            _add_metrics(per_allele.setdefault(allele, {}) if several else entry, ann)
        entry["allele_metrics"] = [
            {"allele": allele, **metrics} for allele, metrics in per_allele.items() if metrics
        ]
        for ann in ordered:
            a: dict = {"source": ann.source, "magnitude": ann.magnitude}
            if ann.significance:
                a["significance"] = ann.significance
            if ann.review_status:
                a["reviewStatus"] = ann.review_status
                a["reviewStars"] = _review_stars(ann.review_status)
            # Label and value pairs for the detail panel: ClinVar's own
            # classification and the split (ADR-0045), a GWAS row's effect;
            # always present, as the panel script reads it.
            a["clinvar"] = [list(pair) for pair in detail_fields(ann)]
            effect = effect_words(ann)
            if effect:
                a["clinvar"].append(["Effect", effect])
            if ann.condition:
                a["condition"] = ann.condition
            if ann.description:
                a["description"] = ann.description
            if ann.references:
                a["references"] = ann.references
            if ann.attribution:
                a["attribution"] = ann.attribution
            a["repute"] = _get_repute(ann)
            if ann.category:
                a["category"] = ann.category
            entry["annotations"].append(a)
        variant_data.append(entry)

    json_str = json.dumps(variant_data, ensure_ascii=False)
    return json_str.replace("</", "<\\/")


# ---------------------------------------------------------------------------
# Row rendering
# ---------------------------------------------------------------------------


def _build_search_text(group: list[Annotation]) -> str:
    """Concatenate all searchable fields from all annotations in a group."""
    parts: list[str] = []
    display = max(group, key=lambda a: a.magnitude)
    parts.extend(
        [
            annotation_label(display),
            display.gene,
            display.genotype_match,
            display.zygosity,
            most_severe([a.plausibility for a in group]),
        ]
    )
    for ann in group:
        parts.extend(
            [
                ann.source,
                ann.significance,
                ann.condition,
                ann.description,
                ann.attribution,
                ann.gene,
                ann.am_class,
            ]
        )
    return " ".join(p for p in parts if p).lower()


def _summary_text(group: list[Annotation]) -> str:
    """Build the summary cell content."""
    display = max(group, key=lambda a: a.magnitude)
    source_label = _escape(display.attribution or display.source)
    extra = len(group) - 1
    prefix = f"[{source_label} +{extra}]" if extra > 0 else f"[{source_label}]"
    text = _escape(display.condition or display.description or "")
    return f"{prefix} {text}"


# ADR-0041: hover text for the warning icon beside a flagged genotype.
_PLAUSIBILITY_TITLES: dict[str, str] = {
    IMPLAUSIBLE: (
        "Implausible call: homozygous for an allele that is very rare in gnomAD. "
        "Often a genotyping artifact; confirm before relying on it."
    ),
    SUSPICIOUS: (
        "Suspicious call: gnomAD reports no carriers of this allele. Confirm before relying on it."
    ),
}


# ADR-0050: why a call's plausibility was not assessed.
_NOT_ASSESSED_TEXT: dict[str, str] = {
    NO_FREQUENCY: "gnomAD has no frequency for this allele",
    GNOMAD_NOT_RUN: "gnomAD was not used in this run",
    CHROMOSOME: "no rule for this call's copy count on this chromosome",
    REASON_HALF_CALL: "one allele is missing, so the copy count is unknown",
    NO_ALLELE: "the annotation names no carried allele to look up",
}


def _not_assessed_reason(group: list[Annotation]) -> str:
    """The display reason of the group's first not-assessed annotation."""
    for a in group:
        if a.plausibility == NOT_ASSESSED:
            return _NOT_ASSESSED_TEXT.get(a.plausibility_reason, a.plausibility_reason)
    return ""


def _plausibility_icon(plausibility: str, reason: str = "") -> str:
    if plausibility == NOT_ASSESSED:
        title = f"Plausibility not assessed: {reason}" if reason else "Plausibility not assessed"
        return (
            f'<span class="plaus plaus-not_assessed" title="{_escape(title)}"'
            f' aria-label="{_escape(title)}">○</span>'
        )
    warning = _PLAUSIBILITY_TITLES.get(plausibility)
    if warning is None:
        return ""
    return (
        f'<span class="plaus plaus-{plausibility}" title="{_escape(warning)}"'
        f' aria-label="{_escape(warning)}">⚠</span>'
    )


def _row_html(group: list[Annotation], row_id: int) -> str:
    """Render a single variant group as an HTML table row."""
    display = max(group, key=lambda a: a.magnitude)
    repute = _get_repute(display)
    mag = max(a.magnitude for a in group)
    search_text = _build_search_text(group)
    plausibility = most_severe([a.plausibility for a in group])

    return (
        f'<tr data-row-id="{row_id}"'
        f' data-magnitude="{mag:.1f}"'
        f' data-gene="{_escape(display.gene)}"'
        f' data-genotype="{_escape(display.genotype_match)}"'
        f' data-repute="{_escape(repute)}"'
        f' data-plausibility="{_escape(plausibility)}"'
        f' data-search-text="{_escape(search_text)}">'
        f'<td class="mag-cell">'
        f'<span class="badge badge-{repute}">{int(mag)}</span></td>'
        f'<td class="gene-cell">{_escape(display.gene) or "—"}</td>'
        f'<td class="gt-cell">{_escape(display.genotype_match)}'
        f"{_plausibility_icon(plausibility, _not_assessed_reason(group))}</td>"
        f'<td class="repute-cell">'
        f'<span class="pill pill-{repute}">{repute.capitalize()}</span></td>'
        f'<td class="sum-cell">{_summary_text(group)}</td>'
        "</tr>"
    )


# ---------------------------------------------------------------------------
# Inline JavaScript
# ---------------------------------------------------------------------------

_SCRIPT = """\
<script>
document.addEventListener("DOMContentLoaded", function() {
  var variantData = JSON.parse(
    document.getElementById("variant-data").textContent
  );
  var searchInput = document.getElementById("search");
  var tbody = document.querySelector("#variant-table tbody");
  if (!tbody) return;
  var panel = document.getElementById("detail-panel");
  var backdrop = document.getElementById("backdrop");
  var panelBody = document.getElementById("panel-body");
  var panelPos = document.getElementById("panel-position");
  var activeFilter = "all";
  var selectedIndex = -1;

  /* --- Filter counts --- */
  var allRows = Array.from(tbody.querySelectorAll("tr"));
  document.getElementById("count-all").textContent = allRows.length;
  document.getElementById("count-bad").textContent =
    allRows.filter(function(r){ return r.dataset.repute === "bad"; }).length;
  document.getElementById("count-good").textContent =
    allRows.filter(function(r){ return r.dataset.repute === "good"; }).length;
  document.getElementById("count-neutral").textContent =
    allRows.filter(function(r){ return r.dataset.repute === "neutral"; }).length;

  /* --- Search + filter --- */
  function getVisibleRows() {
    return Array.from(tbody.querySelectorAll("tr")).filter(
      function(r) { return r.style.display !== "none"; }
    );
  }

  function applyFilters() {
    var term = searchInput.value.toLowerCase();
    allRows.forEach(function(row) {
      var matchesSearch = !term || row.dataset.searchText.indexOf(term) !== -1;
      var matchesFilter =
        activeFilter === "all" || row.dataset.repute === activeFilter;
      row.style.display = matchesSearch && matchesFilter ? "" : "none";
    });
    selectedIndex = -1;
    allRows.forEach(function(r) { r.classList.remove("selected"); });
  }

  searchInput.addEventListener("input", applyFilters);

  document.querySelectorAll(".filter-btn").forEach(function(btn) {
    btn.addEventListener("click", function() {
      document.querySelectorAll(".filter-btn").forEach(function(b) {
        b.classList.remove("active");
      });
      btn.classList.add("active");
      activeFilter = btn.dataset.filter;
      applyFilters();
    });
  });

  /* --- Sort --- */
  document.querySelectorAll("th.sortable").forEach(function(th) {
    th.addEventListener("click", function() {
      var key = th.dataset.sort;
      var rows = Array.from(tbody.querySelectorAll("tr"));

      var wasAsc = th.classList.contains("sort-asc");
      document.querySelectorAll("th.sortable").forEach(function(h) {
        h.classList.remove("sort-asc", "sort-desc", "sort-active");
        var arrow = h.querySelector(".sort-arrow");
        if (arrow) arrow.textContent = "";
      });

      var dir = wasAsc ? -1 : 1;
      th.classList.add("sort-active", dir === 1 ? "sort-asc" : "sort-desc");
      var arrow = th.querySelector(".sort-arrow");
      if (arrow) arrow.textContent = dir === 1 ? "\\u25B2" : "\\u25BC";

      var reputeOrder = {"bad": 0, "neutral": 1, "good": 2};
      rows.sort(function(a, b) {
        var av = a.dataset[key];
        var bv = b.dataset[key];
        if (key === "magnitude") return dir * (Number(av) - Number(bv));
        if (key === "repute") {
          var ao = reputeOrder[av] !== undefined ? reputeOrder[av] : 1;
          var bo = reputeOrder[bv] !== undefined ? reputeOrder[bv] : 1;
          return dir * (ao - bo);
        }
        return dir * av.localeCompare(bv);
      });
      rows.forEach(function(row) { tbody.appendChild(row); });
    });
  });

  /* --- Sidebar --- */
  function populatePanel(row) {
    var idx = Number(row.dataset.rowId);
    var v = variantData[idx];
    if (!v) return;

    var visible = getVisibleRows();
    var pos = visible.indexOf(row) + 1;
    panelPos.textContent = pos + " of " + visible.length;

    var h = '<h2 class="panel-gene">' + esc(v.gene || "\\u2014") + "</h2>";
    h += '<div class="panel-rsid">' + esc(v.rsid) + "</div>";

    var dispAnn = v.annotations[0];
    var repute = dispAnn ? dispAnn.repute : "neutral";
    var mag = dispAnn ? Math.floor(dispAnn.magnitude) : 0;
    h += '<div class="panel-meta">';
    h += '<span class="badge badge-' + repute + '">' + mag + "</span>";
    h += '<span class="pill pill-' + repute + '">' +
         repute.charAt(0).toUpperCase() + repute.slice(1) + "</span>";
    h += "</div>";

    h += fieldRow("Genotype", v.genotype);
    h += fieldRow("Zygosity", v.zygosity);
    if (v.plausibility === "implausible")
      h += fieldRow("Plausibility",
        "\u26a0 Implausible: homozygous for an allele that is very rare " +
        "in gnomAD; often a genotyping artifact");
    else if (v.plausibility === "suspicious")
      h += fieldRow("Plausibility",
        "\u26a0 Suspicious: gnomAD reports no carriers of this allele");
    else if (v.plausibility === "not_assessed")
      h += fieldRow("Plausibility", "Not assessed: " + v.plausibility_reason);

    var hasMetrics = v.allele_frequency != null ||
                     v.am_pathogenicity != null || v.cadd_phred != null;
    if (hasMetrics || v.allele_metrics.length) {
      h += '<div class="source-header">VARIANT METRICS</div>';
      if (hasMetrics) h += metricRows(v, "");
      // More than one allele: each value is shown with the allele it is for.
      for (var mi = 0; mi < v.allele_metrics.length; mi++) {
        var m = v.allele_metrics[mi];
        h += metricRows(m, m.allele + ": ");
      }
    }

    var groups = {};
    var groupOrder = [];
    for (var i = 0; i < v.annotations.length; i++) {
      var a = v.annotations[i];
      var src = a.source;
      if (!groups[src]) { groups[src] = []; groupOrder.push(src); }
      groups[src].push(a);
    }
    for (var gi = 0; gi < groupOrder.length; gi++) {
      var src = groupOrder[gi];
      var anns = groups[src];
      h += '<div class="source-header">' + esc(src) + "</div>";
      for (var ai = 0; ai < anns.length; ai++) {
        var a = anns[ai];
        h += '<div class="ann-entry">';
        h += fieldRow("Magnitude", a.magnitude.toFixed(1));
        if (a.significance)
          h += fieldRow("Significance", a.significance);
        if (a.reviewStars)
          h += fieldRow("Review Status", a.reviewStars);
        for (var ci = 0; ci < a.clinvar.length; ci++)
          h += fieldRow(a.clinvar[ci][0], a.clinvar[ci][1]);
        if (a.condition) h += fieldRow("Condition", a.condition);
        if (a.description)
          h += fieldRow("Description", a.description);
        if (a.attribution)
          h += fieldRow("Attribution", a.attribution);
        if (a.references && a.references.length)
          h += fieldRow("References", a.references.join(" "));
        h += "</div>";
      }
    }
    panelBody.innerHTML = h;
  }

  function metricRows(m, prefix) {
    var rows = "";
    if (m.allele_frequency != null)
      rows += fieldRow("Frequency",
        prefix + (m.allele_frequency * 100).toFixed(1) + "%");
    if (m.am_pathogenicity != null) {
      var amText = m.am_pathogenicity.toFixed(3);
      if (m.am_class) amText += " (" + m.am_class + ")";
      rows += fieldRow("AlphaMissense", prefix + amText);
    }
    if (m.cadd_phred != null) {
      var cs = m.cadd_phred;
      var tier = cs >= 30 ? "top 0.1%"
               : cs >= 20 ? "top 1%"
               : cs >= 10 ? "top 10%" : "";
      var label = cs.toFixed(1) +
        (tier ? " (" + tier + " most deleterious)" : "");
      rows += fieldRow("CADD PHRED", prefix + label);
    }
    return rows;
  }

  function fieldRow(label, value) {
    return '<div class="field-row"><div class="field-label">' +
           esc(label) + '</div><div class="field-value">' +
           esc(String(value)) + "</div></div>";
  }

  function esc(s) {
    var d = document.createElement("div");
    d.appendChild(document.createTextNode(s));
    return d.innerHTML;
  }

  function openPanel(row) {
    populatePanel(row);
    panel.classList.add("open");
    backdrop.classList.add("open");
  }

  function closePanel() {
    panel.classList.remove("open");
    backdrop.classList.remove("open");
  }

  function selectRow(row) {
    allRows.forEach(function(r) { r.classList.remove("selected"); });
    row.classList.add("selected");
    row.scrollIntoView({ block: "nearest" });
    if (panel.classList.contains("open")) {
      populatePanel(row);
    }
  }

  /* Row click */
  tbody.addEventListener("click", function(e) {
    var row = e.target.closest("tr");
    if (!row) return;
    var rows = getVisibleRows();
    selectedIndex = rows.indexOf(row);
    selectRow(row);
    openPanel(row);
  });

  /* Close */
  document.getElementById("panel-close").addEventListener("click", closePanel);
  backdrop.addEventListener("click", closePanel);

  /* Prev / Next */
  document.getElementById("panel-prev").addEventListener("click", function() {
    var rows = getVisibleRows();
    if (selectedIndex > 0) {
      selectedIndex--;
      selectRow(rows[selectedIndex]);
    }
  });
  document.getElementById("panel-next").addEventListener("click", function() {
    var rows = getVisibleRows();
    if (selectedIndex < rows.length - 1) {
      selectedIndex++;
      selectRow(rows[selectedIndex]);
    }
  });

  /* Keyboard */
  document.addEventListener("keydown", function(e) {
    if (document.activeElement === searchInput) return;
    var rows = getVisibleRows();
    if (!rows.length) return;

    if (e.key === "Escape") { closePanel(); return; }
    if (e.key === "ArrowDown") {
      e.preventDefault();
      selectedIndex = Math.min(selectedIndex + 1, rows.length - 1);
      selectRow(rows[selectedIndex]);
    }
    if (e.key === "ArrowUp") {
      e.preventDefault();
      selectedIndex = Math.max(selectedIndex - 1, 0);
      selectRow(rows[selectedIndex]);
    }
    if (e.key === "Enter" && selectedIndex >= 0) {
      openPanel(rows[selectedIndex]);
    }
  });

  /* --- Theme toggle --- */
  var toggleBtn = document.getElementById("theme-toggle");
  function setTheme(theme) {
    document.documentElement.setAttribute("data-theme", theme);
    toggleBtn.textContent = theme === "dark" ? "\\u2600" : "\\u263E";
  }
  toggleBtn.addEventListener("click", function() {
    var current = document.documentElement.getAttribute("data-theme");
    if (!current) {
      var prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
      setTheme(prefersDark ? "light" : "dark");
    } else {
      setTheme(current === "dark" ? "light" : "dark");
    }
  });
});
</script>
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def render_html(
    result: AnalysisResult,
    *,
    output_path: Path,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    rsids: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
    title: str = "Allelix Genotype Report",
    diff: DiffResult | None = None,
    high_value_no_calls: list[str] | None = None,
) -> int:
    """Write a self-contained HTML report. Returns the number of annotations rendered."""
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        rsids=rsids,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    # Use display_name (e.g. "ClinPGx") in the subtitle, not the lowercase
    # internal identifier (e.g. "pharmgkb"). The terminal renderer already
    # does this via `a.display_name` — keep HTML consistent.
    from allelix.annotators import get_annotator_class

    def _display(name: str) -> str:
        cls = get_annotator_class(name)
        return cls.display_name if cls is not None else name

    annotators_str = ", ".join(
        f"{_display(name)} ({version or 'unknown'})" for name, version in result.annotators_used
    )

    build_warn = ""
    diag = result.build_diagnostics
    if diag is not None:
        # The same wording as the CLI (ADR-0044); the report shows the label
        # mismatch and the tie, the notices that change how to read it.
        for notice in build_notices(diag):
            if notice.kind in (MISMATCH, INCONCLUSIVE):
                lead, _, rest = notice.text.partition(": ")
                build_warn += (
                    f'<div class="notice-warn"><strong>{_escape(lead)}.</strong> '
                    f"{_escape(rest)}</div>"
                )
    build_card = _escape(result.build)
    if diag is not None:
        build_card += f'<div class="build-source">{_escape(short_source(diag))}</div>'

    diff_banner = ""
    if diff is not None:
        from allelix.reports.diff import summarize_diff

        diff_banner = (
            f'<div class="notice"><strong>Diff: </strong>{_escape(summarize_diff(diff))}</div>'
        )

    floor_note = ""
    if source_min_magnitudes and min_magnitude > 0:
        ran = {name for name, _version in result.annotators_used}
        lower_floors = {
            src: floor
            for src, floor in source_min_magnitudes.items()
            # A GWAS floor never lowers the global one (ADR-0024).
            if floor < min_magnitude and src in ran and src != "gwas"
        }
        if lower_floors:
            parts = ", ".join(
                f"{src.upper()} (mag &ge; {floor:.1f})"
                for src, floor in sorted(lower_floors.items())
            )
            floor_note = (
                '<div class="notice">'
                f"<strong>Source-specific thresholds.</strong> "
                f"The global minimum magnitude is {min_magnitude:.1f}, "
                f"but lower thresholds apply to: {parts}. "
                "Some rows below the global threshold may appear from these sources."
                "</div>"
            )

    # Group annotations by variant
    sorted_groups = _group_annotations(filtered)

    # Count reputes for summary cards
    def _display_repute(g: list[Annotation]) -> str:
        return _get_repute(max(g, key=lambda a: a.magnitude))

    bad_count = sum(1 for g in sorted_groups if _display_repute(g) == "bad")
    good_count = sum(1 for g in sorted_groups if _display_repute(g) == "good")

    if sorted_groups:
        rows_html = "\n".join(_row_html(g, i) for i, g in enumerate(sorted_groups))
        variant_json = _build_variant_data(sorted_groups)

        body = (
            '<table id="variant-table">'
            "<thead><tr>"
            '<th data-sort="magnitude" class="sortable sort-active sort-desc">'
            'Mag<span class="sort-arrow">&#x25BC;</span></th>'
            '<th data-sort="gene" class="sortable">'
            'Gene<span class="sort-arrow"></span></th>'
            '<th data-sort="genotype" class="sortable">'
            'GT<span class="sort-arrow"></span></th>'
            '<th data-sort="repute" class="sortable">'
            'Repute<span class="sort-arrow"></span></th>'
            "<th>Summary</th>"
            "</tr></thead>"
            f"<tbody>{rows_html}</tbody></table>"
            f'<script id="variant-data" type="application/json">{variant_json}</script>'
        )

        sidebar = (
            '<div class="sidebar-backdrop" id="backdrop"></div>'
            '<div class="detail-panel" id="detail-panel">'
            '<div class="panel-header">'
            '<button class="panel-close" id="panel-close" aria-label="Close">&times;</button>'
            '<div class="panel-nav">'
            '<button id="panel-prev" aria-label="Previous">&lsaquo; Prev</button>'
            '<span id="panel-position"></span>'
            '<button id="panel-next" aria-label="Next">Next &rsaquo;</button>'
            "</div></div>"
            '<div class="panel-body" id="panel-body"></div>'
            "</div>"
        )
    else:
        body = '<div class="empty">No annotations matched the current filters.</div>'
        sidebar = ""
        variant_json = ""

    controls = (
        '<div class="controls">'
        '<input type="text" id="search" placeholder="Search variants..." '
        'aria-label="Search variants">'
        '<div class="filters">'
        '<button class="filter-btn active" data-filter="all">'
        'All <span id="count-all"></span></button>'
        '<button class="filter-btn" data-filter="bad">'
        'Bad <span id="count-bad"></span></button>'
        '<button class="filter-btn" data-filter="good">'
        'Good <span id="count-good"></span></button>'
        '<button class="filter-btn" data-filter="neutral">'
        'Neutral <span id="count-neutral"></span></button>'
        "</div>"
        '<button class="theme-toggle" id="theme-toggle" '
        'aria-label="Toggle dark mode" title="Toggle dark/light mode">'
        "&#x263E;</button>"
        "</div>"
    )

    input_quality = vcf_quality_lines(result.vcf_quality) + array_quality_lines(
        result.array_quality, result.probe_classes
    )
    common = result.clinvar_common_allele
    # An imputed input's caution (confirm by direct genotyping) stays in view.
    imputed = is_imputed(result.input_method)
    counts = [
        _vcf_quality_banner(input_quality),
        _clinvar_banner(result.clinvar_scored_down, result.clinvar_unrecognised, common),
        _strand_banner(result.strand_abstentions),
    ]
    summary_line = _details_summary(
        unchecked=input_left_unchecked(
            result.vcf_quality, result.array_quality, result.probe_classes
        )
        or bool(strand_abstention_lines(result.strand_abstentions)),
        scored_down=result.clinvar_scored_down + common,
        not_shown=result.clinvar_unrecognised,
    )
    details = _technical_details(
        [
            '<div class="notice"><strong>Informational only.</strong> '
            f"{_escape(REGULATORY_NOTICE)}</div>",
            _sources_off_line(result.sources),
            "" if imputed else _input_method_banner(result.input_method),
            *counts,
        ],
        summary=summary_line,
    )
    summary_cards = "\n".join(
        f'<div class="card{css}"><div class="label">{label}</div>'
        f'<div class="value">{value}</div></div>'
        for label, value, css in [
            ("Sample", _escape(result.sample_id) or "(unknown)", ""),
            ("Format", _escape(result.parser_display_name), ""),
            ("Build", build_card, ""),
            ("Variants", f"{len(sorted_groups):,}", ""),
            ("Bad", str(bad_count), " card-bad"),
            ("Good", str(good_count), " card-good"),
            ("Total Annotations", f"{len(filtered):,}", ""),
        ]
    )

    document = (
        "<!DOCTYPE html>"
        "<html lang='en'><head><meta charset='utf-8'>"
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,'
        "%3Csvg viewBox='0 0 32 32' fill='none'"
        " xmlns='http://www.w3.org/2000/svg'%3E"
        "%3Ccircle cx='16' cy='16' r='14' stroke='%236366f1'"
        " stroke-width='1.5' opacity='0.3'/%3E"
        "%3Cpath d='M10 6C10 6 22 12 22 16C22 20 10 26 10 26'"
        " stroke='%236366f1' stroke-width='2'"
        " stroke-linecap='round' fill='none'/%3E"
        "%3Cpath d='M22 6C22 6 10 12 10 16C10 20 22 26 22 26'"
        " stroke='%236366f1' stroke-width='2'"
        " stroke-linecap='round' fill='none' opacity='0.4'/%3E"
        "%3Ccircle cx='10' cy='10' r='2' fill='%236366f1'/%3E"
        "%3Ccircle cx='22' cy='10' r='2' fill='%236366f1'"
        " opacity='0.4'/%3E"
        "%3Ccircle cx='16' cy='16' r='2' fill='%236366f1'"
        " opacity='0.7'/%3E"
        "%3Ccircle cx='22' cy='22' r='2' fill='%236366f1'/%3E"
        "%3Ccircle cx='10' cy='22' r='2' fill='%236366f1'"
        " opacity='0.4'/%3E"
        '%3C/svg%3E">'
        f"<title>{_escape(title)}</title>"
        f"<style>{_CSS}</style>"
        "</head><body>"
        f"<h1>{_escape(title)}</h1>"
        f'<p class="subtitle">Source: <code>{_escape(result.file_path.name)}</code> · '
        f"Annotators: {_escape(annotators_str)} · "
        f"Generated {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}</p>"
        f'<div class="notice disclaimer">{_SHORT_DISCLAIMER}</div>'
        f"{build_warn}"
        f"{_report_warnings_banner(result.report_warnings)}"
        f"{_sources_failed_banner(result.sources, result.report_warnings)}"
        f"{_hv_nocall_banner(high_value_no_calls)}"
        f"{_input_method_banner(result.input_method) if imputed else ''}"
        f"{details}"
        f"{_panel_coverage_section(result.panel_coverage(filtered), result.parser_name)}"
        f"{_education_section(result.input_method)}"
        f"{_MAGNITUDE_LEGEND}"
        f"{diff_banner}"
        f'<div class="summary">{summary_cards}</div>'
        f"{floor_note}"
        f"{controls}"
        f"{body}"
        f"{sidebar}"
        f"{_SCRIPT}"
        f"<footer>Generated by Allelix v{_escape(__version__)} — "
        "<a href='https://github.com/allelix/allelix'>github.com/allelix/allelix</a>. "
        "All variant classifications attributed to their source databases."
        f"{_license_attributions(result.annotators_used)}</footer>"
        "</body></html>"
    )
    atomic_write_text(output_path, document, input_path=result.file_path)
    return len(filtered)

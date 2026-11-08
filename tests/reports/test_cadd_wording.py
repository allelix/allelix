# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The HTML report frames a CADD score as a prediction, not a verdict (2.6.2).

"10.4 (top 10% most deleterious)" next to a ClinVar protective finding
marked Good read as "bad for you". The panel says "predicted deleterious",
and the legend says that CADD predicts how disruptive a change is, not
whether it is good or bad for you. Synthetic data only.
"""

from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING

import pytest

from allelix.models import Annotation
from allelix.reports.html import render_html
from tests.reports.test_html import _result

if TYPE_CHECKING:
    from pathlib import Path

_LEGEND = (
    "CADD predicts how disruptive a change is, not whether its effect is good or bad "
    "for you: ClinVar and the other sources say that, and a protective variant can "
    "score high."
)


def _protective(cadd_phred: float) -> Annotation:
    return Annotation(
        source="clinvar",
        rsid="rs3999508001",
        significance="clinvar_protective",
        category="clinical",
        magnitude=2.0,
        description="ClinVar classifies this allele as protective",
        attribution="ClinVar",
        genotype_match="AG",
        gene="GENE508",
        condition="Synthetic condition, protection against",
        cadd_phred=cadd_phred,
    )


def _rendered(tmp_path: Path, cadd_phred: float) -> tuple[str, dict]:
    out = tmp_path / "report.html"
    render_html(_result([_protective(cadd_phred)]), output_path=out)
    body = out.read_text(encoding="utf-8")
    start = body.index(">", body.index('id="variant-data"')) + 1
    (entry,) = json.loads(body[start : body.index("</script>", start)].replace("<\\/", "</"))
    return body, entry


def test_a_protective_finding_shows_its_cadd_score_as_a_prediction(tmp_path):
    body, entry = _rendered(tmp_path, 10.4)

    assert entry["annotations"][0]["repute"] == "good"
    assert entry["cadd_label"] == "10.4 (top 10% predicted deleterious)"
    assert 'fieldRow("CADD PHRED", prefix + m.cadd_label)' in body
    assert "most deleterious" not in body


def test_the_legend_says_cadd_does_not_judge_good_or_bad(tmp_path):
    body, _ = _rendered(tmp_path, 10.4)
    text = " ".join(re.sub(r"<[^>]+>", " ", body).split())

    assert "CADD (predicted deleteriousness)" in text
    assert _LEGEND in text
    for tier in ("Top 0.1%", "Top 1%", "Top 10%"):
        assert f"{tier} predicted deleterious" in text
    assert "Below the top 10%" in text


@pytest.mark.parametrize(
    ("phred", "label"),
    [
        (38.0, "38.0 (top 0.1% predicted deleterious)"),
        (30.0, "30.0 (top 0.1% predicted deleterious)"),
        (20.0, "20.0 (top 1% predicted deleterious)"),
        (10.0, "10.0 (top 10% predicted deleterious)"),
        (9.9, "9.9"),
        (0.1, "0.1"),
    ],
)
def test_the_tiers_stay_10_20_30(phred, label):
    from allelix.reports.html import _cadd_label

    assert _cadd_label(phred) == label

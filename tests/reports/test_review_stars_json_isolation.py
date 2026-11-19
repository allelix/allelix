"""HTML review stars must leave the shared findings and their JSON untouched."""

import json

import pytest

import allelix.reports.html as html_report
import allelix.reports.json_report as json_report
from tests.reports.test_review_stars_json_unchanged import (
    _FixedClock,
    _result,
    _review_stars_270,
)


@pytest.mark.parametrize("min_magnitude", [0.0, 5.0])
def test_html_stars_preserve_raw_statuses_and_all_json_bytes(tmp_path, monkeypatch, min_magnitude):
    monkeypatch.setattr(json_report, "datetime", _FixedClock)
    fixed_stars = html_report._review_stars
    outputs = {}
    pages = {}
    for name, stars in (("2.7.0", _review_stars_270), ("fixed", fixed_stars)):
        monkeypatch.setattr(html_report, "_review_stars", stars)
        result = _result()
        expected = {a.rsid: a.review_status for a in result.annotations}
        before = tmp_path / f"{name}-before.json"
        after = tmp_path / f"{name}-after.json"
        json_report.render_json(result, output_path=before, min_magnitude=min_magnitude)
        rows = json.loads(before.read_text())["annotations"]
        assert {a["rsid"]: a["review_status"] for a in rows} == expected
        page = tmp_path / f"{name}.html"
        html_report.render_html(result, output_path=page)
        json_report.render_json(result, output_path=after, min_magnitude=min_magnitude)
        assert after.read_bytes() == before.read_bytes(), "HTML changed the shared result's JSON"
        outputs[name] = after.read_bytes()
        pages[name] = page.read_text()
    assert outputs["2.7.0"] == outputs["fixed"]
    assert pages["2.7.0"] != pages["fixed"]

"""The no-table controller exits before dereferencing its absent data node."""

import re

import pytest
from test_html_site_identity import ann, result

from allelix.reports.html import render_html


@pytest.mark.parametrize("shape", ["empty", "below_floor"])
def test_no_table_page_guards_missing_data_before_access(tmp_path, shape):
    output = tmp_path / "empty.html"
    render_html(
        result([] if shape == "empty" else [ann(magnitude=0)]), output_path=output, min_magnitude=5
    )
    text = output.read_text()
    assert 'id="variant-table"' not in text and 'id="variant-data"' not in text
    script = re.findall(r"<script>(.*?)</script>", text, re.S)[0]
    lookup = re.search(
        r'var\s+(\w+)\s*=\s*document.querySelector\("#variant-table tbody"\)', script
    )
    assert lookup, "The existing no-table controller guard must remain"
    guard = re.search(r"if\s*\(!" + lookup[1] + r"\)\s*(?:\{\s*)?return\b", script)
    read = script.index('document.getElementById("variant-data")')
    assert guard and guard.start() < read, (
        "An empty page must return before accessing missing variant-data"
    )


def test_nonempty_page_keeps_data_and_controller(tmp_path):
    output = tmp_path / "one.html"
    render_html(result([ann()]), output_path=output)
    text = output.read_text()
    assert 'id="variant-table"' in text and 'id="variant-data"' in text
    assert len(re.findall(r"<script\b", text)) == 2

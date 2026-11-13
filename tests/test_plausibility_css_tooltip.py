"""Warning explanations use immediate CSS hover, preserving accessible text."""

import re
from html.parser import HTMLParser

import pytest

from allelix.reports import html
from tests.reports.test_html import _ann, _result


class Icons(HTMLParser):
    def __init__(self, page):
        super().__init__(convert_charrefs=True)
        self.icons = []
        self.feed(page)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "plaus" in attrs.get("class", "").split():
            self.icons.append(attrs)


def _page(tmp_path, flag):
    path = tmp_path / "report.html"
    html.render_html(_result([_ann(plausibility=flag)]), output_path=path)
    return path.read_text()


def _rules(page):
    css = "\n".join(re.findall(r"<style[^>]*>(.*?)</style>", page, re.S))
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    return [
        (selector.strip(), body.strip())
        for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css)
    ]


def _tooltip_rules(page):
    rules = [
        (selector, body)
        for selector, body in _rules(page)
        if ".plaus" in selector and ("::after" in selector or "::before" in selector)
    ]
    assert rules, "The warning needs a CSS tooltip; a native title waits for the browser"
    return rules


@pytest.mark.parametrize("flag", ("implausible", "suspicious"))
@pytest.mark.parametrize("surface", ("icon", "page"))
def test_warning_has_aria_explanation_but_no_native_delayed_title(tmp_path, flag, surface):
    page = html._plausibility_icon(flag) if surface == "icon" else _page(tmp_path, flag)
    icons = Icons(page).icons
    assert len(icons) == 1
    assert icons[0]["aria-label"] == html._PLAUSIBILITY_TITLES[flag]
    assert "title" not in icons[0], "Do not retain a second, delayed browser tooltip"
    assert "⚠" in page


@pytest.mark.parametrize("flag", ("implausible", "suspicious"))
def test_css_tooltip_uses_the_accessible_explanation_and_help_cursor(tmp_path, flag):
    page = _page(tmp_path, flag)
    rules = _tooltip_rules(page)
    assert any(re.search(r"content\s*:\s*attr\(\s*aria-label\s*\)", body) for _, body in rules), (
        "Show the same explanation as assistive technology"
    )
    assert any(
        re.search(r"cursor\s*:\s*help\b", body)
        for selector, body in _rules(page)
        if ".plaus" in selector
    )


@pytest.mark.parametrize("flag", ("implausible", "suspicious"))
def test_hover_reveals_the_tooltip_without_timed_visibility(tmp_path, flag):
    rules = _tooltip_rules(_page(tmp_path, flag))
    hover = [body for selector, body in rules if ":hover" in selector]
    assert hover and any(
        re.search(
            r"(?:visibility\s*:\s*visible|opacity\s*:\s*1\b|display\s*:\s*(?:block|inline))", body
        )
        for body in hover
    ), "CSS hover must reveal the explanation"
    assert any(
        re.search(r"(?:visibility\s*:\s*hidden|opacity\s*:\s*0\b|display\s*:\s*none)", body)
        for selector, body in rules
        if ":hover" not in selector
    ), "Hide it until hover"
    for _, body in rules:
        for prop, value in re.findall(r"([\w-]+)\s*:\s*([^;]+)", body):
            if prop.startswith(("transition", "animation")):
                times = re.findall(r"(-?[\d.]+)\s*(?:ms|s)\b", value)
                assert all(float(t) == 0 for t in times), "No tooltip timer or fade delay"


@pytest.mark.parametrize("flag", ("implausible", "suspicious"))
def test_explanation_attribute_round_trips_markup_characters(tmp_path, monkeypatch, flag):
    text = "Warning \"quoted\" & <literal> 'single' \\ tail"
    monkeypatch.setitem(html._PLAUSIBILITY_TITLES, flag, text)
    page = _page(tmp_path, flag)
    icons = Icons(page).icons
    assert icons[0]["aria-label"] == text
    assert "<literal>" not in page
    assert "title" not in icons[0]


@pytest.mark.parametrize("flag", ("not_assessed", "plausible"))
def test_no_new_tooltip_or_mark_for_unflagged_calls(tmp_path, flag):
    assert Icons(_page(tmp_path, flag)).icons == []
    assert html._plausibility_icon(flag) == ""

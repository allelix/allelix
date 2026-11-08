# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Option B: bounded default disclosures, lossless full-detail reports."""

from __future__ import annotations

import contextlib
import json
from html.parser import HTMLParser

import pytest
from click.testing import CliRunner
from rich.console import Console
from test_annotation_parity_and_named_disclosures import _clinvar, _identity
from test_named_unchecked_calls import CONFIG, _chip

from allelix.cli import _helpers as helpers
from allelix.cli import main

BUILDS = ("GRCh37", "GRCh38")
COMMANDS = ("analyze", "methylation", "pharmacogenomics")
GROUPS = (
    "site(s) at rsIDs ClinVar classifies Pathogenic or Likely pathogenic with a call "
    "written as D or I, a deletion or insertion code Allelix cannot check against its "
    "sources. Not assessed: neither a finding nor a normal result.",
    "site(s) at rsIDs with no Pathogenic or Likely pathogenic row of their own with a call "
    "written as D or I, at or near a Pathogenic or Likely pathogenic variant. Allelix "
    "cannot check a deletion or insertion code against its sources.",
    "site(s) at vendor probes without an rsID with a call written as D or I, at or near a "
    "Pathogenic or Likely pathogenic variant. ClinVar does not match a probe to a variant, "
    "so ClinVar did not check them.",
    "site(s) at vendor probes without an rsID with a base call, at a Pathogenic or Likely "
    "pathogenic position. ClinVar does not match a probe to the variant there, so ClinVar "
    "did not check them.",
    "site(s) at a Pathogenic or Likely pathogenic position with a call whose alleles no "
    "source knows for its rsID, so ClinVar did not check them against the variant there.",
    "site(s) at a Pathogenic or Likely pathogenic position with a call naming a base that "
    "is not one of its own rsID's alleles, so ClinVar did not check them against that "
    "variant.",
)
IDS = tuple(f"rs39999749{i:02d}" for i in range(10))
COUNTED = (IDS[0], IDS[1], "i3999974902", "i3999974903", IDS[4], IDS[5])
INLINE = ("rs6025", IDS[7], IDS[8])
HINT = "Run again with --show-unchecked to list every site; a JSON report lists them all."


@pytest.fixture(autouse=True)
def fresh_console(monkeypatch):
    # Replacing the Console avoids leaking its width property into later tests.
    monkeypatch.setattr(helpers, "console", Console(width=2000, color_system=None, emoji=False))


def _fixture(tmp_path, build, *, only_inline=False, reverse=False):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG)
    rows, identities, calls = [], [], []

    def row(rsid, pos, ref, alt, sig="Pathogenic", gene="SYNTH"):
        rows.append((rsid, "1", pos, ref, alt, gene, sig, sig, 1, 1))

    def ident(rsid, pos, ref="G", alt="A"):
        identities.append((rsid, "1", pos, ref, alt, pos, ref, alt, 0))

    for i in range(10):
        pos = 2000000 + i * 100
        rsid = IDS[i]
        other = f"rs39999759{i:02d}"
        if i == 0:
            row(rsid, pos, "G", "GA")
            calls.append((rsid, pos, "DI"))
        elif i == 1:
            row(rsid, pos, "G", "GA", "Benign", "OWN")
            row(other, pos, "GTT", "G", gene="COVERS")
            calls.append((rsid, pos + 1, "DI"))
        elif i == 2:
            row(other, pos, "GTT", "G")
            calls.append((COUNTED[2], pos + 1, "DI"))
        elif i == 3:
            row(other, pos, "G", "A")
            calls.append((COUNTED[3], pos, "AG"))
        elif i == 4:
            row(other, pos, "G", "C")
            calls.append((rsid, pos, "CC"))
        elif i in (5, 7):
            row(rsid, pos, "G", "A", "Benign" if i == 5 else "Pathogenic")
            row(other, pos, "G", "C")
            ident(rsid, pos)
            calls.append((rsid, pos, "CC"))
        elif i == 6:
            calls.append(("rs6025", pos, "DI"))
        elif i == 8:
            row(rsid, pos, "G", "GA")
            row(other, pos, "G", "C")
            calls.append((rsid, pos, "CC"))
        else:
            row(rsid, pos, "G", "A")
            ident(rsid, pos)
            calls.append((rsid, pos, "AG"))
    for b in BUILDS:
        _clinvar(data, b, rows)
    _identity(data, identities)
    if only_inline:
        calls = [c for c in calls if c[0] in (*INLINE, IDS[9])]
    if reverse:
        calls.reverse()
    return data, _chip(tmp_path, calls, build)


def _run(tmp_path, data, source, build, command="analyze", surface="json", show=False):
    output = tmp_path / f"{command}-{surface}-{show}.{surface}"
    args = [command, str(source), "--build", build.lower(), "--data-dir", str(data), "--no-update"]
    if surface != "terminal":
        args += ["--output", str(output)]
    if show:
        args.append("--show-unchecked")
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    value = output.read_text() if surface != "terminal" else result.output
    return (json.loads(value) if surface == "json" else value), result.output


def _flat(text):
    return " ".join(text.split())


class Banner(HTMLParser):
    """Read only the unchecked banner, retaining top-level and detail list items."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.active = False
        self.depth = 0
        self.details = 0
        self.li = None
        self.top = []
        self.all = []
        self.open_details = []
        self.text = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "div" and attrs.get("class") == "notice-warn" and not self.active:
            self.active, self.depth = True, 1
            self.top, self.all, self.text = [], [], []
        elif self.active and tag == "div":
            self.depth += 1
        if not self.active:
            return
        if tag == "details":
            self.details += 1
            self.open_details.append("open" in attrs)
        if tag == "li":
            self.li = []

    def handle_data(self, data):
        if self.active:
            self.text.append(data)
            if self.li is not None:
                self.li.append(data)

    def handle_endtag(self, tag):
        if not self.active:
            return
        if tag == "li" and self.li is not None:
            (self.all if self.details else self.top).append("".join(self.li))
            self.li = None
        if tag == "details":
            self.details -= 1
        if tag == "div":
            self.depth -= 1
            if not self.depth:
                self.active = False
                if "High-value sites not checked." in "".join(self.text):
                    # Stop before an unrelated later notice can replace the result.
                    raise StopIteration


def _banner(raw):
    p = Banner()
    with contextlib.suppress(StopIteration):
        p.feed(raw)
    return p


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", COMMANDS)
def test_default_terminal_counts_six_groups_and_keeps_inline_sites(tmp_path, build, command):
    data, source = _fixture(tmp_path, build)
    report, output = _run(tmp_path, data, source, build, command)
    sites = report["high_value_unchecked"]
    assert {s["rsid"] for s in sites} == set(COUNTED + INLINE), sites
    flat = _flat(output)
    assert (
        "9 high-value or Pathogenic site(s) could not be checked: 3 named here,"
        " 6 counted by reason." in flat
    )
    for site in sites:
        assert (_flat(site["message"]) in flat) == (site["rsid"] in INLINE)
    offsets = [flat.index("1 " + sentence) for sentence in GROUPS]
    assert offsets == sorted(offsets)
    assert HINT in flat


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", COMMANDS)
def test_show_unchecked_preserves_every_original_line(tmp_path, build, command):
    data, source = _fixture(tmp_path, build)
    report, output = _run(tmp_path, data, source, build, command, show=True)
    flat = _flat(output)
    assert "9 high-value or Pathogenic site(s) could not be checked:" in flat
    for site in report["high_value_unchecked"]:
        assert _flat(site["message"]) in flat
    assert "counted by reason" not in flat
    assert HINT not in flat


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", COMMANDS)
def test_json_and_warning_are_unchanged_by_the_terminal_flag(tmp_path, build, command):
    data, source = _fixture(tmp_path, build)
    before, _ = _run(tmp_path, data, source, build, command)
    after, _ = _run(tmp_path, data, source, build, command, show=True)
    before.pop("generated_at", None)
    after.pop("generated_at", None)
    assert before == after
    assert before["schema_version"] == "13"
    assert all(
        set(s) == {"rsid", "gene", "name", "reason", "call", "message"}
        for s in before["high_value_unchecked"]
    )
    warning = next(w for w in before["report_warnings"] if w["code"] == "di_calls_not_checked")
    assert warning["count"] == 4
    for s in before["high_value_unchecked"]:
        if s["reason"] == "di_call":
            assert s["message"] in warning["message"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("command", COMMANDS)
def test_html_collapses_the_complete_list_without_losing_order(tmp_path, build, command):
    data, source = _fixture(tmp_path, build)
    report, _ = _run(tmp_path, data, source, build, command)
    raw, _ = _run(tmp_path, data, source, build, command, "html")
    b = _banner(raw)
    sites = report["high_value_unchecked"]
    expected = [s["message"] for s in sites if s["rsid"] in INLINE] + ["1 " + s for s in GROUPS]
    assert b.top == expected
    assert b.all == [s["message"] for s in sites]
    assert b.open_details == [False]
    assert 'class="unchecked-all"' in raw
    assert "All 9 sites not checked" in raw
    assert raw.index('class="unchecked-all"') < raw.index("<table")


@pytest.mark.parametrize("build", BUILDS)
def test_html_flag_does_not_expand_or_change_the_banner(tmp_path, build):
    data, source = _fixture(tmp_path, build)
    before, _ = _run(tmp_path, data, source, build, surface="html")
    after, _ = _run(tmp_path, data, source, build, surface="html", show=True)
    a, b = _banner(before), _banner(after)
    assert (a.top, a.all, a.open_details, a.text) == (b.top, b.all, b.open_details, b.text)
    assert a.open_details == [False]


@pytest.mark.parametrize("build", BUILDS)
def test_no_counted_sites_keep_the_released_terminal_and_html_block(tmp_path, build):
    data, source = _fixture(tmp_path, build, only_inline=True)
    report, output = _run(tmp_path, data, source, build)
    sites = report["high_value_unchecked"]
    assert {s["rsid"] for s in sites} == set(INLINE)
    assert "3 high-value or Pathogenic site(s) could not be checked:" in _flat(output)
    assert HINT not in output and "counted by reason" not in output
    old_block = "Warning: 3 high-value or Pathogenic site(s) could not be checked:\n"
    old_block += "".join(f"  ⚠ {site['message']}\n" for site in sites)
    assert old_block in output
    for s in sites:
        assert s["message"] in output
    raw, _ = _run(tmp_path, data, source, build, surface="html")
    b = _banner(raw)
    assert b.top == [s["message"] for s in sites]
    assert not b.all and not b.open_details


@pytest.mark.parametrize("build", BUILDS)
def test_input_order_does_not_change_rollup_or_json_site_order(tmp_path, build):
    data, source = _fixture(tmp_path, build)
    before, output1 = _run(tmp_path, data, source, build)
    lines = source.read_text().splitlines(True)
    source.write_text("".join(lines[:3] + list(reversed(lines[3:]))))
    after, output2 = _run(tmp_path, data, source, build)
    assert before["high_value_unchecked"] == after["high_value_unchecked"]
    for sentence in GROUPS:
        assert "1 " + sentence in _flat(output1)
        assert "1 " + sentence in _flat(output2)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Eligibility, lossless serialized data and bounded display for option B."""

from __future__ import annotations

import contextlib
import hashlib
import html
import json
import sqlite3

import pytest
from rich.console import Console
from test_annotation_parity_and_named_disclosures import _clinvar, _identity, _vcf
from test_named_unchecked_calls import CONFIG, OTHER, PLP, _chip, _list_data, _rsid_list
from test_unchecked_rollup_contract import (
    BUILDS,
    GROUPS,
    HINT,
    IDS,
    _banner,
    _fixture,
    _flat,
    _run,
)

from allelix.cli import _helpers as helpers
from allelix.reports.high_value import load_high_value_snps


@pytest.fixture(autouse=True)
def fresh_console(monkeypatch):
    monkeypatch.setattr(helpers, "console", Console(width=2000, color_system=None, emoji=False))


def _digest(report):
    stable = {
        "sites": report["high_value_unchecked"],
        "di_warning": [
            w for w in report["report_warnings"] if w["code"] == "di_calls_not_checked"
        ],
    }
    return hashlib.sha256(
        json.dumps(stable, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@pytest.mark.parametrize("build", BUILDS)
def test_full_json_disclosures_equal_released_261(tmp_path, build):
    data, source = _fixture(tmp_path, build)
    report, _ = _run(tmp_path, data, source, build)
    # Captured using only released 8983c5c, not the implementation branch.
    assert _digest(report) == "47f8871becf55e5dd301d82cb48fa0875c9c5297be0cf5b8fd710af6cc3391f1"


@pytest.mark.parametrize("build", BUILDS)
def test_own_plp_is_decided_on_the_input_build(tmp_path, build):
    data, source = _fixture(tmp_path, build)
    with contextlib.closing(sqlite3.connect(data / "clinvar.GRCh38.sqlite")) as conn, conn:
        conn.execute(
            "UPDATE clinvar_variants SET clinical_significance='Benign', "
            "aggregate='Benign' WHERE rsid=?",
            (IDS[7],),
        )
    report, output = _run(tmp_path, data, source, build)
    site = next(s for s in report["high_value_unchecked"] if s["rsid"] == IDS[7])
    assert site["reason"] == "foreign_alleles"
    flat = _flat(output)
    assert (_flat(site["message"]) in flat) == (build == "GRCh37")
    assert ("1 " if build == "GRCh37" else "2 ") + GROUPS[5] in flat
    assert (
        "3 named here, 6 counted by reason."
        if build == "GRCh37"
        else "2 named here, 7 counted by reason."
    ) in flat


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "index", [0, 1, 4, 5], ids=["own-plp-di", "near-di", "unknown", "foreign"]
)
def test_user_high_value_file_keeps_every_reason_inline_and_html_escapes(
    tmp_path, monkeypatch, build, index
):
    data, source = _fixture(tmp_path, build)
    path = tmp_path / "high-value.yaml"
    path.write_text(f'- rsid: {IDS[index]}\n  gene: "TEST<&>"\n  name: "Name<&>"\n')
    merged = load_high_value_snps([path])
    monkeypatch.setattr(helpers, "load_high_value_snps", lambda: merged)
    report, output = _run(tmp_path, data, source, build)
    site = next(s for s in report["high_value_unchecked"] if s["rsid"] == IDS[index])
    assert _flat(site["message"]) in _flat(output)
    assert "4 named here, 5 counted by reason." in _flat(output)
    raw, _ = _run(tmp_path, data, source, build, surface="html")
    b = _banner(raw)
    assert site["message"] in b.top
    assert b.all == [s["message"] for s in report["high_value_unchecked"]]
    assert "Name&lt;&amp;&gt;" in raw
    assert "Name<&>" not in raw


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("availability", ["off", "absent"])
def test_no_clinvar_claim_from_disabled_or_missing_cache(tmp_path, build, availability):
    data, source = _fixture(tmp_path, build)
    if availability == "off":
        p = data / "config.toml"
        p.write_text(p.read_text().replace("clinvar = true", "clinvar = false"))
    else:
        for p in data.glob("clinvar.*.sqlite"):
            p.unlink()
    report, output = _run(tmp_path, data, source, build)
    assert [s["rsid"] for s in report["high_value_unchecked"]] == ["rs6025"]
    assert "1 high-value or Pathogenic site(s) could not be checked:" in _flat(output)
    assert HINT not in output
    assert not any(sentence in _flat(output) for sentence in GROUPS)


@pytest.mark.parametrize("build", BUILDS)
def test_thousand_vendor_probes_do_not_expand_default_banner(tmp_path, build):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG)
    for b in BUILDS:
        _clinvar(
            data, b, [(IDS[0], "1", 2000000, "G", "A", "SYNTH", "Pathogenic", "Pathogenic", 1, 1)]
        )
    _identity(data, [])
    probes = [(f"i3999949{i:04d}", 2000000, "AG") for i in range(1001)]
    source = _chip(tmp_path, probes, build)
    report, output = _run(tmp_path, data, source, build)
    assert len(report["high_value_unchecked"]) == 1001
    assert report["input"]["probe_classes"]["i_probes_at_plp_loci"] == 1001
    flat = _flat(output)
    assert (
        "1,001 high-value or Pathogenic site(s) could not be checked: 0 named "
        "here, 1,001 counted by reason." in flat
    )
    assert "1,001 " + GROUPS[3] in flat
    assert not any(p[0] in output for p in probes)
    assert flat.count(GROUPS[3]) == 1
    raw, _ = _run(tmp_path, data, source, build, surface="html")
    b = _banner(raw)
    assert b.top == ["1,001 " + GROUPS[3]]
    assert len(b.all) == 1001 and b.open_details == [False]
    assert "All 1,001 sites not checked" in raw


@pytest.mark.parametrize("build", BUILDS)
def test_no_count_html_is_byte_identical_to_released_banner(tmp_path, build):
    data, source = _fixture(tmp_path, build, only_inline=True)
    report, _ = _run(tmp_path, data, source, build)
    raw, _ = _run(tmp_path, data, source, build, surface="html")
    items = "".join(
        f"<li>{html.escape(s['message'])}</li>" for s in report["high_value_unchecked"]
    )
    old = (
        '<div class="notice-warn"><strong>High-value sites not '
        "checked.</strong> These clinically important sites could not be fully "
        "checked; each line says what was not checked:<ul>" + items + "</ul></div>"
    )
    assert old in raw


def test_unplaced_own_plp_stays_inline_without_rollup(tmp_path):
    data = _list_data(tmp_path)
    source = _rsid_list(tmp_path, [(PLP, "AT"), (OTHER, "CT"), ("rs3999980709", "AG")])
    # Lists are GRCh38 only; explicit GRCh38 is accepted by analyze.
    report, output = _run(tmp_path, data, source, "GRCh38")
    sites = report["high_value_unchecked"]
    assert len(sites) == 1 and sites[0]["reason"] == "unplaced"
    assert sites[0]["message"] in output
    assert HINT not in output


@pytest.mark.parametrize("build", BUILDS)
def test_identity_conflict_own_plp_stays_inline(tmp_path, build):
    data, source = _fixture(tmp_path, build, only_inline=True)
    text = source.read_text().replace(f"{IDS[9]}\t1\t2000900\tAG", f"{IDS[9]}\t1\t2000950\tAG")
    source.write_text(text)
    report, output = _run(tmp_path, data, source, build)
    site = next(s for s in report["high_value_unchecked"] if s["rsid"] == IDS[9])
    assert site["reason"] == "identity_locus_conflict"
    assert site["message"] in output
    assert HINT not in output


@pytest.mark.parametrize("build", BUILDS)
def test_absent_reference_carried_high_value_site_stays_inline(tmp_path, build):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG)
    rsid, pos = "rs776746", 2000000
    _identity(data, [(rsid, "7", pos, "C", "T", pos, "C", "T", 0)])
    _clinvar(data, build, [(rsid, "7", pos, "C", "C", "CYP3A5", "Pathogenic", "Pathogenic", 1, 1)])
    source = _vcf(tmp_path / "calls.vcf", [("rs3999974999", "1", 1000, "A", "G", "0/1", ".")])
    report, output = _run(tmp_path, data, source, build)
    (site,) = report["high_value_unchecked"]
    assert site["reason"] == "reference_carried_not_in_file"
    assert site["message"] in output
    assert HINT not in output

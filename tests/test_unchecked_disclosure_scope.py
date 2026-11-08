# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Name counted vendor probes and avoid claiming disabled sources ran.

The new section describes checks actually made, not just possible source
capabilities. Base-valued vendor probes already have a ClinVar P/LP-locus
count and need the same named disclosure as D/I probes. Neither disclosure
may imply a genotype was checked by a source that was switched off or absent.
"""

import pytest
from test_named_unchecked_calls import (
    F5,
    OTHER,
    PLP,
    POS,
    _chip,
    _data,
    _list_data,
    _rsid_list,
    _run,
)


@pytest.mark.parametrize("case", ["unplaced", "foreign"])
@pytest.mark.parametrize("availability", ["disabled", "missing"])
@pytest.mark.parametrize("surface", ["json", "terminal", "html"])
def test_notice_does_not_claim_unavailable_sources_read_the_call(
    tmp_path, case, availability, surface
):
    if case == "unplaced":
        data = _list_data(tmp_path)
        source = _rsid_list(tmp_path, [(F5, "AC"), (OTHER, "CT"), ("rs3999980709", "AG")])
        build = None
    else:
        data = _data(
            tmp_path,
            [(PLP, "1", POS, "G", "C", "GENEP", "Pathogenic", "Pathogenic", 1, 1)],
            [(F5, "1", POS, "G", "A", POS, "G", "A", 0)],
        )
        source = _chip(tmp_path, [(F5, POS, "CC")])
        build = "GRCh37"
    if availability == "missing":
        config = data / "config.toml"
        text = config.read_text()
        for name in ("pharmgkb", "snpedia", "gwas"):
            text = text.replace(f"{name} = false", f"{name} = true")
        config.write_text(text)
    report, _ = _run(tmp_path, data, source, build, surface)
    if surface == "json":
        states = {s["name"]: s for s in report["sources"]}
        for name in ("pharmgkb", "snpedia", "gwas"):
            assert states[name]["state"] == "skipped"
            assert (states[name]["reason"] == "disabled") == (availability == "disabled")
        entries = [e for e in report["high_value_unchecked"] if e["rsid"] == F5]
        assert len(entries) == 1
        text = entries[0]["message"]
    else:
        text = report
        assert F5 in text
    text = " ".join(text.split())
    assert "still did" not in text, text
    assert "other sources read it" not in text, text


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["json", "terminal", "html"])
def test_counted_base_valued_vendor_probe_is_named(tmp_path, build, surface):
    data = _data(
        tmp_path,
        [(PLP, "1", POS, "G", "A", "GENEP", "Pathogenic", "Pathogenic", 1, 1)],
        [(PLP, "1", POS, "G", "A", POS, "G", "A", 0)],
    )
    probe, outside = "i3999980713", "i3999980714"
    source = _chip(tmp_path, [(probe, POS, "AG"), (outside, POS + 90, "AG")], build)
    report, _ = _run(tmp_path, data, source, build, surface, show_unchecked=surface == "terminal")
    if surface == "json":
        assert report["input"]["probe_classes"]["i_probes_at_plp_loci"] == 1
        entries = report.get("high_value_unchecked") or []
        assert {e["rsid"] for e in entries} == {probe}
        (entry,) = entries
        assert entry["gene"] == "GENEP"
        assert entry["call"] == "A/G"
        assert not [a for a in report["annotations"] if a["source"] == "clinvar"]
    else:
        assert probe in report
        assert outside not in report

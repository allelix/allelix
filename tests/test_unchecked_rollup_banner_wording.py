"""The HTML roll-up banner agrees in number with its counts: one site is singular, and
a banner that names no site says so instead of "0 are named here and the rest"."""

from __future__ import annotations

from allelix.reports.html import _hv_unchecked_banner

SITE = "SYNTH (i4000001): called A/G at a vendor probe without an rsID."
COUNT = "1 site(s) at vendor probes without an rsID with a base call."


def test_one_counted_site_is_singular_and_none_is_named():
    banner = _hv_unchecked_banner([SITE], [], [COUNT])
    assert "1 clinically important site could not be fully checked;" in banner
    assert "none is named here and it is counted by reason:" in banner
    assert "<summary>The 1 site not checked</summary>" in banner
    assert "0 are named" not in banner
    assert "the rest" not in banner


def test_several_sites_none_named_are_all_counted():
    banner = _hv_unchecked_banner([SITE, SITE.replace("1)", "2)")], [], ["2 site(s) ..."])
    assert "2 clinically important sites could not be fully checked;" in banner
    assert "none is named here and all are counted by reason:" in banner
    assert "<summary>All 2 sites not checked</summary>" in banner


def test_one_named_and_the_rest_counted():
    sites = [SITE, SITE.replace("1)", "2)"), SITE.replace("1)", "3)")]
    banner = _hv_unchecked_banner(sites, [sites[0]], ["2 site(s) ..."])
    assert "3 clinically important sites could not be fully checked;" in banner
    assert "1 is named here and the rest are counted by reason:" in banner


def test_several_named_keep_the_numeral():
    sites = [SITE.replace("1)", f"{i})") for i in range(1, 5)]
    banner = _hv_unchecked_banner(sites, sites[:2], ["2 site(s) ..."])
    assert "4 clinically important sites could not be fully checked;" in banner
    assert "2 are named here and the rest are counted by reason:" in banner

"""High-value placement naming uses the withheld site, not its raw record IDs."""

import contextlib
import sqlite3

import pytest
from rc3_verification_helpers import OTHER, POS, analyze, caches, warning, write_vcf

from allelix.annotators.clinvar import clinvar_db_filename

HV = "rs4680"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("style", ["dot", "positional", "written"])
@pytest.mark.parametrize("identity_on", [True, False])
def test_withheld_contradiction_is_named_at_the_high_value_placement(
    tmp_path, build, style, identity_on
):
    pairs = [(HV, "A"), (HV, "C")]
    # Keep ClinVar live but give it no P/LP row at the tested site. This proves
    # high-value placement naming independently of the carried-ALT P/LP route.
    data = caches(tmp_path / "data", [(OTHER, "T")], pairs if identity_on else [], [], [])
    for b in POS:
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(b))) as conn, conn:
            conn.execute("UPDATE clinvar_variants SET position=position+100000")
    ids = (
        (HV, HV)
        if style == "written"
        else ((".", ".") if style == "dot" else (f"1:{POS[build]}:G:A", f"1:{POS[build]}:G:C"))
    )
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [dict(id=ids[0], alt="A", gt="1/1"), dict(id=ids[1], alt="C", gt="0/1")],
    )
    report, console = analyze(data, source, tmp_path / "report.json", build)
    entries = report.get("high_value_unchecked", [])
    if identity_on or style == "written":
        assert len(entries) == 1, "The held site's high-value identity must survive resolution"
        assert entries[0]["rsid"] == HV and entries[0]["reason"] == "split_site"
        assert HV in console
        assert "contradict" in entries[0]["message"].lower()
    else:
        assert not entries, "Without an ID or a placement, do not invent a high-value identity"
    assert warning(report, "split_site_calls")[0]["count"] == 2
    assert not report["annotations"]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_another_locus_does_not_borrow_high_value_identity(tmp_path, build):
    data = caches(tmp_path / "data", [], [(HV, "A"), (HV, "C")], [], [])
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [
            dict(id=".", pos=POS[build] + 100, alt="A", gt="1/1"),
            dict(id=".", pos=POS[build] + 100, alt="C", gt="0/1"),
        ],
    )
    report, _ = analyze(data, source, tmp_path / "report.json", build)
    assert not report.get("high_value_unchecked", [])
    assert warning(report, "split_site_calls")[0]["count"] == 2

"""A withheld split site named by a ClinVar P/LP row that has no rsID says why (A-11).

The P/LP row is evidence: with no rsID of its own it is named by its locus,
and its gene stays in the evidence, never as the records' own gene. Records
that name no rsID are listed under the locus with that evidence; records
that name another rsID are listed under it with the same evidence.
"""

import pytest
from split_site_contract_helpers import BUILDS, OTHER, POS, data_dir, record, run, write_vcf


def _split_named(report):
    return [s for s in report.get("high_value_unchecked", []) if s["reason"] == "split_site"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("rsid", (".", OTHER), ids=("unnamed", "other_rsid"))
def test_plp_row_without_rsid_is_named_by_its_locus_with_its_gene(tmp_path, build, rsid):
    # ClinVar's P/LP rows at the site (G>A, G>T, G>C, gene SYNGENE) have no rsID.
    data = data_dir(tmp_path / "data", identity=False, rsid="")
    rows = [record("A", "1/1", rsid=rsid), record("C", "0/1", rsid=rsid)]
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    report, text = run(data, source, tmp_path / "out.json", build, extra=("--show-unchecked",))
    assert not report["annotations"]
    locus = f"1:{POS[build]}"
    named = _split_named(report)
    assert [s["rsid"] for s in named] == [OTHER if rsid == OTHER else locus]
    message = named[0]["message"]
    # The evidence names the row by its locus and keeps its gene.
    assert f"carry: {locus} (SYNGENE)." in message, message
    assert ":  (" not in message and ": (" not in message, message
    # Evidence, not the records' own gene.
    assert named[0]["gene"] == ""
    assert f"{locus} (SYNGENE)" in " ".join(text.split())

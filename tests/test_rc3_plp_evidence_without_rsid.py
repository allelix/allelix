"""A P/LP row with no rsID still has meaningful locus and gene evidence."""

import pytest
from rc3_verification_helpers import FOREIGN, OTHER, POS, analyze, caches, warning, write_vcf


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("record_id", [".", FOREIGN])
@pytest.mark.parametrize("row_id", ["", OTHER])
def test_carried_plp_evidence_is_named_even_without_an_rsid(tmp_path, build, record_id, row_id):
    data = caches(tmp_path / "data", [(row_id, "A"), (row_id, "C")], [], [], [])
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [dict(id=record_id, alt="A", gt="1/1"), dict(id=record_id, alt="C", gt="0/1")],
    )
    report, console = analyze(data, source, tmp_path / "report.json", build)
    entries = report.get("high_value_unchecked", [])
    assert len(entries) == 1
    entry = entries[0]
    assert entry["reason"] == "split_site"
    assert entry["rsid"], "The named site must never have a blank label"
    text = " ".join(entry["message"].split())
    assert "VERIFYGENE" in text, "The row's gene must survive even when it has no rsID"
    if row_id:
        assert row_id in text
    else:
        assert f"1:{POS[build]}" in text or any(str(n) in text for n in (88101, 88102)), (
            "Evidence must identify the row by locus or ClinVar allele ID, not an empty rsID"
        )
    assert "carry: (" not in text and "carry:  (" not in entry["message"]
    assert "VERIFYGENE" in console, "This named site must also be visible in the default console"
    assert not report["annotations"]
    assert warning(report, "split_site_calls")[0]["count"] == 2


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_an_uncarried_plp_alt_does_not_lend_its_gene_or_a_named_site(tmp_path, build):
    data = caches(tmp_path / "data", [("", "T")], [], [], [])
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [dict(id=FOREIGN, alt="A", gt="1/1"), dict(id=FOREIGN, alt="C", gt="0/1")],
    )
    report, console = analyze(data, source, tmp_path / "report.json", build)
    assert not report.get("high_value_unchecked", [])
    assert "VERIFYGENE" not in console
    assert warning(report, "split_site_calls")[0]["count"] == 2

"""Contradicting passing records are withheld, counted, and visibly named."""

import pytest
from split_site_contract_helpers import BUILDS, RSID, by_source, data_dir, record, run, write_vcf

CONTRADICTIONS = (
    (record("A", "1/1"), record("T", "0/1")),
    (record("A", "0/1"), record("T", "1/1")),
    (record("A", "1/0"), record("T", "1/1")),
    (record("A", "1"), record("T", "0/1")),
    (record("A", "1/0"), record("T", "0/1"), record("C", "0/1")),
    (record("A", "0|1", ps="41"), record("T", "0|1", ps="41")),
    (record("A", "1|0", ps="41"), record("T", "1|0", ps="41")),
    (record("A", "1|.", ps="41"), record("T", "1|.", ps="41")),
    (record("A", "1/."), record("T", "1/1")),
)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("rows", CONTRADICTIONS)
def test_contradicting_site_is_counted_and_named_at_a_pathogenic_rsid(
    tmp_path, build, reverse, rows
):
    # Keep GWAS A unambiguous so an accidental half-call leak would be observable.
    data = data_dir(tmp_path / "data", alts="AC")
    rows = [dict(r, alt={"T": "C", "C": "T"}.get(r["alt"], r["alt"])) for r in rows]
    source = write_vcf(tmp_path / "input.vcf", build, rows[::-1] if reverse else rows)
    report, text = run(data, source, tmp_path / "out.json", build, extra=("--show-unchecked",))
    assert all(not found for found in by_source(report).values())
    warning = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    assert len(warning) == 1 and warning[0]["count"] == len(rows)
    named = report.get("high_value_unchecked", [])
    assert any(r.get("rsid") == RSID for r in named), "A count cannot identify the withheld site"
    assert RSID in text and "SYNGENE" in text, (
        "The named line must be visible, not only serialized"
    )
    assert "ClinVar and the GWAS Catalog assessed the called allele" not in text


@pytest.mark.parametrize("build", BUILDS)
def test_high_value_rsid_is_named_even_without_a_pathogenic_cache_row(tmp_path, build):
    rid = "rs6025"
    data = data_dir(tmp_path / "data", rsid=rid, pathogenic=False, alts="AT")
    rows = [record("A", "1/1", rsid=rid), record("T", "0/1", rsid=rid)]
    report, text = run(
        data,
        write_vcf(tmp_path / "input.vcf", build, rows),
        tmp_path / "out.json",
        build,
        extra=("--show-unchecked",),
    )
    assert not report["annotations"]
    assert any(r.get("rsid") == rid for r in report.get("high_value_unchecked", []))
    assert rid in text


@pytest.mark.parametrize("build", BUILDS)
def test_contradiction_is_named_in_html_as_well_as_json(tmp_path, build):
    data = data_dir(tmp_path / "data", alts="AT")
    source = write_vcf(tmp_path / "input.vcf", build, CONTRADICTIONS[0])
    report, _ = run(data, source, tmp_path / "out.json", build)
    # Render the public result through the CLI as the owner receives it.
    from click.testing import CliRunner

    from allelix.cli import main

    output = tmp_path / "out.html"
    res = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            build.lower(),
            "--output",
            str(output),
        ],
    )
    assert res.exit_code == 0, (res.output, res.exception)
    html = output.read_text()
    assert any(r.get("rsid") == RSID for r in report.get("high_value_unchecked", []))
    # The rsID occurs in the rendered disclosure, not merely in embedded findings (empty).
    import html as html_module
    import re

    visible = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", html, flags=re.S | re.I)
    visible = html_module.unescape(re.sub(r"<[^>]*>", " ", visible))
    assert RSID in visible and "SYNGENE" in visible

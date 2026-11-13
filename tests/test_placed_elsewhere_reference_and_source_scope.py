"""Reference-carried exact matches are disclosed; other sources are not implied."""

import contextlib
import json
import sqlite3
from dataclasses import replace

import pytest
from click.testing import CliRunner
from test_clinvar_own_placement_disclosure import (
    BUILDS,
    HIGH,
    ORDINARY,
    POS,
    TARGET,
    _data,
    _run,
    _variant,
    _vcf,
)

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.cli import main
from allelix.models import NO_CALL_MARKER


def _reference_data(data, build, rid, *, own_plp=False):
    _data(data, build, rid, own_plp=own_plp)
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn, conn:
        conn.execute("UPDATE clinvar_variants SET alt=ref WHERE rsid=?", (TARGET,))


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", (False, True))
@pytest.mark.parametrize(
    "alleles,carried",
    [(("G", "A"), True), (("G", "G"), True), (("G", NO_CALL_MARKER), True), (("A", "A"), False)],
    ids=["heterozygous", "hom-reference", "half-reference", "reference-not-carried"],
)
def test_reference_carried_row_is_an_exact_match_for_the_placement_count(
    tmp_path, build, batch, alleles, carried
):
    data = tmp_path / "data"
    _reference_data(data, build, ORDINARY)
    call = replace(_variant(build, ORDINARY), allele1=alleles[0], allele2=alleles[1])

    def run(annotator, variant):
        return list(annotator.batch_annotate([variant])) if batch else annotator.annotate(variant)

    with ClinVarAnnotator(data, builds=(build,)) as control:
        rows = run(control, replace(call, rsid=TARGET))
        assert bool(rows) is carried, "Establish that the reference base is the matched ALT"
        assert all(r.alt == "G" and r.reference_carries_alt for r in rows)
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        assert not run(annotator, call), "The existing placement guard must stay"
        assert annotator.run_counts()["placed_elsewhere_calls"] == int(carried)
        assert ((ORDINARY, "placed_elsewhere") in annotator.unchecked()) is carried


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("kind", ("high_value", "own_plp", "ordinary"))
def test_cli_reference_carried_exact_match_is_counted_and_named_where_important(
    tmp_path, build, explicit, kind
):
    rid = HIGH if kind == "high_value" else ORDINARY
    data = tmp_path / "data"
    _reference_data(data, build, rid, own_plp=kind == "own_plp")
    path = _vcf(tmp_path / "input.vcf", [(rid, "1", POS, "G", "A", "0/1", ".")])
    length = 249250621 if build == "GRCh37" else 248956422
    path.write_text(
        path.read_text().replace(
            "##fileformat=VCFv4.2\n",
            f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n",
        )
    )
    output = tmp_path / "report.json"
    args = ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(output)]
    if explicit:
        args += ["--build", build.lower()]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    report = json.loads(output.read_text())
    assert report["input"]["build"] == build
    assert not report["annotations"]
    warnings = [
        w for w in report["report_warnings"] if w["code"] == "clinvar_rsid_placed_elsewhere"
    ]
    assert len(warnings) == 1 and warnings[0]["count"] == 1
    named = report.get("high_value_unchecked", [])
    assert [r["rsid"] for r in named] == ([] if kind == "ordinary" else [rid])
    if kind == "own_plp":
        assert named[0]["gene"] == "OWNGENE"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("surface", ("warning", "named", "terminal", "html"))
def test_placement_disclosure_does_not_claim_disabled_sources_read_the_call(
    tmp_path, build, surface
):
    if surface in ("warning", "named"):
        report = _run(tmp_path, build, True, HIGH)
        others = [s for s in report["sources"] if s["name"] in {"snpedia", "pharmgkb", "gwas"}]
        assert len(others) == 3 and all(s["state"] == "skipped" for s in others)
        assert not report["annotations"]
        if surface == "warning":
            text = next(
                w["message"]
                for w in report["report_warnings"]
                if w["code"] == "clinvar_rsid_placed_elsewhere"
            )
        else:
            text = report["high_value_unchecked"][0]["message"]
    else:
        text = _run(tmp_path, build, True, HIGH, surface=surface)
    text = " ".join(text.split())
    assert "ClinVar" in text and "elsewhere" in text
    assert "still read the call" not in text, (
        "Disclose ClinVar's decision, not that disabled sources ran"
    )

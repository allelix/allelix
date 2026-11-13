"""ClinVar's own placed-elsewhere guard stays, with a counted disclosure."""

import html
import json
import re

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity, _pharm, _vcf

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.cli import main
from allelix.models import Variant

BUILDS = ("GRCh37", "GRCh38")
TARGET = "rs900451001"
ORDINARY = "rs900451002"
HIGH = "rs6025"
POS = 2000001


def _data(path, build, rid, *, own_plp=False, extra_here=False, pharm=False):
    path.mkdir()
    sig = "Pathogenic" if own_plp else "Benign"
    rows = [
        (TARGET, "1", POS, "G", "A", "TARGETGENE", "Pathogenic", "Pathogenic", 1, 1),
        (rid, "2", POS + 100, "C", "T", "OWNGENE", sig, sig, 1, 1),
    ]
    if extra_here:
        rows += [(rid, "1", POS, "G", "A", "OWNGENE", "Pathogenic", "Pathogenic", 1, 1)]
    _clinvar(path, build, rows)
    if pharm:
        _pharm(path, rid, "OWNGENE", ("AG",))
    (path / "config.toml").write_text(
        "[sources]\nclinvar=true\nsnpedia=false\ngwas=false\n"
        f"pharmgkb={str(pharm).lower()}\ngnomad=false\n"
        "alphamissense=false\ncadd=false\nrsid_identity=false\n"
    )


def _variant(build, rid, *, padded=False):
    return Variant(
        rsid=rid,
        chromosome="1",
        position=POS - int(padded),
        allele1="CG" if padded else "G",
        allele2="CA" if padded else "A",
        ref="CG" if padded else "G",
        ref_from_input=True,
        build=build,
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", (False, True))
@pytest.mark.parametrize("padded", (False, True))
def test_single_and_batch_keep_the_guard_and_count_it(tmp_path, build, batch, padded):
    data = tmp_path / "data"
    _data(data, build, ORDINARY)
    a = ClinVarAnnotator(data, builds=(build,))
    try:
        before = a.run_counts()
        wrong = _variant(build, ORDINARY, padded=padded)
        rows = list(a.batch_annotate([wrong])) if batch else a.annotate(wrong)
        assert not rows, "Do not revive the exact allele under a contradictory rsID"
        after = a.run_counts()
        assert sum(after.get(k, 0) - before.get(k, 0) for k in after) == 1, (
            "The unchanged guard must increment a run count"
        )
    finally:
        a.close()


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("kind", ("correct", "unnamed", "alias_unheld", "one_valid_placement"))
def test_nonconflicting_exact_matches_are_unchanged(tmp_path, build, kind):
    data = tmp_path / "data"
    _data(data, build, ORDINARY, extra_here=kind == "one_valid_placement")
    rid = {
        "correct": TARGET,
        "unnamed": "",
        "alias_unheld": "rs900451099",
        "one_valid_placement": ORDINARY,
    }[kind]
    a = ClinVarAnnotator(data, builds=(build,))
    try:
        rows = a.annotate(_variant(build, rid))
        assert rows and all(r.significance == "clinvar_pathogenic" for r in rows)
        assert not any(a.run_counts().values())
    finally:
        a.close()


def _run(
    tmp_path,
    build,
    explicit,
    rid,
    *,
    own_plp=False,
    surface="json",
    pharm=False,
    identity_present=False,
):
    data = tmp_path / "data"
    _data(data, build, rid, own_plp=own_plp, pharm=pharm)
    if identity_present:
        _identity(data, [(rid, "1", POS, "G", "A", POS, "G", "A", 0)])
    path = _vcf(tmp_path / "input.vcf", [(rid, "1", POS, "G", "A", "0/1", ".")])
    length = 249250621 if build == "GRCh37" else 248956422
    path.write_text(
        path.read_text().replace(
            "##fileformat=VCFv4.2\n", f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
        )
    )
    out = tmp_path / f"report.{surface}"
    args = ["analyze", str(path), "--data-dir", str(data), "--no-update", "--min-magnitude", "0"]
    if explicit:
        args += ["--build", build.lower()]
    if surface != "terminal":
        args += ["--output", str(out)]
    result = CliRunner(env={"COLUMNS": "200"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    if surface == "json":
        return json.loads(out.read_text())
    if surface == "terminal":
        return result.output
    text = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", out.read_text(), flags=re.S)
    return html.unescape(re.sub(r"<[^>]+>", " ", text))


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("kind", ("high_value", "own_plp", "ordinary"))
def test_report_counts_the_clinvar_abstention_and_names_important_calls(
    tmp_path, build, explicit, kind
):
    rid = HIGH if kind == "high_value" else ORDINARY
    report = _run(tmp_path, build, explicit, rid, own_plp=kind == "own_plp")
    assert report["input"]["build"] == build
    assert not report["annotations"]
    warnings = [w for w in report["report_warnings"] if w.get("count")]
    assert len(warnings) == 1 and warnings[0]["count"] == 1
    text = warnings[0]["message"]
    assert "ClinVar" in text and "elsewhere" in text.lower()
    assert "identity table" not in text.lower(), "The disabled table supplied no evidence"
    named = report.get("high_value_unchecked", [])
    if kind == "ordinary":
        assert not named
    else:
        assert [r["rsid"] for r in named] == [rid]
        assert "ClinVar" in named[0]["message"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("surface", ("html", "terminal"))
def test_high_value_disclosure_is_visible_not_just_json(tmp_path, build, surface):
    text = _run(tmp_path, build, True, HIGH, surface=surface)
    assert HIGH in text and "ClinVar" in text and "elsewhere" in text.lower()


@pytest.mark.parametrize("build", BUILDS)
def test_clinvar_placement_guard_does_not_silence_other_rsid_sources(tmp_path, build):
    report = _run(tmp_path, build, True, HIGH, pharm=True)
    assert [r["source"] for r in report["annotations"]] == ["pharmgkb"]
    assert any(
        w.get("count") == 1 and "ClinVar" in w["message"] for w in report["report_warnings"]
    )


@pytest.mark.parametrize("build", BUILDS)
def test_an_installed_but_disabled_identity_cache_does_not_hide_the_count(tmp_path, build):
    report = _run(tmp_path, build, True, HIGH, identity_present=True)
    assert not report["annotations"]
    assert any(
        w.get("count") == 1 and "ClinVar" in w["message"] for w in report["report_warnings"]
    )
    assert [r["rsid"] for r in report.get("high_value_unchecked", [])] == [HIGH]

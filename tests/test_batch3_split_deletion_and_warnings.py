"""Cross-combinations for deletion watching, split calls and both disclosures."""

import html
import json
import re

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _vcf
from test_clinvar_own_placement_disclosure import HIGH, POS, _data
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect

from allelix.cli import main
from allelix.parsers.vcf import VcfParser


def _row(pos, alt, gt, *, chrom="1", rid="rs900450900", info="."):
    return rid, chrom, pos, "G", alt, gt, info


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sources", (("clinvar",), ("gwas", "snpedia", "pharmgkb")))
@pytest.mark.parametrize("where", ("before", "after", "between_components"))
@pytest.mark.parametrize("carried", (False, True))
def test_split_buffer_preserves_both_reference_blocking_paths(
    tmp_path, build, sources, where, carried
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    pair = [_row(pos - 5, "A", "1/0"), _row(pos - 5, "T", "0/1")]
    # Between the components, a symbolic record at the SAME position is not
    # a single-base component. It still contributes its downstream END span.
    deletion = _row(
        pos - (5 if where == "between_components" else 10),
        "<DEL>",
        "0/1" if carried else "0/0",
        rid=".",
        info=f"END={pos + 10}",
    )
    rows = {
        "before": [deletion, *pair],
        "after": [*pair, deletion],
        "between_components": [pair[0], deletion, pair[1]],
    }[where]
    rows.append(_row(pos, "A", "0/0", rid=RSID))
    source = _vcf(tmp_path / "input.vcf", rows)
    result = analyze(data, source, build, sources)
    assert {a.source for a in result.annotations} == (set() if carried else set(sources))
    assert not any(w.code == "split_site_calls" for w in result.report_warnings)
    parser = VcfParser()
    parser.combine_split_sites = True
    calls = list(parser.parse(source))
    combined = [v for v in calls if v.position == pos - 5 and v.rsid == "rs900450900"]
    assert len(combined) == 1
    assert set(combined[0].carried_pair) == {"A", "T"}
    assert parser.split_site_calls == 0


@pytest.mark.parametrize("symbol", ("<DEL>", "<DEL:ME>"))
@pytest.mark.parametrize("reverse", (False, True))
def test_a_symbolic_deletion_is_not_a_single_base_split_component(tmp_path, symbol, reverse):
    rows = [_row(100, "A", "1/0"), _row(100, symbol, "0/1", info="END=200")]
    source = _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows)
    parser = VcfParser()
    parser.combine_split_sites = True
    parser.watch_reference_targets({"1": [150]})
    calls = list(parser.parse(source))
    assert [(v.position, set(v.carried_pair)) for v in calls] == [(100, {"G", "A"})]
    assert parser.split_site_calls == 0
    assert parser._deleted_spans == {"1": [(101, 200)]}


@pytest.mark.parametrize("chrom", ("M", "chrM"))
@pytest.mark.parametrize("rcrs", (False, True))
def test_mito_exclusion_still_applies_with_split_buffer_and_deletion_watch(tmp_path, chrom, rcrs):
    rows = [
        _row(100, "A", "1/0", chrom=chrom),
        _row(100, "<DEL>", "0/1", chrom=chrom, info="END=200"),
        _row(100, "T", "0/1", chrom=chrom),
        _row(150, "A", "0/0", chrom=chrom),
    ]
    source = _vcf(tmp_path / "mito.vcf", rows)
    outcomes = []
    for combine in (False, True):
        parser = VcfParser()
        parser.combine_split_sites = combine
        parser.drop_ucsc_mito = True
        parser.rcrs_mito_names = frozenset({chrom}) if rcrs else frozenset()
        parser.watch({("MT", 150)}, set(), rows=True)
        parser.watch_reference_targets({"MT": [150]})
        calls = list(parser.parse(source))
        outcomes.append((parser, calls))
    old, plain = outcomes[0]
    parser, combined = outcomes[1]
    assert parser.ucsc_mito_dropped == old.ucsc_mito_dropped
    assert parser.quality_stats.locus_states == old.quality_stats.locus_states
    assert parser.reference_target_calls == old.reference_target_calls
    assert parser.reference_target_blocked == old.reference_target_blocked
    assert parser.split_site_calls == 0
    if rcrs:
        assert len(plain) == 2 and len(combined) == 1
        assert set(combined[0].carried_pair) == {"A", "T"}
        assert ("MT", 150) in parser.reference_target_blocked
    else:
        assert not plain and not combined
        assert not parser._deleted_spans


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("surface", ("json", "html"))
@pytest.mark.parametrize("shape", ("both", "split_only", "split_and_del", "elsewhere_only"))
def test_both_run_warnings_reach_json_and_html_independently(
    tmp_path, build, explicit, surface, shape
):
    data = tmp_path / "data"
    _data(data, build, HIGH)
    rows = []
    if shape != "elsewhere_only":
        rows += [_row(POS - 10, "A", "0/1"), _row(POS - 10, "T", "0/1")]
    if shape in {"both", "elsewhere_only"}:
        rows.append(_row(POS, "A", "0/1", rid=HIGH))
    if shape == "split_and_del":
        rows.insert(0, _row(POS - 20, "<DEL>", "0/1", info=f"END={POS + 20}"))
    path = _vcf(tmp_path / "input.vcf", rows)
    length = 249250621 if build == "GRCh37" else 248956422
    path.write_text(
        path.read_text().replace(
            "##fileformat=VCFv4.2\n",
            f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n",
        )
    )
    output = tmp_path / f"report.{surface}"
    args = ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(output)]
    if explicit:
        args += ["--build", build.lower()]
    run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    if surface == "json":
        report = json.loads(output.read_text())
        warnings = {w["code"]: w["count"] for w in report["report_warnings"]}
        assert warnings.get("split_site_calls", 0) == (0 if shape == "elsewhere_only" else 2)
        assert warnings.get("clinvar_rsid_placed_elsewhere", 0) == (
            1 if shape in {"both", "elsewhere_only"} else 0
        )
        assert not report["annotations"]
        assert report["input"]["build"] == build
    else:
        text = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", output.read_text(), flags=re.S)
        text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", text)).split())
        assert ("2 call(s) were not assessed" in text) is (shape != "elsewhere_only")
        assert ("1 call(s) carry the exact allele" in text) is (
            shape in {"both", "elsewhere_only"}
        )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ("empty", "structural_only", "malformed"))
def test_truly_unusable_input_remains_refused(tmp_path, build, shape):
    data = tmp_path / "data"
    _data(data, build, HIGH)
    rows = (
        [_row(POS, "<DEL>", "0/1", info=f"END={POS + 20}")] if shape == "structural_only" else []
    )
    source = _vcf(tmp_path / "input.vcf", rows)
    if shape == "malformed":
        source.write_text(source.read_text() + "broken\n")
    output = tmp_path / "report.json"
    run = CliRunner().invoke(
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
    assert run.exit_code == 1, (run.output, run.exception)
    assert not output.exists()
    assert "Error:" in run.output

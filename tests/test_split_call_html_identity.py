"""Combined calls share one display row without pooling unrelated records."""

import json
import re
from dataclasses import replace

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _pharm
from test_html_site_identity import ann, page, pair

from allelix.cli import main

RS = "rs900529201"
BUILDS = ("GRCh37", "GRCh38")
LENGTHS = {"GRCh37": 249250621, "GRCh38": 248956422}


def _data(path, build):
    path.mkdir()
    _clinvar(
        path,
        build,
        [
            (RS, "1", 1000, "G", "A", "GENEA", "Pathogenic", "Pathogenic", 1, 1),
            (RS, "1", 1000, "G", "T", "GENET", "Likely_pathogenic", "Likely_pathogenic", 1, 1),
        ],
    )
    _pharm(path, RS, "PGX", ["AT"])
    (path / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=true\nsnpedia=false\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    return path


def _source(path, build, phase, reverse, merged=False):
    rows = [("A", "1/0"), ("T", "0/1")]
    if phase != "unphased":
        rows = [("A", "1|0"), ("T", "0|1")]
    if reverse:
        rows.reverse()
    if merged:
        rows = [("A,T", "1/2")]
    fmt = "GT:PS" if phase == "explicit_ps" else "GT"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=1,length={LENGTHS[build]}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(
            f"1\t1000\t{RS}\tG\t{alt}\t50\tPASS\t.\t{fmt}\t{gt}"
            + (":77" if phase == "explicit_ps" else "")
            + "\n"
            for alt, gt in rows
        )
    )
    return path


def _run(data, source, build, out, vcf_out=None):
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build.lower(),
        "--min-magnitude",
        "0",
        "--output",
        str(out),
    ]
    if vcf_out:
        args += ["--vcf-out", str(vcf_out)]
    run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    return out.read_text()


def _entries(text):
    match = re.search(
        r'<script id="variant-data" type="application/json">(.*?)</script>', text, re.S
    )
    assert match and "<" not in match[1]
    return json.loads(match[1])


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("phase", ["unphased", "default_ps", "explicit_ps"])
@pytest.mark.parametrize("reverse", [False, True])
def test_one_split_person_call_has_one_html_row_with_all_annotations(
    tmp_path, build, phase, reverse
):
    data = _data(tmp_path / "data", build)
    split = _source(tmp_path / "split.vcf", build, phase, reverse)
    merged = _source(tmp_path / "merged.vcf", build, phase, reverse, merged=True)
    expected = _entries(_run(data, merged, build, tmp_path / "merged.html"))
    text = _run(data, split, build, tmp_path / "split.html")
    actual = _entries(text)
    assert len(expected) == 1 and len(expected[0]["annotations"]) == 3
    assert len(actual) == 1, "Per-allele writer records must not split one person's display row"
    assert actual[0]["genotype"] == expected[0]["genotype"] == "AT"
    assert sorted(json.dumps(a, sort_keys=True) for a in actual[0]["annotations"]) == sorted(
        json.dumps(a, sort_keys=True) for a in expected[0]["annotations"]
    )
    assert len(re.findall(r"<tr\b[^>]*data-row-id=", text)) == 1
    assert actual[0].get("locus", "") == expected[0].get("locus", "")
    assert actual[0]["allele_metrics"] == expected[0]["allele_metrics"]
    for label in ("Variants", "Bad"):
        card = rf'<div class="label">{label}</div>\s*<div class="value">1</div>'
        assert re.search(card, text), "Summary cards must count the combined call once"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", [False, True])
def test_json_parity_and_each_alt_writer_origin_survive_display_fix(tmp_path, build, reverse):
    data = _data(tmp_path / "data", build)
    split = _source(tmp_path / "split.vcf", build, "unphased", reverse)
    merged = _source(tmp_path / "merged.vcf", build, "unphased", reverse, merged=True)
    output = tmp_path / "annotated.vcf"
    actual = json.loads(_run(data, split, build, tmp_path / "split.json", output))
    expected = json.loads(_run(data, merged, build, tmp_path / "merged.json"))
    assert sorted(json.dumps(a, sort_keys=True) for a in actual["annotations"]) == sorted(
        json.dumps(a, sort_keys=True) for a in expected["annotations"]
    )
    rows = [s.split("\t") for s in output.read_text().splitlines() if not s.startswith("#")]
    assert len(rows) == 2
    by_alt = {r[4]: r[7] for r in rows}
    assert "ALLELIX_CLINVAR=" in by_alt["A"] and "ALLELIX_CLINVAR=" in by_alt["T"]
    assert "likely_pathogenic" not in by_alt["A"]
    assert "likely_pathogenic" in by_alt["T"]


@pytest.mark.parametrize("mode", ["vcf", "chip", "same-position-records", "other-chromosome"])
def test_unrelated_same_label_calls_still_have_separate_rows(tmp_path, mode):
    first, second = pair(mode, genotype="AT")
    _, entries, _ = page(tmp_path, [first, second])
    assert len(entries) == 2
    assert {e["gene"] for e in entries} == {"FIRST", "SECOND"}


def test_distinct_same_position_refs_do_not_pool_allele_metrics(tmp_path):
    first = ann(genotype_match="AT", record_key=("1", 1000, "G", "A,T"))
    second = replace(
        first, record_index=1, locus_ref="C", record_key=("1", 1000, "C", "A,T"), cadd_phred=40.0
    )
    _, entries, _ = page(tmp_path, [first, second])
    assert len(entries) == 2
    assert sorted(e["cadd_phred"] for e in entries) == [1.0, 40.0]

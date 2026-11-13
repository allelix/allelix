"""Split-call abstention is not an empty or structural-only input."""

import html
import json
import re

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _vcf

from allelix.cli import main

BUILDS = ("GRCh37", "GRCh38")
POS = 2000001


def _data(path, build):
    # Ready synthetic cache; its row is away from the split site, so findings
    # cannot decide whether the input is considered to contain records.
    _clinvar(
        path,
        build,
        [
            ("rs900450901", "2", POS + 100, "G", "A", "SYNGENE", "Pathogenic", "Pathogenic", 1, 1),
        ],
    )
    (path / "config.toml").write_text(
        "[sources]\nclinvar=true\nsnpedia=false\ngwas=false\npharmgkb=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )


def _row(pos, alt, gt, *, info="."):
    return "rs900450900", "1", pos, "G", alt, gt, info


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("surface", ("json", "html"))
@pytest.mark.parametrize("shape", ("split_only", "split_and_del"))
def test_all_split_calls_abstaining_still_publish_the_report(
    tmp_path, build, explicit, surface, shape
):
    data = tmp_path / "data"
    _data(data, build)
    rows = [_row(POS - 10, "A", "0/1"), _row(POS - 10, "T", "0/1")]
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
        warnings = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
        assert len(warnings) == 1 and warnings[0]["count"] == 2
        assert not report["annotations"]
        assert report["input"]["build"] == build
    else:
        text = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", output.read_text(), flags=re.S)
        text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", text)).split())
        assert "2 call(s) were not assessed" in text


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ("empty", "structural_only", "malformed"))
def test_truly_unusable_input_remains_refused(tmp_path, build, shape):
    data = tmp_path / "data"
    _data(data, build)
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

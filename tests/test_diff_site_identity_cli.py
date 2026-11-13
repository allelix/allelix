"""Public --diff notices/JSON must not pair the same rsID across two loci."""

import contextlib
import json
import re
import sqlite3

import pytest
from click.testing import CliRunner
from test_homref_ref_by_locus import _cv, _vcf, _write_clinvar

from allelix.annotators.clinvar import clinvar_db_filename
from allelix.cli import _helpers, main
from allelix.reports.diff import load_previous_report

RS = "rs900527101"


def invoke(source, data, output, baseline=None):
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--min-magnitude",
        "0",
        "--include-benign",
    ]
    if baseline is not None:
        args += ["--diff", str(baseline)]
    if output is not None:
        args += ["--output", str(output)]
    command = CliRunner().invoke(main, args)
    assert command.exit_code == 0, (command.output, command.exception)
    return output.read_text() if output is not None else command.output


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
@pytest.mark.parametrize("scenario", ["swapped", "moved"])
@pytest.mark.parametrize("reverse", [False, True])
def test_diff_from_a_real_json_baseline_respects_loci(
    tmp_path, monkeypatch, build, surface, scenario, reverse
):
    # Restore the stored width, never write the computed Console.width property.
    monkeypatch.setattr(_helpers.console, "_width", 240)
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=false\nsnpedia=false\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    _write_clinvar(
        data,
        build,
        [
            _cv(RS, "1", 1000, "G", "A", allele_id=1),
            _cv(
                RS,
                "1",
                5000,
                "G",
                "A",
                allele_id=2,
                sig="Benign" if scenario == "swapped" else "Pathogenic",
                agg="Benign" if scenario == "swapped" else "Pathogenic",
            ),
        ],
    )

    def source(name, positions):
        if reverse:
            positions = list(reversed(positions))
        return _vcf(
            tmp_path / name,
            [("1", pos, RS, "G", "A", 50, "PASS", ".", "GT:GQ", "0/1:99") for pos in positions],
            build=build,
        )

    baseline = tmp_path / "baseline.json"
    invoke(source("before.vcf", [1000, 5000] if scenario == "swapped" else [5000]), data, baseline)
    previous = load_previous_report(baseline)["annotations"]
    assert sorted(a["position"] for a in previous) == (
        [1000, 5000] if scenario == "swapped" else [5000]
    )
    if scenario == "swapped":
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn, conn:
            conn.execute(
                "UPDATE clinvar_variants SET clinical_significance = "
                "CASE position WHEN 1000 THEN 'Benign' ELSE 'Pathogenic' END, "
                "aggregate = CASE position WHEN 1000 THEN 'Benign' ELSE 'Pathogenic' END"
            )
    current = source("after.vcf", [1000, 5000] if scenario == "swapped" else [1000])
    output = None if surface == "terminal" else tmp_path / f"current.{surface}"
    text = invoke(current, data, output, baseline)
    if surface == "json":
        payload = json.loads(text)
        diff = payload["diff"]
        if scenario == "swapped":
            assert len(diff["changed"]) == 2, "CLI diff lost the two per-locus changes"
            assert {a["position"] for a in diff["changed"]} == {1000, 5000}
            assert not diff["new"] and not diff["removed"]
        else:
            assert [a["position"] for a in diff["new"]] == [1000]
            assert [a["position"] for a in diff["removed"]] == [5000]
            assert not diff["changed"]
    else:
        visible = re.sub(r"<(script|style)\b[^>]*>.*?</\1>", "", text, flags=re.S | re.I)
        normalized = (
            " ".join(re.sub(r"<[^>]*>", " ", visible).split())
            if surface == "html"
            else " ".join(text.split())
        )
        if scenario == "swapped":
            assert "2 changed" in normalized, "visible diff must announce both changes"
        else:
            assert "1 new" in normalized and "1 removed" in normalized

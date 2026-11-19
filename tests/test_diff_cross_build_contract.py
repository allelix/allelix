"""Different builds do not turn the same named alleles into removals and additions."""

import contextlib
import json
import re
import sqlite3

import pytest
from click.testing import CliRunner
from test_homref_ref_by_locus import _cv, _vcf, _write_clinvar

from allelix.annotators.clinvar import clinvar_db_filename
from allelix.cli import _helpers, main

RS = "rs900615101"
OTHER = "rs900615102"
POSITIONS = {"GRCh37": 200100, "GRCh38": 300100}


def _setup(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=false\nsnpedia=false\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    for build, pos in POSITIONS.items():
        _write_clinvar(
            data,
            build,
            [
                _cv(RS, "1", pos, "G", "A"),
                _cv(RS, "1", pos + 100, "G", "A", allele_id=2),
                _cv(OTHER, "1", pos + 200, "G", "T", allele_id=3),
            ],
        )
    return data


def _source(tmp_path, folder, build, offsets=(0,)):
    directory = tmp_path / folder
    directory.mkdir(exist_ok=True)
    return _vcf(
        directory / "person.vcf",
        [
            (
                "1",
                POSITIONS[build] + off,
                OTHER if off == 200 else RS,
                "G",
                "T" if off == 200 else "A",
                50,
                "PASS",
                ".",
                "GT:GQ",
                "0/1:99",
            )
            for off in offsets
        ],
        build=build,
    )


def _run(source, data, output, build, baseline=None, auto=False):
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        "auto" if auto else build,
        "--min-magnitude",
        "0",
        "--include-benign",
    ]
    if output is not None:
        args += ["--output", str(output)]
    if baseline is not None:
        args += ["--diff", str(baseline)]
    result = CliRunner().invoke(main, args, env={"COLUMNS": "500"})
    assert result.exit_code == 0, (result.output, result.exception)
    return result.output


def _baseline(tmp_path, data, build, offsets=(0,)):
    baseline = tmp_path / "baseline.json"
    _run(_source(tmp_path, "old", build, offsets), data, baseline, build)
    payload = json.loads(baseline.read_text())
    assert len(payload["annotations"]) == len(offsets)
    return baseline


def _notice(text):
    # No new wording/code or report schema is prescribed. One visible notice
    # must identify both builds and explain that they differ.
    lines = [
        line
        for line in text.splitlines()
        if "GRCh37" in line
        and "GRCh38" in line
        and re.search(r"differ|another|cross.build|mismatch", line, re.I)
    ]
    assert len(lines) == 1, text


@pytest.mark.parametrize("old,new", [("GRCh37", "GRCh38"), ("GRCh38", "GRCh37")])
@pytest.mark.parametrize("auto", [False, True])
@pytest.mark.parametrize("schema", ["13", "14"])
def test_same_named_calls_cross_build_have_no_diff(tmp_path, old, new, auto, schema):
    data = _setup(tmp_path)
    baseline = _baseline(tmp_path, data, old, (0, 200))
    previous = json.loads(baseline.read_text())
    previous["schema_version"] = schema
    baseline.write_text(json.dumps(previous))
    out = tmp_path / "current.json"
    _run(_source(tmp_path, "new", new, (0, 200)), data, out, new, baseline, auto)
    payload = json.loads(out.read_text())
    assert len(payload["annotations"]) == 2
    assert all(not payload["diff"][kind] for kind in ("new", "removed", "changed")), payload[
        "diff"
    ]


@pytest.mark.parametrize("old,new", [("GRCh37", "GRCh38"), ("GRCh38", "GRCh37")])
@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_cross_build_comparison_emits_one_build_notice(tmp_path, monkeypatch, old, new, surface):
    monkeypatch.setattr(_helpers.console, "_width", 500)
    data = _setup(tmp_path)
    baseline = _baseline(tmp_path, data, old)
    out = None if surface == "terminal" else tmp_path / f"current.{surface}"
    text = _run(_source(tmp_path, "new", new), data, out, new, baseline)
    _notice(text)


@pytest.mark.parametrize("old,new", [("GRCh37", "GRCh38"), ("GRCh38", "GRCh37")])
@pytest.mark.parametrize("scenario", ["reclassified", "added", "removed", "multiplicity"])
def test_cross_build_pairing_keeps_real_changes(tmp_path, old, new, scenario):
    data = _setup(tmp_path)
    before = (
        (0, 100) if scenario == "multiplicity" else (0, 200) if scenario == "removed" else (0,)
    )
    baseline = _baseline(tmp_path, data, old, before)
    if scenario == "reclassified":
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(new))) as conn, conn:
            conn.execute(
                "UPDATE clinvar_variants SET clinical_significance='Benign', aggregate='Benign'"
            )
    after = (0, 200) if scenario == "added" else (0,)
    out = tmp_path / "current.json"
    _run(_source(tmp_path, "new", new, after), data, out, new, baseline)
    diff = json.loads(out.read_text())["diff"]
    expected = {
        "reclassified": (0, 0, 1),
        "added": (1, 0, 0),
        "removed": (0, 1, 0),
        "multiplicity": (0, 1, 0),
    }[scenario]
    assert tuple(len(diff[k]) for k in ("new", "removed", "changed")) == expected, diff


@pytest.mark.parametrize("build", list(POSITIONS))
@pytest.mark.parametrize("baseline_build", ["same", "absent", "unknown"])
def test_same_or_unknown_build_does_not_erase_known_site_identity(tmp_path, build, baseline_build):
    data = _setup(tmp_path)
    baseline = _baseline(tmp_path, data, build)
    payload = json.loads(baseline.read_text())
    if baseline_build == "absent":
        payload["input"].pop("build")
    elif baseline_build == "unknown":
        payload["input"]["build"] = "unknown"
    baseline.write_text(json.dumps(payload))
    out = tmp_path / "current.json"
    _run(_source(tmp_path, "new", build, (100,)), data, out, build, baseline)
    diff = json.loads(out.read_text())["diff"]
    assert tuple(len(diff[k]) for k in ("new", "removed", "changed")) == (1, 1, 0)

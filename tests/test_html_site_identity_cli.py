"""Outside reviewer A's lost-frequency case through offline analyze, both builds."""

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_homref_ref_by_locus import _cv, _vcf, _write_clinvar, _write_pharm
from test_html_variant_data_script import _assert_page_intact

from allelix.cli import main
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA

RS = "rs900529101"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("pharmgkb", ["on", "off", "missing"])
def test_cli_keeps_each_record_and_its_findings_in_its_own_panel(
    tmp_path, build, reverse, pharmgkb
):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\nsnpedia=false\ngwas=false\n"
        f"pharmgkb={'false' if pharmgkb == 'off' else 'true'}\n"
        f"gnomad={'true' if build == 'GRCh38' else 'false'}\n"
        "alphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    _write_clinvar(
        data,
        build,
        [_cv(RS, "1", p, "G", "A", allele_id=i) for i, p in enumerate((1000, 5000), 1)],
    )
    if pharmgkb != "missing":
        _write_pharm(data, [(RS, "AG", "synthetic-1")])
    if build == "GRCh38":
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn, conn:
            conn.executescript(GNOMAD_SCHEMA)
            conn.executemany(
                "INSERT INTO gnomad_frequencies(chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
                [("1", 1000, "G", "A", RS, 0.4), ("1", 5000, "G", "A", RS, 0.000001)],
            )
            conn.execute(
                "INSERT INTO database_versions(name,source_url,version,downloaded_at,"
                "record_count,local_version_tag) VALUES ('gnomad','test://synthetic',"
                "'synthetic','2001-01-01',2,?)",
                (f"sv:{GNOMAD_SCHEMA_VERSION}",),
            )
    positions = (5000, 1000) if reverse else (1000, 5000)
    source = _vcf(
        tmp_path / "input.vcf",
        [("1", p, RS, "G", "A", 50, "PASS", ".", "GT:GQ", "0/1:99") for p in positions],
        build=build,
    )
    outputs = {}
    for extension in ("json", "html"):
        output = tmp_path / f"report.{extension}"
        invocation = CliRunner().invoke(
            main,
            [
                "analyze",
                str(source),
                "--data-dir",
                str(data),
                "--no-update",
                "--output",
                str(output),
                "--min-magnitude",
                "0",
            ],
        )
        assert invocation.exit_code == 0, (invocation.output, invocation.exception)
        outputs[extension] = output.read_text()
    annotations = json.loads(outputs["json"])["annotations"]
    expected_sources = {"clinvar", "pharmgkb"} if pharmgkb == "on" else {"clinvar"}
    # Prove the fixture produced both loci and the requested sources before HTML.
    assert {(a["position"], a["source"]) for a in annotations} == {
        (p, s) for p in (1000, 5000) for s in expected_sources
    }
    entries = _assert_page_intact(outputs["html"])
    assert len(entries) == 2, "HTML merged two real CLI findings with one rsID/genotype"
    panels = {e["locus"]: e for e in entries}
    assert set(panels) == {"1:1000", "1:5000"}
    for pos in (1000, 5000):
        entry = panels[f"1:{pos}"]
        assert {a["source"] for a in entry["annotations"]} == expected_sources
        assert len(entry["annotations"]) == len(expected_sources)
        if build == "GRCh38":
            assert entry["allele_frequency"] == (0.4 if pos == 1000 else 0.000001)

"""Re-reading an annotated VCF preserves findings, including per-allele records."""

import gzip
import json
from pathlib import Path

import pytest
from click.testing import CliRunner
from split_site_contract_helpers import (
    BUILDS,
    OTHER,
    RSID,
    data_dir,
    record,
    run,
    signature,
    write_vcf,
)

from allelix.cli import main

SPLIT_SHAPES = (
    (record("A", "1/0", rsid="."), record("T", "0/1", rsid=".")),
    (record("A", "1/0", rsid=RSID), record("T", "0/1", rsid=".")),
    (record("A", "1/0", rsid=RSID), record("T", "0/1", rsid=OTHER)),
    (record("A", "1/."), record("T", "./1")),
    (record("A", "1|.", ps="32"), record("T", ".|1", ps="32")),
    (record("A", "0/1"), record("T", "0/1")),
    (record("A", "0/1"), record("T", "0/1", filt="q10")),
    (record("A", "1/1"), record("T", "0/1")),
)


def _public_state(report):
    return (
        signature(report),
        report.get("report_warnings", []),
        report.get("high_value_unchecked", []),
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("auto", (False, True), ids=("explicit", "auto"))
@pytest.mark.parametrize("rows", SPLIT_SHAPES)
def test_split_input_and_annotated_output_have_identical_public_findings(
    tmp_path, build, auto, rows
):
    data = data_dir(tmp_path / "data", identity=True, alts="AT")
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    annotated = tmp_path / "annotated.vcf"
    first, _ = run(
        data,
        source,
        tmp_path / "first.json",
        build,
        auto=auto,
        extra=("--vcf-out", str(annotated)),
    )
    second, _ = run(data, annotated, tmp_path / "second.json", build, auto=auto)
    assert _public_state(second) == _public_state(first), (
        "A report must not change its own meaning on re-read"
    )
    before = [
        line.split("\t") for line in source.read_text().splitlines() if not line.startswith("#")
    ]
    after = [
        line.split("\t") for line in annotated.read_text().splitlines() if not line.startswith("#")
    ]
    assert len(before) == len(after)
    for old, new in zip(before, after, strict=True):
        assert old[:2] == new[:2] and old[3:7] == new[3:7] and old[8:] == new[8:]
    if rows == SPLIT_SHAPES[0]:
        assert len([r for r in first["annotations"] if r["source"] == "clinvar"]) == 2
        for row, alt in zip(after, ("A", "T"), strict=True):
            assert row[4] == alt
            assert "ALLELIX_CLINVAR=" in row[7], "Each allele's finding stays on its own record"


FIXTURES = Path(__file__).parent / "fixtures"
VCFS = sorted(p for p in FIXTURES.rglob("*") if p.name.endswith((".vcf", ".vcf.gz", ".vcf.bgz")))


@pytest.mark.parametrize("fixture", VCFS, ids=lambda p: str(p.relative_to(FIXTURES)))
def test_every_committed_vcf_fixture_keeps_its_findings_after_annotating(
    tmp_path, all_annotators_data_dir, fixture
):
    data = all_annotators_data_dir
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=true\ngwas=true\nsnpedia=false\ngnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=true\n"
    )
    # Both analyses select the same person in multi-sample fixtures.
    opener = gzip.open if fixture.suffix in (".gz", ".bgz") else open
    with opener(fixture, "rt", encoding="utf-8-sig", errors="replace") as stream:
        samples = next(
            line.rstrip().split("\t")[9:] for line in stream if line.startswith("#CHROM\t")
        )
    opts = ["--sample", samples[0]] if samples else []
    # Fixtures with no build signal require the existing explicit override; the dedicated
    # split matrix above separately tests header-driven auto detection on both builds.
    build = (
        "grch37"
        if fixture.name
        in (
            "mock_vcf_non_utf8_header.vcf",
            "mock_vcf_non_utf8_header.vcf.gz",
            "mock_vcf_rsidless.vcf",
            "mock_vcf_mixed_contigs.vcf",
            "mock_vcf_no_signal.vcf",
            "sequencing_com_sv.vcf",
        )
        else "grch38"
    )
    common = [
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build,
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--include-benign",
        *opts,
    ]
    first, second, emitted = (
        tmp_path / n for n in ("fixture-first.json", "fixture-second.json", "fixture-output.vcf")
    )
    one = CliRunner().invoke(
        main, ["analyze", str(fixture), *common, "--output", str(first), "--vcf-out", str(emitted)]
    )
    if fixture.name == "sequencing_com_sv.vcf":
        assert one.exit_code == 1 and "structural" in one.output.lower()
        assert not emitted.exists() and not first.exists()
        return  # Deliberately unsupported input; there is no emitted VCF to re-read.
    assert one.exit_code == 0, (fixture.name, one.output, repr(one.exception))
    two = CliRunner().invoke(main, ["analyze", str(emitted), *common, "--output", str(second)])
    assert two.exit_code == 0, (fixture.name, two.output, repr(two.exception))
    assert _public_state(json.loads(first.read_text())) == _public_state(
        json.loads(second.read_text())
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "ids",
    [
        (".", "."),
        ("22:11924846:G:C", "22:11924847:G:T"),
        ("COSV123", "RCV456"),
        (RSID, "."),
        (".", RSID),
        ("COSV123", RSID),
        (RSID, "22:11924847:G:T"),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_merged_pair_stamps_the_rsid_and_each_original_id_on_both_records(
    tmp_path, build, ids, reverse
):
    data = data_dir(tmp_path / "data")
    rows = [record("A", "1/0", rsid=ids[0]), record("C", "0/1", rsid=ids[1])]
    if reverse:
        rows.reverse()
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    output = tmp_path / "annotated.vcf"
    first, _ = run(data, source, tmp_path / "first.json", build, extra=("--vcf-out", str(output)))
    emitted = [
        line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")
    ]
    assert len(emitted) == 2
    for original, row in zip(rows, emitted, strict=True):
        assert row[2] == RSID, "Both component records must carry the recovered/sole rsID"
        info = dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
        assert info.get("ALLELIX_ORIGINAL_ID") == original["rsid"]
        assert row[4] == original["alt"]
        assert "ALLELIX_CLINVAR" in info
    assert {a["genotype_match"] for a in first["annotations"]} == {"AC"}
    assert {a["source"] for a in first["annotations"]} == {
        "clinvar",
        "gwas",
        "pharmgkb",
        "snpedia",
    }
    second, _ = run(data, output, tmp_path / "second.json", build)
    assert _public_state(first) == _public_state(second)

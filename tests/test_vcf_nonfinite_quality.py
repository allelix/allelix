"""Unusable numeric quality fields behave as absent, not as zero or infinity."""

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import _ready_cache

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.parsers.vcf import VcfParser

BAD = ("inf", "-inf", "+Infinity", "1e999", "-1e999")
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length={length}>\n"
    '##ALT=<ID=NON_REF,Description="Any other allele">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def _write(path, field, value, *, gt="0/1", block=False, other=None, build="GRCh37"):
    values = {"GQ": "99", "DP": "30"}
    values[field] = value
    if other is not None:
        values["DP" if field == "GQ" else "GQ"] = other
    alt, info = ("<NON_REF>", "END=1010") if block else ("G", ".")
    length = 249250621 if build == "GRCh37" else 248956422
    path.write_text(
        HEADER.format(length=length)
        + f"1\t1000\trs900252401\tA\t{alt}\t.\tPASS\t{info}\tGT:GQ:DP\t"
        + f"{gt}:{values['GQ']}:{values['DP']}\n"
    )
    return path


def _read(path, **options):
    parser = VcfParser(**options)
    parser.watch({("1", 1000), ("1", 1005)}, {"rs900252401"}, rows=True)
    parser.watch_reference_targets({"1": [1000]})
    calls = list(parser.parse(path))
    stats = parser.quality_stats
    return (
        [(v.rsid, v.allele1, v.allele2) for v in calls],
        stats.locus_states,
        stats.reference_calls,
        list(stats.reference_rows.values()),
        stats.low_quality_calls,
        stats.half_calls,
        parser.reference_target_calls,
    )


@pytest.mark.parametrize("field", ("GQ", "DP"))
@pytest.mark.parametrize("value", BAD)
@pytest.mark.parametrize(
    "gt,block", (("0/1", False), ("0/0", False), ("0/.", False), ("0/0", True))
)
def test_nonfinite_is_missing_across_calls_watchers_and_reference_targets(
    tmp_path, field, value, gt, block
):
    expected = _read(_write(tmp_path / "missing.vcf", field, ".", gt=gt, block=block))
    observed = _read(_write(tmp_path / "invalid.vcf", field, value, gt=gt, block=block))
    assert observed == expected
    if gt == "0/1":
        assert observed[0] == [("rs900252401", "A", "G")]


@pytest.mark.parametrize("field", ("GQ", "DP"))
@pytest.mark.parametrize("value", (*BAD, "nan", "bad", "."))
def test_other_usable_field_still_rejects_low_confidence(tmp_path, field, value):
    # An unusable GQ cannot override DP=0; an unusable DP cannot override GQ=5.
    other = "0" if field == "GQ" else "5"
    result = _read(_write(tmp_path / "in.vcf", field, value, other=other))
    assert result[0] == [("rs900252401", "-", "-")]
    assert result[4] == 1


@pytest.mark.parametrize("field", ("GQ", "DP"))
@pytest.mark.parametrize("value", ("nan", "bad", ".", "30.9"))
def test_existing_unusable_and_finite_decimal_rules_stay(tmp_path, field, value):
    result = _read(_write(tmp_path / "in.vcf", field, value))
    assert result[0] == [("rs900252401", "A", "G")]
    assert result[4] == 0


@pytest.mark.parametrize("field", ("GQ", "DP"))
@pytest.mark.parametrize("selected", ("SAMPLE", "OTHER"))
def test_only_the_selected_samples_fields_are_read(tmp_path, field, selected):
    path = _write(tmp_path / "two.vcf", field, "inf")
    path.write_text(
        path.read_text().replace("\tSAMPLE\n", "\tSAMPLE\tOTHER\n").rstrip() + "\t0/1:99:30\n"
    )
    result = _read(path, sample=selected)
    assert result[0] == [("rs900252401", "A", "G")]


@pytest.mark.parametrize("field", ("GQ", "DP"))
@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize("explicit", (False, True), ids=("auto", "explicit"))
@pytest.mark.parametrize("command", ("analyze", "stats", "extract", "plink"))
def test_cli_commands_do_not_crash_and_keep_the_call(tmp_path, field, build, explicit, command):
    data = tmp_path / "data"
    _ready_cache(data)
    if build == "GRCh38":
        cache = data / clinvar_db_filename(build)
        (data / clinvar_db_filename("GRCh37")).rename(cache)
        with contextlib.closing(sqlite3.connect(cache)) as conn, conn:
            conn.execute("UPDATE database_versions SET name=?", (clinvar_record_name(build),))
    source = _write(tmp_path / "in.vcf", field, "inf", build=build)
    args = ["export", "plink"] if command == "plink" else [command]
    args += [str(source), "--data-dir", str(data)]
    if explicit:
        args += ["--build", build.lower()]
    output = tmp_path / "report.json"
    if command == "analyze":
        args += ["--no-update", "--output", str(output)]
    elif command == "extract":
        args += ["--snps", "rs900252401"]
    elif command == "plink":
        args += ["--output-prefix", str(tmp_path / "export")]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    if command == "analyze":
        report = json.loads(output.read_text())
        assert report["input"]["build"] == build
        assert [a["rsid"] for a in report["annotations"]] == ["rs900252401"]
    elif command == "extract":
        assert "rs900252401" in result.output and "A/G" in result.output
    elif command == "plink":
        assert "rs900252401" in (tmp_path / "export.bim").read_text()
    else:
        assert build in result.output

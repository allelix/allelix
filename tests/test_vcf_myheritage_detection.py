"""VCF structure wins over a vendor name inside VCF provenance text."""

import contextlib
import gzip
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_output_destinations_and_modes import VCF, _ready_cache

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.parsers import detect_parser

MENTIONS = (
    "##source=bcftools convert MyHeritage_raw_dna_data.csv",
    '##INFO=<ID=ORIGIN,Number=1,Type=String,Description="Converted MyHeritage input">',
    '##SAMPLE=<ID=SAMPLE,Description="MyHeritage sample">',
    "##reference=synthetic-MyHeritage-reference.fa",
)


def _vcf(path, mention, *, build="GRCh37", gvcf=False):
    length = 249250621 if build == "GRCh37" else 248956422
    text = VCF.replace(
        "##fileformat=VCFv4.2\n",
        f"##fileformat=VCFv4.2\n{mention}\n##contig=<ID=1,length={length}>\n",
    )
    if gvcf:
        text = text.replace("\tA\tG\t", "\tA\tG,<NON_REF>\t")
        text += "1\t2000\t.\tA\t<NON_REF>\t.\tPASS\tEND=2100\tGT\t0/0\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("mention", MENTIONS)
@pytest.mark.parametrize("suffix", (".vcf", ".txt", ".vcf.gz"))
@pytest.mark.parametrize("gvcf", (False, True))
def test_header_provenance_cannot_relabel_vcf_as_myheritage(tmp_path, mention, suffix, gvcf):
    source = _vcf(tmp_path / ("MyHeritage-converted" + suffix), mention, gvcf=gvcf)
    parser = detect_parser(source)
    assert parser.name == "vcf"
    calls = list(parser.parse(source))
    assert len(calls) == 1 and calls[0].rsid == "rs900252401"
    assert calls[0].ref == "A" and {calls[0].allele1, calls[0].allele2} == {"A", "G"}


@pytest.mark.parametrize("build", ("GRCh37", "GRCh38"))
@pytest.mark.parametrize("explicit", (False, True), ids=("auto", "explicit"))
@pytest.mark.parametrize("command", ("stats", "extract", "analyze"))
def test_commands_read_the_vcf_rows_and_decide_its_build(tmp_path, build, explicit, command):
    data = tmp_path / "data"
    _ready_cache(data)
    if build == "GRCh38":
        target = data / clinvar_db_filename(build)
        (data / clinvar_db_filename("GRCh37")).rename(target)
        with contextlib.closing(sqlite3.connect(target)) as conn, conn:
            conn.execute("UPDATE database_versions SET name=?", (clinvar_record_name(build),))
    source = _vcf(tmp_path / "converted.vcf", MENTIONS[0], build=build)
    args = [command, str(source), "--data-dir", str(data)]
    if explicit:
        args += ["--build", build.lower()]
    if command == "extract":
        args += ["--snps", "rs900252401"]
    if command == "analyze":
        output = tmp_path / "report.json"
        args += ["--no-update", "--output", str(output)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    if command == "analyze":
        report = json.loads(output.read_text())
        assert report["input"]["format"] == "vcf" and report["input"]["build"] == build
        assert [row["rsid"] for row in report["annotations"]] == ["rs900252401"]
    else:
        assert build in result.output
        assert "VCF" in result.output if command == "stats" else "rs900252401" in result.output


@pytest.mark.parametrize("separator", (",", "\t"))
@pytest.mark.parametrize("compressed", (False, True))
@pytest.mark.parametrize(
    "preamble", ("# MyHeritage DNA\n", '# MyHeritage DNA\n# exported from "sample.vcf"\n')
)
def test_real_myheritage_layout_stays_myheritage(tmp_path, separator, compressed, preamble):
    text = preamble + separator.join(("RSID", "CHROMOSOME", "POSITION", "RESULT")) + "\n"
    text += separator.join(("rs900252401", "1", "1000", "AG")) + "\n"
    path = tmp_path / ("input.csv.gz" if compressed else "input.csv")
    if compressed:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    parser = detect_parser(path)
    assert parser.name == "myheritage"
    calls = list(parser.parse(path))
    assert len(calls) == 1 and calls[0].rsid == "rs900252401"


def test_the_same_csv_without_myheritage_signature_stays_ftdna(tmp_path):
    source = tmp_path / "plain.csv"
    source.write_text("RSID,CHROMOSOME,POSITION,RESULT\nrs900252401,1,1000,AG\n")
    assert detect_parser(source).name == "ftdna"

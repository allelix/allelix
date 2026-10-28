"""VCF POS is numeric: leading zeroes cannot turn reference coverage into a no-call."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import main
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import CLINVAR_SCHEMA


@pytest.mark.parametrize(
    "build,pos,length", [("GRCh37", 1000, 249250621), ("GRCh38", 2000, 248956422)]
)
@pytest.mark.parametrize("command", ["stats", "analyze"])
@pytest.mark.parametrize("padded", [False, True])
def test_explicit_reference_coverage_uses_numeric_pos(
    tmp_path, build, pos, length, command, padded
):
    data = tmp_path / "data"
    data.mkdir()
    write_identity_db(
        [("rs6025", "1", 1000, "A", "G", 2000, "A", "G", 0)],
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://synthetic",
        downloaded_at="2001-01-01",
    )
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.execute(
            "INSERT INTO clinvar_variants (rsid,chromosome,position,ref,alt,"
            "clinical_significance,condition,gene,allele_id) "
            "VALUES ('rs6025','1',?,'A','G','Pathogenic',"
            "'Synthetic condition','F5',900252999)",
            (pos,),
        )
        conn.execute(
            "INSERT INTO database_versions (name,source_url,version,downloaded_at,"
            "record_count,local_version_tag) "
            "VALUES (?,'test://synthetic','synthetic','2001-01-01',1,?)",
            (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        + f"##contig=<ID=1,length={length}>\n"
        + '##ALT=<ID=NON_REF,Description="Any non-reference allele">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + f"1\t{('0' if padded else '') + str(pos)}\trs6025\tA\tG\t.\tPASS\t.\t"
        "GT:GQ:DP\t0/0:50:20\n"
        + "1\t3000\t.\tA\t<NON_REF>\t.\tPASS\tEND=3005\tGT:GQ:DP\t0/0:50:20\n"
    )
    args = [command, str(source), "--data-dir", str(data), "--build", build.lower()]
    output = tmp_path / "report.json"
    if command == "analyze":
        args += ["--no-update", "--output", str(output)]
    result = CliRunner(env={"COLUMNS": "500"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    if command == "analyze":
        warnings = json.loads(output.read_text()).get("high_value_no_calls", [])
        assert not any(w["rsid"] == "rs6025" for w in warnings), warnings
    else:
        assert "rs6025" not in result.output, result.output
